#include "solo/MoleculeFirstResolver.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace molecule_first {
namespace {

using GroupKey = std::pair<std::string, std::string>;
using CorrectionKey = std::tuple<std::string, std::string, std::string>;

std::uint64_t fnv1a(const std::string &value, std::uint64_t seed)
{
    std::uint64_t hash = seed;
    for (unsigned char byte : value) {
        hash ^= static_cast<std::uint64_t>(byte);
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::string stableId(const std::string &prefix, const std::vector<std::string> &parts)
{
    std::string joined;
    for (const std::string &part : parts) {
        joined.append(part);
        joined.push_back('\x1f');
    }
    const std::uint64_t left = fnv1a(joined, 1469598103934665603ULL);
    const std::uint64_t right = fnv1a(joined, 1099511628211ULL);
    std::ostringstream out;
    out << prefix << '_' << std::hex << std::setfill('0')
        << std::setw(16) << left << std::setw(16) << right;
    return out.str();
}

void validateConfig(const Config &config)
{
    if (!std::isfinite(config.temperature) || config.temperature <= 0.0) {
        throw std::invalid_argument("temperature must be finite and positive");
    }
    if (!std::isfinite(config.priorAlpha) || config.priorAlpha <= 0.0) {
        throw std::invalid_argument("prior alpha must be finite and positive");
    }
    if (!std::isfinite(config.priorBeta) || config.priorBeta < 0.0) {
        throw std::invalid_argument("prior beta must be finite and non-negative");
    }
    if (!std::isfinite(config.gateMinPosterior) || config.gateMinPosterior < 0.0
        || config.gateMinPosterior > 1.0 || !std::isfinite(config.gateMinMargin)
        || config.gateMinMargin < 0.0 || config.gateMinMargin > 1.0) {
        throw std::invalid_argument("hard-call gates must lie in [0,1]");
    }
}

std::vector<std::string> candidateNames(const CandidateRead &read)
{
    std::vector<std::string> names;
    names.reserve(read.scores.size());
    for (const CandidateScore &score : read.scores) {
        names.push_back(score.candidate);
    }
    return names;
}

std::vector<std::string> intersection(const std::vector<std::string> &left,
                                      const std::vector<std::string> &right)
{
    std::vector<std::string> output;
    std::set_intersection(left.begin(), left.end(), right.begin(), right.end(),
                          std::back_inserter(output));
    return output;
}

double readScore(const CandidateRead &read, const std::string &candidate)
{
    const auto found = std::lower_bound(
        read.scores.begin(), read.scores.end(), candidate,
        [](const CandidateScore &score, const std::string &name) {
            return score.candidate < name;
        });
    if (found == read.scores.end() || found->candidate != candidate) {
        throw std::logic_error("clique intersection lost candidate score");
    }
    return found->logLikelihood;
}

std::vector<double> normalizeLogWeights(const std::vector<double> &values)
{
    if (values.empty()) {
        throw std::invalid_argument("cannot normalize empty evidence");
    }
    const double maximum = *std::max_element(values.begin(), values.end());
    std::vector<double> output(values.size(), 0.0);
    double total = 0.0;
    for (std::size_t index = 0; index < values.size(); ++index) {
        output[index] = std::exp(values[index] - maximum);
        total += output[index];
    }
    if (!std::isfinite(total) || total <= 0.0) {
        throw std::runtime_error("posterior normalization failed");
    }
    for (double &value : output) {
        value /= total;
    }
    return output;
}

bool hammingOne(const std::string &left, const std::string &right)
{
    if (left.size() != right.size()) {
        return false;
    }
    std::size_t mismatches = 0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        if (left[index] != right[index] && ++mismatches > 1) {
            return false;
        }
    }
    return mismatches == 1;
}

bool nearlyEqual(double left, double right)
{
    const double absoluteTolerance = 1e-12;
    const double relativeTolerance = 1e-10;
    return std::fabs(left - right) <= absoluteTolerance
        + relativeTolerance * std::max(std::fabs(left), std::fabs(right));
}

bool greaterEqualWithinTolerance(double left, double right)
{
    return left > right || nearlyEqual(left, right);
}

} // namespace

std::vector<ReadClique> buildReadCliques(const std::vector<CandidateRead> &inputReads,
                                         const PriorCounts &priorCounts,
                                         const Config &config)
{
    validateConfig(config);
    std::vector<CandidateRead> reads = inputReads;
    std::set<std::string> seenReadIds;
    std::map<GroupKey, std::vector<std::size_t> > grouped;

    for (std::size_t index = 0; index < reads.size(); ++index) {
        CandidateRead &read = reads[index];
        if (read.readId.empty() || read.featureId.empty() || read.rawUmi.empty()) {
            throw std::invalid_argument("read_id, feature_id, and raw_umi are required");
        }
        if (!seenReadIds.insert(read.readId).second) {
            throw std::invalid_argument("duplicate read_id: " + read.readId);
        }
        if (read.scores.empty()) {
            throw std::invalid_argument("read has no candidate: " + read.readId);
        }
        std::sort(read.scores.begin(), read.scores.end(),
                  [](const CandidateScore &left, const CandidateScore &right) {
                      return left.candidate < right.candidate;
                  });
        for (std::size_t scoreIndex = 0; scoreIndex < read.scores.size(); ++scoreIndex) {
            const CandidateScore &score = read.scores[scoreIndex];
            if (score.candidate.empty() || !std::isfinite(score.logLikelihood)) {
                throw std::invalid_argument("candidate and finite likelihood are required");
            }
            if (scoreIndex > 0 && read.scores[scoreIndex - 1].candidate == score.candidate) {
                throw std::invalid_argument("duplicate read/candidate: " + read.readId + "/" + score.candidate);
            }
            if (priorCounts.find(score.candidate) == priorCounts.end()) {
                throw std::invalid_argument("missing exact-read count for candidate: " + score.candidate);
            }
        }
        grouped[GroupKey(read.featureId, read.rawUmi)].push_back(index);
    }

    struct Partition {
        std::vector<std::string> candidates;
        std::vector<std::size_t> members;
    };

    std::vector<ReadClique> output;
    for (const auto &group : grouped) {
        std::vector<std::size_t> order = group.second;
        std::sort(order.begin(), order.end(), [&](std::size_t left, std::size_t right) {
            if (reads[left].readId != reads[right].readId) {
                return reads[left].readId < reads[right].readId;
            }
            return candidateNames(reads[left]) < candidateNames(reads[right]);
        });

        std::vector<Partition> partitions;
        for (std::size_t readIndex : order) {
            const std::vector<std::string> names = candidateNames(reads[readIndex]);
            bool haveChoice = false;
            std::vector<std::string> bestIntersection;
            std::size_t bestPartition = 0;
            for (std::size_t partitionIndex = 0; partitionIndex < partitions.size(); ++partitionIndex) {
                std::vector<std::string> overlap = intersection(partitions[partitionIndex].candidates, names);
                if (overlap.empty()) {
                    continue;
                }
                if (!haveChoice || overlap < bestIntersection
                    || (overlap == bestIntersection && partitionIndex < bestPartition)) {
                    haveChoice = true;
                    bestIntersection.swap(overlap);
                    bestPartition = partitionIndex;
                }
            }
            if (haveChoice) {
                partitions[bestPartition].candidates = bestIntersection;
                partitions[bestPartition].members.push_back(readIndex);
            } else {
                Partition partition;
                partition.candidates = names;
                partition.members.push_back(readIndex);
                partitions.push_back(partition);
            }
        }

        for (const Partition &partition : partitions) {
            ReadClique clique;
            clique.featureId = group.first.first;
            clique.rawUmi = group.first.second;
            clique.candidates = partition.candidates;
            for (std::size_t readIndex : partition.members) {
                clique.memberReadIds.push_back(reads[readIndex].readId);
            }
            std::sort(clique.memberReadIds.begin(), clique.memberReadIds.end());

            for (const std::string &candidate : clique.candidates) {
                double likelihood = 0.0;
                for (std::size_t readIndex : partition.members) {
                    likelihood += readScore(reads[readIndex], candidate);
                }
                const double prior = std::log(static_cast<double>(priorCounts.at(candidate))
                                              + config.priorAlpha);
                clique.logLikelihoodSums.push_back(likelihood);
                clique.logReadPriors.push_back(prior);
                clique.logEvidence.push_back(likelihood / config.temperature
                                             + config.priorBeta * prior);
            }
            clique.posterior = normalizeLogWeights(clique.logEvidence);
            std::vector<std::string> idParts;
            idParts.push_back(clique.featureId);
            idParts.push_back(clique.rawUmi);
            idParts.insert(idParts.end(), clique.memberReadIds.begin(), clique.memberReadIds.end());
            idParts.insert(idParts.end(), clique.candidates.begin(), clique.candidates.end());
            clique.cliqueId = stableId("clq", idParts);
            output.push_back(clique);
        }
    }
    std::sort(output.begin(), output.end(), [](const ReadClique &left, const ReadClique &right) {
        return left.cliqueId < right.cliqueId;
    });
    return output;
}

UmiCorrections correctedUmis(const std::vector<ReadClique> &cliques,
                             const std::string &umiMode)
{
    if (umiMode != "exact" && umiMode != "1mm_cr") {
        throw std::invalid_argument("UMI mode must be exact or 1mm_cr");
    }
    std::map<GroupKey, std::map<std::string, double> > support;
    for (const ReadClique &clique : cliques) {
        for (std::size_t index = 0; index < clique.candidates.size(); ++index) {
            support[GroupKey(clique.featureId, clique.candidates[index])][clique.rawUmi]
                += clique.posterior[index];
        }
    }

    UmiCorrections output;
    for (const auto &group : support) {
        std::vector<std::string> ordered;
        for (const auto &item : group.second) {
            ordered.push_back(item.first);
        }
        std::sort(ordered.begin(), ordered.end(), [&](const std::string &left, const std::string &right) {
            const double leftSupport = group.second.at(left);
            const double rightSupport = group.second.at(right);
            return !nearlyEqual(leftSupport, rightSupport) ? leftSupport > rightSupport : left < right;
        });

        std::map<std::string, std::string> roots;
        for (std::size_t index = 0; index < ordered.size(); ++index) {
            const std::string &umi = ordered[index];
            std::string root = umi;
            if (umiMode == "1mm_cr") {
                std::vector<std::string> eligible;
                for (std::size_t parentIndex = 0; parentIndex < index; ++parentIndex) {
                    const std::string &parent = ordered[parentIndex];
                    if (hammingOne(umi, parent)
                        && greaterEqualWithinTolerance(
                            group.second.at(parent), 2.0 * group.second.at(umi) - 1.0)) {
                        eligible.push_back(parent);
                    }
                }
                if (!eligible.empty()) {
                    std::sort(eligible.begin(), eligible.end(), [&](const std::string &left,
                                                                   const std::string &right) {
                        const double leftSupport = group.second.at(left);
                        const double rightSupport = group.second.at(right);
                        return !nearlyEqual(leftSupport, rightSupport) ? leftSupport > rightSupport : left < right;
                    });
                    root = roots.at(eligible.front());
                }
            }
            roots[umi] = root;
            output[CorrectionKey(group.first.first, group.first.second, umi)] = root;
        }
    }
    return output;
}

std::vector<Occupancy> weightedOccupancies(const std::vector<ReadClique> &cliques,
                                           const std::string &umiMode)
{
    const UmiCorrections corrections = correctedUmis(cliques, umiMode);
    using OccupancyKey = std::tuple<std::string, std::string, std::string>;
    std::map<OccupancyKey, std::vector<std::pair<std::string, double> > > grouped;
    for (const ReadClique &clique : cliques) {
        for (std::size_t index = 0; index < clique.candidates.size(); ++index) {
            const std::string &candidate = clique.candidates[index];
            const std::string corrected = corrections.at(
                CorrectionKey(clique.featureId, candidate, clique.rawUmi));
            grouped[OccupancyKey(clique.featureId, corrected, candidate)].push_back(
                std::make_pair(clique.cliqueId, clique.posterior[index]));
        }
    }

    std::vector<Occupancy> output;
    for (auto &group : grouped) {
        Occupancy row;
        row.umiMode = umiMode;
        row.featureId = std::get<0>(group.first);
        row.correctedUmi = std::get<1>(group.first);
        row.candidate = std::get<2>(group.first);
        double absent = 1.0;
        for (const auto &value : group.second) {
            absent *= 1.0 - value.second;
            row.readCliqueIds.push_back(value.first);
        }
        row.expectedCount = 1.0 - absent;
        std::sort(row.readCliqueIds.begin(), row.readCliqueIds.end());
        output.push_back(row);
    }
    return output;
}

std::pair<std::string, double> topCandidate(const ReadClique &clique)
{
    if (clique.candidates.empty() || clique.posterior.size() != clique.candidates.size()) {
        throw std::invalid_argument("invalid clique posterior");
    }
    std::size_t best = 0;
    for (std::size_t index = 1; index < clique.candidates.size(); ++index) {
        if ((!nearlyEqual(clique.posterior[index], clique.posterior[best])
             && clique.posterior[index] > clique.posterior[best])
            || (nearlyEqual(clique.posterior[index], clique.posterior[best])
                && clique.candidates[index] < clique.candidates[best])) {
            best = index;
        }
    }
    return std::make_pair(clique.candidates[best], clique.posterior[best]);
}

std::vector<HardCall> gatedHardCalls(const std::vector<ReadClique> &cliques,
                                     const Config &config)
{
    validateConfig(config);
    std::vector<HardCall> output;
    for (const ReadClique &clique : cliques) {
        const std::pair<std::string, double> top = topCandidate(clique);
        double second = 0.0;
        for (std::size_t index = 0; index < clique.candidates.size(); ++index) {
            if (clique.candidates[index] != top.first && clique.posterior[index] > second) {
                second = clique.posterior[index];
            }
        }
        HardCall call;
        call.cliqueId = clique.cliqueId;
        call.posterior = top.second;
        call.margin = nearlyEqual(top.second, second) ? 0.0 : top.second - second;
        call.assigned = call.posterior >= config.gateMinPosterior
            && call.margin >= config.gateMinMargin;
        call.candidate = call.assigned ? top.first : "";
        call.reason = call.assigned ? "posterior_and_margin" : "gate_failed";
        output.push_back(call);
    }
    return output;
}

std::vector<Molecule> policyMolecules(const std::vector<ReadClique> &cliques,
                                      const std::string &umiMode,
                                      const std::string &product,
                                      const std::vector<HardCall> &hardCalls)
{
    if (product != "strict" && product != "hard" && product != "gated_hard") {
        throw std::invalid_argument("unknown molecule product: " + product);
    }
    const UmiCorrections corrections = correctedUmis(cliques, umiMode);
    std::map<std::string, HardCall> callByClique;
    for (const HardCall &call : hardCalls) {
        callByClique[call.cliqueId] = call;
    }

    using MoleculeKey = std::tuple<std::string, std::string, std::string>;
    std::map<MoleculeKey, std::vector<const ReadClique *> > grouped;
    for (const ReadClique &clique : cliques) {
        std::string candidate;
        if (product == "strict") {
            if (clique.candidates.size() != 1) {
                continue;
            }
            candidate = clique.candidates.front();
        } else if (product == "hard") {
            candidate = topCandidate(clique).first;
        } else {
            const auto call = callByClique.find(clique.cliqueId);
            if (call == callByClique.end() || !call->second.assigned) {
                continue;
            }
            candidate = call->second.candidate;
        }
        const std::string corrected = corrections.at(
            CorrectionKey(clique.featureId, candidate, clique.rawUmi));
        grouped[MoleculeKey(clique.featureId, corrected, candidate)].push_back(&clique);
    }

    std::vector<Molecule> output;
    for (const auto &group : grouped) {
        Molecule molecule;
        molecule.umiMode = umiMode;
        molecule.product = product;
        molecule.featureId = std::get<0>(group.first);
        molecule.correctedUmi = std::get<1>(group.first);
        molecule.candidate = std::get<2>(group.first);
        std::set<std::string> members;
        for (const ReadClique *clique : group.second) {
            members.insert(clique->memberReadIds.begin(), clique->memberReadIds.end());
            molecule.readCliqueIds.push_back(clique->cliqueId);
        }
        molecule.memberReadIds.assign(members.begin(), members.end());
        std::sort(molecule.readCliqueIds.begin(), molecule.readCliqueIds.end());
        std::vector<std::string> idParts;
        idParts.push_back(umiMode);
        idParts.push_back(product);
        idParts.push_back(molecule.featureId);
        idParts.push_back(molecule.correctedUmi);
        idParts.push_back(molecule.candidate);
        idParts.insert(idParts.end(), molecule.memberReadIds.begin(), molecule.memberReadIds.end());
        molecule.moleculeId = stableId("mol", idParts);
        output.push_back(molecule);
    }
    return output;
}

} // namespace molecule_first
