#ifndef STAR_MOLECULE_FIRST_RESOLVER_H
#define STAR_MOLECULE_FIRST_RESOLVER_H

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace molecule_first {

struct CandidateScore {
    std::string candidate;
    double logLikelihood = 0.0;
};

struct CandidateRead {
    std::string readId;
    std::string featureId;
    std::string rawUmi;
    std::vector<CandidateScore> scores;
};

struct Config {
    double temperature = 1.0;
    double priorAlpha = 1.0;
    double priorBeta = 1.0;
    double gateMinPosterior = 0.95;
    double gateMinMargin = 0.90;
};

struct ReadClique {
    std::string cliqueId;
    std::string featureId;
    std::string rawUmi;
    std::vector<std::string> memberReadIds;
    std::vector<std::string> candidates;
    std::vector<double> logLikelihoodSums;
    std::vector<double> logReadPriors;
    std::vector<double> logEvidence;
    std::vector<double> posterior;
};

struct HardCall {
    std::string cliqueId;
    bool assigned = false;
    std::string candidate;
    double posterior = 0.0;
    double margin = 0.0;
    std::string reason;
};

struct Occupancy {
    std::string umiMode;
    std::string featureId;
    std::string correctedUmi;
    std::string candidate;
    double expectedCount = 0.0;
    std::vector<std::string> readCliqueIds;
};

struct Molecule {
    std::string umiMode;
    std::string product;
    std::string moleculeId;
    std::string featureId;
    std::string correctedUmi;
    std::string candidate;
    std::vector<std::string> memberReadIds;
    std::vector<std::string> readCliqueIds;
};

struct GexReconciliationStats {
    std::uint64_t groups = 0;
    std::uint64_t accepted = 0;
    std::uint64_t correctedCountTies = 0;
    std::uint64_t originalUmiDominanceRejected = 0;
    double inputExpectedMass = 0.0;
    double outputExpectedMass = 0.0;
};

using PriorCounts = std::map<std::string, std::uint64_t>;
using UmiCorrections = std::map<std::tuple<std::string, std::string, std::string>, std::string>;

std::vector<ReadClique> buildReadCliques(const std::vector<CandidateRead> &reads,
                                         const PriorCounts &priorCounts,
                                         const Config &config);

UmiCorrections correctedUmis(const std::vector<ReadClique> &cliques,
                             const std::string &umiMode);

std::vector<Occupancy> weightedOccupancies(const std::vector<ReadClique> &cliques,
                                           const std::string &umiMode);

std::vector<HardCall> gatedHardCalls(const std::vector<ReadClique> &cliques,
                                     const Config &config);

std::vector<Molecule> policyMolecules(const std::vector<ReadClique> &cliques,
                                      const std::string &umiMode,
                                      const std::string &product,
                                      const std::vector<HardCall> &hardCalls);

// GEX-only variants. Candidate-specific UMI correction occurs first, followed
// by the same MultiGeneUMI_CR helper used by standard Solo collapse. These are
// opt-in so the existing Flex resolver behavior is unchanged.
std::vector<Molecule> gexPolicyMolecules(
    const std::vector<ReadClique> &cliques, const std::string &umiMode,
    const std::string &product, const std::vector<HardCall> &hardCalls,
    GexReconciliationStats *stats = nullptr);

std::vector<Occupancy> gexWeightedOccupancies(
    const std::vector<ReadClique> &cliques, const std::string &umiMode,
    GexReconciliationStats *stats = nullptr);

std::pair<std::string, double> topCandidate(const ReadClique &clique);

} // namespace molecule_first

#endif
