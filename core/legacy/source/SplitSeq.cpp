#include "SplitSeq.h"
#include "BarcodeEditDistance.h"
#include "BarcodeAlignmentLikelihood.h"
#include <cstdio>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace splitseq {
namespace {
Region region(const std::string& value) {
    const auto separator = value.find('_');
    if (separator == std::string::npos || separator == 0 || separator + 1 == value.size() ||
        value.find_first_not_of("0123456789", 0) != separator ||
        value.find_first_not_of("0123456789", separator + 1) != std::string::npos)
        throw std::runtime_error("positions must be zero-based offset_length: " + value);
    const auto offset = std::stoull(value.substr(0, separator));
    const auto length = std::stoull(value.substr(separator + 1));
    if (offset > 1000000 || length == 0 || length > 31)
        throw std::runtime_error("invalid split-seq position: " + value);
    return {static_cast<unsigned>(offset), static_cast<unsigned>(length)};
}

bool pack(const std::string& sequence, size_t start, size_t length, uint64_t& value) {
    value = 0;
    for (size_t i = start; i < start + length; ++i) {
        uint64_t digit;
        switch (sequence[i]) {
            case 'A': digit = 0; break;
            case 'C': digit = 1; break;
            case 'G': digit = 2; break;
            case 'T': digit = 3; break;
            default: return false;
        }
        value = (value << 2) | digit;
    }
    return true;
}
bool pack(const std::string& sequence, uint64_t& value) { return pack(sequence, 0, sequence.size(), value); }

std::string unpack(uint64_t value, unsigned length) {
    std::string result(length, 'A');
    for (unsigned i = length; i > 0; --i) { result[i-1] = "ACGT"[value & 3]; value >>= 2; }
    return result;
}
} // namespace

Decoder::Decoder(const std::vector<std::string>& positions, const std::string& umiPosition,
                 const std::vector<std::string>& whitelists, Policy policy) : policy_(policy) {
    if (positions.size() < 2 || positions.size() > 4 || whitelists.size() != positions.size())
        throw std::runtime_error("requires 2-4 CB positions and one allowlist per position, in cell-ID round order");
    umi_ = region(umiPosition);
    if (umi_.length > 16) throw std::runtime_error("split-seq UMI length must be 1-16 bases");
    minimumReadLength_ = umi_.offset + umi_.length;
    std::vector<Region> occupied{umi_};
    for (const auto& position : positions) {
        const auto r = region(position);
        for (const auto& other : occupied)
            if (r.offset < other.offset + other.length && other.offset < r.offset + r.length)
                throw std::runtime_error("split-seq CB/UMI positions must not overlap");
        occupied.push_back(r);
        regions_.push_back(r);
        cellLength_ += r.length;
        minimumReadLength_ = std::max(minimumReadLength_, r.offset + r.length);
    }
    if (cellLength_ > 31) throw std::runtime_error("total split-seq cell barcode length must be <=31 bases");

    for (size_t i = 0; i < whitelists.size(); ++i) {
        std::ifstream input(whitelists[i]);
        if (!input) throw std::runtime_error("cannot open split-seq allowlist: " + whitelists[i]);
        Tier tier;
        std::string line;
        size_t lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            std::istringstream row(line);
            std::string observed, canonical, extra;
            if (!(row >> observed)) continue;
            if (!(row >> canonical)) canonical = observed;
            uint64_t key, value;
            if ((row >> extra) || observed.size() != regions_[i].length || canonical.size() != regions_[i].length ||
                !pack(observed, key) || !pack(canonical, value))
                throw std::runtime_error("invalid split-seq allowlist row " + whitelists[i] + ":" + std::to_string(lineNumber) +
                                         " (expected ACGT sequence, optionally followed by same-length canonical sequence)");
            if (!tier.exact.emplace(key, value).second)
                throw std::runtime_error("duplicate observed split-seq barcode: " + observed);
        }
        if (input.bad() || tier.exact.empty()) throw std::runtime_error("empty or unreadable split-seq allowlist: " + whitelists[i]);
        for (const auto& entry : tier.exact) {
            if (policy_ == Policy::Edit2)
                tier.targets.emplace_back(entry.second, unpack(entry.first, regions_[i].length));
            for (unsigned base = 0; base < regions_[i].length; ++base) {
                for (uint64_t change = 1; change < 4; ++change) {
                    const uint64_t neighbor = entry.first ^ (change << (2 * base));
                    auto inserted = tier.neighbors.emplace(neighbor, entry.first);
                    // Two observed entries one substitution away are only a
                    // conflict when they stand for different canonical
                    // sequences; aliases of one sequence resolve to it.
                    if (!inserted.second && inserted.first->second != UINT64_MAX &&
                        tier.exact.at(inserted.first->second) != entry.second)
                        inserted.first->second = UINT64_MAX;
                }
            }
        }
        std::sort(tier.targets.begin(), tier.targets.end());
        tiers_.push_back(std::move(tier));
    }
}

Decoded Decoder::decode(const std::string& read) const {
    if (policy_ == Policy::Edit2) return decodeEdit2(read);
    Decoded result;
    if (read.size() < minimumReadLength_) return result;
    result.umi = read.substr(umi_.offset, umi_.length);
    for (size_t i = 0; i < regions_.size(); ++i) {
        const auto& r = regions_[i];
        result.rawCell.append(read, r.offset, r.length);
        uint64_t observed;
        auto found = tiers_[i].exact.end();
        if (!pack(read, r.offset, r.length, observed)) {
            // One non-ACGT base is one substitution, as in ordinary barcode
            // matching: the round is corrected when the entries reachable by
            // filling that base stand for a single canonical sequence.
            size_t unknown = std::string::npos;
            for (size_t p = r.offset; p < r.offset + r.length; ++p) {
                if (read[p] == 'A' || read[p] == 'C' || read[p] == 'G' || read[p] == 'T') continue;
                if (unknown != std::string::npos) { result.status = Status::InvalidBarcode; return result; }
                unknown = p;
            }
            std::string filled = read.substr(r.offset, r.length);
            for (char base : {'A', 'C', 'G', 'T'}) {
                filled[unknown - r.offset] = base;
                pack(filled, observed);
                const auto candidate = tiers_[i].exact.find(observed);
                if (candidate == tiers_[i].exact.end()) continue;
                if (found != tiers_[i].exact.end() && found->second != candidate->second) {
                    result.status = Status::Ambiguous; return result;
                }
                found = candidate;
            }
            if (found == tiers_[i].exact.end()) { result.status = Status::NoMatch; return result; }
            result.correctedRounds |= 1U << i;
        } else if ((found = tiers_[i].exact.find(observed)) == tiers_[i].exact.end()) {
            const auto neighbor = tiers_[i].neighbors.find(observed);
            if (neighbor == tiers_[i].neighbors.end()) { result.status = Status::NoMatch; return result; }
            if (neighbor->second == UINT64_MAX) { result.status = Status::Ambiguous; return result; }
            found = tiers_[i].exact.find(neighbor->second);
            result.correctedRounds |= 1U << i;
        }
        // Apply the reference's observed->canonical replacement once, after correction.
        result.cell = (result.cell << (2 * r.length)) | found->second;
        result.canonicalCell += unpack(found->second, r.length);
    }
    result.status = Status::Accepted;
    return result;
}

RoundMatch Decoder::matchEdit2(const std::string& read, size_t round) const {
    const Region r = regions_[round];
    const auto& tier = tiers_[round];
    uint64_t word;
    RoundMatch match;
    if (pack(read, r.offset, r.length, word)) {
        const auto exact = tier.exact.find(word);
        if (exact != tier.exact.end()) {
            match.distance = 0;
            match.canonical = exact->second;
            return match;
        }
        // At equal lengths the only possible single edit is a substitution.
        const auto neighbor = tier.neighbors.find(word);
        if (neighbor != tier.neighbors.end() && neighbor->second != UINT64_MAX) {
            match.distance = 1;
            match.canonical = tier.exact.at(neighbor->second);
            return match;
        }
    }
    const barcode_edit::Myers64Pattern query(read.data() + r.offset, r.length);
    for (const auto& target : tier.targets) {
        const unsigned dist = query.distance(target.second.data(), r.length);
        if (dist > 2 || dist > match.distance) continue;
        if (dist < match.distance) {
            match.distance = dist;
            match.canonical = target.first;
            match.alternatives.clear();
        } else if (target.first != match.canonical &&
                   (match.alternatives.empty() || match.alternatives.back() != target.first)) {
            match.alternatives.push_back(target.first);
        }
    }
    return match;
}

Decoded Decoder::decodeEdit2(const std::string& read) const {
    Decoded result;
    if (read.size() < minimumReadLength_) return result;
    result.umi = read.substr(umi_.offset, umi_.length);
    bool ambiguous = false;
    for (size_t i = 0; i < regions_.size(); ++i) {
        const auto r = regions_[i];
        result.rawCell.append(read, r.offset, r.length);
        auto& match = result.matches[i];
        match = matchEdit2(read, i);
        result.totalDistance += match.distance;
        if (match.distance > 2 || result.totalDistance > 6) {
            result.status = Status::NoMatch;
            return result;
        }
        if (match.distance) result.correctedRounds |= 1U << i;
        ambiguous |= !match.alternatives.empty();
    }
    // The Cartesian product of per-round minima gives every minimum-total tuple.
    // Preserve it in factored form. Never emit an arbitrary cell for ties.
    if (ambiguous) {
        result.status = Status::Ambiguous;
        return result;
    }
    for (size_t i = 0; i < regions_.size(); ++i) {
        const auto length = regions_[i].length;
        result.cell = (result.cell << (2 * length)) | result.matches[i].canonical;
        result.canonicalCell += unpack(result.matches[i].canonical, length);
    }
    result.status = Status::Accepted;
    return result;
}

void countExact(const Decoded& result, ExactCounts& counts, ExactCounts* prefixCounts) {
    uint64_t raw;
    if (result.status == Status::Accepted && result.totalDistance == 0 && pack(result.rawCell, raw)) {
        ++counts[raw];
        if (prefixCounts) ++(*prefixCounts)[raw];
    }
}

static uint64_t uniquePrefixWinner(const std::map<uint64_t, std::array<uint64_t, 2>>& support) {
    uint64_t bestCount = 0, winner = UINT64_MAX;
    for (const auto& cell : support) {
        const uint64_t count = cell.second[1];
        if (count > bestCount) {
            bestCount = count;
            winner = cell.first;
        } else if (count == bestCount) {
            winner = UINT64_MAX;
        }
    }
    return winner;
}

Posterior Decoder::resolveTie(Decoded& result, const std::string& barcodeQual,
                              const ExactCounts& counts, double cutoff, const ExactCounts* prefixCounts) const {
    if (result.status != Status::Ambiguous || barcodeQual.size() != cellLength_)
        throw std::logic_error("invalid split-seq tie or barcode qualities");
    struct Oligo { uint64_t raw, canonical; double logLikelihood; };
    std::vector<std::vector<Oligo>> options(regions_.size());
    size_t offset = 0;
    for (size_t i = 0; i < regions_.size(); ++i) {
        const unsigned len = regions_[i].length;
        const auto word = result.rawCell.substr(offset, len);
        const auto qual = barcodeQual.substr(offset, len);
        offset += len;
        const barcode_edit::Myers64Pattern query(word.data(), len);
        double maximum = -INFINITY;
        for (const auto& target : tiers_[i].targets) {
            if (query.distance(target.second.data(), len) != result.matches[i].distance) continue;
            uint64_t raw;
            pack(target.second, raw);
            const double likelihood = barcode_quality::alignmentLogLikelihood(word, qual, target.second);
            options[i].push_back({raw, target.first, likelihood});
            maximum = std::max(maximum, likelihood);
        }
        // Subtract the same constant from every candidate in a round. The
        // best joint likelihood is now one, even at very low base qualities.
        for (auto& option : options[i]) option.logLikelihood -= maximum;
    }
    Posterior posterior;
    std::function<void(size_t, uint64_t, uint64_t, double)> visit =
        [&](size_t round, uint64_t raw, uint64_t canonical, double logLikelihood) {
            if (round == options.size()) {
                const auto it = counts.find(raw);
                const uint64_t count = it == counts.end() ? 0 : it->second;
                posterior.cells[canonical] += (count + 1.0) * std::exp(logLikelihood);
                if (prefixCounts) {
                    auto& support = posterior.exactSupport[canonical];
                    support[0] += count;
                    const auto prefix = prefixCounts->find(raw);
                    if (prefix != prefixCounts->end()) support[1] += prefix->second;
                }
                return;
            }
            const unsigned bits = 2 * regions_[round].length;
            for (const auto& option : options[round])
                visit(round + 1, (raw << bits) | option.raw, (canonical << bits) | option.canonical,
                      logLikelihood + option.logLikelihood);
        };
    visit(0, 0, 0, 0.0);
    double sum = 0;
    for (const auto& cell : posterior.cells) sum += cell.second;
    uint64_t winner = 0;
    double runnerUp = 0;
    for (auto& cell : posterior.cells) {
        cell.second /= sum;
        if (cell.second > posterior.best) {
            runnerUp = posterior.best;
            posterior.best = cell.second;
            winner = cell.first;
        } else {
            runnerUp = std::max(runnerUp, cell.second);
        }
    }
    // A permissive cutoff still needs a unique winner after floating-point normalization.
    if (posterior.best < cutoff || posterior.best - runnerUp <= 1e-12) {
        winner = uniquePrefixWinner(posterior.exactSupport);
        if (winner == UINT64_MAX) return posterior;
        posterior.frequencyRescued = true;
    }
    result.cell = winner;
    result.canonicalCell = unpack(winner, cellLength_);
    result.status = Status::Accepted;
    return posterior;
}

TieWriter::TieWriter(const std::string& path, const std::string& pendingPath) : pendingPath_(pendingPath) {
    out_.exceptions(std::ios::badbit | std::ios::failbit);
    out_.open(path);
    out_ << "read_name\tread_file_index\traw_barcode\traw_umi\ttotal_edits\tround_edits\tcanonical_candidates"
            "\tbest_posterior\tassigned_barcode\tcandidate_posteriors\tassignment_method\tassigned_posterior\tcandidate_exact_counts\n";
}

uint64_t TieWriter::defer(const std::string& name, unsigned fileIndex, const std::string& read,
                          const std::string& cbQual, const std::string& umiQual, int umiStatus) {
    std::lock_guard<std::mutex> guard(mutex_);
    if (!pending_.is_open()) {
        pending_.exceptions(std::ios::badbit | std::ios::failbit);
        pending_.open(pendingPath_, std::ios::in | std::ios::out | std::ios::trunc);
    }
    pending_ << fileIndex << '\t' << read << '\t' << cbQual << '\t' << umiQual << '\t'
             << umiStatus << '\t' << name << '\n';
    return pendingBit | rows_++;
}

void TieWriter::finish(const Decoder& decoder, const ExactCounts& counts, double cutoff, const Decision& record,
                       const ExactCounts* prefixCounts) {
    if (rows_) {
        pending_.seekg(0);
        decisions_.reserve(rows_);
        std::string name, read, cbQual, umiQual;
        unsigned fileIndex;
        int umiStatus;
        for (uint64_t i = 0; i < rows_; ++i) {
            pending_ >> fileIndex >> read >> cbQual >> umiQual >> umiStatus;
            pending_.get(); // tab before the FASTQ identifier, which can contain spaces
            std::getline(pending_, name);
            auto result = decoder.decode(read);
            const auto posterior = decoder.resolveTie(result, cbQual, counts, cutoff, prefixCounts);
            decisions_.push_back(result.status == Status::Accepted ? result.cell : UINT64_MAX);
            write(name, fileIndex, result, decoder.regions(), posterior);
            record(result, umiStatus, cbQual, umiQual);
        }
        pending_.close();
        std::remove(pendingPath_.c_str());
    }
    out_.close();
}

void TieWriter::write(const std::string& readName, unsigned fileIndex, const Decoded& result,
                      const std::vector<Region>& regions, const Posterior& posterior) {
    std::ostringstream row;
    row << readName << '\t' << fileIndex << '\t' << result.rawCell << '\t' << result.umi << '\t' << result.totalDistance;
    for (size_t i = 0; i < regions.size(); ++i)
        row << (i == 0 ? '\t' : ',') << result.matches[i].distance;
    for (size_t i = 0; i < regions.size(); ++i) {
        const auto& match = result.matches[i];
        row << (i == 0 ? '\t' : '|') << unpack(match.canonical, regions[i].length);
        for (auto candidate : match.alternatives) row << ',' << unpack(candidate, regions[i].length);
    }
    row << std::setprecision(17) << '\t' << posterior.best << '\t'
        << (result.status == Status::Accepted ? result.canonicalCell : "-") << '\t';
    for (auto it = posterior.cells.begin(); it != posterior.cells.end(); ++it) {
        if (it != posterior.cells.begin()) row << ',';
        row << unpack(it->first, result.rawCell.size()) << '=' << it->second;
    }
    if (result.status == Status::Accepted) {
        row << (posterior.frequencyRescued ? "\tprefix_frequency\t" : "\tbayesian\t")
            << posterior.cells.at(result.cell);
    } else {
        row << "\tunresolved\t-";
    }
    row << '\t';
    if (posterior.exactSupport.empty()) row << '-';
    for (auto it = posterior.exactSupport.begin(); it != posterior.exactSupport.end(); ++it) {
        if (it != posterior.exactSupport.begin()) row << ',';
        row << unpack(it->first, result.rawCell.size()) << '=' << it->second[0] << '/' << it->second[1];
    }
    out_ << row.str() << '\n';
}

void writeMetrics(std::ostream& out, const std::array<uint64_t, 13>& stats) {
    const char* names[] = {"reads_total", "barcode_exact", "barcode_corrected", "read_too_short",
        "barcode_non_acgt", "barcode_no_match", "barcode_ambiguous", "umi_non_acgt",
        "umi_homopolymer", "round1_corrected", "round2_corrected", "round3_corrected", "round4_corrected"};
    out << "metric\tvalue\n";
    for (size_t i = 0; i < stats.size(); ++i) out << names[i] << '\t' << stats[i] << '\n';
}
} // namespace splitseq
