#include "PipSeq.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace pipseq {
namespace {
const unsigned lengths[] = {8, 6, 6, 8};
const unsigned starts[] = {0, 11, 20, 31};
int base(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}
bool pack(const std::string& s, size_t start, size_t length, uint64_t& value) {
    value = 0;
    if (start + length > s.size()) return false;
    for (size_t i = start; i < start + length; ++i) {
        const int b = base(s[i]);
        if (b < 0) return false;
        value = (value << 2) | static_cast<unsigned>(b);
    }
    return true;
}
unsigned distance(uint64_t a, uint64_t b, unsigned length) {
    uint64_t d = a ^ b;
    unsigned n = 0;
    for (unsigned i = 0; i < length; ++i, d >>= 2) n += (d & 3) != 0;
    return n;
}
}

Decoder::Decoder(const std::vector<std::string>& paths, MultiMatch multiMatch) {
    if (paths.size() != 4) throw std::runtime_error("PIP requires four tier allowlists in physical read order (8/6/6/8 bases)");
    for (unsigned t = 0; t < 4; ++t) {
        std::ifstream in(paths[t]);
        if (!in) throw std::runtime_error("Cannot open PIP tier allowlist: " + paths[t]);
        std::string line;
        std::vector<uint64_t> packed;
        std::unordered_set<uint64_t> seen;
        while (std::getline(in, line)) {
            std::istringstream row(line);
            std::string sequence, extra;
            if (!(row >> sequence)) continue;
            uint64_t value;
            if ((row >> extra) || sequence.size() != lengths[t] || !pack(sequence, 0, sequence.size(), value))
                throw std::runtime_error("PIP tier allowlist must contain one uppercase A/C/G/T sequence per line of length " + std::to_string(lengths[t]) + ": " + paths[t]);
            if (!seen.insert(value).second) throw std::runtime_error("Duplicate PIP tier sequence: " + sequence);
            tiers_[t].push_back(sequence);
            packed.push_back(value);
        }
        if (!in.eof() || packed.empty()) throw std::runtime_error("Unreadable or empty PIP tier allowlist: " + paths[t]);
        // Small independent tier tables; never materialize the Cartesian product.
        lookup_[t].resize(1u << (2 * lengths[t]));
        for (size_t q = 0; q < lookup_[t].size(); ++q) {
            Match& m = lookup_[t][q];
            int exact = -1, firstNeighbor = -1;
            unsigned neighbors = 0;
            for (size_t i = 0; i < packed.size(); ++i) {
                const unsigned d = distance(q, packed[i], lengths[t]);
                m.best = std::min(m.best, d);
                if (d == 0) exact = i;
                else if (d == 1 && neighbors++ == 0) firstNeighbor = i;
            }
            // An exact entry always wins: a listed sequence is never rewritten
            // to a listed neighbor. Otherwise one substitution is corrected.
            if (exact >= 0) m.entry = exact;
            else if (neighbors == 1 || (neighbors > 1 && multiMatch == MultiMatch::FirstListed)) m.entry = firstNeighbor;
            else if (neighbors > 1) m.ambiguous = true;
        }
    }
}

Decoded Decoder::decode(const std::string& read) const {
    Decoded result;
    unsigned best = std::numeric_limits<unsigned>::max();
    bool tied = false;
    std::array<Match, 4> selected;
    for (unsigned phase = 0; phase < 4; ++phase) {
        std::array<Match, 4> matches;
        unsigned score = 0;
        bool valid = true;
        for (unsigned t = 0; t < 4; ++t) {
            uint64_t value;
            if (!pack(read, phase + starts[t], lengths[t], value)) { valid = false; break; }
            matches[t] = lookup_[t][value];
            score += matches[t].best;
        }
        if (!valid) continue;
        if (score < best) {
            best = score; result.phase = phase; selected = matches; tied = false;
        } else if (score == best) tied = true;
    }
    if (result.phase < 0) return result;
    if (tied) { result.status = DecodeStatus::PhaseTie; return result; }
    // A tier farther than one substitution from every entry is a non-match,
    // whatever the other tiers are; otherwise an unresolved tier is ambiguous.
    result.status = DecodeStatus::NoMatch;
    for (unsigned t = 0; t < 4; ++t) if (selected[t].best > 1) return result;
    result.status = DecodeStatus::Ambiguous;
    for (unsigned t = 0; t < 4; ++t) if (selected[t].entry < 0) return result;
    for (unsigned t = 0; t < 4; ++t) {
        uint64_t value;
        const std::string& sequence = tiers_[t][selected[t].entry];
        pack(sequence, 0, lengths[t], value);
        result.cell = (result.cell << (2 * lengths[t])) | value;
        result.mismatches += selected[t].best;
    }
    result.status = DecodeStatus::Accepted;
    uint64_t bin;
    if (pack(read, result.phase + 39, 3, bin)) result.bin = bin;
    return result;
}

std::string Decoder::cellSequence(uint64_t cell) {
    std::string s(28, 'A');
    for (int i = 27; i >= 0; --i, cell >>= 2) s[i] = "ACGT"[cell & 3];
    return s;
}

namespace {
bool fragmentLess(const Fragment& a, const Fragment& b) {
    if (a.cell != b.cell) return a.cell < b.cell;
    if (a.gene != b.gene) return a.gene < b.gene;
    return a.key < b.key;
}
bool sameFragment(const Fragment& a, const Fragment& b) {
    return a.cell == b.cell && a.gene == b.gene && a.key == b.key;
}
const size_t compactAfter = size_t(1) << 20; // unsorted records before a compaction
}

std::string Decoder::rawCell(const std::string& read, int phase) {
    if (phase < 0 || read.size() < size_t(phase) + 39) return "-";
    std::string s;
    for (unsigned t = 0; t < 4; ++t) s.append(read, phase + starts[t], lengths[t]);
    return s;
}
std::string Decoder::rawBin(const std::string& read, int phase) {
    if (phase < 0 || read.size() < size_t(phase) + 42) return "-";
    return read.substr(phase + 39, 3);
}

void Fragments::add(uint64_t cell, uint32_t gene, unsigned bin, uint64_t coordinate, bool reverse) {
    if (bin >= 64 || coordinate > (std::numeric_limits<uint64_t>::max() >> 7))
        throw std::runtime_error("PIP fragment key out of range");
    records_.push_back({cell, (coordinate << 7) | (uint64_t(reverse) << 6) | bin, gene, 1});
    // Keep memory near the number of distinct fragments, not of reads.
    if (records_.size() - sorted_ >= std::max(compactAfter, sorted_)) finalize();
}
void Fragments::collapse() {
    size_t out = 0;
    for (size_t i = 0; i < records_.size(); ++i) {
        if (out && sameFragment(records_[out - 1], records_[i])) {
            const uint64_t reads = uint64_t(records_[out - 1].reads) + records_[i].reads;
            if (reads > std::numeric_limits<uint32_t>::max()) throw std::runtime_error("PIP fragment read count exceeds 32 bits");
            records_[out - 1].reads = reads;
        } else records_[out++] = records_[i];
    }
    records_.resize(out);
    sorted_ = out;
}
void Fragments::finalize() {
    if (sorted_ == records_.size()) return;
    std::sort(records_.begin() + sorted_, records_.end(), fragmentLess);
    std::inplace_merge(records_.begin(), records_.begin() + sorted_, records_.end(), fragmentLess);
    collapse();
}
void Fragments::absorb(Fragments& other) {
    finalize();
    other.finalize();
    if (records_.empty()) records_.swap(other.records_);
    else {
        const size_t middle = records_.size();
        records_.insert(records_.end(), other.records_.begin(), other.records_.end());
        std::inplace_merge(records_.begin(), records_.begin() + middle, records_.end(), fragmentLess);
    }
    other.clear();
    collapse();
}
void Fragments::clear() {
    std::vector<Fragment>().swap(records_);
    sorted_ = 0;
}

double expectedMolecules(unsigned bins) {
    if (bins > 64) throw std::runtime_error("PIP occupied bin count exceeds 64");
    // Expected number of uniform draws to first observe B distinct bins.
    // This is the coupon-collector stopping-time mean, not log inversion.
    double mean = 0;
    for (unsigned i = 0; i < bins; ++i) mean += 64.0 / (64 - i);
    return mean;
}
Calibration calibrate(const Fragments& fragments, double suppliedIPM) {
    Calibration c;
    if (!std::isfinite(suppliedIPM) || (suppliedIPM != 0 && suppliedIPM < 1))
        throw std::runtime_error("PIP IPM must be 0 (estimate) or a finite value >= 1");
    if (suppliedIPM != 0) { c.raw = c.effective = suppliedIPM; c.source = "supplied"; return c; }
    long double sum = 0;
    bool needsIPM = false;
    fragments.forEachGroup([&](const GroupSummary& g) {
        needsIPM |= g.bins > 10;
        if (g.bins >= 5 && g.bins <= 32) {
            sum += g.fragments / expectedMolecules(g.bins);
            ++c.groups;
        }
    });
    if (!c.groups) {
        if (needsIPM) throw std::runtime_error("No PIP cell/gene groups with 5-32 bins for IPM calibration; supply --soloPipseqIPM >= 1 or process more reads");
        c.source = "unused_low_occupancy";
    } else {
        c.raw = sum / c.groups;
        // Enforce the documented bound: transcripts cannot exceed fragments.
        c.effective = std::max(1.0, c.raw);
        c.source = "estimated";
    }
    return c;
}
uint64_t transcripts(unsigned bins, uint64_t fragments, double ipm) {
    if (!std::isfinite(ipm) || ipm < 1) throw std::runtime_error("Invalid PIP IPM");
    if (bins <= 10) return bins;
    return std::max<uint64_t>(bins, std::min<uint64_t>(fragments, std::floor(fragments / ipm)));
}

const char* const kitAdapters5p = "AGAGTGAATGGG,TCAACGCAGAGT";
const char* const kitAdapters3p = "polyA,polyG";

uint32_t clip5pAdapter(const char* read, uint32_t length, const char* adapters, uint32_t adaptersLength,
                       double maxMismatchFraction, uint32_t minOverlap) {
    int bestScore = 0;
    uint32_t bestClip = 0;
    for (uint32_t first = 0; first < adaptersLength;) {
        // One adapter: a run of A/C/G/T codes.
        while (first < adaptersLength && (adapters[first] < 0 || adapters[first] > 3)) ++first;
        uint32_t adapterLength = 0;
        while (first + adapterLength < adaptersLength && adapters[first + adapterLength] >= 0 && adapters[first + adapterLength] <= 3) ++adapterLength;
        const char* adapter = adapters + first;
        first += adapterLength;
        if (adapterLength == 0 || length == 0) continue;
        const uint32_t overlapMin = std::max<uint32_t>(1, std::min(minOverlap, adapterLength));
        // start: read position of the adapter's first base; negative when the
        // adapter hangs off the read start. It never hangs off the read end.
        for (int64_t start = -int64_t(adapterLength - overlapMin); start + adapterLength <= length; ++start) {
            const uint32_t skip = start < 0 ? uint32_t(-start) : 0;
            const uint32_t overlap = adapterLength - skip;
            const uint32_t allowed = uint32_t(maxMismatchFraction * overlap + 1e-9);
            uint32_t mismatches = 0;
            for (uint32_t i = skip; i < adapterLength && mismatches <= allowed; ++i)
                mismatches += read[start + i] != adapter[i];
            if (mismatches > allowed) continue;
            const int score = int(overlap) - 2 * int(mismatches);
            const uint32_t clip = uint32_t(start + adapterLength);
            if (score > bestScore || (score == bestScore && clip > bestClip)) { bestScore = score; bestClip = clip; }
        }
    }
    return bestClip;
}

uint32_t clip3pHomopolymer(const char* read, uint32_t length, char base, double maxMismatchFraction, uint32_t window) {
    if (window == 0 || length < window) return 0;
    const uint32_t allowed = uint32_t(maxMismatchFraction * window + 1e-9);
    uint32_t other = 0;
    for (uint32_t i = 0; i < window; ++i) other += read[i] != base;
    uint32_t best = length, fewest = allowed + 1;
    for (uint32_t start = 0;; ++start) {
        if (other < fewest) { best = start; fewest = other; if (other == 0) break; }
        if (start + window >= length) break;
        other += (read[start + window] != base);
        other -= (read[start] != base);
    }
    return length - best;
}
} // namespace pipseq
