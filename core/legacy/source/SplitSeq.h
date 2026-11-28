#ifndef STAR_SPLITSEQ_H
#define STAR_SPLITSEQ_H

#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <map>
#include <mutex>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace splitseq {

struct Region {
    unsigned offset, length;
};

enum class Status { Accepted, ShortRead, InvalidBarcode, NoMatch, Ambiguous };
enum class Policy { Hamming1, Edit2 };

struct RoundMatch {
    unsigned distance = 3; // beyond the two-edit radius
    uint64_t canonical = 0;
    std::vector<uint64_t> alternatives; // other canonical wells at the same minimum
};

struct Decoded {
    Status status = Status::ShortRead;
    uint64_t cell = 0;
    unsigned correctedRounds = 0;
    std::string rawCell, canonicalCell, umi;
    unsigned totalDistance = 0;
    std::array<RoundMatch, 4> matches;
};

using ExactCounts = std::unordered_map<uint64_t, uint64_t>;
constexpr uint64_t frequencyPrefixReads = 1000000;
struct Posterior {
    std::map<uint64_t, double> cells;
    std::map<uint64_t, std::array<uint64_t, 2>> exactSupport; // full input, first million
    double best = 0;
    bool frequencyRescued = false;
};

// Only exact raw oligo tuples enter the prior, before UMI or alignment filters.
void countExact(const Decoded& result, ExactCounts& counts, ExactCounts* prefixCounts = nullptr);

// Fixed coordinates and reference data only; no kit-specific pipeline behavior.
// Immutable after construction, shared by all mapping threads.
class Decoder {
public:
    Decoder(const std::vector<std::string>& positions, const std::string& umiPosition,
            const std::vector<std::string>& whitelists, Policy policy = Policy::Hamming1);
    Decoded decode(const std::string& read) const;
    Posterior resolveTie(Decoded& result, const std::string& barcodeQual,
                         const ExactCounts& counts, double cutoff, const ExactCounts* prefixCounts = nullptr) const;
    const std::vector<Region>& regions() const { return regions_; }
    Region umiRegion() const { return umi_; }
    unsigned cellLength() const { return cellLength_; }
    unsigned minimumReadLength() const { return minimumReadLength_; }
    size_t tierSize(size_t tier) const { return tiers_.at(tier).exact.size(); }

private:
    struct Tier {
        std::unordered_map<uint64_t, uint64_t> exact; // observed -> canonical
        std::unordered_map<uint64_t, uint64_t> neighbors; // variant -> observed, or ambiguous
        std::vector<std::pair<uint64_t, std::string>> targets; // sorted by canonical well
    };
    RoundMatch matchEdit2(const std::string& read, size_t round) const;
    Decoded decodeEdit2(const std::string& read) const;
    Policy policy_;
    std::vector<Region> regions_;
    Region umi_;
    std::vector<Tier> tiers_;
    unsigned cellLength_ = 0, minimumReadLength_ = 0;
};

// One file per mapping run, shared by its mapping threads.
class TieWriter {
public:
    TieWriter(const std::string& path, const std::string& pendingPath);
    uint64_t defer(const std::string& name, unsigned fileIndex, const std::string& read,
                   const std::string& cbQual, const std::string& umiQual, int umiStatus);
    using Decision = std::function<void(const Decoded&, int, const std::string&, const std::string&)>;
    void finish(const Decoder& decoder, const ExactCounts& counts, double cutoff, const Decision& record,
                const ExactCounts* prefixCounts = nullptr);

    // At most 31 bases in a physical cell key. Bit 62 marks a pending row in the Solo spool.
    static constexpr uint64_t pendingBit = UINT64_C(1) << 62;
    uint64_t cellFor(uint64_t key) const {
        return key & pendingBit ? decisions_.at(key ^ pendingBit) : key;
    }
private:
    void write(const std::string& name, unsigned fileIndex, const Decoded& result,
               const std::vector<Region>& regions, const Posterior& posterior);
    std::ofstream out_;
    std::fstream pending_;
    std::string pendingPath_;
    uint64_t rows_ = 0;
    std::vector<uint64_t> decisions_;
    std::mutex mutex_;
};

// Splitseq.metrics.tsv from the summed per-thread counters of
// SoloReadBarcode::splitStats.
void writeMetrics(std::ostream& out, const std::array<uint64_t, 13>& stats);

} // namespace splitseq
#endif
