#ifndef STAR_SPLITSEQ_H
#define STAR_SPLITSEQ_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace splitseq {

struct Region {
    unsigned offset, length;
};

enum class Status { Accepted, ShortRead, InvalidBarcode, NoMatch, Ambiguous };

struct Decoded {
    Status status = Status::ShortRead;
    uint64_t cell = 0;
    unsigned correctedRounds = 0;
    std::string rawCell, canonicalCell, umi;
};

// Fixed coordinates and reference data only; no kit-specific pipeline behavior.
// Immutable after construction, shared by all mapping threads.
class Decoder {
public:
    Decoder(const std::vector<std::string>& positions, const std::string& umiPosition,
            const std::vector<std::string>& whitelists);
    Decoded decode(const std::string& read) const;
    const std::vector<Region>& regions() const { return regions_; }
    Region umiRegion() const { return umi_; }
    unsigned cellLength() const { return cellLength_; }
    unsigned minimumReadLength() const { return minimumReadLength_; }
    size_t tierSize(size_t tier) const { return tiers_.at(tier).exact.size(); }

private:
    struct Tier {
        std::unordered_map<uint64_t, uint64_t> exact; // observed -> canonical
        std::unordered_map<uint64_t, uint64_t> neighbors; // variant -> observed, or ambiguous
    };
    std::vector<Region> regions_;
    Region umi_;
    std::vector<Tier> tiers_;
    unsigned cellLength_ = 0, minimumReadLength_ = 0;
};

} // namespace splitseq
#endif
