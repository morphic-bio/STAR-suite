#ifndef CODE_CbCorrector
#define CODE_CbCorrector

#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include "htslib/khash.h"

struct CbCandidateRange {
    size_t offset;
    uint32_t count;
};

// Mix all packed bases: whitelist suffix structure must not cluster in the
// low bits used by khash's power-of-two bucket count.
KHASH_INIT(cbCorrectorLookup, uint32_t, uint32_t, 1, __ac_Wang_hash, kh_int_hash_equal)
KHASH_INIT(cbCorrectorRanges, uint32_t, CbCandidateRange, 1, __ac_Wang_hash, kh_int_hash_equal)

// Standard integer types (CbCorrector is self-contained, doesn't need IncludeDefine.h)
// uint32_t and uint8_t are provided by <cstdint>

// CB correction module - self-contained library for Cell Ranger-compatible barcode correction
// Phase 1: exact + 1-hamming rescue, optional N expansion

struct CbMatch {
    uint32_t whitelistIdx;   // 1-based CB index, 0 if no match
    uint8_t hammingDist;     // 0, 1, or >1 (use 255 for >1)
    bool ambiguous;          // true if multiple WL entries share this CB
    std::vector<uint32_t> ambiguousIdx; // optional, only when ambiguous (1-based indices)
    
    CbMatch() : whitelistIdx(0), hammingDist(255), ambiguous(false) {}
};

class CbCorrector {
public:
    // Constructor: initialize with whitelist
    // maxHamming: maximum Hamming distance allowed (default 1, 0 = exact only)
    explicit CbCorrector(const std::vector<std::string> &whitelist, int maxHamming = 1,
                         bool cbqNativeOrder = false);
    CbCorrector(const CbCorrector&) = delete;
    CbCorrector& operator=(const CbCorrector&) = delete;
    CbCorrector(CbCorrector&&) = default;
    CbCorrector& operator=(CbCorrector&&) = default;
    
    // Correct a cell barcode
    // Returns CbMatch with correction result
    CbMatch correct(const std::string &cb) const;

    // Look up an already packed CBQ barcode without materializing ASCII.
    // The input keeps CBQ's native order (the first base is in the low bits).
    // Only exact and unique H1 matches return true; ambiguous matches remain on
    // the quality-aware string path.
    bool correctPackedCbq(uint32_t packedKey, uint32_t &whitelistIdx,
                          uint8_t &hammingDist) const;
    
    // Get whitelist size
    size_t whitelistSize() const { return whitelist_.size(); }
    
    // Get whitelist sequences (for Bayesian resolver)
    const std::vector<std::string>& whitelist() const { return whitelist_; }
    
    // Read-only view into the single candidate array; indices remain 0-based
    // and in whitelist encounter order. Valid for this corrector's lifetime.
    struct CandidateView {
        const uint32_t* data;
        size_t count;
        const uint32_t* begin() const { return data; }
        const uint32_t* end() const { return count ? data + count : data; }
        size_t size() const { return count; }
    };
    template<class Visitor> void forEachAmbiguousVariant(Visitor visit) const {
        const auto* h = ambiguousRanges_.get();
        for (khint_t k = kh_begin(h); k != kh_end(h); ++k) {
            if (kh_exist(h, k)) {
                const auto& range = kh_val(h, k);
                visit(kh_key(h, k), CandidateView{candidateIndices_.data() + range.offset, range.count});
            }
        }
    }
    
    // Get CB length
    size_t getCbLength() const { return cbLength_; }
    
    // Decode packed CB key to string (public helper)
    std::string decodePackedKey(uint32_t packedKey, size_t cbLength) const;
    
private:
    // Whitelist storage (canonical CB strings)
    std::vector<std::string> whitelist_;
    
    struct LookupDeleter {
        void operator()(khash_t(cbCorrectorLookup)* h) const { kh_destroy(cbCorrectorLookup, h); }
    };
    struct RangeDeleter {
        void operator()(khash_t(cbCorrectorRanges)* h) const { kh_destroy(cbCorrectorRanges, h); }
    };
    using LookupTable = std::unique_ptr<khash_t(cbCorrectorLookup), LookupDeleter>;
    // Exact values are 0-based. Variant values are 1-based, with 0 for ambiguity.
    LookupTable exactMap_, variantMap_;
    std::unique_ptr<khash_t(cbCorrectorRanges), RangeDeleter> ambiguousRanges_;
    std::vector<uint32_t> candidateIndices_;

    static const uint32_t* findLookup(const LookupTable& table, uint32_t key) {
        const khint_t k = kh_get(cbCorrectorLookup, table.get(), key);
        return k == kh_end(table.get()) ? nullptr : &kh_val(table.get(), k);
    }
    CandidateView ambiguousCandidates(uint32_t key) const {
        const auto* h = ambiguousRanges_.get();
        const khint_t k = kh_get(cbCorrectorRanges, h, key);
        if (k == kh_end(h)) return CandidateView{nullptr, 0};
        const auto& range = kh_val(h, k);
        return CandidateView{candidateIndices_.data() + range.offset, range.count};
    }

    template<class Visitor> void forEachVariant(Visitor visit) const {
        for (size_t i = 0; i < whitelist_.size(); ++i) {
            PackedCB packed;
            if (!encodeCB(whitelist_[i], packed)) continue;
            for (size_t pos = 0; pos < cbLength_; ++pos) {
                const uint32_t shift = packedShift(pos, cbLength_);
                const uint32_t current = (packed.key >> shift) & 3u;
                for (uint32_t base = 0; base < 4; ++base) {
                    if (base != current)
                        visit((packed.key & ~(3u << shift)) | (base << shift), static_cast<uint32_t>(i));
                }
            }
        }
    }
    
    // CB length (assumed constant for all CBs in whitelist)
    size_t cbLength_;
    
    int maxHamming_;  // Maximum Hamming distance allowed

    // CBQ stores the first base in the least-significant two bits. FASTQ and
    // the legacy STAR packed representation store it at the other end.
    bool cbqNativeOrder_;
    
    // Packed key representation: 32-bit CB + 16-bit N mask
    struct PackedCB {
        uint32_t key;      // 16 bases × 2 bits = 32 bits
        uint16_t nMask;    // Bit set to 1 where N is present (0 = no Ns)
        
        PackedCB() : key(0), nMask(0) {}
        bool hasN() const { return nMask != 0; }
    };
    
    // Helper: encode CB string to packed key + N mask
    // Returns true if encoding successful, false if CB too long or invalid
    bool encodeCB(const std::string &cb, PackedCB &packed) const;

    // Bit shift for a logical sequence position in the configured table order.
    uint32_t packedShift(size_t pos, size_t cbLength) const;
    
    // Helper: decode packed key back to CB string
    std::string decodeCB(const PackedCB &packed, size_t cbLength) const;
    
    // Helper: generate all 1-hamming variants of a packed CB
    void generateVariants(const PackedCB &packed, std::vector<PackedCB> &variants, size_t cbLength) const;
    
    // Helper: expand N positions in a packed CB
    // Returns true if exactly one valid candidate found (early return), false otherwise
    // If ambiguous, corrected will contain the first best hit and ambiguousSeqs will contain up to 5 best sequences
    bool expandN(const PackedCB &packed, size_t cbLength, PackedCB &corrected, int &hammingDist, 
                 std::vector<PackedCB> &ambiguousSeqs) const;
    
    // Helper: check if CB contains N
    bool hasN(const std::string &cb) const;
    
    // Helper: convert 0-based index to 1-based (for API)
    uint32_t to1Based(uint32_t idx) const { return idx + 1; }
};

#endif
