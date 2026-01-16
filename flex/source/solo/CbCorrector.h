#ifndef CODE_CbCorrector
#define CODE_CbCorrector

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

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
    explicit CbCorrector(const std::vector<std::string> &whitelist, int maxHamming = 1);
    
    // Correct a cell barcode
    // Returns CbMatch with correction result
    CbMatch correct(const std::string &cb) const;
    
    // Get whitelist size
    size_t whitelistSize() const { return whitelist_.size(); }
    
    // Get whitelist sequences (for Bayesian resolver)
    const std::vector<std::string>& whitelist() const { return whitelist_; }
    
    // Get ambiguous variants map (for building hash maps)
    const std::unordered_map<uint32_t, std::vector<uint32_t>>& getAmbiguousVariants() const { return ambiguousVariants_; }
    
    // Get CB length
    size_t getCbLength() const { return cbLength_; }
    
    // Decode packed CB key to string (public helper)
    std::string decodePackedKey(uint32_t packedKey, size_t cbLength) const;
    
private:
    // Whitelist storage (canonical CB strings)
    std::vector<std::string> whitelist_;
    
    // Exact lookup: map packed CB key -> 0-based whitelist index
    std::unordered_map<uint32_t, uint32_t> exactMap_;
    
    // Variant hash: map 1-hamming variant key -> whitelist index
    // Value 0 means ambiguous (multiple WL entries generate this variant)
    std::unordered_map<uint32_t, uint32_t> variantMap_;
    
    // Ambiguous variants: map variant key -> list of WL indices (0-based)
    std::unordered_map<uint32_t, std::vector<uint32_t>> ambiguousVariants_;
    
    // CB length (assumed constant for all CBs in whitelist)
    size_t cbLength_;
    
    int maxHamming_;  // Maximum Hamming distance allowed
    
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

