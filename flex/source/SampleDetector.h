#ifndef CODE_SampleDetector
#define CODE_SampleDetector

#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <array>
#include <mutex>

// Global canonical tag table (1-based, index 0 unused)
extern std::vector<std::string> gCanonicalTags;

class ParametersSolo;

// SampleDetector implements FLEX-like sample tag calling from read1
class SampleDetector {
public:
    explicit SampleDetector(const ParametersSolo &p);

    bool loadWhitelist(const std::string &path);
    bool loadProbes(const std::string &path);

    // Returns 1-based sample index on success, 0 if unmatched
    // Respects offsets and flags from ParametersSolo
    uint32_t detectSampleIndex(const uint8_t *seqData, int32_t readLength, bool reverseStrand) const;

    // Convenience helper: parse a BAM record and detect the sample index from its sequence.
    // The record pointer must reference the beginning of a BAM alignment (block_size field first).
    uint32_t detectSampleIndexFromBam(const char *bamRecord, uint32_t bamSize) const;

    // Get label by index (1-based); returns empty if idx invalid
    std::string labelFor(uint32_t sampleIdx) const;

    // Get canonical 8-mer sequence by index; empty if unavailable
    std::string canonicalFor(uint32_t sampleIdx) const;
    std::string canonicalForWhitelistIndex(uint32_t sampleIdx) const;

    bool ready() const { return !canonicalToIndex_.empty(); }

    // Map Solo's sequential sample index to the whitelist's canonical numeric index (BC###)
    uint32_t whitelistIndexForCanonical(uint32_t sampleIdx) const;
    uint32_t sequentialIndexForWhitelist(uint32_t sampleIdx) const;
    std::string labelForWhitelistIndex(uint32_t sampleIdx) const;

    // Token (5-bit) -> sequential index helpers
    static void registerSampleToken(uint8_t token, uint16_t sampleIdx);
    static uint16_t sampleIndexForToken(uint8_t token);

    // Global canonical lookup (1-based index) populated during whitelist load
    static std::string canonicalForIndexStatic(uint32_t sampleIdx);
    static std::string labelForIndexStatic(uint32_t sampleIdx);
    static void setCanonicalTable(const std::vector<std::string>& canon);
    static void setLabelTable(const std::vector<std::string>& labels);

private:
    const ParametersSolo &p_;
    std::vector<std::string> indexToLabel_;
    std::vector<std::string> indexToCanonical_;
    std::unordered_map<std::string,uint32_t> canonicalToIndex_;
    std::unordered_map<std::string,uint32_t> canonicalToWhitelistIndex_;
    std::unordered_map<uint32_t,uint32_t> whitelistIndexToSequential_;
    std::unordered_map<uint32_t,std::string> whitelistIndexToCanonical_;
    std::unordered_map<uint32_t,std::string> whitelistIndexToLabel_;
    std::vector<uint32_t> sampleCodes_;                  // nsamples * 8 encoded sequences (canonical + variants)
    std::vector<uint8_t> variantCountsPerSample_;        // actual number of sequences stored per sample (<=8)

    static std::array<uint16_t, 32> tokenToSampleIdx_;
    static std::mutex tokenLUTMutex_;
    static std::vector<std::string> canonicalByIdx_;
    static std::vector<std::string> labelsByIdx_;

    static inline bool isACGT8(const std::string &s) {
        if (s.size()!=8) return false;
        for (char c: s) { char u = (char)toupper(c); if (u!='A'&&u!='C'&&u!='G'&&u!='T') return false; }
        return true;
    }

    static inline bool isACGTNib(uint32_t nib) {
        return nib == 1u || nib == 2u || nib == 4u || nib == 8u;
    }

    static inline uint32_t complementNib(uint32_t nib) {
        switch (nib) {
            case 1u: return 8u; // A -> T
            case 2u: return 4u; // C -> G
            case 4u: return 2u; // G -> C
            case 8u: return 1u; // T -> A
            default: return 0u;
        }
    }

    static bool encodeStringToCode(const std::string &s, uint32_t &out);
};

#endif
