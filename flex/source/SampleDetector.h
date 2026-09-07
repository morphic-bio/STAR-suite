#ifndef CODE_SampleDetector
#define CODE_SampleDetector

#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <array>
#include <atomic>
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

    // Match a pre-packed 8-base BAM-encoded tag (packed at position 0, length 8) against sample codes.
    // Used when the caller already extracted the 8 bases at sampleProbeOffset.
    // exactOnly restricts the lookup to listed sequences, for nearby offsets;
    // *ambiguous is set when one-mismatch rescue found owners tied at this
    // offset, which callers treat as final rather than searching on.
    uint32_t detectSampleFromPackedTag(const uint8_t *packedTag8, bool exactOnly = false, bool *ambiguous = nullptr) const;
    // Match eight A/C/G/T bases encoded with base 0 in the least-significant
    // two bits (A=0, C=1, G=2, T=3), as stored by CBQ.
    uint32_t detectSampleFromTwoBitTag(uint16_t packedTag8, bool exactOnly = false, bool *ambiguous = nullptr) const;
    // Tag lookup table statistics, for the run log: exact entries, one-mismatch
    // entries accepted, one-mismatch entries rejected as ambiguous.
    void tagTableStats(uint32_t &exact, uint32_t &mismatch1, uint32_t &ambiguous) const;
    uint32_t tagTableMaxOccupancy() const { return tagMaxOccupancy_; }

    // Convenience helper: parse a BAM record and detect the sample index from its sequence.
    // The record pointer must reference the beginning of a BAM alignment (block_size field first).
    uint32_t detectSampleIndexFromBam(const char *bamRecord, uint32_t bamSize) const;

    // Get label by index (1-based); returns empty if idx invalid
    std::string labelFor(uint32_t sampleIdx) const;

    // Get canonical 8-mer sequence by index; empty if unavailable
    std::string canonicalFor(uint32_t sampleIdx) const;
    std::string canonicalForWhitelistIndex(uint32_t sampleIdx) const;

    bool ready() const { return !canonicalToIndex_.empty(); }

    /** Number of sequential samples (1-based indices 1..N) from the whitelist; 0 if none. */
    uint32_t sequentialSampleCount() const { return static_cast<uint32_t>(indexToCanonical_.size()); }

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
    // Split the two-bit-packed tag into 3/3/2-base pieces. Each piece table is
    // keyed by the two unchanged pieces and stores a CSR list of listed-tag
    // ids. A tag within one substitution of a listed tag must occur in at
    // least one of these lists. The complete lookup is small enough to remain
    // cache-resident and does not allocate a 4^8 direct-address table.
    struct TagPieceTable {
        std::vector<uint16_t> offsets;
        std::vector<uint16_t> ids;
    };
    std::array<TagPieceTable, 3> tagPieceTables_;
    std::vector<uint16_t> tagCodes_;
    std::vector<uint16_t> tagOwners_;
    bool tagTableBuilt_ = false;
    uint32_t tagStatExact_ = 0, tagStatMismatch_ = 0, tagStatAmbiguous_ = 0;
    uint32_t tagMaxOccupancy_ = 0;
    void buildTagTable();
    uint32_t lookupIndex(uint16_t index, bool exactOnly, bool *ambiguous) const;
    uint32_t lookupCode(uint32_t code, bool exactOnly, bool *ambiguous) const;
    static inline uint16_t pieceKey(uint16_t index, unsigned table) {
        if (table == 0u) return static_cast<uint16_t>(index >> 6u); // P1,P2
        if (table == 1u) { // P0,P2
            return static_cast<uint16_t>((index & 0x003Fu) | ((index >> 6u) & 0x03C0u));
        }
        return static_cast<uint16_t>(index & 0x0FFFu); // P0,P1
    }
    static inline uint16_t codeToIndex(uint32_t code) {
        // code holds eight BAM nibbles, base 0 in the high nibble; the index holds
        // two bits per base, base 0 in the low bits, matching the CBQ packing.
        uint16_t index = 0;
        for (unsigned i = 0; i < 8; ++i) {
            const uint32_t nib = (code >> (4u * (7u - i))) & 0xFu;
            const uint16_t b = nib == 1u ? 0u : nib == 2u ? 1u : nib == 4u ? 2u : 3u;
            index = static_cast<uint16_t>(index | (b << (2u * i)));
        }
        return index;
    }

    // Token -> sample index. Written once per token (first writer wins, under
    // the mutex); read on every read of every thread, so reads are lock-free.
    static std::array<std::atomic<uint16_t>, 32> tokenToSampleIdx_;
    static std::mutex tokenLUTMutex_;
    static std::once_flag canonicalByIdxOnce_;
    static std::once_flag labelsByIdxOnce_;
    static std::once_flag canonicalTagsOnce_;
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

    uint32_t detectSampleCode(uint32_t code) const;

    static bool encodeStringToCode(const std::string &s, uint32_t &out);
};

#endif
