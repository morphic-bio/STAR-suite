#ifndef CODE_InlineCBCorrection
#define CODE_InlineCBCorrection

#include "IncludeDefine.h"
#include "ParametersSolo.h"
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <mutex>
#include <memory>
#include <limits>
#include <stdexcept>
#include "solo/CbBayesianResolver.h"
#include "hash_shims_cpp_compat.h" // brings in khash and pack helpers

// Hash type for packed CB -> whitelist index (1-based)
KHASH_MAP_INIT_INT64(cbwl, uint32_t)
struct InlineCbCollisionRange { size_t offset; uint32_t count; };
KHASH_INIT(inlineCbRanges, uint64_t, InlineCbCollisionRange, 1, kh_int64_hash_func, kh_int64_hash_equal)

// Inline CB/UB correction module - ports process_features-style correction logic
// to replace Solo's CB/UB data structures with streaming correction during R1 parsing

class InlineCBCorrection {
public:
    // Initialize whitelist hash from pSolo.cbWLstr
    // Uses 2 bits per base encoding (LSB-first)
    static void initializeWhitelist(const ParametersSolo &pSolo);
    
    // Fast-path correction: exact match or single variant
    // Returns: 0 = exact match, 1 = single variant corrected, -1 = no match/ambiguous
    // On success, correctedCB contains the corrected barcode string
    static int fastPathCorrection(const std::string &cbSeq, std::string &correctedCB);
    // Packed equivalent for two-bit inputs. correctedIndex1 is the 1-based
    // whitelist index on success. This performs the same exact/unique-H1
    // lookup without constructing or decoding a barcode string.
    static int fastPathCorrectionPacked(uint64_t packedCB, uint32_t &correctedIndex1);
    
    // Find closest barcodes with single mismatch
    // Returns number of matches found (0, 1, or >1)
    static int findClosestBarcodes(const std::string &cbSeq, std::vector<std::string> &matches);
    
    // Check if sequence contains N and correct if possible (Phase 2)
    // Returns number of valid candidates found (0, 1, or >1)
    // If exactly one candidate found, correctedSeq contains the corrected sequence
    static int checkSequenceAndCorrectForN(const std::string &seq, int maxN, std::string &correctedSeq);
    
    // Clear whitelist hash (for cleanup)
    static void clearWhitelist();
    
private:
    // Encode DNA string to packed uint64 (2 bits per base, LSB-first)
    // For CB: typically 16bp = 32 bits = uint32_t
    // Returns UINT64_MAX on invalid input
    static uint64_t encodeCB(const std::string &cbSeq);
    
    // Decode packed uint64 back to DNA string
    static std::string decodeCB(uint64_t packed, int length);
    
    // Find variant matches at a specific position
    static int findVariantMatch(uint64_t packedCode, int position, std::vector<uint64_t> &variants);
    
    // Whitelist lookups (khash for speed/cache locality)
    // exactHash_: packed CB -> 1-based whitelist index (exact)
    // variantHash_: packed 1MM variant -> 1-based index if unique, 0 if ambiguous
    static khash_t(cbwl) *exactHash_;    // packed CB -> 1-based WL index
    static khash_t(cbwl) *variantHash_;  // packed 1MM variant -> 1-based WL index if unique, 0 if ambiguous
    // Collision tracker: variant -> parent whitelist indices (when ambiguous).
    struct CollisionDeleter {
        void operator()(khash_t(inlineCbRanges)* h) const { kh_destroy(inlineCbRanges, h); }
    };
    static std::unique_ptr<khash_t(inlineCbRanges), CollisionDeleter> variantCollisions_;
    static std::vector<uint32_t> collisionParents_;
    struct CandidateView {
        const uint32_t* data;
        size_t count;
        const uint32_t* begin() const { return data; }
        const uint32_t* end() const { return count ? data + count : data; }
        size_t size() const { return count; }
        bool empty() const { return count == 0; }
    };
    static CandidateView collisionCandidates(uint64_t packed);
    template<class Visitor> static void forEachWhitelistVariant(Visitor visit) {
        for (khiter_t it = kh_begin(exactHash_); it != kh_end(exactHash_); ++it) {
            if (!kh_exist(exactHash_, it)) continue;
            uint64_t packed = kh_key(exactHash_, it);
            uint32_t index = kh_val(exactHash_, it);
            const size_t length = (*wlStrings_)[index - 1].size();
            for (size_t pos = 0; pos < length; ++pos) {
                uint32_t base = (packed >> (pos * 2)) & 3;
                for (uint32_t alternate = 0; alternate < 4; ++alternate) {
                    if (alternate == base) continue;
                    uint64_t variant = (packed & ~(uint64_t{3} << (pos * 2)))
                        | (uint64_t{alternate} << (pos * 2));
                    if (kh_get(cbwl, exactHash_, variant) == kh_end(exactHash_))
                        visit(variant, index);
                }
            }
        }
    }
    static const std::vector<std::string>* wlStrings_;  // pointer to whitelist strings
    static bool whitelistInitialized_;

    // Evidence counters: registry of per-thread WL count vectors (heap-owned, freed on clear)
    static std::vector<std::vector<uint64_t>*> evidenceShards_;
    static std::mutex evidenceMutex_;
    static size_t evidenceSize_;

    // Ambiguous CB capture (Stage 1 stub)
    struct AmbigRecord {
        uint64_t packedVariant;
        std::string cbSeq;
        std::string cbQual;
        std::vector<double> cbLogLikMatch;
        std::vector<double> cbLogLikMismatch;
        uint32_t evidenceReads = 0;
        std::vector<uint32_t> parents; // WL indices (1-based)
        uint32_t count;
    };
    struct AmbigShard {
        std::vector<AmbigRecord> records;
        std::unordered_map<uint64_t, size_t> index; // packedVariant -> records idx
    };
    static std::vector<AmbigShard*> ambigShards_;
    static std::mutex ambigMutex_;
    static uint64_t ambigTotal_;

public:
    static size_t exactMapSize();
    static size_t variantMapSize();
    static size_t variantCollisionSize();
    static size_t variantCollisionMaxFanout();
    static uint64_t parentEvidenceTotal();
    static uint64_t ambigCapturedTotal();
    static size_t ambigUniqueVariants();
    static size_t ambigMaxParents();

    // Utility helpers for inline path
    static uint64_t packCBForLookup(const std::string &cbSeq); // packed or UINT64_MAX on error
    static bool isAmbiguousVariant(uint64_t packedVariant);    // true if variantHash_ marks it ambiguous
    static uint32_t exactIndex(uint64_t packed);                // 1-based WL index or 0
    static uint32_t exactIndex(const std::string &cbSeq);       // helper from string
    static void recordParentEvidence(uint32_t wlIndex1);        // increment evidence for a resolved parent WL index
    static void recordAmbiguousCB(uint64_t packedVariant, const std::string &cbSeq, const std::string &cbQual);
    // Resolve a single ambiguous variant immediately; returns resolved WL idx (1-based) or 0
    static uint32_t resolveAmbiguousVariant(uint64_t packedVariant,
                                            const std::string &cbSeq,
                                            const std::string &cbQual,
                                            const ParametersSolo &pSolo);
    static uint64_t getResolvedAmbigCount();
    static void clearEvidence();
    static void clearAmbiguous();

    // Phase 2: merge ambiguous shards into a single map for downstream resolver
    struct MergedAmbigEntry {
        std::string cbSeq;
        std::string cbQual;
        std::vector<double> cbLogLikMatch;
        std::vector<double> cbLogLikMismatch;
        uint32_t evidenceReads = 0;
        std::vector<uint32_t> parents;
        uint32_t count;
    };
    static void mergeAmbiguousShards(std::unordered_map<uint64_t, MergedAmbigEntry> &out);

    // Optional: resolve merged ambiguous map via Bayesian resolver
    struct AmbigResolveStats {
        uint64_t total = 0;
        uint64_t resolved = 0;
        uint64_t ambiguous = 0;
        uint64_t unresolved = 0;
    };
    static AmbigResolveStats resolveAmbiguousMerged(const std::unordered_map<uint64_t, MergedAmbigEntry> &merged,
                                                    const ParametersSolo &pSolo,
                                                    std::vector<uint32_t> &resolvedIdx);
};

#endif
