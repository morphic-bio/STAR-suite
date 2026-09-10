#ifndef H_FlexHashScreen
#define H_FlexHashScreen

#include <cstdint>
#include <string>
#include "FlexHashCacheStorage.h"
#include <vector>

#include "FlexGdna.h"
#include "FlexProbeHalfIndex.h"

class ParametersSolo;

#include "FlexHashScreenDecision.h"
#include "FlexProbePairKhash.h"

class FlexHashScreenCache {
public:
    struct Record {
        uint64_t seqLo = 0;
        uint64_t seqHi = 0;
        uint32_t resolvedGeneIdx15 = 0;
        uint8_t cacheClass = 0;
        uint8_t negativeCode = 0;
        uint16_t sampleIdx = 0;
        FlexGdnaRegion probeRegion = FlexGdnaUnknown;
    };

    static FlexHashScreenCache& instance();

    /** Load a cache already sorted by (seqHi, seqLo, sampleIdx), as emitted
     *  by the cache writers. Loading does not sort or check record ordering. */
    bool ensureLoaded(const ParametersSolo& pSolo, std::string* errorOut = nullptr);
    FlexHashScreenDecision classifyRead(const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const;
    FlexHashScreenDecision classifyReadH0Only(const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const;
    FlexHashScreenDecision classifyReadH0Offset0(const char* readSeq, uint32_t readLen) const;
    FlexHashScreenDecision classifyReadH0H1Offset0(const char* readSeq, uint32_t readLen) const;
    // Classify a 50-base offset-0 key in CBQ's native LSB-first order. lo
    // contains bases 0..31 and hi bases 32..49. Lookup maps use this order
    // in both input modes and in persisted khash files.
    FlexHashScreenDecision classifyCbqH0Offset0(uint64_t seqLo, uint64_t seqHi) const;
    FlexHashScreenDecision classifyCbqH0H1Offset0(uint64_t seqLo, uint64_t seqHi) const;
    /** Window with exactly one N (nMask has one bit): try the four bases at that
     *  position through the H0/H1 lookup. One gene and no DENY -> Keep; anything
     *  else is a miss (Pass), never a certified negative. Two or more Ns -> Pass. */
    FlexHashScreenDecision classifyCbqH0H1Offset0SingleN(uint64_t seqLo, uint64_t seqHi, uint64_t nMask) const;
    FlexHashScreenDecision classifyReadH0H1Offset0SingleN(const char* readSeq, uint32_t readLen) const;
    FlexHashScreenDecision classifyReadH1X2SeedExtend(const char* readSeq,
                                                      uint32_t readLen) const;
    FlexHashScreenDecision classifyCbqH1X2SeedExtend(uint64_t seqLo,
                                                     uint64_t seqHi,
                                                     uint64_t nMask) const;
    bool hasHalfTables() const { return pair_.ready(); }
    // Fused primary/single-N/extension path reuses half-match parent lists.
    FlexHashScreenDecision classifyReadComplete(const char* read, uint32_t len, bool singleN = true) const;
    FlexHashScreenDecision classifyCbqComplete(uint64_t lo, uint64_t hi, uint64_t nmask, bool singleN = true) const;
    bool hasH1X2() const { return hasH1X2_; }
    bool h1x2ProbeIndexReady() const { return pair_.ready() || probeSeedIndex_.ready(); }
    size_t h1x2ProbeCount() const { return pair_.ready() ? pair_.probeCount() : probeSeedIndex_.probeCount(); }
    size_t h1x2ProbeKeyCount() const {
        return pair_.ready() ? pair_.keys(0) + pair_.keys(1) :
            probeSeedIndex_.leftKeyCount() + probeSeedIndex_.rightKeyCount();
    }
    static uint32_t probeWindowLength();
    size_t recordCount() const { return pair_.ready() ? pair_.recordCount() : storage_.recordCount(); }
    size_t h0RecordCount() const { return pair_.ready() ? pair_.h0Count() : storage_.h0Count(); }
    size_t h1DenyRecordCount() const { return pair_.ready() ? pair_.recordCount()-pair_.h0Count() : storage_.h1Count(); }
    bool hasPersistedTables() const { return pair_.ready() || storage_.persisted(); }
    uint16_t cacheVersion() const { return cacheVersion_; }
    bool hasRegionMetadata() const { return regionMetadataComplete_; }

    /** Pack 50bp ACGT window into (seqLo, seqHi); same encoding as classifyRead. */
    static bool encodeProbeWindow(const char* readSeq, uint32_t offset, uint64_t& seqLo, uint64_t& seqHi);
    /** Write FH01SEQ1 v2/v3 cache (sorted by seqHi, seqLo, sampleIdx). */
    static bool writeHashCacheFile(const std::string& path, std::vector<Record>& records,
                                   std::string* errorOut, bool regionMetadataComplete = false);

private:
    FlexHashScreenCache() = default;

    bool loadFile(const std::string& path, std::string* errorOut);
    bool encodeWindow(const char* readSeq, uint32_t offset, uint64_t& seqLo, uint64_t& seqHi) const;
    bool findRecord(uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx, Record& out, bool h0Only = false) const;
    FlexHashScreenDecision classifyHits(const Record* const* hits, const int8_t* relativeOffsets, size_t nHits, uint16_t runtimeSampleIdx) const;

    using SeqKeyNoSample = FlexProbeKey;
    Record decodeStorageRecord(const FlexProbeRecord& raw) const;

    bool encodeWindowLUT(const char* readSeq, uint32_t offset, uint64_t& seqLo, uint64_t& seqHi) const;
    static SeqKeyNoSample cacheKeyToCbqKey(uint64_t seqLo, uint64_t seqHi);
    static SeqKeyNoSample cbqKeyToCacheKey(uint64_t seqLo, uint64_t seqHi);
    SeqKeyNoSample offset0MapKeyFromCacheKey(uint64_t seqLo, uint64_t seqHi) const;
    FlexHashScreenDecision classifyH0Offset0MapKey(const SeqKeyNoSample& key) const;
    FlexHashScreenDecision classifyH0H1Offset0MapKey(const SeqKeyNoSample& key) const;
    bool buildH1X2ProbeIndex(std::string* errorOut);
    static std::string decodeCacheSequence(uint64_t seqLo, uint64_t seqHi);

    bool initialized_ = false;
    bool enabled_ = false;
    bool regionMetadataComplete_ = false;
    uint16_t cacheVersion_ = 0;
    std::string loadedPath_;
    FlexHashCacheStorage storage_;
    FlexProbePairKhash pair_;
    bool hasH1X2_ = false;
    FlexProbeHalfIndex probeSeedIndex_;
    static uint8_t baseLUT_[256];
    static bool lutInitialized_;
};

const char* flexHashScreenDenyReason(uint8_t negativeCode);

#endif
