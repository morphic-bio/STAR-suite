#ifndef H_FlexHashScreen
#define H_FlexHashScreen

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "FlexGdna.h"

class ParametersSolo;

struct FlexHashScreenDecision {
    enum Action : uint8_t {
        Disabled = 0,
        Pass = 1,
        Keep = 2,
        Deny = 3
    };

    Action action = Disabled;
    uint16_t geneIdx15 = 0;
    uint8_t cacheClass = 0;
    uint8_t negativeCode = 0;
    int8_t offset = 0;
    FlexGdnaRegion probeRegion = FlexGdnaUnknown;
};

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

    bool ensureLoaded(const ParametersSolo& pSolo, std::string* errorOut = nullptr);
    FlexHashScreenDecision classifyRead(const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const;
    FlexHashScreenDecision classifyReadH0Only(const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const;
    FlexHashScreenDecision classifyReadH0Offset0(const char* readSeq, uint32_t readLen) const;
    FlexHashScreenDecision classifyReadH0H1Offset0(const char* readSeq, uint32_t readLen) const;
    // Classify a 50-base offset-0 key in CBQ's native LSB-first order. lo
    // contains bases 0..31 and hi bases 32..49. In CBQ runs the lookup maps
    // are built in this order once, rather than reordering every read.
    FlexHashScreenDecision classifyCbqH0Offset0(uint64_t seqLo, uint64_t seqHi) const;
    FlexHashScreenDecision classifyCbqH0H1Offset0(uint64_t seqLo, uint64_t seqHi) const;
    static uint32_t probeWindowLength();
    size_t recordCount() const { return records_.size(); }
    size_t h0RecordCount() const { return h0Records_.size(); }
    size_t h1DenyRecordCount() const { return h1DenyRecords_.size(); }
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
    bool findRecord(uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx, Record& out) const;
    FlexHashScreenDecision classifyHits(const Record* const* hits, const int8_t* relativeOffsets, size_t nHits, uint16_t runtimeSampleIdx) const;

    bool findRecordInVec(const std::vector<Record>& vec, uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx, Record& out) const;
    void buildTieredVectors();

    struct SeqKeyNoSample {
        uint64_t lo, hi;
        bool operator==(const SeqKeyNoSample& o) const { return lo == o.lo && hi == o.hi; }
    };
    struct SeqKeyNoSampleHash {
        size_t operator()(const SeqKeyNoSample& k) const {
            return k.lo * 0x9E3779B97F4A7C15ULL ^ k.hi * 0xBF58476D1CE4E5B9ULL;
        }
    };
    using H0NoSampleMap = std::unordered_map<SeqKeyNoSample, Record, SeqKeyNoSampleHash>;

    bool encodeWindowLUT(const char* readSeq, uint32_t offset, uint64_t& seqLo, uint64_t& seqHi) const;
    static SeqKeyNoSample cacheKeyToCbqKey(uint64_t seqLo, uint64_t seqHi);
    static SeqKeyNoSample cbqKeyToCacheKey(uint64_t seqLo, uint64_t seqHi);
    SeqKeyNoSample offset0MapKeyFromCacheKey(uint64_t seqLo, uint64_t seqHi) const;
    FlexHashScreenDecision classifyH0Offset0MapKey(const SeqKeyNoSample& key) const;
    FlexHashScreenDecision classifyH0H1Offset0MapKey(const SeqKeyNoSample& key) const;
    void buildH0NoSampleMap();
    void buildH1DenyNoSampleMap();

    bool initialized_ = false;
    bool enabled_ = false;
    bool regionMetadataComplete_ = false;
    uint16_t cacheVersion_ = 0;
    std::string loadedPath_;
    std::vector<Record> records_;
    std::vector<Record> h0Records_;
    std::vector<Record> h1DenyRecords_;
    H0NoSampleMap h0NoSampleMap_;
    H0NoSampleMap h1DenyNoSampleMap_;
    bool offset0MapsUseCbqOrder_ = false;
    static uint8_t baseLUT_[256];
    static bool lutInitialized_;
};

// Binary negative class codes from scripts/flex_h01_pilot.py.
enum FlexHashScreenNegativeCode : uint8_t {
    FlexHashNegNone = 0,
    FlexHashNegProbeAmbig = 1
};

#endif
