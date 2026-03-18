#ifndef H_FlexHashScreen
#define H_FlexHashScreen

#include <cstdint>
#include <string>
#include <vector>

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
    };

    static FlexHashScreenCache& instance();

    bool ensureLoaded(const ParametersSolo& pSolo, std::string* errorOut = nullptr);
    FlexHashScreenDecision classifyRead(const char* readSeq, uint32_t readLen, uint16_t sampleIdx) const;
    size_t recordCount() const { return records_.size(); }

private:
    FlexHashScreenCache() = default;

    bool loadFile(const std::string& path, std::string* errorOut);
    bool encodeWindow(const char* readSeq, uint32_t offset, uint64_t& seqLo, uint64_t& seqHi) const;
    bool findRecord(uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx, Record& out) const;
    FlexHashScreenDecision classifyHits(const Record* const* hits, const int8_t* relativeOffsets, size_t nHits, uint16_t runtimeSampleIdx) const;

    bool initialized_ = false;
    bool enabled_ = false;
    std::string loadedPath_;
    std::vector<Record> records_;
};

// Binary negative class codes from scripts/flex_h01_pilot.py.
enum FlexHashScreenNegativeCode : uint8_t {
    FlexHashNegNone = 0,
    FlexHashNegProbeAmbig = 1
};

#endif
