#ifndef H_SoloReadInfoLoader
#define H_SoloReadInfoLoader

#include <functional>
#include <cstdint>
#include <vector>

class SoloReadFeature;
class SoloReadFlagClass;

struct ReadInfoRecord {
    uint64_t readId;
    uint32_t cbIdx;
    uint32_t umi;
    uint8_t status;
    uint32_t featureId;
    uint32_t readIndex; // (uint32_t)-1 if not present
};

enum class SoloReadInfoMode { Minimal, Counting };

using RecordSink = std::function<void(const ReadInfoRecord&)>;

class SoloReadInfoLoader {
public:
    void load(SoloReadFeature &rf,
              SoloReadInfoMode mode,
              const RecordSink &sink,
              std::vector<uint32_t> &cbReadCountTotal,
              SoloReadFlagClass &readFlagCounts,
              std::vector<uint32_t> &nReadPerCBunique1,
              std::vector<uint32_t> &nReadPerCBmulti1);

    // Convenience for Minimal mode
    void loadMinimal(SoloReadFeature &rf,
                     const RecordSink &sink,
                     std::vector<uint32_t> &cbReadCountTotal);
};

#endif

