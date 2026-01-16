#ifndef H_SoloReadInfoSink
#define H_SoloReadInfoSink

#include <cstdint>
#include <vector>
#include "SoloReadInfoLoader.h" // for ReadInfoRecord

class SoloFeature;
struct ReadInfoRecord;

typedef void (*SinkFinalizeFn)(SoloFeature&);

class ISoloReadInfoSink {
public:
    virtual ~ISoloReadInfoSink() = default;
    virtual void onRecord(SoloFeature &feature, const ReadInfoRecord &rec) = 0;
    virtual void finalize(SoloFeature &feature) = 0;
};

class MinimalSink : public ISoloReadInfoSink {
public:
    void onRecord(SoloFeature &feature, const ReadInfoRecord &rec) override;
    void finalize(SoloFeature &feature) override;
};

class CountingSink : public ISoloReadInfoSink {
public:
    std::vector<std::vector<ReadInfoRecord>> perWL; // buffered per-WL CB records
    void onRecord(SoloFeature &feature, const ReadInfoRecord &rec) override;
    void finalize(SoloFeature &feature) override;
};

#endif
