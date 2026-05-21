#ifndef H_SoloReadInfoSink
#define H_SoloReadInfoSink

#include <cstdint>
#include <vector>
#include <unordered_map>
#include "SoloReadInfoLoader.h" // for ReadInfoRecord

class SoloFeature;
struct ReadInfoRecord;

typedef void (*SinkFinalizeFn)(SoloFeature&);

/** Compact row stored in CountingSink buckets (no readId/status; counting-only). */
struct CountingSinkRow {
    uint32_t featureId;
    uint32_t umi;
    uint32_t readRef; // readIndex, or low bits of readId when readIndex absent
};

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
    // Observed whitelist CB buckets only (not perWL[cbWLsize]).
    std::unordered_map<uint32_t, std::vector<CountingSinkRow>> buckets;
    std::unordered_map<uint32_t, uint32_t> readToCb; // optional conflict guard (debug)
    uint64_t totalRecords = 0;

    void onRecord(SoloFeature &feature, const ReadInfoRecord &rec) override;
    void finalize(SoloFeature &feature) override;
};

#endif
