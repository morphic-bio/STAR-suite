#ifndef H_FlexHashScreenTiered
#define H_FlexHashScreenTiered

#include "FlexHashScreenCommon.h"

// Two-tier cache:
//   h0_   — sample-specific H0 records (cacheClass==0, sampleIdx>0), sorted.
//            53K records, 1.3 MB, L2-hot. Searched first at all 3 offsets
//            with exact-match-only (no sampleIdx=0 fallback). On hit with
//            gene>0 → immediate KEEP. Resolves ~78% of reads.
//
//   rest_ — everything else (H0 globals + H1 + deny), sorted.
//            8.0M records, ~192 MB. Searched with full findRecord (with
//            sampleIdx=0 fallback) + classifyHits. Identical logic to flat.
class TieredCache {
public:
    void init(const std::vector<Record>& allRecords);
    FlexHashScreenDecision classifyRead(const char* readSeq, uint32_t readLen,
                                        uint16_t sampleIdx) const;

    size_t h0Count()      const { return h0Count_; }
    size_t h1Count()      const { return h1Count_; }
    size_t h2Count()      const { return h2Count_; }
    size_t denyCount()    const { return denyCount_; }
    size_t droppedCount() const { return dropped_; }
    size_t crossTierDuplicates() const { return crossTierDups_; }

private:
    bool findExact(const std::vector<Record>& arr,
                   uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx,
                   Record& out) const;

    bool findRecord(const std::vector<Record>& arr,
                    uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx,
                    Record& out) const;

    FlexHashScreenDecision classifyHits(const Record* const* hits,
                                        const int8_t* relativeOffsets,
                                        size_t nHits,
                                        uint16_t runtimeSampleIdx) const;

    std::vector<Record> h0_;     // sample-specific H0 (fast tier)
    std::vector<Record> rest_;   // everything else (slow tier)
    size_t h0Count_ = 0;        // all cacheClass==0 (for stats)
    size_t h1Count_ = 0;
    size_t h2Count_ = 0;
    size_t denyCount_ = 0;
    size_t dropped_ = 0;
    size_t crossTierDups_ = 0;
};

#endif
