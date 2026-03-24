#ifndef H_FlexHashScreenFlat
#define H_FlexHashScreenFlat

#include "FlexHashScreenCommon.h"

class FlatCache {
public:
    void init(const std::vector<Record>& allRecords);
    FlexHashScreenDecision classifyRead(const char* readSeq, uint32_t readLen,
                                        uint16_t sampleIdx) const;
    size_t recordCount() const { return records_.size(); }

private:
    bool findRecord(uint64_t seqLo, uint64_t seqHi, uint16_t sampleIdx,
                    Record& out) const;
    FlexHashScreenDecision classifyHits(const Record* const* hits,
                                        const int8_t* relativeOffsets,
                                        size_t nHits,
                                        uint16_t runtimeSampleIdx) const;
    std::vector<Record> records_;
};

#endif
