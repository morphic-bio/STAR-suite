#ifndef TRIM_QC_H
#define TRIM_QC_H

#include <array>
#include <cstdint>
#include <vector>
#include <cstddef>

struct TrimQcMateStats {
    uint64_t reads = 0;
    std::vector<uint64_t> lengthHist;
    std::vector<uint64_t> posCount;
    std::vector<uint64_t> posQualSum;
    std::vector<std::array<uint64_t, 5>> posBaseCounts; // A,C,G,T,N
    std::array<uint64_t, 101> gcHist{};

    void ensureSize(uint32_t len);
    void addRead(const char* seqNum, const char* qual, uint32_t len,
                 uint32_t qualOffset, uint32_t qualBase);
    void merge(const TrimQcMateStats& other);
};

class TrimQcCollector {
public:
    void init(uint32_t mateCount, uint64_t maxReads, uint32_t qualBase);
    bool enabled() const { return enabled_; }
    uint64_t totalReads() const { return totalReads_; }
    uint32_t qualityBase() const { return qualityBase_; }
    size_t mateCount() const { return mates_.size(); }
    const TrimQcMateStats& mate(size_t i) const { return mates_[i]; }
    void addRead(const char* seqNum, const char* qual, uint32_t len,
                 uint32_t mate, uint32_t qualOffset);
    void merge(const TrimQcCollector& other);

private:
    bool enabled_ = false;
    uint64_t maxReads_ = 0;
    uint64_t totalReads_ = 0;
    uint32_t qualityBase_ = 33;
    std::vector<TrimQcMateStats> mates_;
};

#endif
