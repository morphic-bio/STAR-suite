#include "TrimQc.h"
#include <algorithm>

void TrimQcMateStats::ensureSize(uint32_t len) {
    if (len == 0) {
        return;
    }
    if (lengthHist.size() <= len) {
        lengthHist.resize(len + 1, 0);
    }
    if (posCount.size() < len) {
        size_t oldSize = posCount.size();
        posCount.resize(len, 0);
        posQualSum.resize(len, 0);
        posBaseCounts.resize(len);
        for (size_t i = oldSize; i < posBaseCounts.size(); ++i) {
            posBaseCounts[i].fill(0);
        }
    }
}

void TrimQcMateStats::addRead(const char* seqNum, const char* qual, uint32_t len,
                              uint32_t qualOffset, uint32_t qualBase) {
    if (len == 0) {
        return;
    }
    ensureSize(len);
    ++reads;
    lengthHist[len]++;

    uint32_t gcCount = 0;
    for (uint32_t i = 0; i < len; ++i) {
        uint8_t base = static_cast<uint8_t>(seqNum[i]);
        uint32_t baseIdx = base < 5 ? base : 4;
        posBaseCounts[i][baseIdx]++;

        if (baseIdx == 1 || baseIdx == 2) { // C or G
            ++gcCount;
        }

        uint8_t q = static_cast<uint8_t>(qual[qualOffset + i]);
        int32_t phred = static_cast<int32_t>(q) - static_cast<int32_t>(qualBase);
        if (phred < 0) {
            phred = 0;
        }
        posQualSum[i] += static_cast<uint64_t>(phred);
        posCount[i]++;
    }

    uint32_t gcPct = static_cast<uint32_t>(
        std::min(100.0, std::max(0.0, (100.0 * gcCount) / static_cast<double>(len)))
    );
    gcHist[gcPct]++;
}

void TrimQcMateStats::merge(const TrimQcMateStats& other) {
    reads += other.reads;
    if (other.lengthHist.size() > lengthHist.size()) {
        lengthHist.resize(other.lengthHist.size(), 0);
    }
    for (size_t i = 0; i < other.lengthHist.size(); ++i) {
        lengthHist[i] += other.lengthHist[i];
    }
    if (other.posCount.size() > posCount.size()) {
        size_t oldSize = posCount.size();
        posCount.resize(other.posCount.size(), 0);
        posQualSum.resize(other.posQualSum.size(), 0);
        posBaseCounts.resize(other.posBaseCounts.size());
        for (size_t i = oldSize; i < posBaseCounts.size(); ++i) {
            posBaseCounts[i].fill(0);
        }
    }
    for (size_t i = 0; i < other.posCount.size(); ++i) {
        posCount[i] += other.posCount[i];
        posQualSum[i] += other.posQualSum[i];
        for (size_t b = 0; b < 5; ++b) {
            posBaseCounts[i][b] += other.posBaseCounts[i][b];
        }
    }
    for (size_t i = 0; i < gcHist.size(); ++i) {
        gcHist[i] += other.gcHist[i];
    }
}

void TrimQcCollector::init(uint32_t mateCount, uint64_t maxReads, uint32_t qualBase) {
    enabled_ = true;
    maxReads_ = maxReads;
    totalReads_ = 0;
    qualityBase_ = qualBase;
    mates_.clear();
    mates_.resize(mateCount);
}

void TrimQcCollector::addRead(const char* seqNum, const char* qual, uint32_t len,
                              uint32_t mate, uint32_t qualOffset) {
    if (!enabled_ || mate >= mates_.size()) {
        return;
    }
    if (maxReads_ > 0 && totalReads_ >= maxReads_) {
        return;
    }
    mates_[mate].addRead(seqNum, qual, len, qualOffset, qualityBase_);
    ++totalReads_;
}

void TrimQcCollector::merge(const TrimQcCollector& other) {
    if (!other.enabled_) {
        return;
    }
    if (!enabled_) {
        *this = other;
        return;
    }
    if (other.mateCount() > mates_.size()) {
        mates_.resize(other.mateCount());
    }
    for (size_t i = 0; i < other.mateCount(); ++i) {
        mates_[i].merge(other.mate(i));
    }
    totalReads_ += other.totalReads_;
}
