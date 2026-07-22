#ifndef SLAM_VARIANCE_ANALYSIS_H
#define SLAM_VARIANCE_ANALYSIS_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <array>
#include <tuple>

// Per-position statistics for variance analysis
struct SlamPositionVarianceStats {
    uint64_t readCount = 0;           // Number of reads covering this position
    uint64_t tCount = 0;              // Number of T bases at this position
    uint64_t tcCount = 0;             // Number of T→C conversions at this position
    double qualSum = 0.0;              // Sum of quality scores
    uint64_t qualCount = 0;            // Count of quality scores
    double qualSumSq = 0.0;            // Sum of squares for variance calculation
    double tcRateSum = 0.0;           // Sum of T→C rates (0 or 1 per T base)
    double tcRateSumSq = 0.0;         // Sum of squares of T→C rates for variance
    uint64_t tcRateCount = 0;         // Count of T bases observed

    double meanQual() const {
        return qualCount > 0 ? qualSum / qualCount : 0.0;
    }

    double varianceQual() const {
        if (qualCount < 2) return 0.0;
        double mean = meanQual();
        return (qualSumSq / qualCount) - (mean * mean);
    }

    double stddevQual() const {
        return std::sqrt(varianceQual());
    }

    double meanTcRate() const {
        return tcRateCount > 0 ? tcRateSum / tcRateCount : 0.0;
    }

    double varianceTcRate() const {
        if (tcRateCount < 2) return 0.0;
        double mean = meanTcRate();
        return (tcRateSumSq / tcRateCount) - (mean * mean);
    }

    double stddevTcRate() const {
        return std::sqrt(varianceTcRate());
    }
};

struct SegmentFit {
    double slope = 0.0;
    double intercept = 0.0;
    double sse = 0.0;
};

// Per-mate segmented-regression snapshot (indexes 0=mate1, 1=mate2)
struct SlamVarianceMateSlice {
    int trim5p = 0;
    int trim3p = 0;
    uint32_t kneeBin5p = 0;
    uint32_t kneeBinMid = 0;
    uint32_t kneeBin3p = 0;
    double totalSSE = 0.0;
    SegmentFit seg1{};
    SegmentFit seg2{};
    SegmentFit seg3{};
    SegmentFit seg4{};
    std::vector<double> smoothedCurve;
    bool success = false;
};

// Variance analysis results (mate-aware PE; backward-compatible accessors)
struct SlamVarianceTrimResult {
    std::array<SlamVarianceMateSlice, 2> mates{};
    bool success = false;
    std::string mode;
    uint64_t readsAnalyzed = 0;

    int trim5p() const { return mates[0].trim5p; }
    int trim3p() const { return mates[0].trim3p; }
    SegmentFit segment1(size_t ix = 0) const {
        ix = ix < 2 ? ix : 0;
        return mates[ix].seg1;
    }

    SegmentFit segment2(size_t ix = 0) const {
        ix = ix < 2 ? ix : 0;
        return mates[ix].seg2;
    }
    SegmentFit segment3(size_t ix = 0) const {
        ix = ix < 2 ? ix : 0;
        return mates[ix].seg3;
    }
    SegmentFit segment4(size_t ix = 0) const {
        ix = ix < 2 ? ix : 0;
        return mates[ix].seg4;
    }

    uint32_t kneeBin5p = 0;   // mirrored from mates[0] for legacy QC callers
    uint32_t kneeBinMid = 0;
    uint32_t kneeBin3p = 0;
    double totalSSE = 0.0;
    SegmentFit seg1{};
    SegmentFit seg2{};
    SegmentFit seg3{};
    SegmentFit seg4{};
    std::vector<double> smoothedCurve;

    void syncLegacyFromMate0();
};

class SlamVarianceAnalyzer {
public:
    SlamVarianceAnalyzer(uint32_t maxReads = 100000, uint32_t minReads = 1000,
                         uint32_t smoothWindow = 5, uint32_t minSegLen = 3, uint32_t maxTrim = 15);

    void setSeparateMateHistograms(bool yes) { separateMateHistograms_ = yes; }
    bool separateMateHistograms() const { return separateMateHistograms_; }

    bool recordRead();

    void recordPosition(uint32_t readPos, uint8_t qual, bool isT, bool isTc);
    void recordPositionMate(uint32_t mateLocalPos, uint8_t mateIndex, uint8_t qual, bool isT, bool isTc);

    SlamVarianceTrimResult computeTrimCombined(uint32_t concatReadLength);
    SlamVarianceTrimResult computeTrimPerMate(uint32_t mateLen0, uint32_t mateLen1);
    SlamVarianceTrimResult computeTrimUnified(uint32_t concatLen, uint32_t mateLen0, uint32_t mateLen1, bool separateMates);

    const std::unordered_map<uint32_t, SlamPositionVarianceStats>& getStats() const {
        return positionCombined_;
    }

    const std::unordered_map<uint32_t, SlamPositionVarianceStats>& getStats(size_t mateIndex) const;

    uint64_t readsAnalyzed() const { return readsAnalyzed_; }

    std::tuple<uint64_t, uint64_t, double> computeGlobalTcErrorRate(int trim5p = 0, int trim3p = 0, uint32_t readLength = 0) const;
    std::tuple<uint64_t, uint64_t, double> computeGlobalTcErrorRatePerMate(
        const int trim5p[2],
        const int trim3p[2],
        uint32_t mateLen0,
        uint32_t mateLen1) const;

    void reset();
    void merge(const SlamVarianceAnalyzer& other);

private:
    bool fillTrimSlice(uint32_t readLength,
                       const std::unordered_map<uint32_t, SlamPositionVarianceStats>& posMap,
                       SlamVarianceMateSlice& slice,
                       std::string& mode) const;

    bool separateMateHistograms_ = false;
    std::unordered_map<uint32_t, SlamPositionVarianceStats> positionCombined_;
    std::array<std::unordered_map<uint32_t, SlamPositionVarianceStats>, 2> positionMate_{};
    uint64_t readsAnalyzed_;
    uint32_t maxReads_;
    uint32_t minReads_;
    uint32_t smoothWindow_;
    uint32_t minSegLen_;
    uint32_t maxTrim_;

    static std::vector<double> smoothMedian(const std::vector<double>& values, uint32_t window);
    static void interpolateMissing(std::vector<double>& values);

    static SegmentFit fitSegment(const std::vector<double>& prefixN,
                                 const std::vector<double>& prefixX,
                                  const std::vector<double>& prefixXX,
                                  const std::vector<double>& prefixY,
                                  const std::vector<double>& prefixXY,
                                  const std::vector<double>& prefixYY,
                                  const std::vector<double>& y,
                                  uint32_t start, uint32_t end);

    static std::tuple<uint32_t, uint32_t, uint32_t, double, SegmentFit, SegmentFit, SegmentFit, SegmentFit, uint32_t>
    segmentedRegression(const std::vector<double>& y, uint32_t minSegLen, uint32_t maxTrim);
};

inline void SlamVarianceTrimResult::syncLegacyFromMate0() {
    kneeBin5p = mates[0].kneeBin5p;
    kneeBinMid = mates[0].kneeBinMid;
    kneeBin3p = mates[0].kneeBin3p;
    totalSSE = mates[0].totalSSE;
    seg1 = mates[0].seg1;
    seg2 = mates[0].seg2;
    seg3 = mates[0].seg3;
    seg4 = mates[0].seg4;
    smoothedCurve = mates[0].smoothedCurve;
}

#endif // SLAM_VARIANCE_ANALYSIS_H
