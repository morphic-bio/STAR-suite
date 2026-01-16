#ifndef SLAM_VARIANCE_ANALYSIS_H
#define SLAM_VARIANCE_ANALYSIS_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>

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
    
    // Computed statistics
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
    
    // Variance of T→C rate across reads at this position
    double varianceTcRate() const {
        if (tcRateCount < 2) return 0.0;
        double mean = meanTcRate();
        return (tcRateSumSq / tcRateCount) - (mean * mean);
    }
    
    double stddevTcRate() const {
        return std::sqrt(varianceTcRate());
    }
};

// Segmented regression fit info for a single segment
struct SegmentFit {
    double slope = 0.0;
    double intercept = 0.0;
    double sse = 0.0;
};

// Variance analysis results
struct SlamVarianceTrimResult {
    int trim5p = 0;
    int trim3p = 0;
    bool success = false;
    std::string mode;                  // "auto_segmented", "auto_fallback", "manual"
    uint64_t readsAnalyzed = 0;
    uint32_t kneeBin5p = 0;            // breakpoint b1 (for backward compat)
    uint32_t kneeBinMid = 0;           // breakpoint b2 (4-seg mid)
    uint32_t kneeBin3p = 0;            // breakpoint b2 or b3 (for backward compat)
    double totalSSE = 0.0;             // total SSE of segmented fit
    SegmentFit seg1, seg2, seg3, seg4; // fit info for each segment
    std::vector<double> smoothedCurve; // smoothed stdev curve for QC output
};

// Variance analyzer for auto-trim
class SlamVarianceAnalyzer {
public:
    SlamVarianceAnalyzer(uint32_t maxReads = 100000, uint32_t minReads = 1000,
                         uint32_t smoothWindow = 5, uint32_t minSegLen = 3, uint32_t maxTrim = 15);
    
    // Record a read (call once per read before recording positions)
    bool recordRead();
    
    // Record a position observation (only if readsAnalyzed_ < maxReads_)
    void recordPosition(uint32_t readPos, uint8_t qual, bool isT, bool isTc);
    
    // Compute trim values using segmented regression on smoothed stdev curve
    SlamVarianceTrimResult computeTrim(uint32_t readLength);
    
    // Get per-position statistics
    const std::unordered_map<uint32_t, SlamPositionVarianceStats>& getStats() const {
        return positionStats_;
    }
    
    // Get number of reads analyzed
    uint64_t readsAnalyzed() const { return readsAnalyzed_; }
    
    // Compute global T→C error rate (p_err) from position statistics
    // If trim5p/trim3p are provided (>0), restrict to positions within trimmed window
    // Returns: (t_total, tc_total, p_est)
    std::tuple<uint64_t, uint64_t, double> computeGlobalTcErrorRate(int trim5p = 0, int trim3p = 0, uint32_t readLength = 0) const;
    
    // Reset for new file
    void reset();
    
    // Merge stats from another analyzer (for thread aggregation)
    void merge(const SlamVarianceAnalyzer& other);
    
private:
    std::unordered_map<uint32_t, SlamPositionVarianceStats> positionStats_;
    uint64_t readsAnalyzed_;
    uint32_t maxReads_;
    uint32_t minReads_;
    uint32_t smoothWindow_;   // Median smoothing window (default 5)
    uint32_t minSegLen_;      // Minimum segment length (default 3)
    uint32_t maxTrim_;        // Maximum trim at either end (default 15)
    
    // Smooth a curve using median window
    static std::vector<double> smoothMedian(const std::vector<double>& values, uint32_t window);
    
    // Linear interpolation for missing values (NaN)
    static void interpolateMissing(std::vector<double>& values);
    
    // Fit a linear segment and return slope, intercept, SSE using prefix sums
    static SegmentFit fitSegment(const std::vector<double>& prefixN,
                                 const std::vector<double>& prefixX,
                                  const std::vector<double>& prefixXX,
                                  const std::vector<double>& prefixY,
                                  const std::vector<double>& prefixXY,
                                  const std::vector<double>& prefixYY,
                                  const std::vector<double>& y,
                                  uint32_t start, uint32_t end);
    
    // Segmented regression with 2 breakpoints
    static std::tuple<uint32_t, uint32_t, uint32_t, double, SegmentFit, SegmentFit, SegmentFit, SegmentFit, uint32_t>
    segmentedRegression(const std::vector<double>& y, uint32_t minSegLen, uint32_t maxTrim);
};

#endif // SLAM_VARIANCE_ANALYSIS_H
