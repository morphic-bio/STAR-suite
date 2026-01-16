#ifndef GC_BIAS_H
#define GC_BIAS_H

#include <vector>
#include <string>
#include <array>
#include <cstdint>

// GC Fragment Bias Model
// Tracks GC distribution and computes bias ratios
class GCFragModel {
public:
    static constexpr int GC_BINS = 101;  // 0-100% GC
    
    GCFragModel() : normalized_(false) {
        counts_.fill(0.0);
    }
    
    // Add weighted observation (weight is typically log-probability)
    void inc(int32_t gc_pct, double weight);
    
    // Add weighted observation (linear space)
    void incLinear(int32_t gc_pct, double weight);
    
    // Normalize counts to probabilities
    void normalize();
    
    // Get normalized probability for a GC percentage
    double get(int32_t gc_pct) const;
    
    // Compute bias ratio: observed / expected
    // Returns vector of bias ratios, clamped to [1/maxRatio, maxRatio]
    std::vector<double> computeBiasRatio(const GCFragModel& expected, double maxRatio = 1000.0) const;
    
    // Combine counts from another model (for multi-threaded accumulation)
    void combineCounts(const GCFragModel& other);
    
    // Reset to empty state
    void reset();
    
    // I/O
    bool writeToFile(const std::string& path) const;
    bool loadFromFile(const std::string& path);
    
    // Get raw counts (for debugging)
    const std::array<double, GC_BINS>& getCounts() const { return counts_; }
    
private:
    std::array<double, GC_BINS> counts_;
    bool normalized_;
    
    // Helper to clamp GC percentage to valid range
    int32_t clampGC(int32_t gc_pct) const {
        if (gc_pct < 0) return 0;
        if (gc_pct >= GC_BINS) return GC_BINS - 1;
        return gc_pct;
    }
};

#endif // GC_BIAS_H

