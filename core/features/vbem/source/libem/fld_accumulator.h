#ifndef FLD_ACCUMULATOR_H
#define FLD_ACCUMULATOR_H

#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>
#include <atomic>
#include <mutex>

// Log-space constants
constexpr double LOG_0_FLD = -std::numeric_limits<double>::infinity();
constexpr double LOG_EPSILON_FLD = -1e10;  // Salmon's LOG_EPSILON equivalent

// Fragment Length Distribution Accumulator
// Salmon-compatible implementation: works in log-space, applies kernel at add time
// Based on Salmon's FragmentLengthDistribution.cpp
class FLDAccumulator {
public:
    static constexpr int32_t MAX_FRAG_LEN = 1000;  // Maximum fragment length (Salmon: fragLenDistMax)
    
    // Gaussian prior parameters (Salmon defaults from SalmonDefaults.hpp)
    static constexpr double PRIOR_MEAN = 250.0;   // fragLenPriorMean
    static constexpr double PRIOR_SD = 25.0;      // fragLenPriorSD
    static constexpr double ALPHA = 1.0;          // Prior weight (logged pseudo-counts)
    
    // Binomial kernel parameters (Salmon: AlignmentLibrary.hpp)
    static constexpr size_t KERNEL_N = 4;         // Number of trials
    static constexpr double KERNEL_P = 0.5;       // Success probability
    
    FLDAccumulator();
    
    // Add fragment length observation (with logged mass)
    // This applies the binomial kernel at add time (Salmon-style)
    void add(int32_t frag_len, double mass = 0.0);  // mass is in LOG space (default = log(1) = 0)
    
    // Combine counts from another accumulator (for multi-threaded accumulation)
    // Note: This is an approximation - Salmon uses atomics for true thread-safety
    void combine(const FLDAccumulator& other);
    
    // Get log PMF (probability mass function)
    // Returns log(P(frag_len))
    double pmf(int32_t frag_len) const;
    
    // Get log CMF (cumulative mass function)
    double cmf(int32_t frag_len) const;
    
    // Cache the CMF for faster lookups (call once after burn-in)
    void cacheCMF();
    
    // Legacy API for compatibility (implemented in cpp for proper linkage)
    double getLogProb(int32_t frag_len) const;
    double getLogCMF(int32_t frag_len) const;
    double getProb(int32_t frag_len) const;
    std::vector<double> getPMF() const;
    const std::vector<double>& getSmoothedPMF() const { return cachedPMF_; }
    
    // Get mean fragment length (in real space)
    double getMean() const;
    
    // Get standard deviation
    double getStdDev() const;
    
    // Get total mass (logged)
    double totMass() const { return totMass_; }
    
    // Check if FLD is valid (has observations beyond prior)
    bool isValid() const { return totMass_ > std::log(ALPHA) + 0.1; }
    
    // Get total fragment count (approximation from log-space)
    double getTotalFragments() const { return std::exp(totMass_); }
    
    // Min observed length
    size_t minVal() const { return min_ == MAX_FRAG_LEN ? 1 : min_; }
    
    // Max value
    size_t maxVal() const { return MAX_FRAG_LEN; }
    
    // Reset to initial state (with Gaussian prior)
    void reset();
    
    // Dump PMF for debugging
    void dumpPMF(std::vector<double>& pmfOut, size_t& minV, size_t& maxV) const;
    
    // Convert to string (for debugging)
    std::string toString() const;
    
private:
    // Histogram in LOG space (like Salmon's hist_)
    std::vector<double> hist_;
    
    // Cached CMF/PMF (after cacheCMF() is called)
    std::vector<double> cachedCMF_;
    std::vector<double> cachedPMF_;
    bool haveCachedCMF_;
    mutable std::mutex cacheMutex_;
    
    // Total mass (logged)
    double totMass_;
    
    // Sum of (len * mass) for mean calculation (logged)
    double sum_;
    
    // Minimum observed length
    size_t min_;
    
    // Binomial kernel (precomputed, logged)
    std::vector<double> kernel_;
    
    // Log-add helper (Salmon-style)
    static double logAdd(double a, double b) {
        if (a == LOG_0_FLD) return b;
        if (b == LOG_0_FLD) return a;
        if (a > b) {
            return a + std::log1p(std::exp(b - a));
        } else {
            return b + std::log1p(std::exp(a - b));
        }
    }
    
    // Atomic log-add (for thread-safe accumulation)
    static void incLoopLog(double& target, double val) {
        // Salmon uses atomic<double> with compare-exchange
        // For non-atomic version, just use logAdd
        target = logAdd(target, val);
    }
    
    // Binomial PDF
    static double binomialPDF(size_t n, double p, size_t k);
};

#endif // FLD_ACCUMULATOR_H
