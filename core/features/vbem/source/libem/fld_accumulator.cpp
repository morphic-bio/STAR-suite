/**
 * FLDAccumulator - Salmon-compatible Fragment Length Distribution
 * 
 * This implementation matches Salmon's FragmentLengthDistribution.cpp:
 * - Works in log space
 * - Applies binomial kernel at add() time (not at PMF computation time)
 * - Uses Gaussian prior for initialization
 * 
 * Key difference from previous STAR implementation:
 * - Previous: stored raw counts, applied kernel at getPMF() time
 * - Current: stores log-mass, applies kernel at add() time (like Salmon)
 */

#include "fld_accumulator.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>

// Binomial PDF computation
double FLDAccumulator::binomialPDF(size_t n, double p, size_t k) {
    if (k > n) return 0.0;
    
    // Compute binomial coefficient: n choose k
    double coeff = 1.0;
    for (size_t i = 0; i < k; ++i) {
        coeff = coeff * (n - i) / (i + 1);
    }
    
    // PDF = C(n,k) * p^k * (1-p)^(n-k)
    return coeff * std::pow(p, k) * std::pow(1.0 - p, n - k);
}

FLDAccumulator::FLDAccumulator() 
    : cachedCMF_(MAX_FRAG_LEN + 1, LOG_0_FLD),
      cachedPMF_(MAX_FRAG_LEN + 1, LOG_0_FLD),
      haveCachedCMF_(false),
      totMass_(LOG_0_FLD),
      sum_(LOG_0_FLD),
      min_(MAX_FRAG_LEN) {
    
    // Initialize histogram with Gaussian prior (like Salmon)
    hist_.resize(MAX_FRAG_LEN + 1);
    
    double tot = std::log(ALPHA);  // Prior weight in log space
    
    // Set histogram to Gaussian prior
    if (PRIOR_MEAN > 0.0) {
        // Normal distribution parameters
        double mu = PRIOR_MEAN;
        double sigma = PRIOR_SD;
        
        for (size_t i = 0; i <= MAX_FRAG_LEN; ++i) {
            // Compute Gaussian CDF difference (like Salmon)
            // P(i-0.5 < X < i+0.5) for Normal(mu, sigma)
            double norm_mass = 0.0;
            
            // Using erf for CDF: P(X < x) = 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2))))
            double upper = (i + 0.5 - mu) / (sigma * std::sqrt(2.0));
            double lower = (i - 0.5 - mu) / (sigma * std::sqrt(2.0));
            norm_mass = 0.5 * (std::erf(upper) - std::erf(lower));
            
            double mass = LOG_EPSILON_FLD;
            if (norm_mass > 0.0) {
                mass = tot + std::log(norm_mass);
            }
            hist_[i] = mass;
            
            // Update sum_ (for mean calculation)
            if (i > 0) {
                sum_ = logAdd(sum_, std::log(static_cast<double>(i)) + mass);
            }
            totMass_ = logAdd(totMass_, mass);
        }
    } else {
        // Uniform prior (fallback)
        double uniform_mass = tot - std::log(static_cast<double>(MAX_FRAG_LEN));
        std::fill(hist_.begin(), hist_.end(), uniform_mass);
        hist_[0] = LOG_0_FLD;  // No zero-length fragments
        
        // Compute totMass_ and sum_ for uniform
        totMass_ = tot;
        // sum_ = uniform_mass + log(sum(i for i in 1..MAX_FRAG_LEN))
        sum_ = uniform_mass + std::log(static_cast<double>(MAX_FRAG_LEN * (MAX_FRAG_LEN + 1)) / 2.0);
    }
    
    // Precompute binomial kernel (like Salmon)
    // kernel_[i] = log(PDF(Binomial(KERNEL_N, KERNEL_P), i))
    kernel_.resize(KERNEL_N + 1);
    for (size_t i = 0; i <= KERNEL_N; ++i) {
        double pdf = binomialPDF(KERNEL_N, KERNEL_P, i);
        kernel_[i] = (pdf > 0.0) ? std::log(pdf) : LOG_0_FLD;
    }
}

void FLDAccumulator::add(int32_t frag_len, double mass) {
    // Clamp to valid range
    if (frag_len < 0) return;
    size_t len = static_cast<size_t>(frag_len);
    if (len > MAX_FRAG_LEN) {
        len = MAX_FRAG_LEN;
    }
    
    // Update min
    if (len < min_) {
        min_ = len;
    }
    
    // Apply binomial kernel (Salmon-style)
    // Kernel is centered at len, spans [len - kernel_.size()/2, len + kernel_.size()/2]
    size_t offset = (len >= kernel_.size() / 2) ? (len - kernel_.size() / 2) : 0;
    
    for (size_t i = 0; i < kernel_.size(); ++i) {
        if (offset > 0 && offset < hist_.size()) {
            double kMass = mass + kernel_[i];  // Log-space: mass * kernel[i]
            incLoopLog(hist_[offset], kMass);
            incLoopLog(sum_, std::log(static_cast<double>(offset)) + kMass);
            incLoopLog(totMass_, kMass);
        }
        offset++;
    }
    
    // Invalidate cache
    haveCachedCMF_ = false;
}

void FLDAccumulator::combine(const FLDAccumulator& other) {
    // Combine histograms (log-add each bin)
    for (size_t i = 0; i < hist_.size() && i < other.hist_.size(); ++i) {
        hist_[i] = logAdd(hist_[i], other.hist_[i]);
    }
    totMass_ = logAdd(totMass_, other.totMass_);
    sum_ = logAdd(sum_, other.sum_);
    min_ = std::min(min_, other.min_);
    haveCachedCMF_ = false;
}

double FLDAccumulator::pmf(int32_t frag_len) const {
    if (haveCachedCMF_) {
        size_t len = (frag_len >= 0) ? static_cast<size_t>(frag_len) : 0;
        return (len < cachedPMF_.size()) ? cachedPMF_[len] : cachedPMF_.back();
    } else {
        if (frag_len < 0) return LOG_0_FLD;
        size_t len = static_cast<size_t>(frag_len);
        if (len > MAX_FRAG_LEN) {
            len = MAX_FRAG_LEN;
        }
        return hist_[len] - totMass_;  // Normalize
    }
}

double FLDAccumulator::cmf(int32_t frag_len) const {
    if (haveCachedCMF_) {
        size_t len = (frag_len >= 0) ? static_cast<size_t>(frag_len) : 0;
        return (len < cachedCMF_.size()) ? cachedCMF_[len] : cachedCMF_.back();
    } else {
        if (frag_len < 0) return LOG_0_FLD;
        size_t len = static_cast<size_t>(frag_len);
        if (len > MAX_FRAG_LEN) {
            len = MAX_FRAG_LEN;
        }
        
        double cum = LOG_0_FLD;
        for (size_t i = 0; i <= len; ++i) {
            cum = logAdd(cum, hist_[i]);
        }
        return cum - totMass_;  // Normalize
    }
}

void FLDAccumulator::cacheCMF() {
    std::lock_guard<std::mutex> lock(cacheMutex_);
    
    if (haveCachedCMF_) return;
    
    // Compute PMF (normalized)
    double totMassNorm = LOG_0_FLD;
    cachedPMF_.resize(MAX_FRAG_LEN + 1);
    for (size_t i = 0; i <= MAX_FRAG_LEN; ++i) {
        cachedPMF_[i] = hist_[i] - totMass_;
        totMassNorm = logAdd(totMassNorm, cachedPMF_[i]);
    }
    
    // Renormalize to ensure sum = 1 in log space (sum of exp = 1)
    for (size_t i = 0; i <= MAX_FRAG_LEN; ++i) {
        cachedPMF_[i] -= totMassNorm;
    }
    
    // Compute CMF
    cachedCMF_.resize(MAX_FRAG_LEN + 1);
    double cum = LOG_0_FLD;
    for (size_t i = 0; i <= MAX_FRAG_LEN; ++i) {
        cum = logAdd(cum, cachedPMF_[i]);
        cachedCMF_[i] = cum;
    }
    
    haveCachedCMF_ = true;
}

std::vector<double> FLDAccumulator::getPMF() const {
    // Return PMF in linear space (for compatibility)
    std::vector<double> pmfOut(MAX_FRAG_LEN + 1);
    for (size_t i = 0; i <= MAX_FRAG_LEN; ++i) {
        pmfOut[i] = std::exp(pmf(static_cast<int32_t>(i)));
    }
    return pmfOut;
}

// Legacy API for compatibility
double FLDAccumulator::getLogProb(int32_t frag_len) const {
    return pmf(frag_len);
}

double FLDAccumulator::getLogCMF(int32_t frag_len) const {
    return cmf(frag_len);
}

double FLDAccumulator::getProb(int32_t frag_len) const {
    return std::exp(pmf(frag_len));
}

double FLDAccumulator::getMean() const {
    // mean = sum / totMass (in log space, then exp)
    return std::exp(sum_ - totMass_);
}

double FLDAccumulator::getStdDev() const {
    // Compute standard deviation from log-space histogram
    // stddev = sqrt(E[X^2] - E[X]^2)
    double mean = getMean();
    
    // Compute E[X^2] by iterating over PMF
    double sum_sq = LOG_0_FLD;
    for (size_t i = 0; i <= MAX_FRAG_LEN; ++i) {
        if (hist_[i] > LOG_EPSILON_FLD) {
            // log(i^2 * P(i)) = 2*log(i) + log(P(i))
            double log_i_sq = (i > 0) ? 2.0 * std::log(static_cast<double>(i)) : LOG_0_FLD;
            sum_sq = logAdd(sum_sq, log_i_sq + hist_[i] - totMass_);
        }
    }
    
    double e_x_sq = std::exp(sum_sq);
    double variance = e_x_sq - mean * mean;
    
    return (variance > 0.0) ? std::sqrt(variance) : 0.0;
}

void FLDAccumulator::reset() {
    // Re-initialize with Gaussian prior
    // Note: Cannot use assignment since mutex is not copyable/moveable
    // Instead, manually reset all members
    
    double tot = std::log(ALPHA);  // Prior weight in log space
    
    // Reset with Gaussian prior
    if (PRIOR_MEAN > 0.0) {
        double mu = PRIOR_MEAN;
        double sigma = PRIOR_SD;
        
        totMass_ = LOG_0_FLD;
        sum_ = LOG_0_FLD;
        
        for (size_t i = 0; i <= MAX_FRAG_LEN; ++i) {
            double upper = (i + 0.5 - mu) / (sigma * std::sqrt(2.0));
            double lower = (i - 0.5 - mu) / (sigma * std::sqrt(2.0));
            double norm_mass = 0.5 * (std::erf(upper) - std::erf(lower));
            
            double mass = LOG_EPSILON_FLD;
            if (norm_mass > 0.0) {
                mass = tot + std::log(norm_mass);
            }
            hist_[i] = mass;
            
            if (i > 0) {
                sum_ = logAdd(sum_, std::log(static_cast<double>(i)) + mass);
            }
            totMass_ = logAdd(totMass_, mass);
        }
    } else {
        double uniform_mass = tot - std::log(static_cast<double>(MAX_FRAG_LEN));
        std::fill(hist_.begin(), hist_.end(), uniform_mass);
        hist_[0] = LOG_0_FLD;
        totMass_ = tot;
        sum_ = uniform_mass + std::log(static_cast<double>(MAX_FRAG_LEN * (MAX_FRAG_LEN + 1)) / 2.0);
    }
    
    min_ = MAX_FRAG_LEN;
    haveCachedCMF_ = false;
    std::fill(cachedCMF_.begin(), cachedCMF_.end(), LOG_0_FLD);
    std::fill(cachedPMF_.begin(), cachedPMF_.end(), LOG_0_FLD);
}

void FLDAccumulator::dumpPMF(std::vector<double>& pmfOut, size_t& minV, size_t& maxV) const {
    minV = minVal();
    maxV = maxVal();
    pmfOut.clear();
    pmfOut.reserve(maxV - minV + 1);
    for (size_t i = minV; i <= maxV; ++i) {
        pmfOut.push_back(pmf(static_cast<int32_t>(i)));  // Log PMF
    }
}

std::string FLDAccumulator::toString() const {
    std::stringstream ss;
    for (size_t i = 0; i <= MAX_FRAG_LEN; ++i) {
        ss << std::exp(pmf(static_cast<int32_t>(i)));
        if (i != MAX_FRAG_LEN) {
            ss << '\t';
        }
    }
    ss << '\n';
    return ss.str();
}

