#include "slam_snp_em.h"
#include <cmath>
#include <algorithm>
#include <limits>

// Log-space binomial PMF: log(P(k | n, p))
double log_binom_pmf(uint32_t n, uint32_t k, double p) {
    if (n == 0) {
        return (k == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
    }
    if (k > n) {
        return -std::numeric_limits<double>::infinity();
    }
    if (p <= 0.0) {
        return (k == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
    }
    if (p >= 1.0) {
        return (k == n) ? 0.0 : -std::numeric_limits<double>::infinity();
    }
    
    // log(P(k|n,p)) = log(n choose k) + k*log(p) + (n-k)*log(1-p)
    // Use log-gamma for binomial coefficient: log(n!) - log(k!) - log((n-k)!)
    double log_n_choose_k = 0.0;
    if (k > 0 && k < n) {
        // lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
        log_n_choose_k = std::lgamma(static_cast<double>(n + 1)) -
                         std::lgamma(static_cast<double>(k + 1)) -
                         std::lgamma(static_cast<double>(n - k + 1));
    }
    
    double log_p = std::log(p);
    double log_one_minus_p = std::log1p(-p);
    
    return log_n_choose_k + static_cast<double>(k) * log_p + static_cast<double>(n - k) * log_one_minus_p;
}

// Log-sum-exp for two values
double logsumexp(double a, double b) {
    if (a == -std::numeric_limits<double>::infinity()) {
        return b;
    }
    if (b == -std::numeric_limits<double>::infinity()) {
        return a;
    }
    double max_val = std::max(a, b);
    return max_val + std::log1p(std::exp(std::min(a, b) - max_val));
}

// Log-sum-exp for three values
double logsumexp(double a, double b, double c) {
    return logsumexp(logsumexp(a, b), c);
}

// Log-space binomial tail CDF: log(P[X >= k | n, p])
// Optimized: compute the smaller tail to minimize iterations
// Complexity: O(min(k, n-k+1)) instead of O(n-k+1)
double log_binom_tail_cdf(uint32_t n, uint32_t k, double p) {
    if (n == 0) {
        return (k == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
    }
    if (k > n) {
        return -std::numeric_limits<double>::infinity();
    }
    if (k == 0) {
        return 0.0;  // P[X >= 0] = 1.0
    }
    if (p <= 0.0) {
        return (k == 0) ? 0.0 : -std::numeric_limits<double>::infinity();
    }
    if (p >= 1.0) {
        return 0.0;  // P[X >= k] = 1.0 when p=1
    }
    
    // Optimization: compute whichever tail is smaller
    // For SNP detection, k (mismatches) is typically small, so compute lower tail
    if (k <= n / 2 + 1) {
        // Compute lower tail P[X <= k-1] with only k terms, then use complement
        double log_lower = -std::numeric_limits<double>::infinity();
        for (uint32_t i = 0; i < k; ++i) {
            double log_pmf = log_binom_pmf(n, i, p);
            log_lower = logsumexp(log_lower, log_pmf);
        }
        // P[X >= k] = 1 - P[X <= k-1]
        // Use log1p(-exp(x)) for numerical stability
        if (log_lower >= 0.0) {
            return -std::numeric_limits<double>::infinity();  // Lower tail >= 1
        }
        // log(1 - exp(log_lower))
        double exp_lower = std::exp(log_lower);
        if (exp_lower >= 1.0) {
            return -std::numeric_limits<double>::infinity();
        }
        return std::log1p(-exp_lower);
    } else {
        // Compute upper tail directly with n-k+1 terms
        double log_sum = -std::numeric_limits<double>::infinity();
        for (uint32_t i = k; i <= n; ++i) {
            double log_pmf = log_binom_pmf(n, i, p);
            log_sum = logsumexp(log_sum, log_pmf);
        }
        return log_sum;
    }
}

SlamSnpEM::SlamSnpEM(uint32_t maxIter, double convergeRelLL)
    : maxIter_(maxIter), convergeRelLL_(convergeRelLL) {
}

double SlamSnpEM::initialize_p_ERR(const SnpHistogram& histogram) const {
    // Initialize from sites with f <= 0.01
    uint64_t total_k = 0;
    uint64_t total_n = 0;
    
    for (const auto& kv : histogram) {
        uint32_t packed = kv.first;
        uint32_t n = packed >> 16;
        uint32_t k = packed & 0xFFFF;
        uint64_t weight = kv.second;
        
        if (n > 0 && k <= n) {
            double f = static_cast<double>(k) / static_cast<double>(n);
            if (f <= 0.01) {
                total_k += static_cast<uint64_t>(k) * weight;
                total_n += static_cast<uint64_t>(n) * weight;
            }
        }
    }
    
    if (total_n == 0) {
        return 0.001;  // Fallback default
    }
    
    double p_init = static_cast<double>(total_k) / static_cast<double>(total_n);
    // Clamp to valid range
    if (p_init < p_ERR_MIN_) {
        p_init = p_ERR_MIN_;
    }
    if (p_init > p_ERR_MAX_) {
        p_init = p_ERR_MAX_;
    }
    return p_init;
}

double SlamSnpEM::compute_log_likelihood(const SnpHistogram& histogram) const {
    double ll = 0.0;
    
    for (const auto& kv : histogram) {
        uint32_t packed = kv.first;
        uint32_t n = packed >> 16;
        uint32_t k = packed & 0xFFFF;
        uint64_t weight = kv.second;
        
        if (n == 0) {
            continue;
        }
        
        // Log-likelihood: log(sum_c pi_c * P(k|n, p_c))
        double log_prob_ERR = log_binom_pmf(n, k, p_ERR_);
        double log_prob_HET = log_binom_pmf(n, k, p_HET_);
        double log_prob_HOM = log_binom_pmf(n, k, p_HOM_);
        
        double log_mixture = logsumexp(
            std::log(pi_ERR_) + log_prob_ERR,
            std::log(pi_HET_) + log_prob_HET,
            std::log(pi_HOM_) + log_prob_HOM
        );
        
        ll += log_mixture * static_cast<double>(weight);
    }
    
    return ll;
}

void SlamSnpEM::compute_responsibilities(const SnpHistogram& histogram,
                                         std::unordered_map<uint32_t, double>& r_ERR,
                                         std::unordered_map<uint32_t, double>& r_HET,
                                         std::unordered_map<uint32_t, double>& r_HOM) const {
    r_ERR.clear();
    r_HET.clear();
    r_HOM.clear();
    
    for (const auto& kv : histogram) {
        uint32_t packed = kv.first;
        uint32_t n = packed >> 16;
        uint32_t k = packed & 0xFFFF;
        
        if (n == 0) {
            continue;
        }
        
        double log_prob_ERR = log_binom_pmf(n, k, p_ERR_);
        double log_prob_HET = log_binom_pmf(n, k, p_HET_);
        double log_prob_HOM = log_binom_pmf(n, k, p_HOM_);
        
        double log_weight_ERR = std::log(pi_ERR_) + log_prob_ERR;
        double log_weight_HET = std::log(pi_HET_) + log_prob_HET;
        double log_weight_HOM = std::log(pi_HOM_) + log_prob_HOM;
        
        double log_sum = logsumexp(log_weight_ERR, log_weight_HET, log_weight_HOM);
        
        // Responsibilities in log-space, then normalize
        r_ERR[packed] = std::exp(log_weight_ERR - log_sum);
        r_HET[packed] = std::exp(log_weight_HET - log_sum);
        r_HOM[packed] = std::exp(log_weight_HOM - log_sum);
    }
}

void SlamSnpEM::update_parameters(const SnpHistogram& histogram,
                                   const std::unordered_map<uint32_t, double>& r_ERR,
                                   const std::unordered_map<uint32_t, double>& r_HET,
                                   const std::unordered_map<uint32_t, double>& r_HOM) {
    // Update mixture weights
    double total_weight = 0.0;
    double weight_ERR = 0.0;
    double weight_HET = 0.0;
    double weight_HOM = 0.0;
    
    for (const auto& kv : histogram) {
        uint32_t packed = kv.first;
        uint64_t site_count = kv.second;
        
        double r_err = r_ERR.count(packed) > 0 ? r_ERR.at(packed) : 0.0;
        double r_het = r_HET.count(packed) > 0 ? r_HET.at(packed) : 0.0;
        double r_hom = r_HOM.count(packed) > 0 ? r_HOM.at(packed) : 0.0;
        
        double w = static_cast<double>(site_count);
        total_weight += w;
        weight_ERR += r_err * w;
        weight_HET += r_het * w;
        weight_HOM += r_hom * w;
    }
    
    if (total_weight > 0.0) {
        pi_ERR_ = weight_ERR / total_weight;
        pi_HET_ = weight_HET / total_weight;
        pi_HOM_ = weight_HOM / total_weight;
    }
    
    // Update p_ERR (MLE for error component)
    double sum_k = 0.0;
    double sum_n = 0.0;
    
    for (const auto& kv : histogram) {
        uint32_t packed = kv.first;
        uint32_t n = packed >> 16;
        uint32_t k = packed & 0xFFFF;
        uint64_t site_count = kv.second;
        
        double r_err = r_ERR.count(packed) > 0 ? r_ERR.at(packed) : 0.0;
        double w = static_cast<double>(site_count) * r_err;
        
        sum_k += static_cast<double>(k) * w;
        sum_n += static_cast<double>(n) * w;
    }
    
    if (sum_n > 0.0) {
        p_ERR_ = sum_k / sum_n;
        // Clamp to valid range
        if (p_ERR_ < p_ERR_MIN_) {
            p_ERR_ = p_ERR_MIN_;
        }
        if (p_ERR_ > p_ERR_MAX_) {
            p_ERR_ = p_ERR_MAX_;
        }
    }
}

SlamSnpEMResult SlamSnpEM::fit(const SnpHistogram& histogram) {
    SlamSnpEMResult result;
    
    if (histogram.empty()) {
        result.converged = false;
        return result;
    }
    
    // Initialize parameters
    p_ERR_ = initialize_p_ERR(histogram);
    pi_ERR_ = 0.999;
    pi_HET_ = 0.0009;
    pi_HOM_ = 0.0001;
    
    double prev_ll = compute_log_likelihood(histogram);
    result.iterations = 0;
    
    for (uint32_t iter = 0; iter < maxIter_; ++iter) {
        result.iterations = iter + 1;
        
        // E-step: compute responsibilities
        std::unordered_map<uint32_t, double> r_ERR, r_HET, r_HOM;
        compute_responsibilities(histogram, r_ERR, r_HET, r_HOM);
        
        // M-step: update parameters
        update_parameters(histogram, r_ERR, r_HET, r_HOM);
        
        // Check convergence
        double curr_ll = compute_log_likelihood(histogram);
        double rel_change = std::abs((curr_ll - prev_ll) / (std::abs(prev_ll) + 1e-10));
        
        if (rel_change < convergeRelLL_) {
            result.converged = true;
            result.final_log_likelihood = curr_ll;
            break;
        }
        
        prev_ll = curr_ll;
    }
    
    if (!result.converged) {
        result.final_log_likelihood = compute_log_likelihood(histogram);
    }
    
    result.p_ERR = p_ERR_;
    result.pi_ERR = pi_ERR_;
    result.pi_HET = pi_HET_;
    result.pi_HOM = pi_HOM_;
    
    return result;
}

double SlamSnpEM::posterior(uint32_t n, uint32_t k) const {
    if (n == 0) {
        return 0.0;
    }
    
    double log_prob_ERR = log_binom_pmf(n, k, p_ERR_);
    double log_prob_HET = log_binom_pmf(n, k, p_HET_);
    double log_prob_HOM = log_binom_pmf(n, k, p_HOM_);
    
    double log_weight_ERR = std::log(pi_ERR_) + log_prob_ERR;
    double log_weight_HET = std::log(pi_HET_) + log_prob_HET;
    double log_weight_HOM = std::log(pi_HOM_) + log_prob_HOM;
    
    double log_sum = logsumexp(log_weight_ERR, log_weight_HET, log_weight_HOM);
    
    // Posterior for SNP-like (HET or HOM)
    double log_post_snp = logsumexp(log_weight_HET, log_weight_HOM) - log_sum;
    
    return std::exp(log_post_snp);
}
