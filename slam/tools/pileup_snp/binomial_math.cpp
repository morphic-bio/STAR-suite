// Wrapper to expose binomial math functions from STAR's libem
// This file can be compiled separately and linked

#include <cmath>
#include <limits>
#include <cstdint>
#include <algorithm>

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
    double log_n_choose_k = 0.0;
    if (k > 0 && k < n) {
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

// Log-space binomial tail CDF: log(P[X >= k | n, p])
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
    if (k <= n / 2 + 1) {
        // Compute lower tail P[X <= k-1] with only k terms, then use complement
        double log_lower = -std::numeric_limits<double>::infinity();
        for (uint32_t i = 0; i < k; ++i) {
            double log_pmf = log_binom_pmf(n, i, p);
            log_lower = logsumexp(log_lower, log_pmf);
        }
        // P[X >= k] = 1 - P[X <= k-1]
        if (log_lower >= 0.0) {
            return -std::numeric_limits<double>::infinity();
        }
        double exp_lower = std::exp(log_lower);
        if (exp_lower >= 1.0) {
            return -std::numeric_limits<double>::infinity();
        }
        return std::log1p(-exp_lower);
    } else {
        // Compute upper tail directly
        double log_sum = -std::numeric_limits<double>::infinity();
        for (uint32_t i = k; i <= n; ++i) {
            double log_pmf = log_binom_pmf(n, i, p);
            log_sum = logsumexp(log_sum, log_pmf);
        }
        return log_sum;
    }
}
