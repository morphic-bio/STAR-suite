#ifndef SLAM_SNP_EM_H
#define SLAM_SNP_EM_H

#include <cstdint>
#include <unordered_map>
#include <string>

// Log-space binomial PMF: log(P(k | n, p))
double log_binom_pmf(uint32_t n, uint32_t k, double p);

// Log-space binomial tail CDF: log(P[X >= k | n, p])
double log_binom_tail_cdf(uint32_t n, uint32_t k, double p);

// Histogram key: (n, k) packed as uint32_t = (n << 16) | k
// Value: weight (number of sites with this (n,k) pair)
using SnpHistogram = std::unordered_map<uint32_t, uint64_t>;

struct SlamSnpEMResult {
    double p_ERR = 0.0;           // Learned error rate
    double pi_ERR = 0.0;          // Mixture weight for error component
    double pi_HET = 0.0;          // Mixture weight for heterozygous component
    double pi_HOM = 0.0;          // Mixture weight for homozygous component
    uint32_t iterations = 0;     // Number of EM iterations
    double final_log_likelihood = 0.0;
    bool converged = false;
};

class SlamSnpEM {
public:
    SlamSnpEM(uint32_t maxIter, double convergeRelLL);
    
    // Fit EM model to histogram
    SlamSnpEMResult fit(const SnpHistogram& histogram);
    
    // Compute posterior probability that site is SNP-like (HET or HOM)
    double posterior(uint32_t n, uint32_t k) const;
    
    // Get current model parameters
    double get_p_ERR() const { return p_ERR_; }
    double get_pi_ERR() const { return pi_ERR_; }
    double get_pi_HET() const { return pi_HET_; }
    double get_pi_HOM() const { return pi_HOM_; }
    
private:
    static constexpr double p_HET_ = 0.5;   // Fixed heterozygous rate
    static constexpr double p_HOM_ = 0.95;  // Fixed homozygous rate
    static constexpr double p_ERR_MIN_ = 1e-6;
    static constexpr double p_ERR_MAX_ = 0.05;
    
    uint32_t maxIter_;
    double convergeRelLL_;
    
    double p_ERR_ = 0.0;
    double pi_ERR_ = 0.0;
    double pi_HET_ = 0.0;
    double pi_HOM_ = 0.0;
    
    // Initialize p_ERR from low-mismatch sites
    double initialize_p_ERR(const SnpHistogram& histogram) const;
    
    // Compute log-likelihood
    double compute_log_likelihood(const SnpHistogram& histogram) const;
    
    // E-step: compute responsibilities
    void compute_responsibilities(const SnpHistogram& histogram,
                                   std::unordered_map<uint32_t, double>& r_ERR,
                                   std::unordered_map<uint32_t, double>& r_HET,
                                   std::unordered_map<uint32_t, double>& r_HOM) const;
    
    // M-step: update parameters
    void update_parameters(const SnpHistogram& histogram,
                           const std::unordered_map<uint32_t, double>& r_ERR,
                           const std::unordered_map<uint32_t, double>& r_HET,
                           const std::unordered_map<uint32_t, double>& r_HOM);
};

// Utility: log-space binomial PMF
double log_binom_pmf(uint32_t n, uint32_t k, double p);

// Utility: log-sum-exp for numerical stability
double logsumexp(double a, double b);
double logsumexp(double a, double b, double c);

#endif // SLAM_SNP_EM_H
