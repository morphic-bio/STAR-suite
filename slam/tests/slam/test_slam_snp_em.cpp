/**
 * Unit tests for SLAM SNP EM model (3-component binomial mixture).
 *
 * Tests:
 * - log_binom_pmf edge cases and correctness
 * - logsumexp numerical stability
 * - EM convergence on synthetic data
 * - Posterior computation correctness
 *
 * Compile:
 *   g++ -std=c++11 -O2 -I../../source -I../../source/libem \
 *     test_slam_snp_em.cpp ../../source/libem/slam_snp_em.cpp \
 *     -o test_slam_snp_em -lm
 *
 * Run:
 *   ./test_slam_snp_em
 */

#include "slam_snp_em.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <cassert>

static const double EPSILON = 1e-9;

static bool approxEqual(double a, double b, double tol = 1e-6) {
    if (std::isinf(a) && std::isinf(b)) {
        return (a < 0) == (b < 0);  // Same sign infinity
    }
    return std::abs(a - b) < tol;
}

static int check(bool ok, const std::string& label) {
    if (!ok) {
        std::cerr << "FAIL: " << label << "\n";
        return 1;
    }
    std::cerr << "PASS: " << label << "\n";
    return 0;
}

// Test log_binom_pmf edge cases
int test_log_binom_pmf_edges() {
    int failed = 0;
    
    // n=0, k=0 -> log(1) = 0
    failed += check(approxEqual(log_binom_pmf(0, 0, 0.5), 0.0), "log_binom_pmf(0,0,0.5) = 0");
    
    // n=0, k>0 -> -inf
    failed += check(std::isinf(log_binom_pmf(0, 1, 0.5)) && log_binom_pmf(0, 1, 0.5) < 0,
                    "log_binom_pmf(0,1,0.5) = -inf");
    
    // k > n -> -inf
    failed += check(std::isinf(log_binom_pmf(5, 10, 0.5)) && log_binom_pmf(5, 10, 0.5) < 0,
                    "log_binom_pmf(5,10,0.5) = -inf");
    
    // p=0, k=0 -> log(1) = 0
    failed += check(approxEqual(log_binom_pmf(10, 0, 0.0), 0.0), "log_binom_pmf(10,0,0) = 0");
    
    // p=0, k>0 -> -inf
    failed += check(std::isinf(log_binom_pmf(10, 5, 0.0)) && log_binom_pmf(10, 5, 0.0) < 0,
                    "log_binom_pmf(10,5,0) = -inf");
    
    // p=1, k=n -> log(1) = 0
    failed += check(approxEqual(log_binom_pmf(10, 10, 1.0), 0.0), "log_binom_pmf(10,10,1) = 0");
    
    // p=1, k<n -> -inf
    failed += check(std::isinf(log_binom_pmf(10, 5, 1.0)) && log_binom_pmf(10, 5, 1.0) < 0,
                    "log_binom_pmf(10,5,1) = -inf");
    
    return failed;
}

// Test log_binom_pmf correctness against known values
int test_log_binom_pmf_values() {
    int failed = 0;
    
    // B(10, 5, 0.5) = C(10,5) * 0.5^10 = 252 * (1/1024) = 252/1024
    // log(252/1024) = log(252) - log(1024) ≈ -1.404
    double expected = std::log(252.0 / 1024.0);
    failed += check(approxEqual(log_binom_pmf(10, 5, 0.5), expected, 1e-4),
                    "log_binom_pmf(10,5,0.5) correctness");
    
    // B(20, 1, 0.1) = C(20,1) * 0.1 * 0.9^19 = 20 * 0.1 * 0.9^19
    double expected2 = std::log(20.0) + std::log(0.1) + 19 * std::log(0.9);
    failed += check(approxEqual(log_binom_pmf(20, 1, 0.1), expected2, 1e-4),
                    "log_binom_pmf(20,1,0.1) correctness");
    
    // Verify probabilities sum to ~1 (in regular space)
    double sum = 0.0;
    for (uint32_t k = 0; k <= 10; ++k) {
        sum += std::exp(log_binom_pmf(10, k, 0.3));
    }
    failed += check(approxEqual(sum, 1.0, 1e-6), "B(10,k,0.3) sums to 1");
    
    return failed;
}

// Test logsumexp numerical stability
int test_logsumexp() {
    int failed = 0;
    
    // Basic case
    double a = std::log(2.0);
    double b = std::log(3.0);
    double expected = std::log(5.0);
    failed += check(approxEqual(logsumexp(a, b), expected), "logsumexp(log2,log3) = log5");
    
    // Large difference (should not overflow)
    double large = 1000.0;
    double small = -1000.0;
    failed += check(approxEqual(logsumexp(large, small), large, 1e-10),
                    "logsumexp(1000,-1000) ≈ 1000");
    
    // Both -inf
    double neg_inf = -std::numeric_limits<double>::infinity();
    // logsumexp(-inf, -inf) should be -inf
    double result = logsumexp(neg_inf, neg_inf);
    failed += check(std::isinf(result) && result < 0, "logsumexp(-inf,-inf) = -inf");
    
    // One -inf
    failed += check(approxEqual(logsumexp(0.0, neg_inf), 0.0), "logsumexp(0,-inf) = 0");
    
    // Three-argument version
    double c = std::log(5.0);
    double expected3 = std::log(10.0);  // 2 + 3 + 5 = 10
    failed += check(approxEqual(logsumexp(a, b, c), expected3), "logsumexp(log2,log3,log5) = log10");
    
    return failed;
}

// Test EM convergence on synthetic error-only data
int test_em_error_only() {
    int failed = 0;
    
    // Create histogram with only error-like sites (low mismatch fraction)
    SnpHistogram hist;
    // Sites with n=50, k=0 (no mismatches)
    hist[(50 << 16) | 0] = 10000;
    // Sites with n=50, k=1 (very low mismatch)
    hist[(50 << 16) | 1] = 500;
    // Sites with n=50, k=2
    hist[(50 << 16) | 2] = 50;
    
    SlamSnpEM em(100, 1e-8);
    SlamSnpEMResult result = em.fit(hist);
    
    failed += check(result.converged, "error_only: converged");
    failed += check(result.p_ERR > 0 && result.p_ERR < 0.05, "error_only: p_ERR reasonable");
    failed += check(result.pi_ERR > 0.95, "error_only: pi_ERR dominant");
    
    // Posterior for error-like site should be low
    double post_error = em.posterior(50, 1);
    failed += check(post_error < 0.1, "error_only: posterior(50,1) < 0.1");
    
    return failed;
}

// Test EM convergence on mixed data (error + SNPs)
int test_em_mixed() {
    int failed = 0;
    
    // Create histogram with mixed sites
    SnpHistogram hist;
    
    // Error-like sites (low mismatch)
    hist[(100 << 16) | 0] = 8000;
    hist[(100 << 16) | 1] = 1000;
    hist[(100 << 16) | 2] = 200;
    
    // Heterozygous SNP-like sites (~50% mismatch)
    hist[(100 << 16) | 45] = 30;
    hist[(100 << 16) | 50] = 50;
    hist[(100 << 16) | 55] = 30;
    
    // Homozygous SNP-like sites (~95% mismatch)
    hist[(100 << 16) | 93] = 10;
    hist[(100 << 16) | 95] = 20;
    hist[(100 << 16) | 97] = 10;
    
    SlamSnpEM em(100, 1e-8);
    SlamSnpEMResult result = em.fit(hist);
    
    failed += check(result.converged, "mixed: converged");
    failed += check(result.p_ERR > 0 && result.p_ERR < 0.05, "mixed: p_ERR reasonable");
    failed += check(result.pi_ERR > 0.8, "mixed: pi_ERR still dominant");
    failed += check(result.pi_HET > 0.001, "mixed: pi_HET > 0");
    failed += check(result.pi_HOM > 0.001, "mixed: pi_HOM > 0");
    
    // Check posteriors
    double post_het = em.posterior(100, 50);
    double post_hom = em.posterior(100, 95);
    double post_err = em.posterior(100, 1);
    
    failed += check(post_het > 0.9, "mixed: posterior(100,50) > 0.9 (het)");
    failed += check(post_hom > 0.9, "mixed: posterior(100,95) > 0.9 (hom)");
    failed += check(post_err < 0.1, "mixed: posterior(100,1) < 0.1 (error)");
    
    return failed;
}

// Test posterior boundary conditions
int test_posterior_boundaries() {
    int failed = 0;
    
    SlamSnpEM em(50, 1e-7);
    
    // Fit on some data first
    SnpHistogram hist;
    hist[(50 << 16) | 0] = 1000;
    hist[(50 << 16) | 1] = 100;
    em.fit(hist);
    
    // n=0 should return 0
    failed += check(em.posterior(0, 0) == 0.0, "posterior(0,0) = 0");
    
    // Posterior should be in [0,1]
    double p1 = em.posterior(100, 50);
    double p2 = em.posterior(100, 1);
    double p3 = em.posterior(100, 95);
    
    failed += check(p1 >= 0.0 && p1 <= 1.0, "posterior(100,50) in [0,1]");
    failed += check(p2 >= 0.0 && p2 <= 1.0, "posterior(100,1) in [0,1]");
    failed += check(p3 >= 0.0 && p3 <= 1.0, "posterior(100,95) in [0,1]");
    
    return failed;
}

// Test empty histogram handling
int test_empty_histogram() {
    int failed = 0;
    
    SlamSnpEM em(50, 1e-7);
    SnpHistogram empty;
    
    SlamSnpEMResult result = em.fit(empty);
    
    failed += check(!result.converged, "empty: not converged");
    failed += check(result.iterations == 0, "empty: 0 iterations");
    
    return failed;
}

int main() {
    int failed = 0;
    
    std::cerr << "=== log_binom_pmf edge cases ===\n";
    failed += test_log_binom_pmf_edges();
    
    std::cerr << "\n=== log_binom_pmf values ===\n";
    failed += test_log_binom_pmf_values();
    
    std::cerr << "\n=== logsumexp ===\n";
    failed += test_logsumexp();
    
    std::cerr << "\n=== EM error-only data ===\n";
    failed += test_em_error_only();
    
    std::cerr << "\n=== EM mixed data ===\n";
    failed += test_em_mixed();
    
    std::cerr << "\n=== posterior boundaries ===\n";
    failed += test_posterior_boundaries();
    
    std::cerr << "\n=== empty histogram ===\n";
    failed += test_empty_histogram();
    
    std::cerr << "\n";
    if (failed == 0) {
        std::cout << "ALL TESTS PASSED\n";
        return 0;
    } else {
        std::cout << failed << " TEST(S) FAILED\n";
        return 1;
    }
}
