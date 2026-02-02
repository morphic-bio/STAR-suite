/**
 * Basic unit checks for SlamVbOverdispSolver behavior.
 *
 * Compile:
 *   g++ -std=c++11 -I../../slam/source/libem -o test_slam_vb_overdisp \
 *     tests/slam/test_slam_vb_overdisp.cpp ../../slam/source/libem/slam_vb_overdisp.cpp -lm
 *
 * Run:
 *   ./test_slam_vb_overdisp
 */

#include <cmath>
#include <iostream>
#include <map>

#include "slam_vb_overdisp.h"

static double log_binom_pmf(uint16_t n, uint8_t k, double p) {
    if (p <= 0.0 || p >= 1.0) return -1e30;
    double nn = static_cast<double>(n);
    double kk = static_cast<double>(k);
    double log_coeff = std::lgamma(nn + 1.0) - std::lgamma(kk + 1.0) - std::lgamma(nn - kk + 1.0);
    return log_coeff + kk * std::log(p) + (nn - kk) * std::log(1.0 - p);
}

static double log_beta_binom_pmf(uint16_t n, uint8_t k, double p, double phi) {
    if (phi <= 0.0 || p <= 0.0 || p >= 1.0) return -1e30;
    double nn = static_cast<double>(n);
    double kk = static_cast<double>(k);
    double alpha = p * phi;
    double beta = (1.0 - p) * phi;
    double log_coeff = std::lgamma(nn + 1.0) - std::lgamma(kk + 1.0) - std::lgamma(nn - kk + 1.0);
    double log_beta_num = std::lgamma(kk + alpha) + std::lgamma(nn - kk + beta) - std::lgamma(nn + alpha + beta);
    double log_beta_den = std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);
    return log_coeff + log_beta_num - log_beta_den;
}

static std::map<uint16_t, double> make_hist(uint16_t n, double pi, double p_err, double p_conv, double total_reads) {
    std::map<uint16_t, double> hist;
    for (uint8_t k = 0; k <= n; ++k) {
        double p_old = std::exp(log_binom_pmf(n, k, p_err));
        double p_new = std::exp(log_binom_pmf(n, k, p_conv));
        double p_mix = (1.0 - pi) * p_old + pi * p_new;
        double count = total_reads * p_mix;
        if (count > 0.0) {
            uint16_t key = static_cast<uint16_t>((n << 8) | k);
            hist[key] = count;
        }
    }
    return hist;
}

static std::map<uint16_t, double> make_hist_bb(uint16_t n, double pi, double p_err, double p_conv, double phi, double total_reads) {
    std::map<uint16_t, double> hist;
    for (uint8_t k = 0; k <= n; ++k) {
        double p_old = std::exp(log_beta_binom_pmf(n, k, p_err, phi));
        double p_new = std::exp(log_beta_binom_pmf(n, k, p_conv, phi));
        double p_mix = (1.0 - pi) * p_old + pi * p_new;
        double count = total_reads * p_mix;
        if (count > 0.0) {
            uint16_t key = static_cast<uint16_t>((n << 8) | k);
            hist[key] = count;
        }
    }
    return hist;
}

static bool approx(double a, double b, double tol) {
    return std::fabs(a - b) <= tol;
}

int main() {
    const uint16_t n = 20;
    const double total_reads = 10000.0;
    const double p_err = 0.001;
    const double p_conv = 0.05;

    // Test 1: recover pi ~ 0.3 under low dispersion
    {
        double pi_true = 0.30;
        std::map<uint16_t, double> hist = make_hist(n, pi_true, p_err, p_conv, total_reads);
        SlamVbOverdispSolver solver(p_err, p_conv, 1000.0, 1.0, 1.0);
        VbOverdispResult res = solver.solve(hist);
        if (!approx(res.ntr_map, pi_true, 0.05)) {
            std::cerr << "FAIL: pi recovery expected ~0.30, got " << res.ntr_map << "\n";
            return 1;
        }
    }

    // Test 2: monotonicity (pi=0.6 > pi=0.3)
    {
        std::map<uint16_t, double> hist_low = make_hist(n, 0.30, p_err, p_conv, total_reads);
        std::map<uint16_t, double> hist_high = make_hist(n, 0.60, p_err, p_conv, total_reads);
        SlamVbOverdispSolver solver(p_err, p_conv, 1000.0, 1.0, 1.0);
        VbOverdispResult res_low = solver.solve(hist_low);
        VbOverdispResult res_high = solver.solve(hist_high);
        if (!(res_high.ntr_map > res_low.ntr_map)) {
            std::cerr << "FAIL: monotonicity failed (high <= low)\n";
            return 1;
        }
    }

    // Test 3: dispersion effect (match phi improves recovery)
    {
        double pi_true = 0.30;
        double phi_true = 5.0;
        std::map<uint16_t, double> hist = make_hist_bb(n, pi_true, p_err, p_conv, phi_true, total_reads);
        SlamVbOverdispSolver solver_match(p_err, p_conv, phi_true, 1.0, 1.0);
        SlamVbOverdispSolver solver_mismatch(p_err, p_conv, 1000.0, 1.0, 1.0);
        VbOverdispResult res_match = solver_match.solve(hist);
        VbOverdispResult res_mismatch = solver_mismatch.solve(hist);
        double err_match = std::fabs(res_match.ntr_map - pi_true);
        double err_mismatch = std::fabs(res_mismatch.ntr_map - pi_true);
        if (!(err_match <= err_mismatch + 1e-6)) {
            std::cerr << "FAIL: dispersion match did not improve recovery\n";
            return 1;
        }
    }

    // Test 4: stronger prior lowers pi
    {
        std::map<uint16_t, double> hist = make_hist(n, 0.30, p_err, p_conv, total_reads);
        SlamVbOverdispSolver solver_weak(p_err, p_conv, 1000.0, 1.0, 1.0);
        SlamVbOverdispSolver solver_strong(p_err, p_conv, 1000.0, 1.0, 99.0);
        VbOverdispResult res_weak = solver_weak.solve(hist);
        VbOverdispResult res_strong = solver_strong.solve(hist);
        if (!(res_strong.ntr_map < res_weak.ntr_map)) {
            std::cerr << "FAIL: strong prior did not lower pi\n";
            return 1;
        }
    }

    std::cout << "PASS: SlamVbOverdispSolver basic checks\n";
    return 0;
}
