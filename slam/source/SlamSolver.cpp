#include "SlamSolver.h"

#include <cmath>
#include <limits>

double SlamSolver::log_binom_pmf(uint16_t n, uint8_t k, double p) const {
    if (p <= 0.0 || p >= 1.0) {
        return -std::numeric_limits<double>::infinity();
    }
    double nn = static_cast<double>(n);
    double kk = static_cast<double>(k);
    double log_coeff = std::lgamma(nn + 1.0) - std::lgamma(kk + 1.0) - std::lgamma(nn - kk + 1.0);
    return log_coeff + kk * std::log(p) + (nn - kk) * std::log(1.0 - p);
}

double SlamSolver::calc_log_likelihood(const MismatchHistogram& data, double pi) const {
    double ll = 0.0;
    for (const auto& entry : data) {
        uint16_t n = entry.first >> 8;
        uint8_t k = static_cast<uint8_t>(entry.first & 0xFFu);
        double count = entry.second;

        double log_old = std::log(1.0 - pi) + log_binom_pmf(n, k, p_error_rate_);
        double log_new = std::log(pi) + log_binom_pmf(n, k, p_conversion_rate_);

        double max_log = (log_old > log_new) ? log_old : log_new;
        double sum = std::exp(log_old - max_log) + std::exp(log_new - max_log);
        ll += count * (max_log + std::log(sum));
    }
    return ll;
}

SlamResult SlamSolver::solve(const MismatchHistogram& gene_data) const {
    SlamResult result;
    if (gene_data.empty()) {
        return result;
    }

    double total = 0.0;
    for (const auto& entry : gene_data) {
        total += entry.second;
    }
    if (total <= 0.0) {
        return result;
    }

    double pi = 0.1;
    double prev_ll = -std::numeric_limits<double>::infinity();
    const int max_iters = 1000;
    const double tol = 1e-6;

    for (int iter = 0; iter < max_iters; ++iter) {
        double num = 0.0;
        double ll = 0.0;
        for (const auto& entry : gene_data) {
            uint16_t n = entry.first >> 8;
            uint8_t k = static_cast<uint8_t>(entry.first & 0xFFu);
            double count = entry.second;

            double log_old = std::log(1.0 - pi) + log_binom_pmf(n, k, p_error_rate_);
            double log_new = std::log(pi) + log_binom_pmf(n, k, p_conversion_rate_);
            double max_log = (log_old > log_new) ? log_old : log_new;
            double sum = std::exp(log_old - max_log) + std::exp(log_new - max_log);
            double gamma = std::exp(log_new - max_log) / sum;
            num += count * gamma;
            ll += count * (max_log + std::log(sum));
        }

        double pi_new = num / total;
        if (std::abs(pi_new - pi) < tol && std::abs(ll - prev_ll) < tol) {
            pi = pi_new;
            result.converged = true;
            result.log_likelihood = ll;
            break;
        }
        pi = pi_new;
        prev_ll = ll;
        result.log_likelihood = ll;
    }

    result.ntr = pi;
    result.sigma = std::sqrt(pi * (1.0 - pi) / total);
    return result;
}
