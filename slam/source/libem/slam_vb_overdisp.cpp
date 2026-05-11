#include "slam_vb_overdisp.h"

#include <cmath>
#include <limits>

static inline double clamp_prob(double p) {
    const double eps = 1e-12;
    if (p < eps) return eps;
    if (p > 1.0 - eps) return 1.0 - eps;
    return p;
}

double SlamVbOverdispSolver::log_beta_binom_pmf(uint16_t n, uint16_t tc, double p, double phi) const {
    if (tc > n || phi <= 0.0) {
        return -std::numeric_limits<double>::infinity();
    }
    p = clamp_prob(p);
    double nn = static_cast<double>(n);
    double kk = static_cast<double>(tc);
    double alpha = p * phi;
    double beta = (1.0 - p) * phi;

    // log [ C(n,k) * Beta(k+alpha, n-k+beta) / Beta(alpha,beta) ]
    double log_coeff = std::lgamma(nn + 1.0) - std::lgamma(kk + 1.0) - std::lgamma(nn - kk + 1.0);
    double log_beta_num = std::lgamma(kk + alpha) + std::lgamma(nn - kk + beta) - std::lgamma(nn + alpha + beta);
    double log_beta_den = std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);
    return log_coeff + log_beta_num - log_beta_den;
}

double SlamVbOverdispSolver::calc_log_likelihood(const MismatchHistogram& data, double pi) const {
    double ll = 0.0;
    for (const auto& entry : data) {
        uint16_t n = slamKeyNT(entry.first);
        uint16_t tc = slamKeyTC(entry.first);
        double count = entry.second;

        double log_old = std::log(1.0 - pi) + log_beta_binom_pmf(n, tc, p_error_rate_, dispersion_phi_);
        double log_new = std::log(pi) + log_beta_binom_pmf(n, tc, p_conversion_rate_, dispersion_phi_);

        double max_log = (log_old > log_new) ? log_old : log_new;
        double sum = std::exp(log_old - max_log) + std::exp(log_new - max_log);
        ll += count * (max_log + std::log(sum));
    }
    return ll;
}

VbOverdispResult SlamVbOverdispSolver::solve(const MismatchHistogram& gene_data) const {
    VbOverdispResult result;
    result.p_err_used = p_error_rate_;
    result.p_conv_used = p_conversion_rate_;
    result.dispersion_phi = dispersion_phi_;

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

    const double tol = 1e-6;
    const int max_iters = 2000;
    double pi = 0.05;
    double prev_ll = -std::numeric_limits<double>::infinity();

    // Ensure prior is well-posed for MAP update
    double prior_a = (prior_alpha_ > 0.0) ? prior_alpha_ : 1.0;
    double prior_b = (prior_beta_ > 0.0) ? prior_beta_ : 1.0;

    for (int iter = 0; iter < max_iters; ++iter) {
        double sum_gamma = 0.0;
        double ll = 0.0;

        for (const auto& entry : gene_data) {
            uint16_t n = slamKeyNT(entry.first);
            uint16_t tc = slamKeyTC(entry.first);
            double count = entry.second;

            double log_old = std::log(1.0 - pi) + log_beta_binom_pmf(n, tc, p_error_rate_, dispersion_phi_);
            double log_new = std::log(pi) + log_beta_binom_pmf(n, tc, p_conversion_rate_, dispersion_phi_);
            double max_log = (log_old > log_new) ? log_old : log_new;
            double sum = std::exp(log_old - max_log) + std::exp(log_new - max_log);
            double gamma = std::exp(log_new - max_log) / sum;
            sum_gamma += count * gamma;
            ll += count * (max_log + std::log(sum));
        }

        // MAP update for pi with Beta(prior_a, prior_b)
        double denom = prior_a + prior_b - 2.0 + total;
        double numer = prior_a - 1.0 + sum_gamma;
        double pi_new;
        if (denom > 0.0) {
            pi_new = numer / denom;
        } else {
            pi_new = (prior_a + sum_gamma) / (prior_a + prior_b + total);
        }
        pi_new = clamp_prob(pi_new);

        if (std::abs(pi_new - pi) < tol && std::abs(ll - prev_ll) < tol) {
            pi = pi_new;
            result.converged = true;
            result.log_likelihood = ll;
            result.iters = iter + 1;
            break;
        }
        pi = pi_new;
        prev_ll = ll;
        result.log_likelihood = ll;
        result.iters = iter + 1;
    }

    // Posterior mean for pi under Beta prior
    double sum_gamma = 0.0;
    for (const auto& entry : gene_data) {
        uint16_t n = slamKeyNT(entry.first);
        uint16_t tc = slamKeyTC(entry.first);
        double count = entry.second;
        double log_old = std::log(1.0 - pi) + log_beta_binom_pmf(n, tc, p_error_rate_, dispersion_phi_);
        double log_new = std::log(pi) + log_beta_binom_pmf(n, tc, p_conversion_rate_, dispersion_phi_);
        double max_log = (log_old > log_new) ? log_old : log_new;
        double sum = std::exp(log_old - max_log) + std::exp(log_new - max_log);
        double gamma = std::exp(log_new - max_log) / sum;
        sum_gamma += count * gamma;
    }
    result.ntr_map = pi;
    result.ntr_mean = (prior_a + sum_gamma) / (prior_a + prior_b + total);

    return result;
}
