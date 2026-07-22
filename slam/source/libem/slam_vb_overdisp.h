#ifndef SLAM_VB_OVERDISP_H
#define SLAM_VB_OVERDISP_H

#include "SlamSolver.h"

#include <cstdint>
#include <vector>

struct VbOverdispResult {
    double ntr_map = 0.0;         // MAP estimate for fraction new
    double ntr_mean = 0.0;        // posterior mean for fraction new
    double log_likelihood = 0.0;
    double p_err_used = 0.0;
    double p_conv_used = 0.0;
    double dispersion_phi = 0.0;
    int iters = 0;
    bool converged = false;
};

class SlamVbOverdispSolver {
public:
    SlamVbOverdispSolver(double error_rate, double conversion_rate, double phi,
                         double prior_alpha, double prior_beta)
        : p_error_rate_(error_rate),
          p_conversion_rate_(conversion_rate),
          dispersion_phi_(phi),
          prior_alpha_(prior_alpha),
          prior_beta_(prior_beta) {}

    VbOverdispResult solve(const MismatchHistogram& gene_data) const;

private:
    double p_error_rate_;
    double p_conversion_rate_;
    double dispersion_phi_;
    double prior_alpha_;
    double prior_beta_;

    double log_beta_binom_pmf(uint16_t n, uint16_t tc, double p, double phi) const;
    double calc_log_likelihood(const MismatchHistogram& data, double pi) const;
};

#endif
