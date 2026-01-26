/**
 * @file gmm.c
 * @brief Implementation of 1D 2-component GMM with tied covariance
 * 
 * This implements CR9-compatible GMM feature calling using the EM algorithm.
 */

#include "../include/gmm.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#define PI 3.14159265358979323846
#define LOG_2PI 1.8378770664093453  /* log(2*pi) */
#define MIN_VARIANCE 1e-6

/* Simple linear congruential generator for reproducibility */
typedef struct {
    unsigned int state;
} rng_state;

static void rng_init(rng_state *rng, unsigned int seed) {
    rng->state = seed;
}

static double rng_uniform(rng_state *rng) {
    /* LCG parameters from Numerical Recipes */
    rng->state = 1664525 * rng->state + 1013904223;
    return (double)(rng->state & 0x7FFFFFFF) / (double)0x7FFFFFFF;
}

/* Compute log of Gaussian PDF */
static double log_gaussian(double x, double mean, double variance) {
    double diff = x - mean;
    return -0.5 * (LOG_2PI + log(variance) + diff * diff / variance);
}

/* Log-sum-exp for numerical stability */
static double log_sum_exp(double a, double b) {
    if (a > b) {
        return a + log(1.0 + exp(b - a));
    } else {
        return b + log(1.0 + exp(a - b));
    }
}

/* Single EM run */
static gmm_result gmm_em_single(const double *data, int n, double init_mean0, 
                                 double init_mean1, double init_var, 
                                 int max_iter, double tol) {
    gmm_result result;
    memset(&result, 0, sizeof(result));
    
    if (n < 2) {
        result.converged = 0;
        return result;
    }
    
    /* Initialize parameters */
    double mean0 = init_mean0;
    double mean1 = init_mean1;
    double var = init_var > MIN_VARIANCE ? init_var : MIN_VARIANCE;
    double weight0 = 0.5;
    double weight1 = 0.5;
    
    /* Allocate responsibilities */
    double *resp = malloc(n * sizeof(double));  /* responsibility for component 1 */
    if (!resp) {
        result.converged = 0;
        return result;
    }
    
    double prev_ll = -DBL_MAX;
    int iter;
    
    for (iter = 0; iter < max_iter; iter++) {
        /* E-step: compute responsibilities */
        double ll = 0.0;
        for (int i = 0; i < n; i++) {
            double log_p0 = log(weight0) + log_gaussian(data[i], mean0, var);
            double log_p1 = log(weight1) + log_gaussian(data[i], mean1, var);
            double log_sum = log_sum_exp(log_p0, log_p1);
            resp[i] = exp(log_p1 - log_sum);
            ll += log_sum;
        }
        
        /* Check convergence */
        if (iter > 0 && fabs(ll - prev_ll) < tol) {
            result.converged = 1;
            break;
        }
        prev_ll = ll;
        
        /* M-step: update parameters */
        double sum_resp0 = 0.0, sum_resp1 = 0.0;
        double weighted_sum0 = 0.0, weighted_sum1 = 0.0;
        
        for (int i = 0; i < n; i++) {
            double r0 = 1.0 - resp[i];
            double r1 = resp[i];
            sum_resp0 += r0;
            sum_resp1 += r1;
            weighted_sum0 += r0 * data[i];
            weighted_sum1 += r1 * data[i];
        }
        
        /* Update means */
        if (sum_resp0 > 1e-10) {
            mean0 = weighted_sum0 / sum_resp0;
        }
        if (sum_resp1 > 1e-10) {
            mean1 = weighted_sum1 / sum_resp1;
        }
        
        /* Update tied variance */
        double var_sum = 0.0;
        for (int i = 0; i < n; i++) {
            double r0 = 1.0 - resp[i];
            double r1 = resp[i];
            double diff0 = data[i] - mean0;
            double diff1 = data[i] - mean1;
            var_sum += r0 * diff0 * diff0 + r1 * diff1 * diff1;
        }
        var = var_sum / n;
        if (var < MIN_VARIANCE) var = MIN_VARIANCE;
        
        /* Update weights */
        weight0 = sum_resp0 / n;
        weight1 = sum_resp1 / n;
        
        /* Ensure weights are valid */
        if (weight0 < 1e-10) weight0 = 1e-10;
        if (weight1 < 1e-10) weight1 = 1e-10;
        double weight_sum = weight0 + weight1;
        weight0 /= weight_sum;
        weight1 /= weight_sum;
    }
    
    result.n_iter = iter;
    result.log_likelihood = prev_ll;
    
    /* Order components so mean_low < mean_high */
    if (mean0 <= mean1) {
        result.mean_low = mean0;
        result.mean_high = mean1;
        result.weight_low = weight0;
        result.weight_high = weight1;
        result.positive_component = 1;
    } else {
        result.mean_low = mean1;
        result.mean_high = mean0;
        result.weight_low = weight1;
        result.weight_high = weight0;
        result.positive_component = 0;
    }
    result.variance = var;
    
    free(resp);
    return result;
}

/* Compute data statistics for initialization */
static void compute_stats(const double *data, int n, double *min_val, double *max_val,
                          double *mean, double *var) {
    *min_val = data[0];
    *max_val = data[0];
    double sum = 0.0;
    
    for (int i = 0; i < n; i++) {
        if (data[i] < *min_val) *min_val = data[i];
        if (data[i] > *max_val) *max_val = data[i];
        sum += data[i];
    }
    *mean = sum / n;
    
    double var_sum = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = data[i] - *mean;
        var_sum += diff * diff;
    }
    *var = var_sum / n;
    if (*var < MIN_VARIANCE) *var = MIN_VARIANCE;
}

gmm_result gmm_fit_1d(const double *data, int n, int n_init, int max_iter,
                      double tol, unsigned int random_state) {
    gmm_result best_result;
    memset(&best_result, 0, sizeof(best_result));
    best_result.log_likelihood = -DBL_MAX;
    
    if (n < 2) {
        return best_result;
    }
    
    /* Set defaults */
    if (n_init <= 0) n_init = 10;
    if (max_iter <= 0) max_iter = 100;
    if (tol <= 0) tol = 1e-3;
    
    /* Compute data statistics */
    double min_val, max_val, data_mean, data_var;
    compute_stats(data, n, &min_val, &max_val, &data_mean, &data_var);
    
    /* Initialize RNG */
    rng_state rng;
    rng_init(&rng, random_state);
    
    /* Run multiple initializations */
    for (int init = 0; init < n_init; init++) {
        /* Generate random initial means */
        double range = max_val - min_val;
        double init_mean0, init_mean1;
        
        if (init == 0) {
            /* First init: use quantile-like initialization */
            init_mean0 = min_val + 0.25 * range;
            init_mean1 = min_val + 0.75 * range;
        } else {
            /* Random initialization */
            init_mean0 = min_val + rng_uniform(&rng) * range;
            init_mean1 = min_val + rng_uniform(&rng) * range;
        }
        
        /* Run EM */
        gmm_result result = gmm_em_single(data, n, init_mean0, init_mean1,
                                          data_var, max_iter, tol);
        
        /* Keep best result */
        if (result.log_likelihood > best_result.log_likelihood) {
            best_result = result;
        }
    }
    
    return best_result;
}

gmm_assignments* gmm_predict(const double *data, int n, const gmm_result *result) {
    if (!data || n <= 0 || !result) return NULL;
    
    gmm_assignments *assignments = malloc(sizeof(gmm_assignments));
    if (!assignments) return NULL;
    
    assignments->labels = malloc(n * sizeof(int));
    assignments->responsibilities = malloc(n * sizeof(double));
    assignments->n = n;
    
    if (!assignments->labels || !assignments->responsibilities) {
        gmm_free_assignments(assignments);
        return NULL;
    }
    
    double mean0 = result->mean_low;
    double mean1 = result->mean_high;
    double var = result->variance;
    double weight0 = result->weight_low;
    double weight1 = result->weight_high;
    
    for (int i = 0; i < n; i++) {
        double log_p0 = log(weight0) + log_gaussian(data[i], mean0, var);
        double log_p1 = log(weight1) + log_gaussian(data[i], mean1, var);
        double log_sum = log_sum_exp(log_p0, log_p1);
        double resp1 = exp(log_p1 - log_sum);
        
        assignments->responsibilities[i] = resp1;
        assignments->labels[i] = (resp1 >= 0.5) ? 1 : 0;
    }
    
    return assignments;
}

void gmm_free_assignments(gmm_assignments *assignments) {
    if (assignments) {
        free(assignments->labels);
        free(assignments->responsibilities);
        free(assignments);
    }
}

int gmm_call_feature(const int *counts, int n_cells, int min_umi_threshold,
                     int *positive_calls, int *umi_threshold) {
    if (!counts || n_cells <= 0 || !positive_calls) return -1;
    
    /* Initialize outputs */
    memset(positive_calls, 0, n_cells * sizeof(int));
    if (umi_threshold) *umi_threshold = 0;
    
    /* Check if all counts are zero */
    int max_count = 0;
    for (int i = 0; i < n_cells; i++) {
        if (counts[i] > max_count) max_count = counts[i];
    }
    
    if (max_count == 0 || n_cells < 2) {
        /* No positive calls possible */
        return 0;
    }
    
    /* Transform counts: x = log10(1 + counts) */
    double *transformed = malloc(n_cells * sizeof(double));
    if (!transformed) return -1;
    
    for (int i = 0; i < n_cells; i++) {
        transformed[i] = log10(1.0 + counts[i]);
    }
    
    /* Fit GMM */
    gmm_result gmm = gmm_fit_1d(transformed, n_cells, 10, 100, 1e-3, 0);
    
    /* Get assignments */
    gmm_assignments *assignments = gmm_predict(transformed, n_cells, &gmm);
    free(transformed);
    
    if (!assignments) {
        return -1;
    }
    
    /* Make calls: positive if assigned to high component AND count >= threshold */
    int min_positive_umi = INT_MAX;
    for (int i = 0; i < n_cells; i++) {
        if (assignments->labels[i] == 1 && counts[i] >= min_umi_threshold) {
            positive_calls[i] = 1;
            if (counts[i] < min_positive_umi) {
                min_positive_umi = counts[i];
            }
        }
    }
    
    if (umi_threshold) {
        *umi_threshold = (min_positive_umi == INT_MAX) ? 0 : min_positive_umi;
    }
    
    gmm_free_assignments(assignments);
    return 0;
}
