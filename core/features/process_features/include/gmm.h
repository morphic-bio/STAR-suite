/**
 * @file gmm.h
 * @brief 1D 2-component Gaussian Mixture Model with tied covariance
 * 
 * This module implements a GMM compatible with CR9's feature calling:
 * - 2 components with tied (shared) variance
 * - EM algorithm with deterministic initialization
 * - n_init=10 restarts to match sklearn's default
 * - random_state=0 for reproducibility
 */

#ifndef GMM_H
#define GMM_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* GMM fit result */
typedef struct {
    double mean_low;         /* Mean of lower component */
    double mean_high;        /* Mean of higher component */
    double variance;         /* Shared (tied) variance */
    double weight_low;       /* Mixing weight for lower component */
    double weight_high;      /* Mixing weight for higher component */
    int positive_component;  /* 0 or 1 - component with higher mean */
    double log_likelihood;   /* Final log-likelihood */
    int converged;           /* 1 if converged, 0 otherwise */
    int n_iter;              /* Number of EM iterations */
} gmm_result;

/* Component assignment for each data point */
typedef struct {
    int *labels;             /* Component assignment (0 or 1) for each point */
    double *responsibilities; /* Posterior probability of positive component */
    int n;                   /* Number of data points */
} gmm_assignments;

/**
 * Fit a 2-component GMM with tied covariance to 1D data.
 * 
 * @param data Array of data points
 * @param n Number of data points
 * @param n_init Number of random initializations (default: 10)
 * @param max_iter Maximum EM iterations per init (default: 100)
 * @param tol Convergence tolerance for log-likelihood (default: 1e-3)
 * @param random_state Seed for deterministic initialization (default: 0)
 * @return GMM fit result
 */
gmm_result gmm_fit_1d(const double *data, int n, int n_init, int max_iter, 
                      double tol, unsigned int random_state);

/**
 * Predict component assignments for data points.
 * 
 * @param data Array of data points
 * @param n Number of data points  
 * @param result GMM fit result from gmm_fit_1d
 * @return Component assignments (caller must free with gmm_free_assignments)
 */
gmm_assignments* gmm_predict(const double *data, int n, const gmm_result *result);

/**
 * Free GMM assignments structure.
 */
void gmm_free_assignments(gmm_assignments *assignments);

/**
 * CR9-style feature calling using GMM.
 * 
 * For a single feature:
 * 1. Transform counts: x = log10(1 + counts)
 * 2. Fit 2-component GMM with tied covariance
 * 3. Positive component = component with higher mean
 * 4. Cell is positive if:
 *    - counts >= min_umi_threshold
 *    - assigned to positive component
 *
 * @param counts UMI counts for this feature across all cells
 * @param n_cells Number of cells
 * @param min_umi_threshold Minimum UMI count for a call (default: 3)
 * @param positive_calls Output: 1 if positive, 0 if negative (must be pre-allocated)
 * @param umi_threshold Output: detected UMI threshold (min UMI among positive cells)
 * @return 0 on success, -1 on error
 */
int gmm_call_feature(const int *counts, int n_cells, int min_umi_threshold,
                     int *positive_calls, int *umi_threshold);

#ifdef __cplusplus
}
#endif

#endif /* GMM_H */
