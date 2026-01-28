/**
 * @file nbem.h
 * @brief NB-EM feature calling with cell-specific background (SCEPTRE-style)
 *
 * This module implements negative binomial mixture model feature calling:
 * - Cell-specific background mean via Poisson IRLS on covariates
 * - Per-feature EM to fit mixture proportion and effect size
 * - Posterior-based calling with fallback threshold
 *
 * Model for feature f, cell i with count g_i:
 *   Background: mu0_i = exp(beta · x_i)  where x_i are covariates
 *   Perturbed:  mu1_i = mu0_i * exp(delta_f)
 *   Mixture:    g_i ~ (1-pi_f)*NB(mu0_i, phi) + pi_f*NB(mu1_i, phi)
 *   Call:       positive if P(perturbed|g_i) >= prob_threshold
 */

#ifndef NBEM_H
#define NBEM_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Small Matrix Linear Algebra (Householder QR)
 * ============================================================================ */

/**
 * Solve a linear system Ax = b using Householder QR decomposition.
 * Suitable for small dense systems (n <= 10).
 *
 * @param A     Matrix in column-major order (n x n), MODIFIED in place
 * @param b     Right-hand side vector (n), MODIFIED to contain solution
 * @param n     System dimension
 * @return      0 on success, -1 if singular or numerical failure
 */
int qr_solve(double *A, double *b, int n);

/**
 * Solve least squares min ||Ax - b||_2 using Householder QR.
 * For overdetermined systems (m >= n).
 *
 * @param A     Matrix in column-major order (m x n), MODIFIED in place
 * @param b     Right-hand side vector (m), MODIFIED to contain solution in first n elements
 * @param m     Number of rows
 * @param n     Number of columns (n <= m)
 * @return      0 on success, -1 if rank deficient or numerical failure
 */
int qr_lstsq(double *A, double *b, int m, int n);

/* ============================================================================
 * Poisson IRLS for Background Model
 * ============================================================================ */

/**
 * Fit Poisson GLM via IRLS: E[y] = exp(X * beta)
 *
 * @param X         Design matrix (n_cells x n_covars), column-major
 * @param y         Response vector (counts), length n_cells
 * @param n_cells   Number of observations
 * @param n_covars  Number of covariates (including intercept if present)
 * @param beta      Output: fitted coefficients, length n_covars
 * @param max_iter  Maximum IRLS iterations (default: 25)
 * @param tol       Convergence tolerance on deviance (default: 1e-6)
 * @return          0 on success, -1 on failure (non-convergence, singular)
 */
int poisson_irls(const double *X, const int *y, int n_cells, int n_covars,
                 double *beta, int max_iter, double tol);

/**
 * Compute predicted means from Poisson GLM: mu_i = exp(X_i * beta)
 *
 * @param X         Design matrix (n_cells x n_covars), column-major
 * @param beta      Coefficients, length n_covars
 * @param n_cells   Number of cells
 * @param n_covars  Number of covariates
 * @param mu        Output: predicted means, length n_cells
 */
void poisson_predict(const double *X, const double *beta, int n_cells, int n_covars,
                     double *mu);

/* ============================================================================
 * Negative Binomial Distribution
 * ============================================================================ */

/**
 * Log probability mass function of NB(mu, phi).
 * Parameterization: mean = mu, variance = mu + mu^2/phi
 * (phi is the "size" or "dispersion" parameter; larger phi = less overdispersion)
 *
 * @param k     Count (non-negative integer)
 * @param mu    Mean (positive)
 * @param phi   Dispersion parameter (positive)
 * @return      log P(X = k)
 */
double nb_logpmf(int k, double mu, double phi);

/**
 * Estimate global dispersion (phi) using method of moments.
 * phi = mean^2 / (var - mean), clamped to [phi_min, phi_max]
 *
 * @param counts    Count array
 * @param n         Length of array
 * @param phi_min   Minimum phi (default: 0.01)
 * @param phi_max   Maximum phi (default: 1000.0)
 * @return          Estimated phi
 */
double nb_estimate_phi_mom(const int *counts, int n, double phi_min, double phi_max);

/**
 * Streaming phi accumulator for memory-efficient global phi estimation.
 * Uses Welford's online algorithm for mean and variance.
 */
typedef struct {
    int64_t n;          /* Number of observations */
    double mean;        /* Running mean */
    double M2;          /* Sum of squared differences from mean */
} nbem_phi_accum;

/**
 * Initialize streaming phi accumulator.
 */
void nbem_phi_accum_init(nbem_phi_accum *acc);

/**
 * Add a batch of counts to the streaming phi accumulator.
 *
 * @param acc       Accumulator state
 * @param counts    Count array for one feature
 * @param n         Length of array
 */
void nbem_phi_accum_add(nbem_phi_accum *acc, const int *counts, int n);

/**
 * Finalize streaming phi estimate.
 *
 * @param acc       Accumulator state
 * @param phi_min   Minimum phi (default: 0.01)
 * @param phi_max   Maximum phi (default: 1000.0)
 * @return          Estimated phi
 */
double nbem_phi_accum_finish(const nbem_phi_accum *acc, double phi_min, double phi_max);

/* ============================================================================
 * NB-EM Per-Feature Fitting
 * ============================================================================ */

/** Result of NB-EM fitting for a single feature */
typedef struct {
    double pi;              /* Mixture proportion (fraction perturbed) */
    double delta;           /* Log fold-change: mu1 = mu0 * exp(delta) */
    double phi;             /* Dispersion used (may be global) */
    double log_likelihood;  /* Final log-likelihood */
    int converged;          /* 1 if converged, 0 otherwise */
    int n_iter;             /* Number of EM iterations */
    int used_fallback;      /* 1 if EM failed and fallback was used */
    int failure_reason;     /* 0=none, >0 indicates failure reason */
} nbem_fit_result;

/** NB-EM failure reasons (for debug instrumentation) */
typedef enum {
    NBEM_FAIL_NONE = 0,
    NBEM_FAIL_INVALID_INPUT = 1,
    NBEM_FAIL_ALLOC = 2,
    NBEM_FAIL_NAN_INF = 3,
    NBEM_FAIL_LL_DECREASE = 4,
    NBEM_FAIL_DEGENERATE_PI = 5,
    NBEM_FAIL_MAX_ITER = 6
} nbem_failure_reason;

/**
 * Fit NB mixture model for a single feature using EM.
 *
 * @param counts        UMI counts for this feature, length n_cells
 * @param mu0           Background means (from Poisson GLM), length n_cells
 * @param n_cells       Number of cells
 * @param phi           Dispersion parameter (fixed)
 * @param max_iter      Maximum EM iterations (default: 100)
 * @param tol           Convergence tolerance on log-likelihood (default: 1e-4)
 * @param pi_init       Initial mixing proportion (default: 0.1)
 * @param delta_init    Initial log fold-change (default: 1.0)
 * @return              Fit result
 */
nbem_fit_result nbem_fit_feature(const int *counts, const double *mu0, int n_cells,
                                  double phi, int use_poisson, int max_iter, double tol,
                                  double pi_init, double delta_init);

/**
 * Fit NB mixture model with per-iteration debug logging.
 *
 * @param counts        UMI counts for this feature, length n_cells
 * @param mu0           Background means (from Poisson GLM), length n_cells
 * @param n_cells       Number of cells
 * @param phi           Dispersion parameter (fixed)
 * @param max_iter      Maximum EM iterations (default: 100)
 * @param tol           Convergence tolerance on log-likelihood (default: 1e-4)
 * @param pi_init       Initial mixing proportion (default: 0.1)
 * @param delta_init    Initial log fold-change (default: 1.0)
 * @param feature_name  Feature name for logging (optional)
 * @param debug_fp      File handle to write per-iteration logs (optional)
 * @return              Fit result
 */
nbem_fit_result nbem_fit_feature_debug(const int *counts, const double *mu0, int n_cells,
                                       double phi, int max_iter, double tol,
                                       double pi_init, double delta_init,
                                       int use_poisson, const char *feature_name, FILE *debug_fp);

/**
 * Fit a Poisson mixture (with fixed mu0) to seed NB-EM parameters.
 *
 * @param counts     Feature counts per cell
 * @param mu0        Background mean per cell
 * @param n_cells    Number of cells
 * @param max_iter   Max iterations for Poisson EM (default: 20)
 * @param tol        Convergence tolerance (default: 1e-4)
 * @param pi_init    Initial pi
 * @param delta_init Initial delta (log fold change)
 * @param out_pi     Output pi seed
 * @param out_delta  Output delta seed
 * @return 0 on success, -1 on failure
 */

/**
 * Compute posterior probability of perturbation for each cell.
 *
 * @param counts        UMI counts, length n_cells
 * @param mu0           Background means, length n_cells
 * @param n_cells       Number of cells
 * @param pi            Mixing proportion
 * @param delta         Log fold-change
 * @param phi           Dispersion
 * @param posteriors    Output: P(perturbed|count), length n_cells
 */
void nbem_posteriors(const int *counts, const double *mu0, int n_cells,
                     double pi, double delta, double phi, int use_poisson, double *posteriors);

/**
 * Make calls based on posteriors or fallback threshold.
 *
 * @param posteriors        Posterior probabilities, length n_cells
 * @param counts            UMI counts, length n_cells
 * @param n_cells           Number of cells
 * @param prob_threshold    Posterior threshold for call (default: 0.8)
 * @param backup_threshold  UMI threshold for fallback (default: 5)
 * @param use_fallback      If true, use backup_threshold instead of posteriors
 * @param calls             Output: 1 if positive, 0 if negative, length n_cells
 * @return                  Number of positive calls
 */
int nbem_make_calls(const double *posteriors, const int *counts, int n_cells,
                    double prob_threshold, int backup_threshold, int use_fallback,
                    int *calls);

/* ============================================================================
 * Covariate Computation
 * ============================================================================ */

/** Per-cell covariates for background model */
typedef struct {
    double *grna_log_n_umis;      /* log(1 + total gRNA UMIs) */
    double *grna_log_n_nonzero;   /* log(1 + number of detected features) */
    double *gex_log_n_umis;       /* log(1 + GEX UMIs), NULL if not provided */
    double *gex_log_n_nonzero;    /* log(1 + GEX genes), NULL if not provided */
    int n_cells;
    int has_gex;                  /* 1 if GEX covariates present */
} nbem_covariates;

/**
 * Compute gRNA covariates from feature count matrix.
 *
 * @param feature_counts    Per-feature count arrays [num_features][n_cells]
 * @param num_features      Number of features
 * @param n_cells           Number of cells
 * @return                  Covariates structure (caller must free with nbem_free_covariates)
 */
nbem_covariates* nbem_compute_grna_covariates(int **feature_counts, int num_features, int n_cells);

/**
 * Load GEX covariates from TSV file and merge with existing covariates.
 * TSV format: barcode,response_n_umis,response_n_nonzero
 *
 * @param covariates        Existing covariates (modified in place)
 * @param tsv_path          Path to covariate TSV
 * @param barcodes          Barcode strings for matching, length n_cells
 * @return                  0 on success, -1 on failure
 */
int nbem_load_gex_covariates(nbem_covariates *covariates, const char *tsv_path,
                              char **barcodes);

/**
 * Build design matrix from covariates.
 * Columns: [intercept, grna_log_n_umis, grna_log_n_nonzero, (gex if present)]
 *
 * @param covariates        Covariate structure
 * @param X                 Output: design matrix (n_cells x n_covars), column-major
 * @return                  Number of covariates (columns)
 */
int nbem_build_design_matrix(const nbem_covariates *covariates, double *X);

/**
 * Free covariates structure.
 */
void nbem_free_covariates(nbem_covariates *covariates);

/* ============================================================================
 * MOI Estimation
 * ============================================================================ */

/** MOI classification */
typedef enum {
    NBEM_MOI_LOW = 0,
    NBEM_MOI_HIGH = 1
} nbem_moi_class;

/** MOI estimation result */
typedef struct {
    nbem_moi_class classification;
    double p0;              /* Fraction of cells with 0 guides */
    double lambda;          /* Estimated Poisson rate: -log(p0) */
    double p_multi_exp;     /* Expected multi-guide fraction under Poisson */
    double p_multi_obs;     /* Observed multi-guide fraction */
    int n_cells_0;          /* Cells with 0 guides */
    int n_cells_1;          /* Cells with 1 guide */
    int n_cells_multi;      /* Cells with >=2 guides */
} nbem_moi_result;

/**
 * Estimate MOI from guide presence.
 *
 * @param feature_counts    Per-feature count arrays [num_features][n_cells]
 * @param num_features      Number of features
 * @param n_cells           Number of cells
 * @param min_umi           Minimum UMI to count as present (default: 1)
 * @param pmulti_threshold  Threshold for high MOI classification (default: 0.05)
 * @return                  MOI estimation result
 */
nbem_moi_result nbem_estimate_moi(int **feature_counts, int num_features, int n_cells,
                                   int min_umi, double pmulti_threshold);

#ifdef __cplusplus
}
#endif

#endif /* NBEM_H */
