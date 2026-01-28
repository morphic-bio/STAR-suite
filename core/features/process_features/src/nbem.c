/**
 * @file nbem.c
 * @brief NB-EM feature calling implementation (SCEPTRE-style)
 */

#include "../include/nbem.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

/* ============================================================================
 * Constants
 * ============================================================================ */

#define NBEM_EPS 1e-10
#define NBEM_LOG_EPS -23.025850929940457  /* log(1e-10) */

/* ============================================================================
 * Householder QR Decomposition
 * ============================================================================ */

/**
 * Compute the 2-norm of a vector segment.
 */
static double vec_norm2(const double *v, int n, int stride) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double val = v[i * stride];
        sum += val * val;
    }
    return sqrt(sum);
}

/**
 * Apply Householder reflection H = I - 2*v*v'/||v||^2 to a matrix column.
 * Modifies A in place. v is the Householder vector.
 */
static void apply_householder(double *A, int m, int n, int col, const double *v, double v_norm_sq) {
    if (v_norm_sq < NBEM_EPS) return;
    
    double scale = 2.0 / v_norm_sq;
    
    /* Apply to each column j >= col */
    for (int j = col; j < n; j++) {
        /* Compute v' * A(:,j) for rows >= col */
        double dot = 0.0;
        for (int i = col; i < m; i++) {
            dot += v[i - col] * A[j * m + i];
        }
        
        /* A(:,j) -= (2 * dot / ||v||^2) * v */
        double coef = scale * dot;
        for (int i = col; i < m; i++) {
            A[j * m + i] -= coef * v[i - col];
        }
    }
}

/**
 * Apply Householder reflection to a vector (RHS in least squares).
 */
static void apply_householder_vec(double *b, int m, int col, const double *v, double v_norm_sq) {
    if (v_norm_sq < NBEM_EPS) return;
    
    double scale = 2.0 / v_norm_sq;
    double dot = 0.0;
    for (int i = col; i < m; i++) {
        dot += v[i - col] * b[i];
    }
    double coef = scale * dot;
    for (int i = col; i < m; i++) {
        b[i] -= coef * v[i - col];
    }
}

int qr_solve(double *A, double *b, int n) {
    return qr_lstsq(A, b, n, n);
}

int qr_lstsq(double *A, double *b, int m, int n) {
    if (m < n || n <= 0) return -1;
    
    /* Allocate Householder vector */
    double *v = (double *)malloc(m * sizeof(double));
    if (!v) return -1;
    
    /* QR decomposition with Householder reflections */
    for (int k = 0; k < n; k++) {
        /* Extract column k from row k onwards */
        int len = m - k;
        for (int i = 0; i < len; i++) {
            v[i] = A[k * m + k + i];
        }
        
        /* Compute Householder vector: v = x - ||x|| * e1 (or +||x|| for stability) */
        double norm_x = vec_norm2(v, len, 1);
        if (norm_x < NBEM_EPS) {
            /* Column is essentially zero - rank deficient */
            free(v);
            return -1;
        }
        
        /* Choose sign to avoid cancellation */
        double sign = (v[0] >= 0) ? 1.0 : -1.0;
        v[0] += sign * norm_x;
        
        /* Compute ||v||^2 */
        double v_norm_sq = 0.0;
        for (int i = 0; i < len; i++) {
            v_norm_sq += v[i] * v[i];
        }
        
        /* Apply reflection to A and b */
        apply_householder(A, m, n, k, v, v_norm_sq);
        apply_householder_vec(b, m, k, v, v_norm_sq);
    }
    
    free(v);
    
    /* Back substitution: R * x = Q' * b (stored in b[0:n]) */
    for (int i = n - 1; i >= 0; i--) {
        double diag = A[i * m + i];
        if (fabs(diag) < NBEM_EPS) {
            return -1;  /* Singular */
        }
        
        double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[j * m + i] * b[j];
        }
        b[i] = sum / diag;
    }
    
    return 0;
}

/* ============================================================================
 * Poisson IRLS
 * ============================================================================ */

int poisson_irls(const double *X, const int *y, int n_cells, int n_covars,
                 double *beta, int max_iter, double tol) {
    if (!X || !y || !beta || n_cells <= 0 || n_covars <= 0) return -1;
    if (max_iter <= 0) max_iter = 25;
    if (tol <= 0) tol = 1e-6;
    
    /* Allocate working arrays */
    double *mu = (double *)malloc(n_cells * sizeof(double));
    double *z = (double *)malloc(n_cells * sizeof(double));
    double *w = (double *)malloc(n_cells * sizeof(double));
    double *XtWX = (double *)malloc(n_covars * n_covars * sizeof(double));
    double *XtWz = (double *)malloc(n_covars * sizeof(double));
    double *X_work = (double *)malloc(n_cells * n_covars * sizeof(double));
    
    if (!mu || !z || !w || !XtWX || !XtWz || !X_work) {
        free(mu); free(z); free(w); free(XtWX); free(XtWz); free(X_work);
        return -1;
    }
    
    /* Initialize beta to zero (intercept-only gives mu = 1) */
    memset(beta, 0, n_covars * sizeof(double));
    
    /* Compute initial mean: use data mean for intercept */
    double y_mean = 0.0;
    for (int i = 0; i < n_cells; i++) {
        y_mean += y[i];
    }
    y_mean /= n_cells;
    if (y_mean < NBEM_EPS) y_mean = NBEM_EPS;
    beta[0] = log(y_mean);  /* Assumes first column is intercept */
    
    double prev_deviance = DBL_MAX;
    int iter;
    
    for (iter = 0; iter < max_iter; iter++) {
        /* Compute mu = exp(X * beta) */
        for (int i = 0; i < n_cells; i++) {
            double eta = 0.0;
            for (int j = 0; j < n_covars; j++) {
                eta += X[j * n_cells + i] * beta[j];
            }
            /* Clamp eta to avoid overflow */
            if (eta > 20.0) eta = 20.0;
            if (eta < -20.0) eta = -20.0;
            mu[i] = exp(eta);
            if (mu[i] < NBEM_EPS) mu[i] = NBEM_EPS;
        }
        
        /* Compute deviance */
        double deviance = 0.0;
        for (int i = 0; i < n_cells; i++) {
            if (y[i] > 0) {
                deviance += 2.0 * (y[i] * log((double)y[i] / mu[i]) - (y[i] - mu[i]));
            } else {
                deviance += 2.0 * mu[i];
            }
        }
        
        /* Check convergence */
        if (fabs(deviance - prev_deviance) < tol * (fabs(prev_deviance) + 1.0)) {
            break;
        }
        prev_deviance = deviance;
        
        /* Compute working response z and weights w */
        for (int i = 0; i < n_cells; i++) {
            double eta = log(mu[i]);
            z[i] = eta + (y[i] - mu[i]) / mu[i];
            w[i] = mu[i];  /* Poisson variance = mean */
            if (w[i] < NBEM_EPS) w[i] = NBEM_EPS;
        }
        
        /* Form weighted normal equations: X'WX * beta = X'Wz */
        /* XtWX[j,k] = sum_i X[i,j] * w[i] * X[i,k] */
        memset(XtWX, 0, n_covars * n_covars * sizeof(double));
        memset(XtWz, 0, n_covars * sizeof(double));
        
        for (int j = 0; j < n_covars; j++) {
            for (int k = j; k < n_covars; k++) {
                double sum = 0.0;
                for (int i = 0; i < n_cells; i++) {
                    sum += X[j * n_cells + i] * w[i] * X[k * n_cells + i];
                }
                XtWX[k * n_covars + j] = sum;
                XtWX[j * n_covars + k] = sum;  /* Symmetric */
            }
            
            double sumz = 0.0;
            for (int i = 0; i < n_cells; i++) {
                sumz += X[j * n_cells + i] * w[i] * z[i];
            }
            XtWz[j] = sumz;
        }
        
        /* Solve normal equations */
        memcpy(X_work, XtWX, n_covars * n_covars * sizeof(double));
        if (qr_solve(X_work, XtWz, n_covars) != 0) {
            free(mu); free(z); free(w); free(XtWX); free(XtWz); free(X_work);
            return -1;
        }
        
        memcpy(beta, XtWz, n_covars * sizeof(double));
    }
    
    free(mu); free(z); free(w); free(XtWX); free(XtWz); free(X_work);
    return 0;
}

void poisson_predict(const double *X, const double *beta, int n_cells, int n_covars,
                     double *mu) {
    for (int i = 0; i < n_cells; i++) {
        double eta = 0.0;
        for (int j = 0; j < n_covars; j++) {
            eta += X[j * n_cells + i] * beta[j];
        }
        if (eta > 20.0) eta = 20.0;
        if (eta < -20.0) eta = -20.0;
        mu[i] = exp(eta);
        if (mu[i] < NBEM_EPS) mu[i] = NBEM_EPS;
    }
}

/* ============================================================================
 * Negative Binomial Distribution
 * ============================================================================ */

/**
 * Log gamma function (simple approximation for positive integers).
 * For k >= 1, lgamma(k+1) = log(k!)
 */
static double log_factorial(int k) {
    if (k <= 1) return 0.0;
    return lgamma((double)(k + 1));
}

static double poisson_logpmf(int k, double mu) {
    if (mu <= 0.0) mu = NBEM_EPS;
    return k * log(mu) - mu - log_factorial(k);
}

double nb_logpmf(int k, double mu, double phi) {
    if (mu <= 0.0) mu = NBEM_EPS;
    if (phi <= 0.0) phi = NBEM_EPS;
    
    /* NB parameterization: mean = mu, var = mu + mu^2/phi
     * P(X=k) = Gamma(k+phi)/(k!*Gamma(phi)) * (phi/(phi+mu))^phi * (mu/(phi+mu))^k
     */
    double r = phi;  /* "size" parameter */
    double p = phi / (phi + mu);  /* "prob" parameter */
    
    double log_p = lgamma(k + r) - log_factorial(k) - lgamma(r)
                   + r * log(p) + k * log(1.0 - p);
    
    return log_p;
}

double nb_estimate_phi_mom(const int *counts, int n, double phi_min, double phi_max) {
    if (n <= 1) return phi_min;
    
    /* Compute mean and variance */
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        sum += counts[i];
        sum_sq += (double)counts[i] * counts[i];
    }
    
    double mean = sum / n;
    double var = (sum_sq - sum * sum / n) / (n - 1);  /* Bessel's correction */
    
    if (mean < NBEM_EPS) return phi_max;  /* All zeros, no overdispersion */
    if (var <= mean) return phi_max;  /* Underdispersed, use high phi */
    
    /* Method of moments: var = mean + mean^2/phi => phi = mean^2 / (var - mean) */
    double phi = mean * mean / (var - mean);
    
    if (phi < phi_min) phi = phi_min;
    if (phi > phi_max) phi = phi_max;
    
    return phi;
}

/* Streaming phi accumulator using Welford's online algorithm */

void nbem_phi_accum_init(nbem_phi_accum *acc) {
    if (!acc) return;
    acc->n = 0;
    acc->mean = 0.0;
    acc->M2 = 0.0;
}

void nbem_phi_accum_add(nbem_phi_accum *acc, const int *counts, int n) {
    if (!acc || !counts || n <= 0) return;
    
    /* Welford's online algorithm for mean and variance */
    for (int i = 0; i < n; i++) {
        double x = (double)counts[i];
        acc->n++;
        double delta = x - acc->mean;
        acc->mean += delta / acc->n;
        double delta2 = x - acc->mean;
        acc->M2 += delta * delta2;
    }
}

double nbem_phi_accum_finish(const nbem_phi_accum *acc, double phi_min, double phi_max) {
    if (!acc || acc->n <= 1) return phi_min;
    
    double mean = acc->mean;
    double var = acc->M2 / (acc->n - 1);  /* Bessel's correction */
    
    if (mean < NBEM_EPS) return phi_max;  /* All zeros, no overdispersion */
    if (var <= mean) return phi_max;  /* Underdispersed, use high phi */
    
    /* Method of moments: var = mean + mean^2/phi => phi = mean^2 / (var - mean) */
    double phi = mean * mean / (var - mean);
    
    if (phi < phi_min) phi = phi_min;
    if (phi > phi_max) phi = phi_max;
    
    return phi;
}

/* ============================================================================
 * NB-EM Per-Feature Fitting
 * ============================================================================ */

/**
 * Log-sum-exp for numerical stability.
 */
static double log_sum_exp(double a, double b) {
    if (a > b) {
        return a + log(1.0 + exp(b - a));
    } else {
        return b + log(1.0 + exp(a - b));
    }
}


nbem_fit_result nbem_fit_feature(const int *counts, const double *mu0, int n_cells,
                                  double phi, int use_poisson, int max_iter, double tol,
                                  double pi_init, double delta_init) {
    nbem_fit_result result;
    memset(&result, 0, sizeof(result));
    
    if (!counts || !mu0 || n_cells <= 0) {
        result.failure_reason = NBEM_FAIL_INVALID_INPUT;
        result.used_fallback = 1;
        return result;
    }
    
    /* Set defaults */
    if (max_iter <= 0) max_iter = 100;
    if (tol <= 0) tol = 1e-4;
    if (pi_init <= 0 || pi_init >= 1) pi_init = 0.1;
    if (delta_init <= 0) delta_init = 1.0;
    if (phi <= 0) phi = 1.0;
    
    /* Initialize parameters */
    double pi = pi_init;
    double delta = delta_init;
    result.phi = phi;
    
    /* Allocate responsibilities */
    double *gamma = (double *)malloc(n_cells * sizeof(double));
    if (!gamma) {
        result.failure_reason = NBEM_FAIL_ALLOC;
        result.used_fallback = 1;
        return result;
    }
    
    double prev_ll = -DBL_MAX;
    int iter;
    
    for (iter = 0; iter < max_iter; iter++) {
        /* E-step: compute responsibilities gamma_i = P(perturbed | count_i) */
        double ll = 0.0;
        double sum_gamma = 0.0;
        
        for (int i = 0; i < n_cells; i++) {
            double mu1 = mu0[i] * exp(delta);
            
            double log_p0 = log(1.0 - pi + NBEM_EPS) +
                (use_poisson ? poisson_logpmf(counts[i], mu0[i]) : nb_logpmf(counts[i], mu0[i], phi));
            double log_p1 = log(pi + NBEM_EPS) +
                (use_poisson ? poisson_logpmf(counts[i], mu1) : nb_logpmf(counts[i], mu1, phi));
            double log_sum = log_sum_exp(log_p0, log_p1);
            
            gamma[i] = exp(log_p1 - log_sum);
            ll += log_sum;
            sum_gamma += gamma[i];
            
            /* Check for NaN */
            if (isnan(gamma[i]) || isinf(gamma[i])) {
                gamma[i] = 0.5;
            }
        }
        
        /* Check for numerical issues */
        if (isnan(ll) || isinf(ll)) {
            result.converged = 0;
            result.failure_reason = NBEM_FAIL_NAN_INF;
            result.used_fallback = 1;
            break;
        }
        
        /* Check convergence */
        if (iter > 0 && fabs(ll - prev_ll) < tol * (fabs(prev_ll) + 1.0)) {
            result.converged = 1;
            result.log_likelihood = ll;
            break;
        }

        /* Check for likelihood decrease (numerical instability) */
        if (iter > 0 && ll < prev_ll - tol) {
            /* Likelihood decreased: treat as failure and trigger fallback */
            result.converged = 0;
            result.failure_reason = NBEM_FAIL_LL_DECREASE;
            result.used_fallback = 1;
            break;
        }
        
        prev_ll = ll;
        result.log_likelihood = ll;
        
        /* M-step: update pi and delta */
        /* pi = mean(gamma) */
        double new_pi = sum_gamma / n_cells;
        if (new_pi < NBEM_EPS) new_pi = NBEM_EPS;
        if (new_pi > 1.0 - NBEM_EPS) new_pi = 1.0 - NBEM_EPS;
        
        /* delta: weighted mean log-ratio estimate */
        double sum_weighted_ratio = 0.0;
        double sum_weight = 0.0;
        for (int i = 0; i < n_cells; i++) {
            if (counts[i] > 0 && mu0[i] > NBEM_EPS) {
                sum_weighted_ratio += gamma[i] * (double)counts[i] / mu0[i];
                sum_weight += gamma[i];
            }
        }
        
        double new_delta;
        if (sum_weight > NBEM_EPS) {
            double ratio = sum_weighted_ratio / sum_weight;
            if (ratio < NBEM_EPS) ratio = NBEM_EPS;
            new_delta = log(ratio);
        } else {
            new_delta = delta;  /* Keep previous */
        }
        
        /* Clamp delta to reasonable range */
        if (new_delta < -5.0) new_delta = -5.0;
        if (new_delta > 10.0) new_delta = 10.0;

        /* Guarded delta update: back off if LL decreases */
        double step = 1.0;
        double ll_probe = ll;
        double pi_probe = new_pi;
        double delta_probe = new_delta;
        int improved = 0;
        for (int ls_iter = 0; ls_iter < 6; ls_iter++) {
            double try_delta = delta + step * (new_delta - delta);
            double try_ll = 0.0;
            for (int i = 0; i < n_cells; i++) {
                double mu1 = mu0[i] * exp(try_delta);
                double log_p0 = log(1.0 - pi_probe + NBEM_EPS) +
                    (use_poisson ? poisson_logpmf(counts[i], mu0[i]) : nb_logpmf(counts[i], mu0[i], phi));
                double log_p1 = log(pi_probe + NBEM_EPS) +
                    (use_poisson ? poisson_logpmf(counts[i], mu1) : nb_logpmf(counts[i], mu1, phi));
                try_ll += log_sum_exp(log_p0, log_p1);
            }
            if (try_ll >= ll) {
                ll_probe = try_ll;
                delta_probe = try_delta;
                improved = 1;
                break;
            }
            step *= 0.5;
        }

        if (!improved) {
            result.converged = 0;
            result.failure_reason = NBEM_FAIL_LL_DECREASE;
            result.used_fallback = 1;
            break;
        }
        
        /* Check for degenerate solution */
        if (new_pi < 1e-6 || new_pi > 1.0 - 1e-6) {
            result.converged = 0;
            result.failure_reason = NBEM_FAIL_DEGENERATE_PI;
            result.used_fallback = 1;
            break;
        }
        
        pi = new_pi;
        delta = delta_probe;
        result.log_likelihood = ll_probe;
    }
    
    result.pi = pi;
    result.delta = delta;
    result.n_iter = iter;
    
    /* If didn't converge and no fallback yet, mark as needing fallback */
    if (!result.converged) {
        if (result.failure_reason == NBEM_FAIL_NONE) {
            result.failure_reason = NBEM_FAIL_MAX_ITER;
        }
        result.used_fallback = 1;
    }
    
    free(gamma);
    return result;
}

nbem_fit_result nbem_fit_feature_debug(const int *counts, const double *mu0, int n_cells,
                                       double phi, int max_iter, double tol,
                                       double pi_init, double delta_init,
                                       int use_poisson, const char *feature_name, FILE *debug_fp) {
    nbem_fit_result result;
    memset(&result, 0, sizeof(result));

    if (debug_fp) {
        fprintf(debug_fp, "feature,iter,ll,pi,delta,sum_gamma,min_gamma,max_gamma,mu0_mean,mu0_min,mu0_max,mu1_mean,mu1_min,mu1_max,converged,used_fallback,failure_reason\n");
    }

    if (!counts || !mu0 || n_cells <= 0) {
        result.failure_reason = NBEM_FAIL_INVALID_INPUT;
        result.used_fallback = 1;
        return result;
    }

    if (max_iter <= 0) max_iter = 100;
    if (tol <= 0) tol = 1e-4;
    if (pi_init <= 0 || pi_init >= 1) pi_init = 0.1;
    if (delta_init <= 0) delta_init = 1.0;
    if (phi <= 0) phi = 1.0;

    double pi = pi_init;
    double delta = delta_init;
    result.phi = phi;

    double *gamma = (double *)malloc(n_cells * sizeof(double));
    if (!gamma) {
        result.failure_reason = NBEM_FAIL_ALLOC;
        result.used_fallback = 1;
        return result;
    }

    double prev_ll = -DBL_MAX;
    int iter;

    for (iter = 0; iter < max_iter; iter++) {
        double ll = 0.0;
        double sum_gamma = 0.0;
        double min_gamma = 1.0;
        double max_gamma = 0.0;
        double mu0_sum = 0.0;
        double mu0_min = DBL_MAX;
        double mu0_max = 0.0;
        double mu1_sum = 0.0;
        double mu1_min = DBL_MAX;
        double mu1_max = 0.0;

        for (int i = 0; i < n_cells; i++) {
            double mu1 = mu0[i] * exp(delta);
            double log_p0 = log(1.0 - pi + NBEM_EPS) +
                (use_poisson ? poisson_logpmf(counts[i], mu0[i]) : nb_logpmf(counts[i], mu0[i], phi));
            double log_p1 = log(pi + NBEM_EPS) +
                (use_poisson ? poisson_logpmf(counts[i], mu1) : nb_logpmf(counts[i], mu1, phi));
            double log_sum = log_sum_exp(log_p0, log_p1);

            gamma[i] = exp(log_p1 - log_sum);
            ll += log_sum;
            sum_gamma += gamma[i];
            if (gamma[i] < min_gamma) min_gamma = gamma[i];
            if (gamma[i] > max_gamma) max_gamma = gamma[i];
            mu0_sum += mu0[i];
            if (mu0[i] < mu0_min) mu0_min = mu0[i];
            if (mu0[i] > mu0_max) mu0_max = mu0[i];
            mu1_sum += mu1;
            if (mu1 < mu1_min) mu1_min = mu1;
            if (mu1 > mu1_max) mu1_max = mu1;

            if (isnan(gamma[i]) || isinf(gamma[i])) {
                gamma[i] = 0.5;
            }
        }

        if (debug_fp) {
            double mu0_mean = mu0_sum / n_cells;
            double mu1_mean = mu1_sum / n_cells;
            fprintf(debug_fp, "%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d\n",
                    feature_name ? feature_name : "Unknown",
                    iter, ll, pi, delta, sum_gamma, min_gamma, max_gamma,
                    mu0_mean, mu0_min, mu0_max,
                    mu1_mean, mu1_min, mu1_max,
                    result.converged, result.used_fallback, result.failure_reason);
        }

        if (isnan(ll) || isinf(ll)) {
            result.converged = 0;
            result.failure_reason = NBEM_FAIL_NAN_INF;
            result.used_fallback = 1;
            break;
        }

        if (iter > 0 && fabs(ll - prev_ll) < tol * (fabs(prev_ll) + 1.0)) {
            result.converged = 1;
            result.log_likelihood = ll;
            break;
        }

        prev_ll = ll;
        result.log_likelihood = ll;

        double new_pi = sum_gamma / n_cells;
        if (new_pi < NBEM_EPS) new_pi = NBEM_EPS;
        if (new_pi > 1.0 - NBEM_EPS) new_pi = 1.0 - NBEM_EPS;

        double sum_weighted_ratio = 0.0;
        double sum_weight = 0.0;
        for (int i = 0; i < n_cells; i++) {
            if (counts[i] > 0 && mu0[i] > NBEM_EPS) {
                sum_weighted_ratio += gamma[i] * (double)counts[i] / mu0[i];
                sum_weight += gamma[i];
            }
        }

        double new_delta;
        if (sum_weight > NBEM_EPS) {
            double ratio = sum_weighted_ratio / sum_weight;
            if (ratio < NBEM_EPS) ratio = NBEM_EPS;
            new_delta = log(ratio);
        } else {
            new_delta = delta;
        }

        if (new_delta < -5.0) new_delta = -5.0;
        if (new_delta > 10.0) new_delta = 10.0;

        if (new_pi < 1e-6 || new_pi > 1.0 - 1e-6) {
            result.converged = 0;
            result.failure_reason = NBEM_FAIL_DEGENERATE_PI;
            result.used_fallback = 1;
            break;
        }

        pi = new_pi;
        delta = new_delta;
    }

    result.pi = pi;
    result.delta = delta;
    result.n_iter = iter;

    if (!result.converged) {
        if (result.failure_reason == NBEM_FAIL_NONE) {
            result.failure_reason = NBEM_FAIL_MAX_ITER;
        }
        result.used_fallback = 1;
    }

    free(gamma);
    return result;
}

void nbem_posteriors(const int *counts, const double *mu0, int n_cells,
                     double pi, double delta, double phi, int use_poisson, double *posteriors) {
    if (!counts || !mu0 || !posteriors) return;
    
    for (int i = 0; i < n_cells; i++) {
        double mu1 = mu0[i] * exp(delta);
        
        double log_p0 = log(1.0 - pi + NBEM_EPS) +
            (use_poisson ? poisson_logpmf(counts[i], mu0[i]) : nb_logpmf(counts[i], mu0[i], phi));
        double log_p1 = log(pi + NBEM_EPS) +
            (use_poisson ? poisson_logpmf(counts[i], mu1) : nb_logpmf(counts[i], mu1, phi));
        double log_sum = log_sum_exp(log_p0, log_p1);
        
        posteriors[i] = exp(log_p1 - log_sum);
        
        if (isnan(posteriors[i]) || isinf(posteriors[i])) {
            posteriors[i] = 0.0;
        }
    }
}

int nbem_make_calls(const double *posteriors, const int *counts, int n_cells,
                    double prob_threshold, int backup_threshold, int use_fallback,
                    int *calls) {
    if (!calls) return 0;
    
    if (prob_threshold <= 0) prob_threshold = 0.8;
    if (backup_threshold <= 0) backup_threshold = 5;
    
    int n_positive = 0;
    
    for (int i = 0; i < n_cells; i++) {
        if (use_fallback) {
            calls[i] = (counts[i] >= backup_threshold) ? 1 : 0;
        } else {
            calls[i] = (posteriors && posteriors[i] >= prob_threshold) ? 1 : 0;
        }
        n_positive += calls[i];
    }
    
    return n_positive;
}

/* ============================================================================
 * Covariate Computation
 * ============================================================================ */

nbem_covariates* nbem_compute_grna_covariates(int **feature_counts, int num_features, int n_cells) {
    if (!feature_counts || num_features <= 0 || n_cells <= 0) return NULL;
    
    nbem_covariates *cov = (nbem_covariates *)calloc(1, sizeof(nbem_covariates));
    if (!cov) return NULL;
    
    cov->n_cells = n_cells;
    cov->has_gex = 0;
    
    cov->grna_log_n_umis = (double *)malloc(n_cells * sizeof(double));
    cov->grna_log_n_nonzero = (double *)malloc(n_cells * sizeof(double));
    
    if (!cov->grna_log_n_umis || !cov->grna_log_n_nonzero) {
        nbem_free_covariates(cov);
        return NULL;
    }
    
    for (int i = 0; i < n_cells; i++) {
        int total_umis = 0;
        int n_nonzero = 0;
        
        for (int f = 0; f < num_features; f++) {
            int count = feature_counts[f][i];
            total_umis += count;
            if (count > 0) n_nonzero++;
        }
        
        cov->grna_log_n_umis[i] = log(1.0 + total_umis);
        cov->grna_log_n_nonzero[i] = log(1.0 + n_nonzero);
    }
    
    return cov;
}

int nbem_load_gex_covariates(nbem_covariates *covariates, const char *tsv_path,
                              char **barcodes) {
    if (!covariates || !tsv_path || !barcodes) return -1;
    
    FILE *fp = fopen(tsv_path, "r");
    if (!fp) {
        fprintf(stderr, "Warning: could not open covariate file %s, using gRNA-only\n", tsv_path);
        return -1;
    }
    
    /* Allocate GEX covariate arrays */
    covariates->gex_log_n_umis = (double *)calloc(covariates->n_cells, sizeof(double));
    covariates->gex_log_n_nonzero = (double *)calloc(covariates->n_cells, sizeof(double));
    
    if (!covariates->gex_log_n_umis || !covariates->gex_log_n_nonzero) {
        fclose(fp);
        return -1;
    }
    
    /* Build barcode hash for fast lookup */
    /* For simplicity, use linear search (OK for typical cell counts) */
    /* TODO: use khash for large datasets */
    
    char line[4096];
    int header_skipped = 0;
    int n_matched = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        /* Skip header */
        if (!header_skipped) {
            if (strstr(line, "barcode") || strstr(line, "response")) {
                header_skipped = 1;
                continue;
            }
            header_skipped = 1;  /* No header, process this line */
        }
        
        /* Parse: barcode,response_n_umis,response_n_nonzero */
        char barcode[256];
        int n_umis, n_nonzero;
        
        if (sscanf(line, "%255[^,\t],%d,%d", barcode, &n_umis, &n_nonzero) >= 3 ||
            sscanf(line, "%255[^,\t]\t%d\t%d", barcode, &n_umis, &n_nonzero) >= 3) {
            
            /* Find barcode in cell list */
            for (int i = 0; i < covariates->n_cells; i++) {
                if (strcmp(barcodes[i], barcode) == 0) {
                    covariates->gex_log_n_umis[i] = log(1.0 + n_umis);
                    covariates->gex_log_n_nonzero[i] = log(1.0 + n_nonzero);
                    n_matched++;
                    break;
                }
            }
        }
    }
    
    fclose(fp);
    
    if (n_matched == 0) {
        fprintf(stderr, "Warning: no barcodes matched in covariate file\n");
        free(covariates->gex_log_n_umis);
        free(covariates->gex_log_n_nonzero);
        covariates->gex_log_n_umis = NULL;
        covariates->gex_log_n_nonzero = NULL;
        return -1;
    }
    
    covariates->has_gex = 1;
    fprintf(stderr, "Loaded GEX covariates for %d/%d cells\n", n_matched, covariates->n_cells);
    
    return 0;
}

int nbem_build_design_matrix(const nbem_covariates *covariates, double *X) {
    if (!covariates || !X) return 0;
    
    int n = covariates->n_cells;
    int p = covariates->has_gex ? 5 : 3;  /* intercept + 2 or 4 covariates */
    
    /* Column 0: intercept */
    for (int i = 0; i < n; i++) {
        X[0 * n + i] = 1.0;
    }
    
    /* Column 1: grna_log_n_umis */
    for (int i = 0; i < n; i++) {
        X[1 * n + i] = covariates->grna_log_n_umis[i];
    }
    
    /* Column 2: grna_log_n_nonzero */
    for (int i = 0; i < n; i++) {
        X[2 * n + i] = covariates->grna_log_n_nonzero[i];
    }
    
    if (covariates->has_gex) {
        /* Column 3: gex_log_n_umis */
        for (int i = 0; i < n; i++) {
            X[3 * n + i] = covariates->gex_log_n_umis[i];
        }
        
        /* Column 4: gex_log_n_nonzero */
        for (int i = 0; i < n; i++) {
            X[4 * n + i] = covariates->gex_log_n_nonzero[i];
        }
    }
    
    return p;
}

void nbem_free_covariates(nbem_covariates *covariates) {
    if (!covariates) return;
    free(covariates->grna_log_n_umis);
    free(covariates->grna_log_n_nonzero);
    free(covariates->gex_log_n_umis);
    free(covariates->gex_log_n_nonzero);
    free(covariates);
}

/* ============================================================================
 * MOI Estimation
 * ============================================================================ */

nbem_moi_result nbem_estimate_moi(int **feature_counts, int num_features, int n_cells,
                                   int min_umi, double pmulti_threshold) {
    nbem_moi_result result;
    memset(&result, 0, sizeof(result));
    
    if (!feature_counts || num_features <= 0 || n_cells <= 0) {
        result.classification = NBEM_MOI_LOW;
        return result;
    }
    
    if (min_umi <= 0) min_umi = 1;
    if (pmulti_threshold <= 0) pmulti_threshold = 0.05;
    
    /* Count cells by number of guides present */
    int *n_guides_per_cell = (int *)calloc(n_cells, sizeof(int));
    if (!n_guides_per_cell) {
        result.classification = NBEM_MOI_LOW;
        return result;
    }
    
    for (int i = 0; i < n_cells; i++) {
        for (int f = 0; f < num_features; f++) {
            if (feature_counts[f][i] >= min_umi) {
                n_guides_per_cell[i]++;
            }
        }
    }
    
    /* Count cells with 0, 1, >=2 guides */
    for (int i = 0; i < n_cells; i++) {
        if (n_guides_per_cell[i] == 0) {
            result.n_cells_0++;
        } else if (n_guides_per_cell[i] == 1) {
            result.n_cells_1++;
        } else {
            result.n_cells_multi++;
        }
    }
    
    free(n_guides_per_cell);
    
    /* Compute statistics */
    result.p0 = (double)result.n_cells_0 / n_cells;
    result.p_multi_obs = (double)result.n_cells_multi / n_cells;
    
    /* Poisson fit */
    if (result.p0 > 0 && result.p0 < 1.0) {
        result.lambda = -log(result.p0);
        /* P(X >= 2) = 1 - P(X=0) - P(X=1) = 1 - exp(-lambda) - lambda*exp(-lambda) */
        result.p_multi_exp = 1.0 - (1.0 + result.lambda) * exp(-result.lambda);
    } else if (result.p0 == 0) {
        result.lambda = 10.0;  /* High MOI */
        result.p_multi_exp = 1.0;
    } else {
        result.lambda = 0.0;
        result.p_multi_exp = 0.0;
    }
    
    /* Classify */
    if (result.p_multi_obs >= pmulti_threshold || result.p_multi_exp >= pmulti_threshold) {
        result.classification = NBEM_MOI_HIGH;
    } else {
        result.classification = NBEM_MOI_LOW;
    }
    
    return result;
}
