/**
 * @file call_features.h
 * @brief Dominant feature calling for perturbation assignment
 * 
 * This module implements deterministic feature calling to assign
 * a dominant perturbation per cell (barcode) based on UMI counts.
 * 
 * Algorithm:
 * 1. Read per-cell per-feature UMI counts
 * 2. For each cell:
 *    - If no features pass threshold: feature_call = None
 *    - If one feature: assign it
 *    - If multiple: pick dominant if:
 *      - max_count >= min_deduped_counts
 *      - max_count / total >= dominance_fraction
 *      - max_count - second_max >= dominance_margin
 *    - Else mark as ambiguous
 */

#ifndef CALL_FEATURES_H
#define CALL_FEATURES_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Configuration
 * ============================================================================ */

typedef struct cf_config {
    int min_deduped_counts;      /* Minimum UMI count for a feature to be considered (default: 2) */
    double dominance_fraction;   /* Required fraction of total for dominant (default: 0.8) */
    int dominance_margin;        /* Required count margin over second-best (default: 1) */
    int include_ambiguous;       /* Include ambiguous calls in output (default: 1) */
} cf_config;

/**
 * Create a config with default values.
 */
cf_config* cf_config_create(void);

/**
 * Destroy a config.
 */
void cf_config_destroy(cf_config *config);

/* ============================================================================
 * Call Results
 * ============================================================================ */

/* Call categories */
typedef enum {
    CF_CALL_NONE = 0,         /* No features passed threshold */
    CF_CALL_ASSIGNED = 1,     /* Single dominant feature assigned */
    CF_CALL_AMBIGUOUS = 2,    /* Multiple features, no clear dominant */
} cf_call_type;

/* Per-cell call result */
typedef struct cf_cell_call {
    char *barcode;              /* Cell barcode string */
    cf_call_type call_type;     /* Type of call */
    int feature_index;          /* Assigned feature index (-1 if none/ambiguous) */
    char *feature_name;         /* Assigned feature name (NULL if none/ambiguous) */
    int umi_count;              /* UMI count for assigned feature */
    int total_umi_count;        /* Total UMI count across all features */
    int num_features;           /* Number of features with counts > 0 */
    int second_feature_index;   /* Second-best feature index (-1 if none) */
    int second_umi_count;       /* Second-best feature UMI count */
} cf_cell_call;

/* Collection of all cell calls */
typedef struct cf_call_results {
    cf_cell_call *calls;        /* Array of per-cell calls */
    int num_calls;              /* Number of cells with calls */
    int num_assigned;           /* Number of assigned calls */
    int num_ambiguous;          /* Number of ambiguous calls */
    int num_none;               /* Number of cells with no features */
    char **feature_names;       /* Feature names array (borrowed) */
    int num_features;           /* Number of features */
} cf_call_results;

/**
 * Free call results.
 */
void cf_free_results(cf_call_results *results);

/* ============================================================================
 * MEX Matrix Reader
 * ============================================================================ */

/* Sparse matrix entry */
typedef struct cf_matrix_entry {
    uint32_t row;               /* Feature index (0-based) */
    uint32_t col;               /* Barcode index (0-based) */
    uint32_t value;             /* UMI count */
} cf_matrix_entry;

/* Sparse matrix in COO format */
typedef struct cf_sparse_matrix {
    int num_rows;               /* Number of features */
    int num_cols;               /* Number of barcodes */
    int num_entries;            /* Number of non-zero entries */
    cf_matrix_entry *entries;   /* Array of entries */
    char **row_names;           /* Feature names */
    char **col_names;           /* Barcode names */
    char *row_names_storage;    /* Storage for feature names */
    char *col_names_storage;    /* Storage for barcode names */
} cf_sparse_matrix;

/**
 * Load a sparse matrix from MEX directory (matrix.mtx, barcodes.txt, features.txt).
 * @param mex_dir Directory containing MEX files.
 * @return Sparse matrix, or NULL on failure.
 */
cf_sparse_matrix* cf_load_mex(const char *mex_dir);

/**
 * Free a sparse matrix.
 */
void cf_free_matrix(cf_sparse_matrix *matrix);

/* ============================================================================
 * Feature Calling API
 * ============================================================================ */

/**
 * Call dominant features for all cells in a matrix.
 * @param matrix Sparse matrix of feature x barcode counts.
 * @param config Calling configuration.
 * @return Call results, or NULL on failure.
 */
cf_call_results* cf_call_features(const cf_sparse_matrix *matrix, const cf_config *config);

/**
 * Write call results to a CSV file.
 * @param results Call results to write.
 * @param output_path Path to output CSV file.
 * @return 0 on success, -1 on failure.
 */
int cf_write_calls_csv(const cf_call_results *results, const char *output_path);

/**
 * Write call summary statistics to a file.
 * @param results Call results.
 * @param output_path Path to output summary file.
 * @return 0 on success, -1 on failure.
 */
int cf_write_summary(const cf_call_results *results, const char *output_path);

/**
 * Convenience function: load MEX, call features, write outputs.
 * @param mex_dir Directory containing MEX files.
 * @param output_dir Directory for output files.
 * @param config Calling configuration (NULL for defaults).
 * @return 0 on success, -1 on failure.
 */
int cf_process_mex_dir(const char *mex_dir, const char *output_dir, const cf_config *config);

/* ============================================================================
 * CR9-Compatible GMM Calling
 * ============================================================================ */

/* GMM call configuration */
typedef struct cf_gmm_config {
    int min_umi_threshold;       /* Minimum UMI count for a call (default: 3) */
    int n_init;                  /* Number of GMM initializations (default: 10) */
} cf_gmm_config;

/* Per-feature GMM results */
typedef struct cf_feature_gmm_result {
    int feature_index;
    char *feature_name;
    int *positive_cells;         /* Array of cell indices that are positive */
    int num_positive;            /* Number of positive cells */
    int umi_threshold;           /* Detected UMI threshold */
    int total_umis;              /* Total UMIs for this feature */
} cf_feature_gmm_result;

/* GMM call results (CR9-compatible) */
typedef struct cf_gmm_results {
    /* Per-cell results */
    char **cell_barcodes;        /* Array of cell barcodes */
    char **feature_calls;        /* Feature call string per cell ("|" joined if multiple) */
    int *num_features;           /* Number of features called per cell */
    int *num_umis;               /* UMI counts for called features */
    int num_cells;
    
    /* Per-feature results */
    cf_feature_gmm_result *feature_results;
    int num_features_total;
    char **feature_names;
    
    /* Summary counts */
    int cells_no_molecules;      /* Cells with no UMIs */
    int cells_no_call;           /* Cells with UMIs but no confident call */
    int cells_1_feature;         /* Cells with exactly 1 feature called */
    int cells_multi_feature;     /* Cells with >1 features called */
} cf_gmm_results;

/**
 * Create GMM config with defaults.
 */
cf_gmm_config* cf_gmm_config_create(void);

/**
 * Destroy GMM config.
 */
void cf_gmm_config_destroy(cf_gmm_config *config);

/**
 * Call features using CR9-style GMM algorithm.
 * @param matrix Sparse matrix of feature x barcode counts.
 * @param config GMM configuration (NULL for defaults).
 * @return GMM call results, or NULL on failure.
 */
cf_gmm_results* cf_call_features_gmm(const cf_sparse_matrix *matrix, const cf_gmm_config *config);

/**
 * Free GMM results.
 */
void cf_free_gmm_results(cf_gmm_results *results);

/**
 * Write CR9-compatible protospacer_calls_per_cell.csv.
 * @param results GMM call results.
 * @param output_path Path to output CSV file.
 * @return 0 on success, -1 on failure.
 */
int cf_write_protospacer_calls_per_cell(const cf_gmm_results *results, const char *output_path);

/**
 * Write CR9-compatible protospacer_calls_summary.csv.
 * @param results GMM call results.
 * @param output_path Path to output CSV file.
 * @return 0 on success, -1 on failure.
 */
int cf_write_protospacer_calls_summary(const cf_gmm_results *results, const char *output_path);

/**
 * Write CR9-compatible protospacer_umi_thresholds.csv.
 * @param results GMM call results.
 * @param output_path Path to output CSV file.
 * @return 0 on success, -1 on failure.
 */
int cf_write_protospacer_umi_thresholds(const cf_gmm_results *results, const char *output_path);

/**
 * Write CR9-compatible protospacer_umi_thresholds.json.
 * @param results GMM call results.
 * @param output_path Path to output JSON file.
 * @return 0 on success, -1 on failure.
 */
int cf_write_protospacer_umi_thresholds_json(const cf_gmm_results *results, const char *output_path);

/**
 * Convenience function: load MEX, call features with GMM, write CR9-compatible outputs.
 * @param mex_dir Directory containing MEX files.
 * @param output_dir Directory for output files.
 * @param config GMM configuration (NULL for defaults).
 * @return 0 on success, -1 on failure.
 */
int cf_process_mex_dir_gmm(const char *mex_dir, const char *output_dir, const cf_gmm_config *config);

/* ============================================================================
 * NB-EM Feature Calling (SCEPTRE-style)
 * ============================================================================ */

/** MOI mode for NB-EM */
typedef enum {
    CF_MOI_LOW = 0,     /* Low MOI: single assignment per cell */
    CF_MOI_HIGH = 1,    /* High MOI: allow multiple assignments per cell */
    CF_MOI_AUTO = 2     /* Auto: determine from Poisson fit */
} cf_moi_mode;

/** NB-EM call configuration */
typedef struct cf_nbem_config {
    /* MOI settings */
    cf_moi_mode moi_mode;           /* MOI mode (default: AUTO) */
    int moi_min_umi;                /* Min UMI for guide presence in MOI calc (default: 1) */
    double moi_pmulti_threshold;    /* Threshold for high MOI classification (default: 0.05) */
    
    /* EM settings */
    double prob_threshold;          /* Posterior threshold for call (default: 0.8) */
    int backup_threshold;           /* UMI threshold for fallback (default: 5) */
    int max_iter;                   /* Max EM iterations (default: 100) */
    double tol;                     /* Convergence tolerance (default: 1e-4) */
    int use_poisson;                /* Use Poisson mixture instead of NB (default: 0) */
    int sceptre_parity;             /* SCEPTRE-like behavior (Poisson + per-guide calls, no fallback) */
    
    /* Dispersion */
    double global_phi;              /* Fixed global phi, or 0 for auto (default: 0 = auto) */
    double phi_min;                 /* Minimum phi for MoM estimate (default: 0.01) */
    double phi_max;                 /* Maximum phi for MoM estimate (default: 1000.0) */
    
    /* Covariates */
    const char *covariate_tsv;      /* Optional TSV with GEX covariates (default: NULL) */

    /* Debug instrumentation */
    int debug;                      /* Write per-feature debug CSV (default: 0) */
    int debug_max_features;         /* Limit debug rows (0 = all) */
    const char *debug_csv;          /* Debug CSV path (default: NULL) */
    const char *debug_feature;      /* Feature name to log per-iteration EM */
    const char *debug_iter_csv;     /* Per-iteration debug CSV path */
} cf_nbem_config;

/** Per-feature NB-EM result */
typedef struct cf_feature_nbem_result {
    int feature_index;
    const char *feature_name;
    double pi;                      /* Mixing proportion */
    double delta;                   /* Log fold-change */
    double phi;                     /* Dispersion used */
    double log_likelihood;
    int converged;
    int n_iter;
    int used_fallback;
    int num_positive;               /* Number of positive cells */
    int total_umis;                 /* Total UMIs for this feature */
} cf_feature_nbem_result;

/** MOI estimation result */
typedef struct cf_moi_result {
    cf_moi_mode classification;     /* LOW or HIGH */
    double p0;                      /* Fraction with 0 guides */
    double lambda;                  /* Poisson rate estimate */
    double p_multi_exp;             /* Expected multi-guide fraction */
    double p_multi_obs;             /* Observed multi-guide fraction */
    int n_cells_0;
    int n_cells_1;
    int n_cells_multi;
} cf_moi_result;

/** NB-EM call results */
typedef struct cf_nbem_results {
    /* Per-cell results */
    char **cell_barcodes;           /* Array of cell barcodes (borrowed) */
    char **feature_calls;           /* Feature call string per cell ("|" joined if multiple) */
    int *num_features;              /* Number of features called per cell */
    int *num_umis;                  /* UMI counts for called features */
    int num_cells;
    
    /* Per-cell tracking for LOW-MOI disambiguation (no dense posteriors matrix) */
    double *best_post;              /* Best posterior per cell */
    double *second_post;            /* Second-best posterior per cell */
    int *best_feature;              /* Feature index of best posterior (-1 if none) */
    int *second_feature;            /* Feature index of second-best posterior (-1 if none) */
    
    /* Per-feature results */
    cf_feature_nbem_result *feature_results;
    int num_features_total;
    char **feature_names;           /* Feature names array (borrowed) */
    
    /* Summary counts */
    int cells_no_molecules;
    int cells_no_call;
    int cells_1_feature;
    int cells_multi_feature;
    int num_em_failures;            /* Number of features where EM failed (used fallback) */
    
    /* MOI result */
    cf_moi_result moi;
    
    /* Global parameters */
    double global_phi;              /* Global dispersion used */
    
    /* Posteriors output file (written during processing, not stored in memory) */
    const char *posteriors_mtx_path; /* Path to posteriors MTX file (borrowed) */
    int posteriors_nnz;             /* Number of non-zero posteriors written */
} cf_nbem_results;

/**
 * Create NB-EM config with defaults.
 */
cf_nbem_config* cf_nbem_config_create(void);

/**
 * Destroy NB-EM config.
 */
void cf_nbem_config_destroy(cf_nbem_config *config);

/**
 * Call features using NB-EM algorithm (SCEPTRE-style).
 * @param matrix Sparse matrix of feature x barcode counts.
 * @param config NB-EM configuration (NULL for defaults).
 * @return NB-EM call results, or NULL on failure.
 */
cf_nbem_results* cf_call_features_nbem(const cf_sparse_matrix *matrix, const cf_nbem_config *config);

/**
 * Free NB-EM results.
 */
void cf_free_nbem_results(cf_nbem_results *results);

/**
 * Write NB-EM protospacer_calls_per_cell.csv (same format as GMM).
 */
int cf_write_nbem_calls_per_cell(const cf_nbem_results *results, const char *output_path);

/**
 * Write NB-EM protospacer_calls_summary.csv (same format as GMM).
 */
int cf_write_nbem_calls_summary(const cf_nbem_results *results, const char *output_path);

/**
 * Write NB-EM feature parameters (pi, delta, phi, etc.).
 */
int cf_write_nbem_feature_params(const cf_nbem_results *results, const char *output_path);

/**
 * Finalize NB-EM cell posteriors MatrixMarket file.
 * Posteriors are streamed during processing; this adds the header.
 * @param temp_path Path to temp file with raw entries
 * @param output_path Final output path
 * @param num_features Number of features (rows)
 * @param num_cells Number of cells (columns)
 * @param nnz Number of non-zero entries written
 * @return 0 on success, -1 on failure
 */
int cf_finalize_nbem_posteriors_mtx(const char *temp_path, const char *output_path,
                                     int num_features, int num_cells, int nnz);

/**
 * Write NB-EM summary JSON (MOI, parameters, etc.).
 */
int cf_write_nbem_summary_json(const cf_nbem_results *results, const cf_nbem_config *config,
                                const char *output_path);

/**
 * Convenience function: load MEX, call features with NB-EM, write outputs.
 * @param mex_dir Directory containing MEX files.
 * @param output_dir Directory for output files.
 * @param config NB-EM configuration (NULL for defaults).
 * @return 0 on success, -1 on failure.
 */
int cf_process_mex_dir_nbem(const char *mex_dir, const char *output_dir, const cf_nbem_config *config);

#ifdef __cplusplus
}
#endif

#endif /* CALL_FEATURES_H */
