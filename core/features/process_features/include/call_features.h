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

#ifdef __cplusplus
}
#endif

#endif /* CALL_FEATURES_H */
