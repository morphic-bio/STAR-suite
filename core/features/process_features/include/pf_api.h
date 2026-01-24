/**
 * @file pf_api.h
 * @brief Public API for process_features library
 * 
 * This header provides a clean C API for feature barcode processing,
 * suitable for integration into STAR-suite or other pipelines.
 * 
 * Usage:
 *   1. Create a config with pf_config_create()
 *   2. Set options on the config
 *   3. Initialize context with pf_init()
 *   4. Load references with pf_load_feature_ref()
 *   5. Process FASTQs with pf_process_fastqs()
 *   6. Clean up with pf_destroy()
 */

#ifndef PF_API_H
#define PF_API_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handles */
typedef struct pf_config pf_config;
typedef struct pf_context pf_context;

/* Error codes */
typedef enum {
    PF_OK = 0,
    PF_ERR_INVALID_ARG = -1,
    PF_ERR_FILE_NOT_FOUND = -2,
    PF_ERR_PARSE_ERROR = -3,
    PF_ERR_OUT_OF_MEMORY = -4,
    PF_ERR_IO_ERROR = -5,
    PF_ERR_NOT_INITIALIZED = -6,
    PF_ERR_ALREADY_INITIALIZED = -7,
} pf_error;

/* Statistics structure returned after processing */
typedef struct {
    size_t total_reads;
    size_t matched_reads;
    size_t unmatched_reads;
    size_t total_barcodes;
    size_t whitelisted_barcodes;
    size_t rescued_barcodes;
    size_t total_features;
    size_t total_deduped_counts;
    size_t total_raw_counts;
    double processing_time_sec;
} pf_stats;

/* ============================================================================
 * Configuration API
 * ============================================================================ */

/**
 * Create a new configuration with default values.
 * @return New config handle, or NULL on failure.
 */
pf_config* pf_config_create(void);

/**
 * Destroy a configuration.
 * @param config Config handle to destroy.
 */
void pf_config_destroy(pf_config *config);

/**
 * Clone a configuration.
 * @param config Config to clone.
 * @return New config handle, or NULL on failure.
 */
pf_config* pf_config_clone(const pf_config *config);

/* Configuration setters */
void pf_config_set_barcode_length(pf_config *config, int length);
void pf_config_set_umi_length(pf_config *config, int length);
void pf_config_set_max_hamming_distance(pf_config *config, int distance);
void pf_config_set_stringency(pf_config *config, int stringency);
void pf_config_set_min_counts(pf_config *config, int min_counts);
void pf_config_set_min_posterior(pf_config *config, double min_posterior);
void pf_config_set_feature_offset(pf_config *config, int offset);
void pf_config_set_barcode_offset(pf_config *config, int offset);
void pf_config_set_max_barcode_mismatches(pf_config *config, int mismatches);
void pf_config_set_max_feature_n(pf_config *config, int max_n);
void pf_config_set_max_barcode_n(pf_config *config, int max_n);
void pf_config_set_threads(pf_config *config, int threads);
void pf_config_set_search_threads(pf_config *config, int threads);
void pf_config_set_consumer_threads(pf_config *config, int threads);
void pf_config_set_debug(pf_config *config, int enable);
void pf_config_set_reverse_complement_whitelist(pf_config *config, int enable);
void pf_config_set_limit_search(pf_config *config, int limit);
void pf_config_set_max_reads(pf_config *config, long long max_reads);
void pf_config_set_translate_nxt(pf_config *config, int enable);
void pf_config_set_use_feature_offset_array(pf_config *config, int enable);
void pf_config_set_use_feature_anchor_search(pf_config *config, int enable);
void pf_config_set_require_feature_anchor_match(pf_config *config, int enable);
void pf_config_set_feature_mode_bootstrap_reads(pf_config *config, int n_reads);

/* EmptyDrops control */
void pf_config_set_skip_emptydrops(pf_config *config, int enable);
void pf_config_set_emptydrops_failure_fatal(pf_config *config, int enable);
void pf_config_set_expected_cells(pf_config *config, int n_cells);
void pf_config_set_emptydrops_use_fdr(pf_config *config, int enable);

/* ============================================================================
 * Context Lifecycle API
 * ============================================================================ */

/**
 * Initialize a processing context with the given configuration.
 * @param config Configuration to use (will be cloned internally).
 * @return New context handle, or NULL on failure.
 */
pf_context* pf_init(const pf_config *config);

/**
 * Destroy a processing context and free all resources.
 * @param ctx Context to destroy.
 */
void pf_destroy(pf_context *ctx);

/**
 * Get the last error message (if any).
 * @param ctx Context handle.
 * @return Error message string, or NULL if no error.
 */
const char* pf_get_error(pf_context *ctx);

/* ============================================================================
 * Reference Loading API
 * ============================================================================ */

/**
 * Load a feature reference CSV file.
 * The CSV must have columns: id, name, sequence (or just name, sequence).
 * @param ctx Context handle.
 * @param feature_csv Path to feature reference CSV.
 * @return PF_OK on success, error code otherwise.
 */
pf_error pf_load_feature_ref(pf_context *ctx, const char *feature_csv);

/**
 * Load a barcode whitelist file (one barcode per line).
 * @param ctx Context handle.
 * @param whitelist_path Path to whitelist file.
 * @return PF_OK on success, error code otherwise.
 */
pf_error pf_load_whitelist(pf_context *ctx, const char *whitelist_path);

/**
 * Load a filtered barcodes file (subset of barcodes to process).
 * @param ctx Context handle.
 * @param filtered_path Path to filtered barcodes file.
 * @return PF_OK on success, error code otherwise.
 */
pf_error pf_load_filtered_barcodes(pf_context *ctx, const char *filtered_path);

/* ============================================================================
 * Processing API
 * ============================================================================ */

/**
 * Process FASTQ files from a directory.
 * @param ctx Context handle.
 * @param fastq_dir Directory containing FASTQ files (R1/R2 pairs).
 * @param output_dir Directory to write output files.
 * @param stats_out Optional pointer to receive processing statistics.
 * @return PF_OK on success, error code otherwise.
 */
pf_error pf_process_fastq_dir(pf_context *ctx, 
                               const char *fastq_dir,
                               const char *output_dir,
                               pf_stats *stats_out);

/**
 * Process explicit FASTQ file lists.
 * @param ctx Context handle.
 * @param barcode_fastqs Array of R1 (barcode) FASTQ paths.
 * @param feature_fastqs Array of R2 (feature) FASTQ paths.
 * @param n_files Number of file pairs.
 * @param output_dir Directory to write output files.
 * @param sample_name Sample name for output subdirectory.
 * @param stats_out Optional pointer to receive processing statistics.
 * @return PF_OK on success, error code otherwise.
 */
pf_error pf_process_fastqs(pf_context *ctx,
                            const char **barcode_fastqs,
                            const char **feature_fastqs,
                            int n_files,
                            const char *output_dir,
                            const char *sample_name,
                            pf_stats *stats_out);

/* ============================================================================
 * Output API
 * ============================================================================ */

/**
 * Get the number of features loaded.
 * @param ctx Context handle.
 * @return Number of features, or 0 if not loaded.
 */
int pf_get_num_features(pf_context *ctx);

/**
 * Get feature name by index.
 * @param ctx Context handle.
 * @param index Feature index (0-based).
 * @return Feature name, or NULL if invalid index.
 */
const char* pf_get_feature_name(pf_context *ctx, int index);

/**
 * Get feature sequence by index.
 * @param ctx Context handle.
 * @param index Feature index (0-based).
 * @return Feature sequence, or NULL if invalid index.
 */
const char* pf_get_feature_sequence(pf_context *ctx, int index);

/* ============================================================================
 * EmptyDrops Filtering API (via libscrna)
 * ============================================================================ */

/**
 * Run EmptyDrops filtering on pre-MEX data (in-memory structures).
 * This is the PIPELINE integration point - called from finalize_processing.
 * 
 * Writes:
 *   - filtered_barcodes.txt (at output_dir)
 *   - EmptyDrops/emptydrops_results.tsv (audit file)
 * 
 * @param umi_counts UMI counts per barcode (length: n_barcodes)
 * @param barcodes Barcode strings (length: n_barcodes)
 * @param n_barcodes Number of barcodes
 * @param features Feature names (length: n_features, or NULL for simple mode)
 * @param n_features Number of features
 * @param sparse_gene_ids Gene IDs for sparse entries (or NULL for simple mode)
 * @param sparse_counts Counts for sparse entries (or NULL for simple mode)
 * @param sparse_cell_index Start index for each cell (or NULL for simple mode)
 * @param n_genes_per_cell Genes per cell (or NULL for simple mode)
 * @param output_dir Directory to write outputs
 * @param n_expected_cells Expected number of cells (0 = auto)
 * @param use_fdr_gate If true, gate tail rescues by FDR instead of raw p-value
 * @param filtered_barcodes_out Output: array of filtered barcode strings (caller frees)
 * @param n_filtered_out Output: number of filtered barcodes
 * @return PF_OK on success, error code otherwise.
 */
pf_error pf_run_emptydrops_premex(
    const uint32_t *umi_counts,
    const char **barcodes,
    uint32_t n_barcodes,
    const char **features,
    uint32_t n_features,
    const uint32_t *sparse_gene_ids,
    const uint32_t *sparse_counts,
    const uint32_t *sparse_cell_index,
    const uint32_t *n_genes_per_cell,
    const char *output_dir,
    int n_expected_cells,
    int use_fdr_gate,
    char ***filtered_barcodes_out,
    uint32_t *n_filtered_out
);

/**
 * Run EmptyDrops filtering on a MEX directory (standalone tool only).
 * NOT for pipeline use - use pf_run_emptydrops_premex() instead.
 * 
 * @param mex_dir Path to MEX directory (containing matrix.mtx, barcodes.tsv, features.tsv)
 * @param output_dir Directory to write outputs
 * @param n_expected_cells Expected number of cells (0 = auto)
 * @param use_fdr_gate If true, gate tail rescues by FDR instead of raw p-value
 * @return PF_OK on success, error code otherwise.
 */
pf_error pf_run_emptydrops_mex(const char *mex_dir, const char *output_dir, int n_expected_cells, int use_fdr_gate);

/* ============================================================================
 * Utility Functions
 * ============================================================================ */

/**
 * Get the library version string.
 * @return Version string (e.g., "1.0.0").
 */
const char* pf_version(void);

/**
 * Initialize global state (called automatically by pf_init).
 * Safe to call multiple times.
 */
void pf_global_init(void);

#ifdef __cplusplus
}
#endif

#endif /* PF_API_H */
