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

#include "pf_split_read.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handles */
typedef struct pf_config pf_config;
typedef struct pf_context pf_context;
typedef struct pf_record_stream pf_record_stream;
typedef struct pf_direct_range_job pf_direct_range_job;

typedef struct {
    const char *data;
    size_t length;
} pf_sequence_view;

/* In-memory read record for non-FASTQ input providers. */
typedef struct {
    const char *barcode_sequence;
    const char *barcode_quality;
    const char *feature_sequence;
    const char *feature_quality;
    const char *feature_sequence2;
    const char *feature_quality2;
} pf_read_record;

/* Borrowed-span read record for native non-FASTQ providers. */
typedef struct {
    pf_sequence_view barcode_sequence;
    pf_sequence_view barcode_quality;
    pf_sequence_view feature_sequence;
    pf_sequence_view feature_quality;
    pf_sequence_view feature_sequence2;
    pf_sequence_view feature_quality2;
} pf_read_record_view;

/* Optional permit hook API for external schedulers */
typedef uint64_t (*pf_permit_acquire_fn)(void *hook_ctx);
typedef void (*pf_permit_release_fn)(
    void *hook_ctx,
    uint64_t wait_ns,
    uint64_t work_units,
    uint64_t work_bytes,
    uint64_t work_ns
);

/* ============================================================================
 * Barcode Namespace Types
 * ============================================================================ */

/**
 * Barcode namespace identifiers.
 * NXT and TRU differ by complementation at 0-based positions 7 and 8.
 * UNION indicates an explicitly opted-in mixed-namespace set (both forms
 * present in the hash via expand_hash_union_namespace).
 */
typedef enum {
    PF_NS_UNKNOWN = 0,
    PF_NS_NXT     = 1,
    PF_NS_TRU     = 2,
    PF_NS_UNION   = 3,
} pf_namespace_t;

const char* pf_namespace_to_string(pf_namespace_t ns);
pf_namespace_t pf_namespace_from_string(const char *s);

/**
 * Translate a barcode between NXT and TRU namespaces in place.
 * Complements 0-based positions 7 and 8.  Idempotent: applying twice
 * yields the original.  No-op if len < 9 or barcode is NULL.
 */
void pf_translate_barcode_inplace(char *barcode, size_t len);

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
    PF_ERR_OFFSET_CONFLICT = -8,       /* Both global offset and per-feature offsets specified */
    PF_ERR_MULTI_OFFSET_DETECTED = -9, /* Multiple offsets detected, need explicit choice */
    PF_ERR_ALLOC = -10,                /* Memory allocation failure */
    PF_ERR_NAMESPACE = -11,            /* Namespace unresolved or ambiguous */
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
void pf_config_set_legacy_cb_rescue(pf_config *config, int enable);
void pf_config_set_feature_offset(pf_config *config, int offset);
void pf_config_set_barcode_offset(pf_config *config, int offset);
void pf_config_set_max_barcode_mismatches(pf_config *config, int mismatches);
void pf_config_set_max_feature_n(pf_config *config, int max_n);
void pf_config_set_max_barcode_n(pf_config *config, int max_n);
void pf_config_set_threads(pf_config *config, int threads);
void pf_config_set_search_threads(pf_config *config, int threads);
void pf_config_set_consumer_threads(pf_config *config, int threads);
void pf_config_set_read_buffer_lines(pf_config *config, int lines);
void pf_config_set_permit_hooks(
    pf_config *config,
    pf_permit_acquire_fn acquire_cb,
    pf_permit_release_fn release_cb,
    void *hook_ctx
);
void pf_config_set_debug(pf_config *config, int enable);
void pf_config_set_reverse_complement_whitelist(pf_config *config, int enable);
void pf_config_set_limit_search(pf_config *config, int limit);
/* 0=in_window_full, 1=in_window_simple (limited-search branch only; search is strictly in-window) */
void pf_config_set_feature_limited_mode(pf_config *config, int mode);
void pf_config_set_feature_limited_fallback(pf_config *config, int mode); /* deprecated alias */
void pf_config_set_max_reads(pf_config *config, long long max_reads);
void pf_config_set_translate_nxt(pf_config *config, int enable);
void pf_config_set_use_feature_offset_array(pf_config *config, int enable);
void pf_config_set_strict_offset_check(pf_config *config, int enable);
void pf_config_set_use_feature_anchor_search(pf_config *config, int enable);
void pf_config_set_require_feature_anchor_match(pf_config *config, int enable);
void pf_config_set_feature_mode_bootstrap_reads(pf_config *config, int n_reads);
void pf_config_set_use_hot_hash(pf_config *config, int enable);
void pf_config_set_skip_heatmaps(pf_config *config, int enable);

void pf_config_set_split_read_layout(pf_config *config, const pf_split_read_layout *layout);
void pf_config_set_split_read_fastq_patterns(pf_config *config,
                                             const char *r1_pattern,
                                             const char *r2_pattern,
                                             const char *r3_pattern);
void pf_config_clear_split_read_layout(pf_config *config);
const pf_split_read_layout *pf_config_get_split_read_layout(const pf_config *config);

/* Prehash memory budget (0 = auto-detect from system memory) */
void pf_config_set_prehash_memory_budget(pf_config *config, unsigned long long budget);

/* EmptyDrops control */
void pf_config_set_skip_emptydrops(pf_config *config, int enable);
void pf_config_set_emptydrops_failure_fatal(pf_config *config, int enable);
void pf_config_set_expected_cells(pf_config *config, int n_cells);
void pf_config_set_emptydrops_use_fdr(pf_config *config, int enable);

/* NXT/TRU auto-detection */
void pf_config_set_autodetect_chemistry(pf_config *config, int enabled);
void pf_config_set_autodetect_chemistry_reads(pf_config *config, int n_reads);
void pf_config_set_autodetect_chemistry_min_hits(pf_config *config, int min_hits);
void pf_config_set_probe_only(pf_config *config, int enabled);
void pf_config_set_skip_qc_outputs(pf_config *config, int enabled);

/* ADT / protein MEX output (assignBarcodes --output-mode adt_mex) */
void pf_config_set_adt_mex_output(pf_config *config, int enable);

/* Hash / HTO / CMO demux (adt_mex extension) */
void pf_config_set_hash_demux_mode(pf_config *config, int mode);
void pf_config_set_hash_feature_selector(pf_config *config, const char *selector);
void pf_config_set_hash_demux_method(pf_config *config, const char *method);
void pf_config_set_library_feature_type(pf_config *config, const char *feature_type);
void pf_config_set_hash_min_total(pf_config *config, int min_total);
void pf_config_set_hash_min_top(pf_config *config, int min_top);
void pf_config_set_hash_min_ratio(pf_config *config, double min_ratio);

/**
 * Allow union (mixed NXT+TRU) whitelists and filtered barcode sets.
 * When enabled, filtered barcodes are expanded at ingress: each barcode's
 * NXT/TRU translation is also inserted into the hash so both namespace
 * forms are available for downstream lookup.  Emits a warning to stderr.
 *
 * filtered_barcode_hash_contains() uses exact-only matching (no runtime
 * translated fallback).  This expansion ensures that union-whitelist
 * workflows (e.g. raw 3M-february-2018.txt containing both NXT and TRU
 * forms) continue to work with exact-only matching.  Without this flag,
 * callers must ensure filtered barcodes are in the same namespace as the
 * assignment output, either by pre-normalizing the file or by setting
 * source_namespace and target_namespace so pf_load_filtered_barcodes()
 * normalizes at ingress.  Using --allow_union_whitelist is the escape
 * hatch for legacy workflows that do
 * not regress union-whitelist workflows.
 *
 * Default: 0 (no expansion at ingress).
 */
void pf_config_set_allow_union_whitelist(pf_config *config, int enable);

/**
 * Set the source (filtered barcode file) and target (assignment output)
 * namespaces.  When both are known single-namespace types (NXT or TRU)
 * and differ, pf_load_filtered_barcodes() normalizes all barcodes to the
 * target namespace at ingress.  Default: PF_NS_UNKNOWN (no normalization).
 */
void pf_config_set_source_namespace(pf_config *config, pf_namespace_t ns);
void pf_config_set_target_namespace(pf_config *config, pf_namespace_t ns);

/**
 * Get detected barcode match mode after processing.
 * Returns "RAW_MATCH", "TRANSLATED_MATCH", "AMBIGUOUS", or "UNKNOWN".
 * Only meaningful after processing completes with autodetect_chemistry enabled.
 */
const char* pf_get_detected_match_mode(pf_context *ctx);

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

pf_error pf_process_split_fastq_dir(pf_context *ctx,
                                    const char *fastq_dir,
                                    const char *output_dir,
                                    pf_stats *stats_out,
                                    pf_split_read_metrics *metrics_out);

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

/**
 * Process in-memory barcode + feature records.
 * This is the adapter surface for native non-FASTQ readers. Records are
 * expected to provide barcode_sequence as CB+UMI and feature_sequence as the
 * feature/protospacer read. feature_sequence2 is optional for dual-orientation
 * feature libraries. All strings must be NUL-terminated for this initial API.
 *
 * @param ctx Context handle.
 * @param records Array of input records.
 * @param n_records Number of records in the array.
 * @param output_dir Directory to write output files.
 * @param sample_name Sample name for output subdirectory.
 * @param stats_out Optional pointer to receive processing statistics.
 * @return PF_OK on success, error code otherwise.
 */
pf_error pf_process_records(pf_context *ctx,
                            const pf_read_record *records,
                            size_t n_records,
                            const char *output_dir,
                            const char *sample_name,
                            pf_stats *stats_out);

/**
 * Begin streaming in-memory barcode + feature records for one sample.
 * The returned stream owns process_features sample state until
 * pf_process_records_end() or pf_process_records_abort() is called.
 */
pf_error pf_process_records_begin(pf_context *ctx,
                                  const char *output_dir,
                                  const char *sample_name,
                                  pf_record_stream **stream_out);

/**
 * Process a batch of NUL-terminated in-memory records through an open stream.
 */
pf_error pf_process_record_batch(pf_record_stream *stream,
                                 const pf_read_record *records,
                                 size_t n_records);

/**
 * Process a batch of borrowed-span records through an open stream.
 * Required sequence spans must have non-NULL data. Missing quality spans are
 * represented by data == NULL and are replaced with default qualities.
 */
pf_error pf_process_record_views(pf_record_stream *stream,
                                 const pf_read_record_view *records,
                                 size_t n_records);

/**
 * Begin direct worker-owned processing for externally range-partitioned input.
 * Each caller thread must use a stable worker_id in [0, nworkers).
 */
pf_error pf_direct_range_begin(pf_context *ctx,
                               const char *output_dir,
                               const char *sample_name,
                               int nworkers,
                               int nreaders,
                               pf_direct_range_job **job_out);

/**
 * Process borrowed-span records on one direct worker. Calls are thread-safe
 * only when different threads use different worker_id values.
 */
pf_error pf_direct_range_process_record_views(pf_direct_range_job *job,
                                              int worker_id,
                                              const pf_read_record_view *records,
                                              size_t n_records);

/**
 * Finish direct worker processing, merge worker-local counts, write outputs,
 * destroy job state, and release the process_features runtime lock.
 */
pf_error pf_direct_range_end(pf_direct_range_job *job,
                             pf_stats *stats_out);

/**
 * Abort direct worker processing without final outputs and release the
 * process_features runtime lock.
 */
void pf_direct_range_abort(pf_direct_range_job *job);

/**
 * Finish a streaming in-memory sample, write outputs, destroy stream state,
 * and release the process_features runtime lock.
 */
pf_error pf_process_records_end(pf_record_stream *stream,
                                pf_stats *stats_out);

/**
 * Abort a streaming in-memory sample without writing final outputs and release
 * the process_features runtime lock.
 */
void pf_process_records_abort(pf_record_stream *stream);

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

/**
 * Query per-feature prehash uniqueness (no_ambiguity) flags.
 * @param ctx Context handle.
 * @param level 1 for le1, 2 for le2.
 * @param index Feature index (0-based).
 * @return 1 if unique at that level, 0 if ambiguous, -1 on error/not built.
 */
int pf_get_feature_no_ambiguity(pf_context *ctx, int level, int index);

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
