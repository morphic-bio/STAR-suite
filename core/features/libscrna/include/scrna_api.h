/**
 * @file scrna_api.h
 * @brief C ABI wrapper for libscrna EmptyDrops functionality
 * 
 * This header provides a C-compatible interface to the EmptyDrops/OrdMag
 * filtering functionality, suitable for use from C code (e.g., process_features).
 */

#ifndef H_scrna_api
#define H_scrna_api

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// ============================================================================
// EmptyDrops Configuration
// ============================================================================

typedef struct {
    // Simple EmptyDrops (OrdMag) parameters
    uint32_t n_expected_cells;       // Expected number of cells (default: 3000)
    double max_percentile;           // Max percentile for robust max (default: 0.99)
    double max_min_ratio;            // Max/min ratio (default: 10.0)
    uint32_t umi_min;                // Minimum UMI threshold (default: 500)
    double umi_min_frac_median;      // Min UMI as fraction of median (default: 0.01)
    uint32_t cand_max_n;             // Maximum candidates (default: 20000)
    uint32_t ind_min;                // Min index for ambient cells (default: 45000)
    uint32_t ind_max;                // Max index for ambient cells (default: 90000)
    
    // EmptyDrops Monte Carlo parameters
    uint32_t sim_n;                  // Number of MC simulations (default: 10000)
    double fdr;                      // FDR threshold (default: 0.01)
    double raw_pvalue_threshold;     // Raw p-value threshold (default: 0.05)
    uint64_t seed;                   // Random seed (default: 1)
    uint32_t lower_testing_bound;    // Lower UMI bound for testing (default: 500)
    uint32_t ambient_umi_max;        // Max UMI for ambient cells (default: 100)
    uint32_t mc_threads;             // Threads for Monte Carlo (0 = single-threaded)
    
    // Occupancy filter (Flex-only, disabled in no-tag/compat mode)
    int disable_occupancy_filter;    // If true, skip occupancy filter (default: 1 for compat mode)
    
    // ED retain window
    uint32_t ed_retain_count;        // Optional cap on retained cells (0 = no cap)
    
    // Gating mode
    int use_fdr_gate;                // If true, use FDR threshold; if false (default), use raw p-value
    
} scrna_ed_config;

// ============================================================================
// EmptyDrops Results
// ============================================================================

// Per-candidate result (Flex-compatible schema)
typedef struct {
    char *barcode;                   // Cell barcode string
    uint32_t cell_index;             // Index in original matrix
    uint32_t umi_count;              // UMI count for this cell
    double p_value;                  // Monte Carlo p-value
    double p_adjusted;               // FDR-adjusted p-value
    int passes_raw_p;                // 1 if passes raw p-value threshold
    int passes_fdr;                  // 1 if passes FDR threshold
    double obs_log_prob;             // Observed log probability (diagnostic)
    int is_simple_cell;              // 1 if cell is a simple cell (auto-pass)
} scrna_ed_candidate;

// Full EmptyDrops result
typedef struct {
    // Filtered barcodes (passing cells)
    char **barcodes;                 // Array of barcode strings
    size_t n_barcodes;               // Number of passing barcodes
    
    // Candidate results (for detailed output)
    scrna_ed_candidate *candidates;  // Array of candidate results
    size_t n_candidates;             // Number of candidates tested
    
    // Summary statistics
    uint32_t n_simple_cells;         // Cells passing simple filter (high UMI)
    uint32_t n_tail_cells;           // Cells tested by ED
    uint32_t n_ed_passers;           // Cells rescued by ED
    uint32_t retain_threshold;       // UMI threshold for simple cells
    uint32_t min_umi;                // Minimum UMI for candidates
    
    // Error handling
    char *error_message;             // Error message (NULL if success)
} scrna_ed_result;

// ============================================================================
// Input Matrix
// ============================================================================

// Simple matrix input (feature counts per cell)
// For single-sample/no-tag mode, this is sufficient
typedef struct {
    uint32_t *umi_counts;            // UMI counts per cell (length: n_cells)
    char **barcodes;                 // Cell barcode strings (length: n_cells)
    char **features;                 // Feature/gene names (length: n_features)
    uint32_t n_cells;                // Number of cells
    uint32_t n_features;             // Number of features
    
    // Sparse matrix data (for ED p-value computation)
    // If NULL, only simple filtering is performed
    uint32_t *sparse_gene_ids;       // Gene IDs for sparse entries
    uint32_t *sparse_counts;         // Counts for sparse entries
    uint32_t *sparse_cell_index;     // Start index for each cell in sparse arrays
    uint32_t *n_genes_per_cell;      // Number of genes per cell
    size_t sparse_nnz;               // Number of non-zero entries
} scrna_matrix_input;

// ============================================================================
// API Functions
// ============================================================================

/**
 * Create EmptyDrops config with default values (matching CR defaults).
 * Caller must free with scrna_ed_config_destroy().
 */
scrna_ed_config* scrna_ed_config_create(void);

/**
 * Destroy EmptyDrops config.
 */
void scrna_ed_config_destroy(scrna_ed_config *config);

/**
 * Run EmptyDrops filtering on a single-sample matrix.
 * 
 * This is the main entry point for no-tag/single-sample mode.
 * If sparse matrix data is NULL, only simple (OrdMag) filtering is performed.
 * 
 * @param input Input matrix data
 * @param config EmptyDrops configuration
 * @param result Output: results (caller must free with scrna_ed_result_free)
 * @return 0 on success, non-zero on error
 */
int scrna_emptydrops_run(
    const scrna_matrix_input *input,
    const scrna_ed_config *config,
    scrna_ed_result *result
);

/**
 * Free EmptyDrops result.
 */
void scrna_ed_result_free(scrna_ed_result *result);

/**
 * Write EmptyDrops outputs to directory.
 * 
 * Outputs:
 *   - filtered_barcodes.txt: passing barcode list
 *   - EmptyDrops/emptydrops_results.tsv: detailed results
 * 
 * @param result EmptyDrops result
 * @param output_dir Output directory path
 * @return 0 on success, non-zero on error
 */
int scrna_emptydrops_write_outputs(
    const scrna_ed_result *result,
    const char *output_dir
);

/**
 * Write filtered barcodes to file.
 * 
 * @param result EmptyDrops result
 * @param filepath Output file path
 * @return 0 on success, non-zero on error
 */
int scrna_write_filtered_barcodes(
    const scrna_ed_result *result,
    const char *filepath
);

#ifdef __cplusplus
}
#endif

#endif // H_scrna_api
