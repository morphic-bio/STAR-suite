/**
 * @file emptydrops_bridge.h
 * @brief Bridge between pf_counts_result and EmptyDrops API
 * 
 * This module converts the in-memory hash-based counts representation
 * to the array-based format required by the EmptyDrops C API.
 */

#ifndef EMPTYDROPS_BRIDGE_H
#define EMPTYDROPS_BRIDGE_H

#include "common.h"
#include "pf_counts.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Arrays prepared for EmptyDrops input
 * 
 * This structure holds the converted array data needed by pf_run_emptydrops_premex().
 * All arrays are owned by this structure and freed by emptydrops_input_free().
 */
typedef struct emptydrops_input {
    char **barcodes;              /* Barcode strings (length: n_barcodes) */
    uint32_t *umi_counts;         /* UMI counts per barcode (length: n_barcodes) */
    uint32_t n_barcodes;          /* Number of barcodes */
    
    char **features;              /* Feature names (length: n_features) */
    uint32_t n_features;          /* Number of features */
    
    /* Sparse matrix data (CSC format) */
    uint32_t *sparse_gene_ids;    /* Gene/feature indices (0-based) */
    uint32_t *sparse_counts;      /* Counts for each sparse entry */
    uint32_t *sparse_cell_index;  /* Start index for each cell (length: n_barcodes + 1) */
    uint32_t *n_genes_per_cell;   /* Number of genes per cell */
    size_t sparse_nnz;            /* Total non-zero entries */
} emptydrops_input;

/**
 * @brief Convert pf_counts_result to EmptyDrops input arrays
 * 
 * Converts the hash-based barcode_to_deduped_hash into the array format
 * required by pf_run_emptydrops_premex(). The returned structure owns all
 * allocated memory.
 * 
 * @param counts Deduped counts from pf_build_deduped_counts()
 * @param feature_names Feature names array (borrowed, not copied)
 * @param translate_nxt If true, apply translate_nxt_inplace() to barcode strings
 * @return Newly allocated emptydrops_input, or NULL on error. Caller must free.
 */
emptydrops_input* emptydrops_input_from_counts(
    pf_counts_result *counts,
    char **feature_names,
    int translate_nxt
);

/**
 * @brief Free EmptyDrops input structure
 * 
 * @param input Structure to free (safe to pass NULL)
 */
void emptydrops_input_free(emptydrops_input *input);

/**
 * @brief Run EmptyDrops on deduped hash and return filtered barcodes hash
 * 
 * This is the main entry point for EmptyDrops filtering. It:
 * 1. Converts barcode_to_deduped_hash to array format
 * 2. Calls pf_run_emptydrops_premex()
 * 3. Writes filtered_barcodes.txt and EmptyDrops/emptydrops_results.tsv
 * 4. Converts the result to a khash_t(strptr)* for use by mex_writer
 * 
 * @param counts Deduped counts from pf_build_deduped_counts()
 * @param feature_names Feature names array
 * @param n_features Number of features
 * @param output_dir Directory to write EmptyDrops outputs
 * @param n_expected_cells Expected cells (0 = auto)
 * @param use_fdr_gate If true, gate tail rescues by FDR instead of raw p-value
 * @param translate_nxt If true, apply NXT translation to barcodes
 * @return New hash of filtered barcodes (caller owns), or NULL on error
 */
khash_t(strptr)* run_emptydrops_on_counts(
    pf_counts_result *counts,
    char **feature_names,
    int n_features,
    const char *output_dir,
    int n_expected_cells,
    int use_fdr_gate,
    int translate_nxt
);

#ifdef __cplusplus
}
#endif

#endif /* EMPTYDROPS_BRIDGE_H */
