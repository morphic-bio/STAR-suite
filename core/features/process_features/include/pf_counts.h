/**
 * @file pf_counts.h
 * @brief Internal structures for deduped counts result
 * 
 * This header defines the pf_counts_result structure used to pass
 * deduped counts between processing stages. This is an internal API
 * not exposed to external consumers of the library.
 */

#ifndef PF_COUNTS_H
#define PF_COUNTS_H

#include "common.h"
#include "khash_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Result structure containing deduped counts for all barcodes
 * 
 * This structure encapsulates the output of the deduplication stage,
 * making it explicit what data flows between assignment and output stages.
 * 
 * Memory ownership:
 * - The caller is responsible for freeing this structure via pf_counts_result_free()
 * - barcode_to_deduped_hash contains nested khash_t(u32u32)* that must be freed
 * - feature_hist array and its vec_u32_t* elements must be freed
 */
typedef struct pf_counts_result {
    /**
     * Hash table mapping barcode keys (uint32_t) to their deduped counts.
     * Key: encoded barcode (uint32_t)
     * Value: khash_t(u32u32)* mapping feature_id -> deduped_count
     */
    khash_t(u32ptr) *barcode_to_deduped_hash;
    
    /**
     * Per-feature histogram of UMI counts.
     * Array of size (n_features + 1), where index 0 is cumulative.
     * Each vec_u32_t* stores frequency at index [count].
     */
    vec_u32_t **feature_hist;
    
    /**
     * Number of features (determines feature_hist array size)
     */
    int n_features;
    
    /**
     * Total deduped counts array per feature (size: n_features + 1)
     * Index 0 unused, feature indices are 1-based.
     */
    int *total_deduped_counts;
    
    /**
     * Total raw/barcoded counts array per feature (size: n_features + 1)
     */
    int *total_barcoded_counts;
    
} pf_counts_result;

/**
 * @brief Initialize a pf_counts_result structure with zeroed fields
 * 
 * @param n_features Number of features to allocate arrays for
 * @return Newly allocated and initialized structure, or NULL on failure
 */
pf_counts_result* pf_counts_result_init(int n_features);

/**
 * @brief Build deduped counts from hashes into a pf_counts_result
 * 
 * This function extracts the deduplication logic from finalize_processing(),
 * building the barcode_to_deduped_hash from the sequence_umi_hash.
 * 
 * @param hashes Data structures containing sequence_umi_hash and filtered_hash
 * @param n_features Number of features
 * @param stringency Stringency parameter for deduplication
 * @param min_counts Minimum counts threshold
 * @return Newly allocated pf_counts_result, or NULL on failure
 */
pf_counts_result* pf_build_deduped_counts(
    data_structures *hashes,
    int n_features,
    uint16_t stringency,
    uint16_t min_counts
);

/**
 * @brief Free all resources associated with a pf_counts_result
 * 
 * Frees:
 * - All nested khash_t(u32u32)* in barcode_to_deduped_hash
 * - The barcode_to_deduped_hash itself
 * - All vec_u32_t* in feature_hist
 * - The feature_hist array
 * - The total_deduped_counts and total_barcoded_counts arrays
 * - The structure itself
 * 
 * @param result Structure to free (safe to pass NULL)
 */
void pf_counts_result_free(pf_counts_result *result);

/**
 * @brief Reset feature histogram data (for reuse between filtered/unfiltered passes)
 * 
 * Clears and frees all vec_u32_t* in feature_hist, leaving the array allocated.
 * 
 * @param result Counts result to reset histograms for
 */
void pf_counts_result_reset_histograms(pf_counts_result *result);

#ifdef __cplusplus
}
#endif

#endif /* PF_COUNTS_H */
