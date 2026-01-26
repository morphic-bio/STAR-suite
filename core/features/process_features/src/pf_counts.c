/**
 * @file pf_counts.c
 * @brief Implementation of pf_counts_result functions
 * 
 * This module provides the deduped counts building logic extracted from
 * finalize_processing(). It serves as an intermediate data structure
 * between assignment and output stages.
 */

#include "../include/pf_counts.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include <stdlib.h>
#include <string.h>

pf_counts_result* pf_counts_result_init(int n_features) {
    pf_counts_result *result = calloc(1, sizeof(pf_counts_result));
    if (!result) return NULL;
    
    result->n_features = n_features;
    
    /* Allocate barcode_to_deduped_hash */
    result->barcode_to_deduped_hash = kh_init(u32ptr);
    if (!result->barcode_to_deduped_hash) {
        free(result);
        return NULL;
    }
    
    /* Allocate feature histogram array (n_features + 1 for 1-based indexing + cumulative at 0) */
    result->feature_hist = calloc(n_features + 1, sizeof(vec_u32_t*));
    if (!result->feature_hist) {
        kh_destroy(u32ptr, result->barcode_to_deduped_hash);
        free(result);
        return NULL;
    }
    
    /* Allocate counts arrays */
    result->total_deduped_counts = calloc(n_features + 1, sizeof(int));
    result->total_barcoded_counts = calloc(n_features + 1, sizeof(int));
    if (!result->total_deduped_counts || !result->total_barcoded_counts) {
        free(result->total_deduped_counts);
        free(result->total_barcoded_counts);
        free(result->feature_hist);
        kh_destroy(u32ptr, result->barcode_to_deduped_hash);
        free(result);
        return NULL;
    }
    
    return result;
}

pf_counts_result* pf_build_deduped_counts(
    data_structures *hashes,
    int n_features,
    uint16_t stringency,
    uint16_t min_counts
) {
    if (!hashes) return NULL;
    
    pf_counts_result *result = pf_counts_result_init(n_features);
    if (!result) return NULL;
    
    /* Call the existing find_deduped_counts to populate barcode_to_deduped_hash */
    find_deduped_counts(hashes, result->barcode_to_deduped_hash, stringency, min_counts);
    
    return result;
}

void pf_counts_result_free(pf_counts_result *result) {
    if (!result) return;
    
    /* Free nested hash tables in barcode_to_deduped_hash */
    if (result->barcode_to_deduped_hash) {
        khint_t k;
        for (k = kh_begin(result->barcode_to_deduped_hash); 
             k != kh_end(result->barcode_to_deduped_hash); ++k) {
            if (kh_exist(result->barcode_to_deduped_hash, k)) {
                khash_t(u32u32) *inner = (khash_t(u32u32)*)kh_val(result->barcode_to_deduped_hash, k);
                if (inner) kh_destroy(u32u32, inner);
            }
        }
        kh_destroy(u32ptr, result->barcode_to_deduped_hash);
    }
    
    /* Free feature histograms */
    if (result->feature_hist) {
        for (int i = 0; i <= result->n_features; i++) {
            if (result->feature_hist[i]) {
                vec_u32_destroy(result->feature_hist[i]);
            }
        }
        free(result->feature_hist);
    }
    
    /* Free count arrays */
    free(result->total_deduped_counts);
    free(result->total_barcoded_counts);
    
    free(result);
}

void pf_counts_result_reset_histograms(pf_counts_result *result) {
    if (!result || !result->feature_hist) return;
    
    for (int i = 0; i <= result->n_features; i++) {
        if (result->feature_hist[i]) {
            vec_u32_destroy(result->feature_hist[i]);
            result->feature_hist[i] = NULL;
        }
    }
    
    /* Also reset the counts arrays */
    if (result->total_deduped_counts) {
        memset(result->total_deduped_counts, 0, (result->n_features + 1) * sizeof(int));
    }
    if (result->total_barcoded_counts) {
        memset(result->total_barcoded_counts, 0, (result->n_features + 1) * sizeof(int));
    }
}
