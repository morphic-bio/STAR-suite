/**
 * @file barcode_filter.h
 * @brief Barcode filtering stage for assignBarcodes pipeline
 * 
 * This module provides the filtering stage, which can use either:
 * 1. An externally provided barcode list (--filtered_barcodes)
 * 2. EmptyDrops cell calling (automatic when no list provided)
 */

#ifndef BARCODE_FILTER_H
#define BARCODE_FILTER_H

#include "common.h"
#include "pf_counts.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Filter status codes returned by pf_filter_barcodes()
 */
typedef enum {
    PF_FILTER_OK = 0,       /* Filter obtained successfully */
    PF_FILTER_SKIPPED = 1,  /* Filtering skipped (--skip-empty-drops) */
    PF_FILTER_FAILED = 2    /* EmptyDrops failed */
} pf_filter_status;

/**
 * @brief Filter configuration options
 */
typedef struct barcode_filter_config {
    const char *output_dir;         /* Output directory for EmptyDrops files */
    char **feature_names;           /* Feature names array (borrowed) */
    int n_features;                 /* Number of features */
    int n_expected_cells;           /* Expected cells for EmptyDrops (0 = auto) */
    int translate_nxt;              /* If true, apply NXT translation to barcodes */
    
    /* EmptyDrops control */
    int skip_emptydrops;            /* 1 = skip EmptyDrops entirely */
    int emptydrops_failure_fatal;   /* 1 = treat ED failure as error (for logging) */
    int emptydrops_use_fdr;         /* 1 = use FDR gate for tail rescue */
} barcode_filter_config;

/**
 * @brief Run barcode filtering on deduped counts
 * 
 * If external_filter is provided, it is returned via filter_out (no copying).
 * If skip_emptydrops is set, returns SKIPPED with filter_out=NULL.
 * Otherwise, EmptyDrops is run and a new filter hash is created.
 * 
 * @param counts Deduped counts from pf_build_deduped_counts()
 * @param config Filter configuration
 * @param external_filter Externally provided filter list (NULL to use EmptyDrops)
 * @param filter_out Output: filter hash (caller owns if newly created, borrowed if external)
 * @return Status code indicating success, skip, or failure
 */
pf_filter_status pf_filter_barcodes(
    pf_counts_result *counts,
    barcode_filter_config *config,
    khash_t(strptr) *external_filter,
    khash_t(strptr) **filter_out
);

/**
 * @brief Check if a filter hash was newly allocated (vs externally provided)
 * 
 * This helper is useful for determining whether to free the filter hash.
 * 
 * @param filter Filter hash to check
 * @param external_filter The external filter that was passed to pf_filter_barcodes()
 * @return 1 if filter was newly allocated, 0 if it's the external filter
 */
static inline int pf_filter_is_owned(khash_t(strptr) *filter, khash_t(strptr) *external_filter) {
    return filter != NULL && filter != external_filter;
}

/**
 * @brief Free a filter hash and its contents
 * 
 * Frees all keys (barcode strings) and the hash itself.
 * Safe to call with NULL.
 * 
 * @param filter Filter hash to free
 */
void pf_filter_free(khash_t(strptr) *filter);

#ifdef __cplusplus
}
#endif

#endif /* BARCODE_FILTER_H */
