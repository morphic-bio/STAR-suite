/**
 * @file barcode_filter.c
 * @brief Barcode filtering stage implementation
 */

#include "../include/barcode_filter.h"
#include "../include/emptydrops_bridge.h"
#include <stdlib.h>
#include <stdio.h>

/* ============================================================================
 * Main Filter Function
 * ============================================================================ */

pf_filter_status pf_filter_barcodes(
    pf_counts_result *counts,
    barcode_filter_config *config,
    khash_t(strptr) *external_filter,
    khash_t(strptr) **filter_out
) {
    if (!counts || !config || !filter_out) {
        fprintf(stderr, "[barcode_filter] Invalid arguments\n");
        if (filter_out) *filter_out = NULL;
        return PF_FILTER_FAILED;
    }
    
    /* If external filter is provided, use it directly */
    if (external_filter) {
        fprintf(stderr, "[barcode_filter] Using externally provided filter (%u barcodes)\n",
                (unsigned int)kh_size(external_filter));
        *filter_out = external_filter;  /* Borrowed reference */
        return PF_FILTER_OK;
    }
    
    /* Skip EmptyDrops if requested */
    if (config->skip_emptydrops) {
        fprintf(stderr, "[barcode_filter] EmptyDrops skipped (--skip-empty-drops)\n");
        *filter_out = NULL;
        return PF_FILTER_SKIPPED;
    }
    
    /* No external filter - run EmptyDrops */
    fprintf(stderr, "[barcode_filter] No external filter provided, running EmptyDrops...\n");
    
    khash_t(strptr) *ed_filter = run_emptydrops_on_counts(
        counts,
        config->feature_names,
        config->n_features,
        config->output_dir,
        config->n_expected_cells,
        config->translate_nxt
    );
    
    if (!ed_filter) {
        fprintf(stderr, "[barcode_filter] EmptyDrops filtering failed\n");
        if (!config->emptydrops_failure_fatal) {
            fprintf(stderr, "[barcode_filter] To provide explicit filtering, use --filtered_barcodes <file>\n");
        }
        *filter_out = NULL;
        return PF_FILTER_FAILED;
    }
    
    fprintf(stderr, "[barcode_filter] EmptyDrops filter created with %u barcodes\n",
            (unsigned int)kh_size(ed_filter));
    
    *filter_out = ed_filter;  /* Caller owns this */
    return PF_FILTER_OK;
}

/* ============================================================================
 * Filter Memory Management
 * ============================================================================ */

void pf_filter_free(khash_t(strptr) *filter) {
    if (!filter) return;
    
    /* Free all keys (barcode strings) */
    khint_t k;
    for (k = kh_begin(filter); k != kh_end(filter); ++k) {
        if (kh_exist(filter, k)) {
            char *key = (char*)kh_key(filter, k);
            free(key);
        }
    }
    
    kh_destroy(strptr, filter);
}
