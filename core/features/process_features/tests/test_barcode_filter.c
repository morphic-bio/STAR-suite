/**
 * @file test_barcode_filter.c
 * @brief Unit tests for barcode_filter and emptydrops_bridge modules
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../include/barcode_filter.h"
#include "../include/emptydrops_bridge.h"
#include "../include/pf_counts.h"
#include "../include/globals.h"
#include "../include/prototypes.h"

/* Test state */
static int tests_run = 0;
static int tests_passed = 0;
static char test_dir[256];

/* ============================================================================
 * Helper functions
 * ============================================================================ */

static void create_test_dir(void) {
    snprintf(test_dir, sizeof(test_dir), "/tmp/test_barcode_filter_%d", getpid());
    mkdir(test_dir, 0755);
}

static void cleanup_test_dir(void) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "rm -rf %s", test_dir);
    int rc = system(cmd);
    (void)rc;
}

/* Initialize globals needed for barcode encoding/decoding */
static void init_globals(void) {
    barcode_length = 16;
    barcode_code_length = 4;
    umi_length = 12;
    umi_code_length = 3;
    number_of_features = 5;
    translate_NXT = 0;
    
    initseq2Code();
    initcode2seq();
}

/* Create a pf_counts_result with test data */
static pf_counts_result* create_test_counts(int n_features, int n_barcodes) {
    pf_counts_result *counts = pf_counts_result_init(n_features);
    if (!counts) return NULL;
    
    /* Add some test barcodes with counts */
    for (int i = 0; i < n_barcodes; i++) {
        /* Create a test barcode key */
        uint32_t barcode_key = 0x10000000 + i * 0x01000000;
        
        /* Create deduped hash for this barcode */
        khash_t(u32u32) *deduped = kh_init(u32u32);
        
        /* Add some feature counts */
        int ret;
        for (int f = 1; f <= n_features && f <= (i % n_features) + 1; f++) {
            khint_t k = kh_put(u32u32, deduped, f, &ret);
            kh_val(deduped, k) = (i + 1) * 10 + f;  /* Count varies by barcode and feature */
        }
        
        /* Add to barcode_to_deduped_hash */
        khint_t k = kh_put(u32ptr, counts->barcode_to_deduped_hash, barcode_key, &ret);
        kh_val(counts->barcode_to_deduped_hash, k) = deduped;
    }
    
    return counts;
}

/* Create test feature names */
static char** create_test_feature_names(int n_features) {
    char **names = calloc(n_features, sizeof(char*));
    if (!names) return NULL;
    
    for (int i = 0; i < n_features; i++) {
        names[i] = malloc(32);
        snprintf(names[i], 32, "Feature_%d", i + 1);
    }
    
    return names;
}

static void free_test_feature_names(char **names, int n_features) {
    if (names) {
        for (int i = 0; i < n_features; i++) {
            free(names[i]);
        }
        free(names);
    }
}

/* ============================================================================
 * Test: emptydrops_input_from_counts basic functionality
 * ============================================================================ */
static int test_emptydrops_input_from_counts(void) {
    int n_features = 3;
    int n_barcodes = 5;
    
    pf_counts_result *counts = create_test_counts(n_features, n_barcodes);
    if (!counts) return 0;
    
    char **feature_names = create_test_feature_names(n_features);
    if (!feature_names) {
        pf_counts_result_free(counts);
        return 0;
    }
    
    /* Convert to EmptyDrops input */
    emptydrops_input *input = emptydrops_input_from_counts(counts, feature_names, 0);
    
    int success = 1;
    
    if (!input) {
        printf("    emptydrops_input_from_counts returned NULL\n");
        success = 0;
    } else {
        /* Check basic properties */
        if (input->n_barcodes != (uint32_t)n_barcodes) {
            printf("    Expected %d barcodes, got %u\n", n_barcodes, input->n_barcodes);
            success = 0;
        }
        
        if (input->n_features != (uint32_t)n_features) {
            printf("    Expected %d features, got %u\n", n_features, input->n_features);
            success = 0;
        }
        
        /* Check that barcodes were decoded */
        if (!input->barcodes || !input->barcodes[0]) {
            printf("    Barcodes not properly decoded\n");
            success = 0;
        }
        
        /* Check sparse arrays */
        if (!input->sparse_gene_ids || !input->sparse_counts || !input->sparse_cell_index) {
            printf("    Sparse arrays not allocated\n");
            success = 0;
        }
        
        /* Verify sparse_cell_index is properly constructed */
        if (input->sparse_cell_index[n_barcodes] != input->sparse_nnz) {
            printf("    sparse_cell_index[n_barcodes] != sparse_nnz\n");
            success = 0;
        }
        
        emptydrops_input_free(input);
    }
    
    free_test_feature_names(feature_names, n_features);
    pf_counts_result_free(counts);
    
    return success;
}

/* ============================================================================
 * Test: emptydrops_input_from_counts with empty input
 * ============================================================================ */
static int test_emptydrops_input_empty(void) {
    int n_features = 3;
    
    /* Create counts with no barcodes */
    pf_counts_result *counts = pf_counts_result_init(n_features);
    if (!counts) return 0;
    
    char **feature_names = create_test_feature_names(n_features);
    if (!feature_names) {
        pf_counts_result_free(counts);
        return 0;
    }
    
    /* Convert to EmptyDrops input - should return NULL for empty input */
    emptydrops_input *input = emptydrops_input_from_counts(counts, feature_names, 0);
    
    int success = (input == NULL);  /* Expected to fail with 0 barcodes */
    
    if (input) {
        printf("    Expected NULL for empty input, got non-NULL\n");
        emptydrops_input_free(input);
    }
    
    free_test_feature_names(feature_names, n_features);
    pf_counts_result_free(counts);
    
    return success;
}

/* ============================================================================
 * Test: pf_filter_barcodes with external filter
 * ============================================================================ */
static int test_filter_with_external(void) {
    int n_features = 3;
    int n_barcodes = 5;
    
    pf_counts_result *counts = create_test_counts(n_features, n_barcodes);
    if (!counts) return 0;
    
    char **feature_names = create_test_feature_names(n_features);
    if (!feature_names) {
        pf_counts_result_free(counts);
        return 0;
    }
    
    /* Create external filter */
    khash_t(strptr) *external_filter = kh_init(strptr);
    int ret;
    khint_t k = kh_put(strptr, external_filter, strdup("AAACCCAAGAAACCAT"), &ret);
    kh_val(external_filter, k) = NULL;
    k = kh_put(strptr, external_filter, strdup("AAACCCAAGAAACCCA"), &ret);
    kh_val(external_filter, k) = NULL;
    
    barcode_filter_config config = {
        .output_dir = test_dir,
        .feature_names = feature_names,
        .n_features = n_features,
        .n_expected_cells = 0,
        .translate_nxt = 0,
        .skip_emptydrops = 0,
        .emptydrops_failure_fatal = 0
    };
    
    /* Call filter with external filter - should return OK with the external filter */
    khash_t(strptr) *result = NULL;
    pf_filter_status status = pf_filter_barcodes(counts, &config, external_filter, &result);
    
    int success = 1;
    
    if (status != PF_FILTER_OK) {
        printf("    Expected PF_FILTER_OK, got %d\n", status);
        success = 0;
    }
    
    if (result != external_filter) {
        printf("    Expected external filter to be returned\n");
        success = 0;
    }
    
    /* Check ownership helper */
    if (pf_filter_is_owned(result, external_filter)) {
        printf("    pf_filter_is_owned should return 0 for external filter\n");
        success = 0;
    }
    
    /* Clean up external filter (we own it) */
    for (k = kh_begin(external_filter); k != kh_end(external_filter); ++k) {
        if (kh_exist(external_filter, k)) {
            free((char*)kh_key(external_filter, k));
        }
    }
    kh_destroy(strptr, external_filter);
    
    free_test_feature_names(feature_names, n_features);
    pf_counts_result_free(counts);
    
    return success;
}

/* ============================================================================
 * Test: pf_filter_free with NULL
 * ============================================================================ */
static int test_filter_free_null(void) {
    /* Should not crash */
    pf_filter_free(NULL);
    return 1;
}

/* ============================================================================
 * Test: EmptyDrops integration (full path without external filter)
 * ============================================================================ */
static int test_emptydrops_integration(void) {
    int n_features = 3;
    int n_barcodes = 10;  /* More barcodes for EmptyDrops to work with */
    
    pf_counts_result *counts = create_test_counts(n_features, n_barcodes);
    if (!counts) return 0;
    
    char **feature_names = create_test_feature_names(n_features);
    if (!feature_names) {
        pf_counts_result_free(counts);
        return 0;
    }
    
    barcode_filter_config config = {
        .output_dir = test_dir,
        .feature_names = feature_names,
        .n_features = n_features,
        .n_expected_cells = 0,  /* Auto-detect */
        .translate_nxt = 0,
        .skip_emptydrops = 0,
        .emptydrops_failure_fatal = 0
    };
    
    /* Call filter WITHOUT external filter - should run EmptyDrops */
    khash_t(strptr) *result = NULL;
    pf_filter_status status = pf_filter_barcodes(counts, &config, NULL, &result);
    
    int success = 1;
    
    /* EmptyDrops may succeed or fail depending on data, but shouldn't crash. */
    if (status == PF_FILTER_OK && result != NULL) {
        /* Verify ownership helper says we own it */
        if (!pf_filter_is_owned(result, NULL)) {
            printf("    pf_filter_is_owned should return 1 for ED result\n");
            success = 0;
        }
        
        /* Verify it's a valid hash */
        if (kh_size(result) == 0) {
            /* EmptyDrops returned empty filter - this is valid if all cells fail */
            printf("    Note: EmptyDrops returned 0 filtered barcodes (expected for synthetic data)\n");
        }
        
        /* Clean up the filter we own */
        pf_filter_free(result);
    } else if (status == PF_FILTER_FAILED) {
        /* EmptyDrops failed - this is acceptable for synthetic data,
         * but we should have gotten a warning message */
        printf("    Note: EmptyDrops returned FAILED (expected for synthetic data without enough cells)\n");
    }
    
    /* Verify EmptyDrops output files were attempted (dir should exist) */
    char ed_dir[512];
    snprintf(ed_dir, sizeof(ed_dir), "%s/EmptyDrops", test_dir);
    struct stat st;
    if (stat(ed_dir, &st) == 0 && S_ISDIR(st.st_mode)) {
        printf("    EmptyDrops output directory created: OK\n");
    }
    
    free_test_feature_names(feature_names, n_features);
    pf_counts_result_free(counts);
    
    return success;
}

/* ============================================================================
 * Test: pf_filter_barcodes with skip_emptydrops flag
 * ============================================================================ */
static int test_filter_skip_emptydrops(void) {
    int n_features = 3;
    int n_barcodes = 5;
    
    pf_counts_result *counts = create_test_counts(n_features, n_barcodes);
    if (!counts) return 0;
    
    char **feature_names = create_test_feature_names(n_features);
    if (!feature_names) {
        pf_counts_result_free(counts);
        return 0;
    }
    
    barcode_filter_config config = {
        .output_dir = test_dir,
        .feature_names = feature_names,
        .n_features = n_features,
        .n_expected_cells = 0,
        .translate_nxt = 0,
        .skip_emptydrops = 1,  /* Skip EmptyDrops */
        .emptydrops_failure_fatal = 0
    };
    
    /* Call filter with skip_emptydrops - should return SKIPPED with NULL filter */
    khash_t(strptr) *result = NULL;
    pf_filter_status status = pf_filter_barcodes(counts, &config, NULL, &result);
    
    int success = 1;
    
    if (status != PF_FILTER_SKIPPED) {
        printf("    Expected PF_FILTER_SKIPPED, got %d\n", status);
        success = 0;
    }
    
    if (result != NULL) {
        printf("    Expected NULL filter when skipping EmptyDrops\n");
        success = 0;
        pf_filter_free(result);
    }
    
    free_test_feature_names(feature_names, n_features);
    pf_counts_result_free(counts);
    
    return success;
}

/* ============================================================================
 * Test: pf_filter_barcodes failure with emptydrops_failure_fatal=0
 * ============================================================================ */
static int test_filter_failure_nonfatal(void) {
    int n_features = 3;
    
    /* Create empty counts to force EmptyDrops failure */
    pf_counts_result *counts = pf_counts_result_init(n_features);
    if (!counts) return 0;
    
    char **feature_names = create_test_feature_names(n_features);
    if (!feature_names) {
        pf_counts_result_free(counts);
        return 0;
    }
    
    barcode_filter_config config = {
        .output_dir = test_dir,
        .feature_names = feature_names,
        .n_features = n_features,
        .n_expected_cells = 0,
        .translate_nxt = 0,
        .skip_emptydrops = 0,
        .emptydrops_failure_fatal = 0  /* Non-fatal mode */
    };
    
    /* Call filter with empty counts - EmptyDrops should fail */
    khash_t(strptr) *result = NULL;
    pf_filter_status status = pf_filter_barcodes(counts, &config, NULL, &result);
    
    int success = 1;
    
    if (status != PF_FILTER_FAILED) {
        printf("    Expected PF_FILTER_FAILED for empty counts, got %d\n", status);
        success = 0;
    }
    
    if (result != NULL) {
        printf("    Expected NULL filter on failure\n");
        success = 0;
        pf_filter_free(result);
    }
    
    free_test_feature_names(feature_names, n_features);
    pf_counts_result_free(counts);
    
    return success;
}

/* ============================================================================
 * Test: pf_filter_barcodes failure with emptydrops_failure_fatal=1
 * ============================================================================ */
static int test_filter_failure_fatal(void) {
    int n_features = 3;
    
    /* Create empty counts to force EmptyDrops failure */
    pf_counts_result *counts = pf_counts_result_init(n_features);
    if (!counts) return 0;
    
    char **feature_names = create_test_feature_names(n_features);
    if (!feature_names) {
        pf_counts_result_free(counts);
        return 0;
    }
    
    barcode_filter_config config = {
        .output_dir = test_dir,
        .feature_names = feature_names,
        .n_features = n_features,
        .n_expected_cells = 0,
        .translate_nxt = 0,
        .skip_emptydrops = 0,
        .emptydrops_failure_fatal = 1  /* Fatal mode */
    };
    
    /* Call filter with empty counts - EmptyDrops should fail */
    khash_t(strptr) *result = NULL;
    pf_filter_status status = pf_filter_barcodes(counts, &config, NULL, &result);
    
    int success = 1;
    
    /* Should still return FAILED - the fatal flag just affects logging,
     * the caller decides what to do with the failure */
    if (status != PF_FILTER_FAILED) {
        printf("    Expected PF_FILTER_FAILED for empty counts, got %d\n", status);
        success = 0;
    }
    
    if (result != NULL) {
        printf("    Expected NULL filter on failure\n");
        success = 0;
        pf_filter_free(result);
    }
    
    free_test_feature_names(feature_names, n_features);
    pf_counts_result_free(counts);
    
    return success;
}

/* ============================================================================
 * Test: run_emptydrops_on_counts n_features validation
 * ============================================================================ */
static int test_nfeatures_mismatch_warning(void) {
    int n_features = 3;
    int n_barcodes = 5;
    
    pf_counts_result *counts = create_test_counts(n_features, n_barcodes);
    if (!counts) return 0;
    
    char **feature_names = create_test_feature_names(n_features);
    if (!feature_names) {
        pf_counts_result_free(counts);
        return 0;
    }
    
    /* Call with mismatched n_features - should warn but not crash */
    int wrong_n_features = n_features + 5;  /* Intentionally wrong */
    
    /* This should log a warning about the mismatch */
    khash_t(strptr) *result = run_emptydrops_on_counts(
        counts,
        feature_names,
        wrong_n_features,  /* Mismatched! */
        test_dir,
        0,
        0
    );
    
    /* May succeed or fail, but shouldn't crash */
    int success = 1;  /* Pass as long as no crash */
    
    if (result) {
        pf_filter_free(result);
    }
    
    free_test_feature_names(feature_names, n_features);
    pf_counts_result_free(counts);
    
    return success;
}

/* ============================================================================
 * Main
 * ============================================================================ */
int main(int argc, char **argv) {
    (void)argc;
    (void)argv;
    
    printf("Running barcode_filter tests...\n\n");
    
    /* Initialize globals */
    init_globals();
    create_test_dir();
    
    /* Run tests */
    printf("  Testing emptydrops_input_from_counts...");
    tests_run++;
    if (test_emptydrops_input_from_counts()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing emptydrops_input_empty...");
    tests_run++;
    if (test_emptydrops_input_empty()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing filter_with_external...");
    tests_run++;
    if (test_filter_with_external()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing filter_free_null...");
    tests_run++;
    if (test_filter_free_null()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing emptydrops_integration...");
    tests_run++;
    if (test_emptydrops_integration()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing nfeatures_mismatch_warning...");
    tests_run++;
    if (test_nfeatures_mismatch_warning()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing filter_skip_emptydrops...");
    tests_run++;
    if (test_filter_skip_emptydrops()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing filter_failure_nonfatal...");
    tests_run++;
    if (test_filter_failure_nonfatal()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing filter_failure_fatal...");
    tests_run++;
    if (test_filter_failure_fatal()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    cleanup_test_dir();
    
    printf("\n%d/%d tests passed\n", tests_passed, tests_run);
    
    return (tests_passed == tests_run) ? 0 : 1;
}
