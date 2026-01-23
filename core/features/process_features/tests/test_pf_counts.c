/**
 * @file test_pf_counts.c
 * @brief Unit tests for pf_counts module
 * 
 * Tests the pf_counts_result structure initialization, cleanup,
 * and basic functionality. Uses a minimal synthetic dataset.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../include/pf_counts.h"
#include "../include/globals.h"
#include "../include/prototypes.h"

/* Test counters */
static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    printf("  Testing %s...", #name); \
    tests_run++; \
    if (test_##name()) { \
        printf(" PASS\n"); \
        tests_passed++; \
    } else { \
        printf(" FAIL\n"); \
    } \
} while(0)

/* ============================================================================
 * Test: pf_counts_result_init
 * ============================================================================ */
static int test_pf_counts_result_init(void) {
    int n_features = 10;
    
    pf_counts_result *result = pf_counts_result_init(n_features);
    if (!result) return 0;
    
    /* Verify fields are initialized */
    if (result->n_features != n_features) {
        pf_counts_result_free(result);
        return 0;
    }
    
    if (!result->barcode_to_deduped_hash) {
        pf_counts_result_free(result);
        return 0;
    }
    
    if (!result->feature_hist) {
        pf_counts_result_free(result);
        return 0;
    }
    
    if (!result->total_deduped_counts || !result->total_barcoded_counts) {
        pf_counts_result_free(result);
        return 0;
    }
    
    /* Verify all histograms are NULL initially */
    for (int i = 0; i <= n_features; i++) {
        if (result->feature_hist[i] != NULL) {
            pf_counts_result_free(result);
            return 0;
        }
    }
    
    /* Verify counts arrays are zeroed */
    for (int i = 0; i <= n_features; i++) {
        if (result->total_deduped_counts[i] != 0 || 
            result->total_barcoded_counts[i] != 0) {
            pf_counts_result_free(result);
            return 0;
        }
    }
    
    pf_counts_result_free(result);
    return 1;
}

/* ============================================================================
 * Test: pf_counts_result_free_null
 * Verify that freeing NULL is safe
 * ============================================================================ */
static int test_pf_counts_result_free_null(void) {
    /* Should not crash */
    pf_counts_result_free(NULL);
    return 1;
}

/* ============================================================================
 * Test: pf_counts_result_reset_histograms
 * ============================================================================ */
static int test_pf_counts_result_reset_histograms(void) {
    int n_features = 5;
    
    pf_counts_result *result = pf_counts_result_init(n_features);
    if (!result) return 0;
    
    /* Manually add some histogram data */
    result->feature_hist[1] = vec_u32_init();
    vec_u32_set(result->feature_hist[1], 0, 10);
    vec_u32_set(result->feature_hist[1], 1, 20);
    
    result->feature_hist[2] = vec_u32_init();
    vec_u32_set(result->feature_hist[2], 0, 5);
    
    result->total_deduped_counts[1] = 100;
    result->total_barcoded_counts[1] = 200;
    
    /* Reset */
    pf_counts_result_reset_histograms(result);
    
    /* Verify histograms are cleared */
    for (int i = 0; i <= n_features; i++) {
        if (result->feature_hist[i] != NULL) {
            pf_counts_result_free(result);
            return 0;
        }
    }
    
    /* Verify counts are reset */
    if (result->total_deduped_counts[1] != 0 ||
        result->total_barcoded_counts[1] != 0) {
        pf_counts_result_free(result);
        return 0;
    }
    
    pf_counts_result_free(result);
    return 1;
}

/* ============================================================================
 * Test: pf_counts_result_with_hash
 * Verify manual insertion into the hash works
 * ============================================================================ */
static int test_pf_counts_result_with_hash(void) {
    int n_features = 3;
    
    pf_counts_result *result = pf_counts_result_init(n_features);
    if (!result) return 0;
    
    /* Manually insert a barcode with feature counts */
    uint32_t barcode_key = 0x12345678;
    khash_t(u32u32) *inner_hash = kh_init(u32u32);
    
    /* Add counts for feature 1 and 2 */
    int ret;
    khint_t k1 = kh_put(u32u32, inner_hash, 1, &ret);
    kh_val(inner_hash, k1) = 10;
    
    khint_t k2 = kh_put(u32u32, inner_hash, 2, &ret);
    kh_val(inner_hash, k2) = 20;
    
    /* Insert into barcode_to_deduped_hash */
    khint_t kb = kh_put(u32ptr, result->barcode_to_deduped_hash, barcode_key, &ret);
    kh_val(result->barcode_to_deduped_hash, kb) = inner_hash;
    
    /* Verify we can retrieve it */
    khint_t kget = kh_get(u32ptr, result->barcode_to_deduped_hash, barcode_key);
    if (kget == kh_end(result->barcode_to_deduped_hash)) {
        pf_counts_result_free(result);
        return 0;
    }
    
    khash_t(u32u32) *retrieved = (khash_t(u32u32)*)kh_val(result->barcode_to_deduped_hash, kget);
    if (kh_size(retrieved) != 2) {
        pf_counts_result_free(result);
        return 0;
    }
    
    /* Free should handle the nested hash */
    pf_counts_result_free(result);
    return 1;
}

/* ============================================================================
 * Main test runner
 * ============================================================================ */
int main(int argc, char **argv) {
    (void)argc;
    (void)argv;
    
    printf("Running pf_counts unit tests...\n\n");
    
    /* Initialize globals that may be needed */
    barcode_length = 16;
    barcode_code_length = 4;
    number_of_features = 10;
    
    TEST(pf_counts_result_init);
    TEST(pf_counts_result_free_null);
    TEST(pf_counts_result_reset_histograms);
    TEST(pf_counts_result_with_hash);
    
    printf("\n%d/%d tests passed\n", tests_passed, tests_run);
    
    return (tests_passed == tests_run) ? 0 : 1;
}
