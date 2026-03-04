/**
 * test_prehash_memory_safety.c
 *
 * Tests that prehash memory budget enforcement works correctly:
 *   (a) low budget forces prehash skip, but assignment still completes
 *   (b) normal budget allows prehash build, correct assignment
 */

#include "../include/pf_api.h"
#include "../include/common.h"
#include "../include/globals.h"
#include "../include/barcode_match.h"
#include "../include/io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define TEST_ASSERT(cond, msg) do { \
    if (!(cond)) { fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); return 1; } \
    else { fprintf(stderr, "  PASS: %s\n", msg); } \
} while(0)

static int write_test_feature_csv(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "name,sequence\n");
    fprintf(f, "feat1,ACGTACGT\n");
    fprintf(f, "feat2,TTTTCCCC\n");
    fprintf(f, "feat3,GGGGAAAA\n");
    fprintf(f, "feat4,ACACACAC\n");
    fprintf(f, "feat5,TGTGTGTG\n");
    fclose(f);
    return 0;
}

static int test_low_budget_skips_prehash(void) {
    fprintf(stderr, "\n--- Test: low budget forces prehash skip ---\n");

    const char *csv = "/tmp/test_prehash_budget_features.csv";
    TEST_ASSERT(write_test_feature_csv(csv) == 0, "write feature CSV");

    initseq2Code();
    initcode2seq();
    initdiff2hamming(diff2Hamming);
    initialize_complement();

    feature_prehash_max_hamming = 2;
    feature_prehash_max_entries = 50000000ULL;
    feature_prehash_memory_budget = 1024; /* 1 KB - way too small for any prehash */

    feature_arrays *features = read_features_file(csv);
    TEST_ASSERT(features != NULL, "features loaded");
    TEST_ASSERT(features->number_of_features == 5, "5 features loaded");
    TEST_ASSERT(features->feature_hamming_le1_enabled == 0, "d1 prehash skipped (low budget)");
    TEST_ASSERT(features->feature_hamming_le2_enabled == 0, "d2 prehash skipped (low budget)");

    /* Exact hash lookup should still work */
    int idx = feature_lookup_seq("ACGTACGT", 8);
    TEST_ASSERT(idx == 1, "exact lookup feat1 works without prehash");
    idx = feature_lookup_seq("TTTTCCCC", 8);
    TEST_ASSERT(idx == 2, "exact lookup feat2 works without prehash");

    free_feature_arrays(features);
    feature_code_hash.h64 = NULL;
    unlink(csv);

    feature_prehash_memory_budget = 0;
    return 0;
}

static int test_normal_budget_builds_prehash(void) {
    fprintf(stderr, "\n--- Test: normal budget allows prehash build ---\n");

    const char *csv = "/tmp/test_prehash_normal_features.csv";
    TEST_ASSERT(write_test_feature_csv(csv) == 0, "write feature CSV");

    initseq2Code();
    initcode2seq();
    initdiff2hamming(diff2Hamming);
    initialize_complement();

    feature_prehash_max_hamming = 2;
    feature_prehash_max_entries = 50000000ULL;
    feature_prehash_memory_budget = 1024ULL * 1024ULL * 1024ULL; /* 1 GB - plenty */

    feature_arrays *features = read_features_file(csv);
    TEST_ASSERT(features != NULL, "features loaded");
    TEST_ASSERT(features->feature_hamming_le1_enabled == 1, "d1 prehash built");
    TEST_ASSERT(features->feature_hamming_le2_enabled == 1, "d2 prehash built");

    /* Exact hash lookup should work */
    int idx = feature_lookup_seq("ACGTACGT", 8);
    TEST_ASSERT(idx == 1, "exact lookup feat1 works with prehash");

    /* Hamming d=1 lookup via prehash should work */
    char query_d1[] = "ACGTACGC"; /* 1 mismatch from ACGTACGT */
    khint_t k = kh_get(stru32, features->feature_hamming_le1_hash, query_d1);
    TEST_ASSERT(k != kh_end(features->feature_hamming_le1_hash), "d1 prehash finds 1-mismatch variant");

    free_feature_arrays(features);
    feature_code_hash.h64 = NULL;
    unlink(csv);

    feature_prehash_memory_budget = 0;
    return 0;
}

static int test_autodetect_budget(void) {
    fprintf(stderr, "\n--- Test: auto-detect memory budget ---\n");
    unsigned long long budget = prehash_detect_memory_budget();
    fprintf(stderr, "  Detected budget: %llu bytes (%.1f GB)\n",
            budget, (double)budget / (1024.0*1024.0*1024.0));
    /* On most Linux systems, this should be > 0 */
    TEST_ASSERT(budget > 0, "auto-detected budget is positive");
    TEST_ASSERT(budget > 1024ULL * 1024ULL, "budget > 1 MB (reasonable minimum)");
    return 0;
}

int main(void) {
    int failures = 0;
    failures += test_low_budget_skips_prehash();
    failures += test_normal_budget_builds_prehash();
    failures += test_autodetect_budget();

    fprintf(stderr, "\n=== %s ===\n", failures ? "SOME TESTS FAILED" : "ALL TESTS PASSED");
    return failures;
}
