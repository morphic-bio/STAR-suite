/**
 * @file test_offset_detection.c
 * @brief Unit tests for pf_api offset auto-detection
 * 
 * Tests:
 *   1. Single offset auto-detect (all features same offset)
 *   2. Multi-offset WARNING (default: warn + proceed with dominant offset)
 *   3. Multi-offset ERROR (strict mode: error on heterogeneous offsets)
 *   4. Conflict error (global offset + per-feature array)
 *   5. Explicit offset bypasses auto-detect
 *   6. No pattern fallback to 0
 *   7. Heterogeneity rule: second-best > 5% of dominant triggers detection
 * 
 * Note: Tests that pass preflight but fail at FASTQ scan will exit() due to
 * organize_fastq_files_by_directory calling exit(). We create minimal dummy
 * FASTQs to get past that check.
 */

#include "../include/pf_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#define TEST_ASSERT(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
        return 1; \
    } \
} while(0)

static char g_temp_dir[256];
static char g_whitelist_path[512];

/* Create a minimal whitelist file */
static int create_whitelist(void) {
    FILE *f = fopen(g_whitelist_path, "w");
    if (!f) return -1;
    fprintf(f, "AAACCCAAGAAACACT\n");
    fprintf(f, "AAACCCAAGAAACCAT\n");
    fprintf(f, "AAACCCAAGAAACCCA\n");
    fclose(f);
    return 0;
}

/* Create minimal dummy FASTQ files so pf_process_fastq_dir doesn't exit() */
static int create_dummy_fastqs(const char *dir) {
    char path[512];
    FILE *f;
    
    /* Create R1 (barcode) file */
    snprintf(path, sizeof(path), "%s/test_R1_001.fastq", dir);
    f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "@read1\nAAACCCAAGAAACACTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n");
    fclose(f);
    
    /* Create R2 (forward) file */
    snprintf(path, sizeof(path), "%s/test_R2_001.fastq", dir);
    f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n");
    fclose(f);
    
    return 0;
}

/* Create a feature CSV with uniform offsets (all pattern has (BC) at position 4) */
static int create_uniform_features(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "name,sequence,pattern\n");
    fprintf(f, "Feature1,ACGTACGT,NNNN(BC)NNNN\n");
    fprintf(f, "Feature2,TGCATGCA,NNNN(BC)NNNN\n");
    fprintf(f, "Feature3,GATCGATC,NNNN(BC)NNNN\n");
    fclose(f);
    return 0;
}

/* Create a feature CSV with heterogeneous offsets */
static int create_hetero_features(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "name,sequence,pattern\n");
    fprintf(f, "Feature1,ACGTACGT,NNNN(BC)NNNN\n");  /* offset 4 */
    fprintf(f, "Feature2,TGCATGCA,NNNN(BC)NNNN\n");  /* offset 4 */
    fprintf(f, "Feature3,GATCGATC,NNNN(BC)NNNN\n");  /* offset 4 */
    fprintf(f, "Feature4,AAAACCCC,NN(BC)NNNNNN\n");  /* offset 2 */
    fprintf(f, "Feature5,GGGGTTTT,NN(BC)NNNNNN\n");  /* offset 2 */
    fclose(f);
    return 0;
}

/* Create a feature CSV without pattern column */
static int create_no_pattern_features(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "name,sequence\n");
    fprintf(f, "Feature1,ACGTACGT\n");
    fprintf(f, "Feature2,TGCATGCA\n");
    fclose(f);
    return 0;
}

/* Create features where second-best offset is BELOW 5% threshold (no heterogeneity) */
static int create_below_threshold_features(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "name,sequence,pattern\n");
    /* 100 features at offset 4 */
    for (int i = 0; i < 100; i++) {
        fprintf(f, "Feature%d,ACGTACGT,NNNN(BC)NNNN\n", i);
    }
    /* 4 features at offset 2 (4% of dominant = below 5% threshold) */
    for (int i = 100; i < 104; i++) {
        fprintf(f, "Feature%d,TGCATGCA,NN(BC)NNNNNN\n", i);
    }
    fclose(f);
    return 0;
}

/* Test 1: Multi-offset WARNING (default behavior: warn + proceed with dominant) */
static int test_multi_offset_warning(void) {
    printf("Test: Multi-offset warning (default)... ");
    fflush(stdout);
    
    char feature_path[512];
    char fastq_dir[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/hetero_features_warn.csv", g_temp_dir);
    snprintf(fastq_dir, sizeof(fastq_dir), "%s/fastqs_warn", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_warn", g_temp_dir);
    
    TEST_ASSERT(create_hetero_features(feature_path) == 0, "create feature file");
    mkdir(fastq_dir, 0755);
    mkdir(output_path, 0755);
    TEST_ASSERT(create_dummy_fastqs(fastq_dir) == 0, "create dummy FASTQs");
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    /* Default: strict_offset_check = 0, so warning + proceed */
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* Default: warn + proceed with dominant offset - should NOT error */
    err = pf_process_fastq_dir(ctx, fastq_dir, output_path, NULL);
    TEST_ASSERT(err != PF_ERR_MULTI_OFFSET_DETECTED, "should NOT get MULTI_OFFSET (warn + proceed)");
    TEST_ASSERT(err != PF_ERR_OFFSET_CONFLICT, "should NOT get OFFSET_CONFLICT");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Test 2: Multi-offset ERROR (strict mode: error on heterogeneous offsets) */
static int test_multi_offset_strict_error(void) {
    printf("Test: Multi-offset error (strict mode)... ");
    fflush(stdout);
    
    char feature_path[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/hetero_features_strict.csv", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_hetero_strict", g_temp_dir);
    
    TEST_ASSERT(create_hetero_features(feature_path) == 0, "create feature file");
    mkdir(output_path, 0755);
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    
    /* Enable strict mode - should error on heterogeneous offsets */
    pf_config_set_strict_offset_check(config, 1);
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* Strict mode: should fail with MULTI_OFFSET_DETECTED in preflight */
    err = pf_process_fastq_dir(ctx, g_temp_dir, output_path, NULL);
    TEST_ASSERT(err == PF_ERR_MULTI_OFFSET_DETECTED, "expected MULTI_OFFSET_DETECTED (strict mode)");
    
    /* Verify error message contains guidance */
    const char *errmsg = pf_get_error(ctx);
    TEST_ASSERT(errmsg != NULL, "error message exists");
    TEST_ASSERT(strstr(errmsg, "Multiple feature offsets") != NULL, "error mentions multiple offsets");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Test 2: Conflict error - both global offset and per-feature array */
static int test_conflict_error(void) {
    printf("Test: Conflict error (global + array)... ");
    fflush(stdout);
    
    char feature_path[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/uniform_features.csv", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_conflict", g_temp_dir);
    
    TEST_ASSERT(create_uniform_features(feature_path) == 0, "create feature file");
    mkdir(output_path, 0755);
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    
    /* Set BOTH global offset AND per-feature array - should conflict */
    pf_config_set_feature_offset(config, 10);
    pf_config_set_use_feature_offset_array(config, 1);
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* Try to process - should fail with OFFSET_CONFLICT in preflight */
    err = pf_process_fastq_dir(ctx, g_temp_dir, output_path, NULL);
    TEST_ASSERT(err == PF_ERR_OFFSET_CONFLICT, "expected OFFSET_CONFLICT");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Test 3: Error codes are defined */
static int test_error_codes(void) {
    printf("Test: Error codes defined... ");
    fflush(stdout);
    
    TEST_ASSERT(PF_ERR_OFFSET_CONFLICT == -8, "PF_ERR_OFFSET_CONFLICT value");
    TEST_ASSERT(PF_ERR_MULTI_OFFSET_DETECTED == -9, "PF_ERR_MULTI_OFFSET_DETECTED value");
    
    printf("PASS\n");
    return 0;
}

/* Test 4: Single offset auto-detect (requires dummy FASTQs) */
static int test_single_offset_auto_detect(void) {
    printf("Test: Single offset auto-detect... ");
    fflush(stdout);
    
    char feature_path[512];
    char fastq_dir[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/uniform_features.csv", g_temp_dir);
    snprintf(fastq_dir, sizeof(fastq_dir), "%s/fastqs_uniform", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_uniform", g_temp_dir);
    
    TEST_ASSERT(create_uniform_features(feature_path) == 0, "create feature file");
    mkdir(fastq_dir, 0755);
    mkdir(output_path, 0755);
    TEST_ASSERT(create_dummy_fastqs(fastq_dir) == 0, "create dummy FASTQs");
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* Process - preflight should pass (auto-detect offset 4), then process sample */
    err = pf_process_fastq_dir(ctx, fastq_dir, output_path, NULL);
    
    /* Should NOT get offset errors - may fail for other reasons but that's OK */
    TEST_ASSERT(err != PF_ERR_MULTI_OFFSET_DETECTED, "should NOT get MULTI_OFFSET");
    TEST_ASSERT(err != PF_ERR_OFFSET_CONFLICT, "should NOT get OFFSET_CONFLICT");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Test 5: Explicit offset bypasses auto-detect (requires dummy FASTQs) */
static int test_explicit_offset_bypass(void) {
    printf("Test: Explicit offset bypasses auto-detect... ");
    fflush(stdout);
    
    char feature_path[512];
    char fastq_dir[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/hetero_features.csv", g_temp_dir);
    snprintf(fastq_dir, sizeof(fastq_dir), "%s/fastqs_explicit", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_explicit", g_temp_dir);
    
    TEST_ASSERT(create_hetero_features(feature_path) == 0, "create feature file");
    mkdir(fastq_dir, 0755);
    mkdir(output_path, 0755);
    TEST_ASSERT(create_dummy_fastqs(fastq_dir) == 0, "create dummy FASTQs");
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    
    /* Set explicit offset - should bypass auto-detect even with hetero features */
    pf_config_set_feature_offset(config, 4);
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* Process - explicit offset should bypass the multi-offset check */
    err = pf_process_fastq_dir(ctx, fastq_dir, output_path, NULL);
    
    /* Should NOT get offset errors - bypassed by explicit setting */
    TEST_ASSERT(err != PF_ERR_MULTI_OFFSET_DETECTED, "should NOT get MULTI_OFFSET");
    TEST_ASSERT(err != PF_ERR_OFFSET_CONFLICT, "should NOT get OFFSET_CONFLICT");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Test 6: Per-feature array flag bypasses auto-detect (requires dummy FASTQs) */
static int test_array_flag_bypass(void) {
    printf("Test: Per-feature array flag bypasses auto-detect... ");
    fflush(stdout);
    
    char feature_path[512];
    char fastq_dir[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/hetero_features.csv", g_temp_dir);
    snprintf(fastq_dir, sizeof(fastq_dir), "%s/fastqs_array", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_array", g_temp_dir);
    
    mkdir(fastq_dir, 0755);
    mkdir(output_path, 0755);
    TEST_ASSERT(create_dummy_fastqs(fastq_dir) == 0, "create dummy FASTQs");
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    
    /* Set per-feature array flag - should bypass auto-detect even with hetero features */
    pf_config_set_use_feature_offset_array(config, 1);
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* Process - array flag should bypass the multi-offset check */
    err = pf_process_fastq_dir(ctx, fastq_dir, output_path, NULL);
    
    /* Should NOT get offset errors - bypassed by array flag */
    TEST_ASSERT(err != PF_ERR_MULTI_OFFSET_DETECTED, "should NOT get MULTI_OFFSET");
    TEST_ASSERT(err != PF_ERR_OFFSET_CONFLICT, "should NOT get OFFSET_CONFLICT");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Test 7: No pattern column falls back to 0 (requires dummy FASTQs) */
static int test_no_pattern_fallback(void) {
    printf("Test: No pattern column falls back to 0... ");
    fflush(stdout);
    
    char feature_path[512];
    char fastq_dir[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/no_pattern_features.csv", g_temp_dir);
    snprintf(fastq_dir, sizeof(fastq_dir), "%s/fastqs_nopattern", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_nopattern", g_temp_dir);
    
    TEST_ASSERT(create_no_pattern_features(feature_path) == 0, "create feature file");
    mkdir(fastq_dir, 0755);
    mkdir(output_path, 0755);
    TEST_ASSERT(create_dummy_fastqs(fastq_dir) == 0, "create dummy FASTQs");
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* Process - should fall back to offset 0 and continue */
    err = pf_process_fastq_dir(ctx, fastq_dir, output_path, NULL);
    
    /* Should NOT get offset errors */
    TEST_ASSERT(err != PF_ERR_MULTI_OFFSET_DETECTED, "should NOT get MULTI_OFFSET");
    TEST_ASSERT(err != PF_ERR_OFFSET_CONFLICT, "should NOT get OFFSET_CONFLICT");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Test 8: Heterogeneity rule - second-best < 5% should NOT trigger warning */
static int test_heterogeneity_below_threshold(void) {
    printf("Test: Heterogeneity rule (below 5%% threshold)... ");
    fflush(stdout);
    
    char feature_path[512];
    char fastq_dir[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/below_threshold_features.csv", g_temp_dir);
    snprintf(fastq_dir, sizeof(fastq_dir), "%s/fastqs_below", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_below", g_temp_dir);
    
    TEST_ASSERT(create_below_threshold_features(feature_path) == 0, "create feature file");
    mkdir(fastq_dir, 0755);
    mkdir(output_path, 0755);
    TEST_ASSERT(create_dummy_fastqs(fastq_dir) == 0, "create dummy FASTQs");
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    
    /* Even with strict mode, below-threshold shouldn't error */
    pf_config_set_strict_offset_check(config, 1);
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* 4% heterogeneity is below 5% threshold - should NOT error even with strict */
    err = pf_process_fastq_dir(ctx, fastq_dir, output_path, NULL);
    TEST_ASSERT(err != PF_ERR_MULTI_OFFSET_DETECTED, "4%% heterogeneity should NOT trigger error");
    TEST_ASSERT(err != PF_ERR_OFFSET_CONFLICT, "should NOT get OFFSET_CONFLICT");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Test 9: Explicit offset 0 bypasses auto-detect (tests feature_offset_explicit flag) */
static int test_explicit_offset_zero(void) {
    printf("Test: Explicit offset 0 bypasses auto-detect... ");
    fflush(stdout);
    
    char feature_path[512];
    char fastq_dir[512];
    char output_path[512];
    snprintf(feature_path, sizeof(feature_path), "%s/hetero_features.csv", g_temp_dir);
    snprintf(fastq_dir, sizeof(fastq_dir), "%s/fastqs_zero", g_temp_dir);
    snprintf(output_path, sizeof(output_path), "%s/out_zero", g_temp_dir);
    
    mkdir(fastq_dir, 0755);
    mkdir(output_path, 0755);
    TEST_ASSERT(create_dummy_fastqs(fastq_dir) == 0, "create dummy FASTQs");
    
    pf_config *config = pf_config_create();
    TEST_ASSERT(config != NULL, "config creation");
    
    /* Explicitly set offset to 0 - should bypass auto-detect even with hetero features */
    pf_config_set_feature_offset(config, 0);
    /* Enable strict mode to ensure it would error if auto-detect ran */
    pf_config_set_strict_offset_check(config, 1);
    
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    TEST_ASSERT(ctx != NULL, "context init");
    
    pf_error err = pf_load_feature_ref(ctx, feature_path);
    TEST_ASSERT(err == PF_OK, "load features");
    
    err = pf_load_whitelist(ctx, g_whitelist_path);
    TEST_ASSERT(err == PF_OK, "load whitelist");
    
    /* Explicit offset 0 should bypass auto-detect - no multi-offset error */
    err = pf_process_fastq_dir(ctx, fastq_dir, output_path, NULL);
    TEST_ASSERT(err != PF_ERR_MULTI_OFFSET_DETECTED, "explicit offset 0 should bypass auto-detect");
    TEST_ASSERT(err != PF_ERR_OFFSET_CONFLICT, "should NOT get OFFSET_CONFLICT");
    
    pf_destroy(ctx);
    printf("PASS\n");
    return 0;
}

/* Cleanup temp directory */
static void cleanup_temp_dir(void) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "rm -rf %s", g_temp_dir);
    (void)system(cmd);  /* Cast to void to silence unused result warning */
}

int main(int argc, char *argv[]) {
    int failures = 0;
    
    /* Create temp directory */
    snprintf(g_temp_dir, sizeof(g_temp_dir), "/tmp/test_offset_detection_%d", getpid());
    mkdir(g_temp_dir, 0755);
    snprintf(g_whitelist_path, sizeof(g_whitelist_path), "%s/whitelist.txt", g_temp_dir);
    
    if (create_whitelist() != 0) {
        fprintf(stderr, "Failed to create whitelist\n");
        cleanup_temp_dir();
        return 1;
    }
    
    /* Initialize pf global state */
    pf_global_init();
    
    printf("\n=== pf_api Offset Detection Tests ===\n\n");
    
    /* Test error code values */
    failures += test_error_codes();
    
    /* Test strict mode errors (fail in preflight before FASTQ scan) */
    failures += test_multi_offset_strict_error();
    failures += test_conflict_error();
    
    /* Test default warn + proceed behavior (needs dummy FASTQs) */
    failures += test_multi_offset_warning();
    
    /* These tests pass preflight and need dummy FASTQs */
    failures += test_single_offset_auto_detect();
    failures += test_explicit_offset_bypass();
    failures += test_array_flag_bypass();
    failures += test_no_pattern_fallback();
    
    /* Test heterogeneity rule (5% threshold) and explicit offset 0 */
    failures += test_heterogeneity_below_threshold();
    failures += test_explicit_offset_zero();
    
    printf("\n=== Results: %d failures ===\n\n", failures);
    
    cleanup_temp_dir();
    return failures > 0 ? 1 : 0;
}
