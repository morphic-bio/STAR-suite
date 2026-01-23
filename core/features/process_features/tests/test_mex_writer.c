/**
 * @file test_mex_writer.c
 * @brief Unit tests for mex_writer module
 * 
 * Tests the mex_writer output functions with synthetic data.
 * Verifies that output files are created with correct format.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <unistd.h>
#include "../include/mex_writer.h"
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

/* Helper: check if file exists */
static int file_exists(const char *path) {
    struct stat st;
    return (stat(path, &st) == 0);
}

/* Helper: read first line of file */
static int read_first_line(const char *path, char *buf, size_t buflen) {
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    if (!fgets(buf, buflen, fp)) {
        fclose(fp);
        return -1;
    }
    fclose(fp);
    return 0;
}

/* Helper: count lines in file */
static int count_lines(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;
    int count = 0;
    char buf[4096];
    while (fgets(buf, sizeof(buf), fp)) count++;
    fclose(fp);
    return count;
}

/* Helper: create test output directory */
static char* create_test_dir(void) {
    static char dir[256];
    snprintf(dir, sizeof(dir), "/tmp/test_mex_writer_%d", getpid());
    mkdir(dir, 0755);
    return dir;
}

/* Helper: remove test directory recursively */
static void cleanup_test_dir(const char *dir) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "rm -rf %s", dir);
    system(cmd);
}

/* Helper: create minimal feature_arrays for testing */
static feature_arrays* create_test_features(int n_features) {
    feature_arrays *features = calloc(1, sizeof(feature_arrays));
    features->number_of_features = n_features;
    features->feature_names = calloc(n_features, sizeof(char*));
    features->feature_names_storage = calloc(n_features * 32, 1);
    
    for (int i = 0; i < n_features; i++) {
        features->feature_names[i] = features->feature_names_storage + i * 32;
        snprintf(features->feature_names[i], 32, "Feature_%d", i + 1);
    }
    
    return features;
}

static void free_test_features(feature_arrays *features) {
    if (features) {
        free(features->feature_names);
        free(features->feature_names_storage);
        free(features);
    }
}

/* Helper: create minimal data_structures for testing */
static data_structures* create_test_hashes(void) {
    data_structures *hashes = calloc(1, sizeof(data_structures));
    hashes->filtered_hash = kh_init(u32ptr);
    hashes->sequence_umi_hash = kh_init(u64ptr);
    hashes->unique_features_match = kh_init(strptr);
    return hashes;
}

static void free_test_hashes(data_structures *hashes) {
    if (hashes) {
        if (hashes->filtered_hash) kh_destroy(u32ptr, hashes->filtered_hash);
        if (hashes->sequence_umi_hash) kh_destroy(u64ptr, hashes->sequence_umi_hash);
        if (hashes->unique_features_match) kh_destroy(strptr, hashes->unique_features_match);
        free(hashes);
    }
}

/* Helper: create test statistics */
static statistics* create_test_stats(void) {
    statistics *stats = calloc(1, sizeof(statistics));
    stats->number_of_reads = 1000;
    stats->total_unmatched_features = 100;
    return stats;
}

/* ============================================================================
 * Test: mex_write_core creates basic files
 * ============================================================================ */
static int test_mex_write_core_creates_files(void) {
    char *test_dir = create_test_dir();
    int success = 0;
    
    /* Setup */
    feature_arrays *features = create_test_features(3);
    data_structures *hashes = create_test_hashes();
    statistics *stats = create_test_stats();
    pf_counts_result *counts = pf_counts_result_init(3);
    
    /* Add a test barcode with counts */
    uint32_t barcode_key = 0x12345678;
    khash_t(u32u32) *inner = kh_init(u32u32);
    int ret;
    khint_t k = kh_put(u32u32, inner, 1, &ret);
    kh_val(inner, k) = 10;
    k = kh_put(u32u32, inner, 2, &ret);
    kh_val(inner, k) = 20;
    
    k = kh_put(u32ptr, counts->barcode_to_deduped_hash, barcode_key, &ret);
    kh_val(counts->barcode_to_deduped_hash, k) = inner;
    
    /* Also add to filtered_hash for iteration */
    feature_counts *fc = calloc(1, sizeof(feature_counts));
    *(uint32_t*)fc->sequence_code = barcode_key;
    fc->counts = kh_init(u32u32);
    k = kh_put(u32ptr, hashes->filtered_hash, barcode_key, &ret);
    kh_val(hashes->filtered_hash, k) = fc;
    
    /* Create config and call mex_write_core */
    mex_writer_config config = {
        .output_dir = test_dir,
        .features = features,
        .hashes = hashes,
        .stats = stats,
        .filtered_barcodes_hash = NULL,
        .min_heatmap_counts = 0
    };
    
    int rc = mex_write_core(&config, counts);
    
    /* Verify files exist */
    char path[512];
    snprintf(path, sizeof(path), "%s/barcodes.txt", test_dir);
    int has_barcodes = file_exists(path);
    
    snprintf(path, sizeof(path), "%s/features.txt", test_dir);
    int has_features = file_exists(path);
    
    snprintf(path, sizeof(path), "%s/matrix.mtx", test_dir);
    int has_matrix = file_exists(path);
    
    success = (rc == 0 && has_barcodes && has_features && has_matrix);
    
    /* Cleanup */
    kh_destroy(u32u32, fc->counts);
    free(fc);
    pf_counts_result_free(counts);
    free_test_hashes(hashes);
    free_test_features(features);
    free(stats);
    cleanup_test_dir(test_dir);
    
    return success;
}

/* ============================================================================
 * Test: matrix.mtx has correct header format
 * ============================================================================ */
static int test_matrix_header_format(void) {
    char *test_dir = create_test_dir();
    int success = 0;
    
    /* Setup */
    feature_arrays *features = create_test_features(3);
    data_structures *hashes = create_test_hashes();
    statistics *stats = create_test_stats();
    pf_counts_result *counts = pf_counts_result_init(3);
    
    mex_writer_config config = {
        .output_dir = test_dir,
        .features = features,
        .hashes = hashes,
        .stats = stats,
        .filtered_barcodes_hash = NULL,
        .min_heatmap_counts = 0
    };
    
    mex_write_core(&config, counts);
    
    /* Check matrix header */
    char path[512];
    snprintf(path, sizeof(path), "%s/matrix.mtx", test_dir);
    char line[1024];
    if (read_first_line(path, line, sizeof(line)) == 0) {
        /* Should start with %%MatrixMarket */
        success = (strncmp(line, "%%MatrixMarket", 14) == 0);
    }
    
    /* Cleanup */
    pf_counts_result_free(counts);
    free_test_hashes(hashes);
    free_test_features(features);
    free(stats);
    cleanup_test_dir(test_dir);
    
    return success;
}

/* ============================================================================
 * Test: features.txt has correct line count
 * ============================================================================ */
static int test_features_line_count(void) {
    char *test_dir = create_test_dir();
    int success = 0;
    
    int n_features = 5;
    feature_arrays *features = create_test_features(n_features);
    data_structures *hashes = create_test_hashes();
    statistics *stats = create_test_stats();
    pf_counts_result *counts = pf_counts_result_init(n_features);
    
    mex_writer_config config = {
        .output_dir = test_dir,
        .features = features,
        .hashes = hashes,
        .stats = stats,
        .filtered_barcodes_hash = NULL,
        .min_heatmap_counts = 0
    };
    
    mex_write_core(&config, counts);
    
    char path[512];
    snprintf(path, sizeof(path), "%s/features.txt", test_dir);
    int lines = count_lines(path);
    success = (lines == n_features);
    
    /* Cleanup */
    pf_counts_result_free(counts);
    free_test_hashes(hashes);
    free_test_features(features);
    free(stats);
    cleanup_test_dir(test_dir);
    
    return success;
}

/* ============================================================================
 * Test: filtered output goes to filtered/ subdirectory
 * ============================================================================ */
static int test_filtered_subdir(void) {
    char *test_dir = create_test_dir();
    int success = 0;
    
    feature_arrays *features = create_test_features(3);
    data_structures *hashes = create_test_hashes();
    statistics *stats = create_test_stats();
    pf_counts_result *counts = pf_counts_result_init(3);
    
    /* Create a filter hash */
    khash_t(strptr) *filter = kh_init(strptr);
    
    mex_writer_config config = {
        .output_dir = test_dir,
        .features = features,
        .hashes = hashes,
        .stats = stats,
        .filtered_barcodes_hash = filter,  /* Non-NULL triggers filtered/ subdir */
        .min_heatmap_counts = 0
    };
    
    mex_write_core(&config, counts);
    
    /* Check that filtered/features.txt exists */
    char path[512];
    snprintf(path, sizeof(path), "%s/filtered/features.txt", test_dir);
    success = file_exists(path);
    
    /* Cleanup */
    kh_destroy(strptr, filter);
    pf_counts_result_free(counts);
    free_test_hashes(hashes);
    free_test_features(features);
    free(stats);
    cleanup_test_dir(test_dir);
    
    return success;
}

/* ============================================================================
 * Test: counts arrays are cleared between passes
 * ============================================================================ */
static int test_counts_cleared_between_passes(void) {
    char *test_dir = create_test_dir();
    int success = 0;
    
    feature_arrays *features = create_test_features(3);
    data_structures *hashes = create_test_hashes();
    statistics *stats = create_test_stats();
    pf_counts_result *counts = pf_counts_result_init(3);
    
    /* Pre-populate counts to verify they get cleared */
    counts->total_deduped_counts[0] = 999;
    counts->total_deduped_counts[1] = 999;
    counts->total_barcoded_counts[0] = 888;
    
    mex_writer_config config = {
        .output_dir = test_dir,
        .features = features,
        .hashes = hashes,
        .stats = stats,
        .filtered_barcodes_hash = NULL,
        .min_heatmap_counts = 0
    };
    
    mex_write_core(&config, counts);
    
    /* After write, counts should be 0 (no actual data in hashes) */
    success = (counts->total_deduped_counts[0] == 0 &&
               counts->total_deduped_counts[1] == 0 &&
               counts->total_barcoded_counts[0] == 0);
    
    /* Cleanup */
    pf_counts_result_free(counts);
    free_test_hashes(hashes);
    free_test_features(features);
    free(stats);
    cleanup_test_dir(test_dir);
    
    return success;
}

/* ============================================================================
 * Main test runner
 * ============================================================================ */
int main(int argc, char **argv) {
    (void)argc;
    (void)argv;
    
    printf("Running mex_writer unit tests...\n\n");
    
    /* Initialize globals */
    barcode_length = 16;
    barcode_code_length = 4;
    umi_length = 12;
    umi_code_length = 3;
    number_of_features = 10;
    translate_NXT = 0;
    min_heatmap = 0;
    
    /* Initialize code tables */
    initseq2Code();
    initcode2seq();
    
    TEST(mex_write_core_creates_files);
    TEST(matrix_header_format);
    TEST(features_line_count);
    TEST(filtered_subdir);
    TEST(counts_cleared_between_passes);
    
    printf("\n%d/%d tests passed\n", tests_passed, tests_run);
    
    return (tests_passed == tests_run) ? 0 : 1;
}
