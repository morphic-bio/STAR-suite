/**
 * @file test_mex_writer_parity.c
 * @brief Parity test: mex_write_all() vs printFeatureCounts()
 * 
 * This test verifies that mex_write_all() produces byte-identical output
 * to the legacy printFeatureCounts() function. It:
 * 1. Creates synthetic test data with known values
 * 2. Calls printFeatureCounts() to generate baseline outputs
 * 3. Calls mex_write_all() to generate new outputs  
 * 4. Compares all output files byte-by-byte
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include "../include/mex_writer.h"
#include "../include/pf_counts.h"
#include "../include/globals.h"
#include "../include/prototypes.h"

/* Test state */
static int tests_run = 0;
static int tests_passed = 0;
static char test_dir_legacy[256];
static char test_dir_new[256];

/* ============================================================================
 * Helper functions
 * ============================================================================ */

static void create_test_dirs(void) {
    snprintf(test_dir_legacy, sizeof(test_dir_legacy), 
             "/tmp/test_mex_parity_legacy_%d", getpid());
    snprintf(test_dir_new, sizeof(test_dir_new), 
             "/tmp/test_mex_parity_new_%d", getpid());
    mkdir(test_dir_legacy, 0755);
    mkdir(test_dir_new, 0755);
}

static void cleanup_test_dirs(void) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "rm -rf %s %s", test_dir_legacy, test_dir_new);
    int rc = system(cmd);
    (void)rc;
}

/* Compare two files byte-by-byte, return 0 if identical */
static int compare_files(const char *path1, const char *path2) {
    FILE *f1 = fopen(path1, "rb");
    FILE *f2 = fopen(path2, "rb");
    
    if (!f1 && !f2) return 0;  /* Both don't exist = equal */
    if (!f1 || !f2) {
        if (f1) fclose(f1);
        if (f2) fclose(f2);
        fprintf(stderr, "    File exists mismatch: %s vs %s\n", path1, path2);
        return -1;
    }
    
    int result = 0;
    int c1, c2;
    long pos = 0;
    
    while (1) {
        c1 = fgetc(f1);
        c2 = fgetc(f2);
        
        if (c1 != c2) {
            fprintf(stderr, "    Diff at byte %ld: 0x%02x vs 0x%02x\n", pos, c1, c2);
            result = -1;
            break;
        }
        
        if (c1 == EOF) break;
        pos++;
    }
    
    fclose(f1);
    fclose(f2);
    return result;
}

/* Compare a specific file in both directories */
static int compare_output_file(const char *filename) {
    char path1[512], path2[512];
    snprintf(path1, sizeof(path1), "%s/%s", test_dir_legacy, filename);
    snprintf(path2, sizeof(path2), "%s/%s", test_dir_new, filename);
    
    int result = compare_files(path1, path2);
    if (result != 0) {
        fprintf(stderr, "    MISMATCH: %s\n", filename);
    }
    return result;
}

/* Create minimal feature_arrays for testing */
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

/* Create data_structures with test data */
static data_structures* create_test_hashes(void) {
    data_structures *hashes = calloc(1, sizeof(data_structures));
    hashes->filtered_hash = kh_init(u32ptr);
    hashes->sequence_umi_hash = kh_init(u64ptr);
    hashes->unique_features_match = kh_init(strptr);
    return hashes;
}

static void free_test_hashes(data_structures *hashes) {
    if (hashes) {
        /* Free feature_counts entries in filtered_hash */
        khint_t k;
        for (k = kh_begin(hashes->filtered_hash); k != kh_end(hashes->filtered_hash); ++k) {
            if (kh_exist(hashes->filtered_hash, k)) {
                feature_counts *fc = (feature_counts*)kh_val(hashes->filtered_hash, k);
                if (fc) {
                    if (fc->counts) kh_destroy(u32u32, fc->counts);
                    free(fc);
                }
            }
        }
        kh_destroy(u32ptr, hashes->filtered_hash);
        kh_destroy(u64ptr, hashes->sequence_umi_hash);
        kh_destroy(strptr, hashes->unique_features_match);
        free(hashes);
    }
}

/* Create test statistics */
static statistics* create_test_stats(void) {
    statistics *stats = calloc(1, sizeof(statistics));
    stats->number_of_reads = 5000;
    stats->total_unmatched_features = 500;
    return stats;
}

/* Add a barcode with feature counts to both hashes and counts_result */
static void add_test_barcode(
    data_structures *hashes,
    pf_counts_result *counts,
    uint32_t barcode_key,
    int *feature_counts_arr,  /* array of counts, one per feature (1-indexed) */
    int n_features
) {
    int ret;
    khint_t k;
    
    /* Add to filtered_hash */
    feature_counts *fc = calloc(1, sizeof(feature_counts));
    *(uint32_t*)fc->sequence_code = barcode_key;
    fc->counts = kh_init(u32u32);
    
    int total = 0;
    for (int i = 1; i <= n_features; i++) {
        if (feature_counts_arr[i] > 0) {
            k = kh_put(u32u32, fc->counts, i, &ret);
            kh_val(fc->counts, k) = feature_counts_arr[i];
            total += feature_counts_arr[i];
        }
    }
    /* Store total at index 0 */
    k = kh_put(u32u32, fc->counts, 0, &ret);
    kh_val(fc->counts, k) = total;
    
    k = kh_put(u32ptr, hashes->filtered_hash, barcode_key, &ret);
    kh_val(hashes->filtered_hash, k) = fc;
    
    /* Add to barcode_to_deduped_hash in counts */
    khash_t(u32u32) *deduped = kh_init(u32u32);
    for (int i = 1; i <= n_features; i++) {
        if (feature_counts_arr[i] > 0) {
            k = kh_put(u32u32, deduped, i, &ret);
            kh_val(deduped, k) = feature_counts_arr[i];
        }
    }
    
    k = kh_put(u32ptr, counts->barcode_to_deduped_hash, barcode_key, &ret);
    kh_val(counts->barcode_to_deduped_hash, k) = deduped;
}

/* ============================================================================
 * Parity Test: Core files (barcodes.txt, features.txt, matrix.mtx)
 * ============================================================================ */
static int test_parity_core_files(void) {
    create_test_dirs();
    int success = 1;
    
    int n_features = 5;
    feature_arrays *features = create_test_features(n_features);
    data_structures *hashes = create_test_hashes();
    statistics *stats = create_test_stats();
    
    /* Create counts result */
    pf_counts_result *counts = pf_counts_result_init(n_features);
    
    /* Add test barcodes with known data */
    /* Barcode 1: features 1=10, 2=5 */
    int fc1[6] = {0, 10, 5, 0, 0, 0};
    add_test_barcode(hashes, counts, 0x11111111, fc1, n_features);
    
    /* Barcode 2: features 2=3, 3=7, 4=2 */
    int fc2[6] = {0, 0, 3, 7, 2, 0};
    add_test_barcode(hashes, counts, 0x22222222, fc2, n_features);
    
    /* Barcode 3: features 1=1, 5=15 */
    int fc3[6] = {0, 1, 0, 0, 0, 15};
    add_test_barcode(hashes, counts, 0x33333333, fc3, n_features);
    
    /* Local arrays for printFeatureCounts */
    int total_deduped[n_features + 1];
    int total_barcoded[n_features + 1];
    vec_u32_t **feature_hist = calloc(n_features + 1, sizeof(vec_u32_t*));
    
    /* --- Run legacy printFeatureCounts --- */
    printFeatureCounts(features, total_deduped, total_barcoded, feature_hist,
                       test_dir_legacy, hashes, stats, 
                       counts->barcode_to_deduped_hash, NULL);
    
    /* Clean up feature_hist for reuse */
    for (int i = 0; i <= n_features; i++) {
        if (feature_hist[i]) {
            vec_u32_destroy(feature_hist[i]);
            feature_hist[i] = NULL;
        }
    }
    
    /* --- Run new mex_write_all --- */
    mex_writer_config config = {
        .output_dir = test_dir_new,
        .features = features,
        .hashes = hashes,
        .stats = stats,
        .filtered_barcodes_hash = NULL,
        .min_heatmap_counts = 0
    };
    
    mex_write_all(&config, counts);
    
    /* --- Compare core output files --- */
    printf("\n    Comparing core files:\n");
    
    if (compare_output_file("barcodes.txt") != 0) success = 0;
    else printf("    barcodes.txt: OK\n");
    
    if (compare_output_file("features.txt") != 0) success = 0;
    else printf("    features.txt: OK\n");
    
    if (compare_output_file("matrix.mtx") != 0) success = 0;
    else printf("    matrix.mtx: OK\n");
    
    if (compare_output_file("stats.txt") != 0) success = 0;
    else printf("    stats.txt: OK\n");
    
    if (compare_output_file("feature_per_cell.csv") != 0) success = 0;
    else printf("    feature_per_cell.csv: OK\n");
    
    /* Cleanup */
    free(feature_hist);
    pf_counts_result_free(counts);
    free_test_hashes(hashes);
    free_test_features(features);
    free(stats);
    cleanup_test_dirs();
    
    return success;
}

/* ============================================================================
 * Parity Test: Filtered output
 * ============================================================================ */
static int test_parity_filtered_output(void) {
    create_test_dirs();
    int success = 1;
    
    int n_features = 3;
    feature_arrays *features = create_test_features(n_features);
    data_structures *hashes = create_test_hashes();
    statistics *stats = create_test_stats();
    pf_counts_result *counts = pf_counts_result_init(n_features);
    
    /* Add test barcodes */
    int fc1[4] = {0, 5, 3, 0};
    add_test_barcode(hashes, counts, 0xAAAAAAAA, fc1, n_features);
    
    int fc2[4] = {0, 0, 2, 8};
    add_test_barcode(hashes, counts, 0xBBBBBBBB, fc2, n_features);
    
    /* Create filter that only includes first barcode */
    khash_t(strptr) *filter = kh_init(strptr);
    char barcode_str[17];
    code2string((unsigned char*)&(uint32_t){0xAAAAAAAA}, barcode_str, barcode_code_length);
    int ret;
    khint_t k = kh_put(strptr, filter, strdup(barcode_str), &ret);
    kh_val(filter, k) = NULL;
    
    /* Local arrays for printFeatureCounts */
    int total_deduped[n_features + 1];
    int total_barcoded[n_features + 1];
    vec_u32_t **feature_hist = calloc(n_features + 1, sizeof(vec_u32_t*));
    
    /* --- Run legacy printFeatureCounts with filter --- */
    printFeatureCounts(features, total_deduped, total_barcoded, feature_hist,
                       test_dir_legacy, hashes, stats,
                       counts->barcode_to_deduped_hash, filter);
    
    /* Clean up feature_hist */
    for (int i = 0; i <= n_features; i++) {
        if (feature_hist[i]) {
            vec_u32_destroy(feature_hist[i]);
            feature_hist[i] = NULL;
        }
    }
    
    /* --- Run new mex_write_all with filter --- */
    mex_writer_config config = {
        .output_dir = test_dir_new,
        .features = features,
        .hashes = hashes,
        .stats = stats,
        .filtered_barcodes_hash = filter,
        .min_heatmap_counts = 0
    };
    
    mex_write_all(&config, counts);
    
    /* --- Compare filtered output files --- */
    printf("\n    Comparing filtered files:\n");
    
    char path1[512], path2[512];
    
    snprintf(path1, sizeof(path1), "%s/filtered/barcodes.txt", test_dir_legacy);
    snprintf(path2, sizeof(path2), "%s/filtered/barcodes.txt", test_dir_new);
    if (compare_files(path1, path2) != 0) success = 0;
    else printf("    filtered/barcodes.txt: OK\n");
    
    snprintf(path1, sizeof(path1), "%s/filtered/features.txt", test_dir_legacy);
    snprintf(path2, sizeof(path2), "%s/filtered/features.txt", test_dir_new);
    if (compare_files(path1, path2) != 0) success = 0;
    else printf("    filtered/features.txt: OK\n");
    
    snprintf(path1, sizeof(path1), "%s/filtered/matrix.mtx", test_dir_legacy);
    snprintf(path2, sizeof(path2), "%s/filtered/matrix.mtx", test_dir_new);
    if (compare_files(path1, path2) != 0) success = 0;
    else printf("    filtered/matrix.mtx: OK\n");
    
    /* Cleanup filter hash */
    for (k = kh_begin(filter); k != kh_end(filter); ++k) {
        if (kh_exist(filter, k)) {
            free((char*)kh_key(filter, k));
        }
    }
    kh_destroy(strptr, filter);
    
    free(feature_hist);
    pf_counts_result_free(counts);
    free_test_hashes(hashes);
    free_test_features(features);
    free(stats);
    cleanup_test_dirs();
    
    return success;
}

/* ============================================================================
 * Main
 * ============================================================================ */
int main(int argc, char **argv) {
    (void)argc;
    (void)argv;
    
    printf("Running mex_writer parity tests (vs printFeatureCounts)...\n\n");
    
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
    
    /* Run parity tests */
    printf("  Testing parity_core_files...");
    tests_run++;
    if (test_parity_core_files()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("  Testing parity_filtered_output...");
    tests_run++;
    if (test_parity_filtered_output()) {
        printf(" PASS\n");
        tests_passed++;
    } else {
        printf(" FAIL\n");
    }
    
    printf("\n%d/%d parity tests passed\n", tests_passed, tests_run);
    
    return (tests_passed == tests_run) ? 0 : 1;
}
