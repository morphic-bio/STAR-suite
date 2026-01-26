/**
 * @file test_pf_api.c
 * @brief Test program for pf_api library functions
 * 
 * Usage: test_pf_api <whitelist> <feature_csv> <fastq_dir> <output_dir>
 */

#include "../include/pf_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_usage(const char *prog) {
    fprintf(stderr, "Usage: %s <whitelist> <feature_csv> <fastq_dir> <output_dir>\n", prog);
    fprintf(stderr, "\nTests the pf_api library by processing FASTQs.\n");
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        print_usage(argv[0]);
        return 1;
    }
    
    const char *whitelist_path = argv[1];
    const char *feature_csv = argv[2];
    const char *fastq_dir = argv[3];
    const char *output_dir = argv[4];
    
    printf("=== process_features API Test ===\n");
    printf("Version: %s\n\n", pf_version());
    
    /* Create config with defaults */
    printf("Creating config...\n");
    pf_config *config = pf_config_create();
    if (!config) {
        fprintf(stderr, "Failed to create config\n");
        return 1;
    }
    
    /* Initialize context */
    printf("Initializing context...\n");
    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);  /* Config is cloned, original can be freed */
    
    if (!ctx) {
        fprintf(stderr, "Failed to initialize context\n");
        return 1;
    }
    
    /* Load feature reference */
    printf("Loading feature reference: %s\n", feature_csv);
    pf_error err = pf_load_feature_ref(ctx, feature_csv);
    if (err != PF_OK) {
        fprintf(stderr, "Failed to load features: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return 1;
    }
    printf("  Loaded %d features\n", pf_get_num_features(ctx));
    
    /* Print first few features */
    int nf = pf_get_num_features(ctx);
    int show = nf < 5 ? nf : 5;
    for (int i = 0; i < show; i++) {
        printf("    [%d] %s: %s\n", i, pf_get_feature_name(ctx, i), pf_get_feature_sequence(ctx, i));
    }
    if (nf > show) {
        printf("    ... and %d more\n", nf - show);
    }
    
    /* Load whitelist */
    printf("\nLoading whitelist: %s\n", whitelist_path);
    err = pf_load_whitelist(ctx, whitelist_path);
    if (err != PF_OK) {
        fprintf(stderr, "Failed to load whitelist: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return 1;
    }
    printf("  Whitelist loaded successfully\n");
    
    /* Process FASTQs */
    printf("\nProcessing FASTQs from: %s\n", fastq_dir);
    printf("Output directory: %s\n", output_dir);
    
    pf_stats stats;
    err = pf_process_fastq_dir(ctx, fastq_dir, output_dir, &stats);
    if (err != PF_OK) {
        fprintf(stderr, "Failed to process FASTQs: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return 1;
    }
    
    /* Print stats */
    printf("\n=== Processing Complete ===\n");
    printf("  Time: %.2f seconds\n", stats.processing_time_sec);
    printf("  Features: %zu\n", stats.total_features);
    
    /* Cleanup */
    pf_destroy(ctx);
    
    printf("\nTest PASSED\n");
    return 0;
}
