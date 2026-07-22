/**
 * @file test_catatac_split_read.c
 * @brief Run CAT-ATAC split-read guide assignment via pf_process_split_fastq_dir.
 *
 * Usage: test_catatac_split_read <fastq_dir> <output_dir> <whitelist> <feature_csv>
 */

#include "../include/pf_api.h"
#include "../include/pf_split_read.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static void catatac_layout(pf_split_read_layout *layout) {
    memset(layout, 0, sizeof(*layout));
    layout->enabled = 1;
    strncpy(layout->barcode_format, "bc:8:23:-", sizeof(layout->barcode_format) - 1);
    layout->barcode_read_idx = 1;
    layout->umi_read_idx = 0;
    layout->umi_start = 0;
    layout->umi_length = 12;
    layout->feature_read_idx = 2;
    layout->capture_read_idx = 0;
    strncpy(layout->capture_sequences,
            "CAAGTTGATAACGGACTAGCC|CAAGTTGTAAACGGACTAGCC",
            sizeof(layout->capture_sequences) - 1);
    layout->capture_max_hamming = 0;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr,
                "Usage: %s <fastq_dir> <output_dir> <whitelist> <feature_csv>\n",
                argv[0]);
        return 1;
    }

    const char *fastq_dir = argv[1];
    const char *output_dir = argv[2];
    const char *whitelist_path = argv[3];
    const char *feature_csv = argv[4];

    pf_config *config = pf_config_create();
    if (!config) {
        return 1;
    }
    pf_config_set_max_hamming_distance(config, 1);
    pf_config_set_limit_search(config, -1);
    pf_config_set_min_counts(config, 0);
    pf_config_set_consumer_threads(config, 1);
    pf_config_set_search_threads(config, 1);
    pf_config_set_skip_heatmaps(config, 1);
    const char *max_reads_env = getenv("CATATAC_GUIDE_MAX_READS");
    if (max_reads_env && max_reads_env[0]) {
        pf_config_set_max_reads(config, atoll(max_reads_env));
    }
    pf_config_set_use_feature_anchor_search(config, 0);
    pf_config_set_require_feature_anchor_match(config, 0);
    pf_split_read_layout layout;
    catatac_layout(&layout);
    pf_config_set_split_read_layout(config, &layout);

    pf_context *ctx = pf_init(config);
    pf_config_destroy(config);
    if (!ctx) {
        return 1;
    }

    if (pf_load_feature_ref(ctx, feature_csv) != PF_OK ||
        pf_load_whitelist(ctx, whitelist_path) != PF_OK) {
        fprintf(stderr, "Failed to load references: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return 1;
    }

    pf_stats stats;
    pf_split_read_metrics metrics;
    pf_error err = pf_process_split_fastq_dir(ctx, fastq_dir, output_dir, &stats, &metrics);
    if (err != PF_OK) {
        fprintf(stderr, "split-read processing failed: %s\n", pf_get_error(ctx));
        pf_destroy(ctx);
        return 1;
    }

    printf("split_read_metrics total_reads=%llu barcode_synth_ok=%llu capture_either=%llu\n",
           metrics.total_reads, metrics.barcode_synth_ok, metrics.capture_either_hits);
    pf_destroy(ctx);
    return 0;
}
