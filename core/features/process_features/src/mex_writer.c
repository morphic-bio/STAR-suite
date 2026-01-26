/**
 * @file mex_writer.c
 * @brief Implementation of MEX format output writer
 * 
 * This module extracts output logic from printFeatureCounts() to provide
 * a modular output stage. All output formatting is preserved exactly to
 * ensure byte-identical results.
 */

#include "../include/mex_writer.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/heatmap.h"
#include "../include/plot_histogram.h"
#include "../include/io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

/* Forward declarations for internal helpers */
static int mex_write_barcodes_matrix(
    const char *output_dir,
    feature_arrays *features,
    data_structures *hashes,
    pf_counts_result *counts,
    khash_t(strptr) *filtered_barcodes_hash,
    int *out_total_barcodes,
    int *out_total_excluded,
    int *out_total_deduped,
    int *out_total_raw
);

static int mex_write_features_file(
    const char *output_dir,
    feature_arrays *features
);

static void update_feature_hist_internal(
    vec_u32_t **feature_hist,
    int n_features,
    uint32_t feature_index,
    uint32_t count
);

static void build_cumulative_feature_hist_internal(
    feature_arrays *features,
    vec_u32_t **feature_hist
);

static int build_and_write_heatmaps(
    const char *output_dir,
    feature_arrays *features,
    pf_counts_result *counts,
    khash_t(strptr) *filtered_barcodes_hash,
    int min_heatmap_counts
);


int mex_write_all(
    const mex_writer_config *config,
    pf_counts_result *counts
) {
    if (!config || !counts || !config->output_dir || !config->features) {
        return -1;
    }
    
    char output_directory[FILENAME_LENGTH];
    strcpy(output_directory, config->output_dir);
    if (config->filtered_barcodes_hash) {
        snprintf(output_directory, FILENAME_LENGTH, "%s/filtered", config->output_dir);
    }
    
    /* Ensure output directory exists */
    mkdir_p(output_directory);
    
    int total_barcodes = 0;
    int total_excluded = 0;
    int total_deduped = 0;
    int total_raw = 0;
    
    /* Write barcodes.txt and matrix.mtx */
    if (mex_write_barcodes_matrix(
            output_directory,
            config->features,
            config->hashes,
            counts,
            config->filtered_barcodes_hash,
            &total_barcodes,
            &total_excluded,
            &total_deduped,
            &total_raw) != 0) {
        return -1;
    }
    
    /* Write features.txt */
    if (mex_write_features_file(output_directory, config->features) != 0) {
        return -1;
    }
    
    /* Write stats.txt */
    if (mex_write_stats(
            output_directory,
            config->stats,
            config->hashes,
            total_raw,
            total_deduped,
            total_barcodes,
            total_excluded) != 0) {
        return -1;
    }
    
    /* Write feature_per_cell.csv */
    if (mex_write_feature_per_cell(
            output_directory,
            counts,
            config->filtered_barcodes_hash) != 0) {
        return -1;
    }
    
    /* Build cumulative histogram and write histogram HTMLs */
    build_cumulative_feature_hist_internal(config->features, counts->feature_hist);
    
    /* Write histograms if data exists */
    vec_u32_t *richness_hist = NULL;
    /* Build richness histogram from counts */
    richness_hist = vec_u32_init();
    khint_t k;
    for (k = kh_begin(counts->barcode_to_deduped_hash); 
         k != kh_end(counts->barcode_to_deduped_hash); ++k) {
        if (!kh_exist(counts->barcode_to_deduped_hash, k)) continue;
        
        uint32_t bcode_key = kh_key(counts->barcode_to_deduped_hash, k);
        khash_t(u32u32) *deduped_hash = (khash_t(u32u32)*)kh_val(counts->barcode_to_deduped_hash, k);
        
        if (config->filtered_barcodes_hash) {
            char barcode[barcode_length + 1];
            code2string((unsigned char *)&bcode_key, barcode, barcode_code_length);
            if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
            khint_t kf = kh_get(strptr, config->filtered_barcodes_hash, barcode);
            if (kf == kh_end(config->filtered_barcodes_hash)) continue;
        }
        
        size_t n_features_in_bc = kh_size(deduped_hash);
        if (n_features_in_bc > 0) {
            if (vec_u32_size(richness_hist) <= n_features_in_bc) {
                vec_u32_set_size(richness_hist, n_features_in_bc + 1);
            }
            vec_u32_inc(richness_hist, n_features_in_bc);
        }
    }
    
    if (richness_hist && vec_u32_size(richness_hist) > 0) {
        plot_simple_histogram(output_directory, "feature_richness_histogram.html", 
                             "Feature Richness per Barcode", "Number of Distinct Features", 
                             "Frequency (Barcodes)", richness_hist);
    }
    vec_u32_destroy(richness_hist);
    
    if (counts->feature_hist[0]) {
        plot_simple_histogram(output_directory, "feature_multiplicity_histogram.html",
                             "Feature Multiplicity", "Deduplicated UMI Count",
                             "Frequency", counts->feature_hist[0]);
    }
    
    /* Generate heatmaps */
    if (mex_write_heatmaps(
            output_directory,
            config->features,
            counts,
            config->filtered_barcodes_hash,
            config->min_heatmap_counts) != 0) {
        /* Heatmap failures are non-fatal - log but continue */
        fprintf(stderr, "Warning: heatmap generation failed\n");
    }
    
    return 0;
}

int mex_write_core(
    const mex_writer_config *config,
    pf_counts_result *counts
) {
    if (!config || !counts || !config->output_dir || !config->features) {
        return -1;
    }
    
    char output_directory[FILENAME_LENGTH];
    strcpy(output_directory, config->output_dir);
    if (config->filtered_barcodes_hash) {
        snprintf(output_directory, FILENAME_LENGTH, "%s/filtered", config->output_dir);
    }
    
    mkdir_p(output_directory);
    
    int total_barcodes = 0, total_excluded = 0, total_deduped = 0, total_raw = 0;
    
    if (mex_write_barcodes_matrix(
            output_directory,
            config->features,
            config->hashes,
            counts,
            config->filtered_barcodes_hash,
            &total_barcodes,
            &total_excluded,
            &total_deduped,
            &total_raw) != 0) {
        return -1;
    }
    
    if (mex_write_features_file(output_directory, config->features) != 0) {
        return -1;
    }
    
    return 0;
}

int mex_write_stats(
    const char *output_dir,
    const statistics *stats,
    const data_structures *hashes,
    int total_raw_counts,
    int total_deduped_counts,
    int total_barcodes,
    int total_excluded_barcodes
) {
    char stats_file[FILENAME_LENGTH];
    snprintf(stats_file, FILENAME_LENGTH, "%s/stats.txt", output_dir);
    
    FILE *statsfp = fopen(stats_file, "w");
    if (!statsfp) {
        fprintf(stderr, "Error opening stats file %s\n", stats_file);
        return -1;
    }
    
    /* Match exact output format from original printFeatureCounts */
    fprintf(stderr, "Total feature counts %d\n", total_raw_counts);
    fprintf(stderr, "Total deduped feature counts %d\n", total_deduped_counts);
    fprintf(stderr, "Total barcodes %d\n", total_barcodes);
    fprintf(stderr, "Total excluded barcodes %d\n", total_excluded_barcodes);
    fprintf(stderr, "Total unique barcode UMIs %d\n", (int)kh_size(hashes->sequence_umi_hash));
    fprintf(stderr, "Total whitelisted barcodes %d\n", (int)kh_size(hashes->filtered_hash));
    fprintf(stderr, "Total feature counts %d total_unmatched_reads %ld\n", 
            total_raw_counts, stats->total_unmatched_features);
    fprintf(stderr, "Percentage reads assigned to barcode %.4f\n", 
            100.0 * (total_raw_counts / (double)(total_raw_counts + stats->total_unmatched_features)));
    
    fprintf(statsfp, "Total feature counts %d\n", total_raw_counts);
    fprintf(statsfp, "Total deduped feature counts %d\n", total_deduped_counts);
    fprintf(statsfp, "Total barcodes %d\n", total_barcodes);
    fprintf(statsfp, "Total excluded barcodes %d\n", total_excluded_barcodes);
    fprintf(statsfp, "Total unique barcode UMIs %d\n", (int)kh_size(hashes->sequence_umi_hash));
    fprintf(statsfp, "Total whitelisted barcodes %d\n", (int)kh_size(hashes->filtered_hash));
    fprintf(statsfp, "Total_unmatched_reads %ld\n", stats->total_unmatched_features);
    fprintf(statsfp, "Percentage reads assigned to barcode %.4f\n",
            100.0 * (total_raw_counts / (double)(total_raw_counts + stats->total_unmatched_features)));
    
    fclose(statsfp);
    return 0;
}

int mex_write_feature_per_cell(
    const char *output_dir,
    pf_counts_result *counts,
    khash_t(strptr) *filtered_barcodes_hash
) {
    char fpc_file[FILENAME_LENGTH];
    snprintf(fpc_file, FILENAME_LENGTH, "%s/feature_per_cell.csv", output_dir);
    
    FILE *fpcfp = fopen(fpc_file, "w");
    if (!fpcfp) {
        fprintf(stderr, "Error opening feature_per_cell.csv\n");
        return -1;
    }
    
    /* Header - match exact format */
    fprintf(fpcfp, "barcode,num_features,top_feature_index,total_deduped_umi\n");
    
    khint_t k;
    for (k = kh_begin(counts->barcode_to_deduped_hash); 
         k != kh_end(counts->barcode_to_deduped_hash); ++k) {
        if (!kh_exist(counts->barcode_to_deduped_hash, k)) continue;
        
        uint32_t barcode_key = kh_key(counts->barcode_to_deduped_hash, k);
        khash_t(u32u32) *features_in_barcode = (khash_t(u32u32)*)kh_val(counts->barcode_to_deduped_hash, k);
        
        char barcode[barcode_length + 1];
        code2string((unsigned char *)&barcode_key, barcode, barcode_code_length);
        if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
        
        if (filtered_barcodes_hash) {
            khint_t kf = kh_get(strptr, filtered_barcodes_hash, barcode);
            if (kf == kh_end(filtered_barcodes_hash)) continue;
        }
        
        size_t n_feats_in_bc = kh_size(features_in_barcode);
        uint32_t total_umi = 0;
        uint32_t top_count = 0;
        uint32_t top_feat = 0;
        int tie = 0;
        
        khint_t kfeat;
        for (kfeat = kh_begin(features_in_barcode); 
             kfeat != kh_end(features_in_barcode); ++kfeat) {
            if (!kh_exist(features_in_barcode, kfeat)) continue;
            uint32_t feat_idx = kh_key(features_in_barcode, kfeat);
            uint32_t feat_cnt = kh_val(features_in_barcode, kfeat);
            total_umi += feat_cnt;
            
            if (feat_cnt > top_count) {
                top_count = feat_cnt;
                top_feat = feat_idx;
                tie = 0;
            } else if (feat_cnt == top_count) {
                tie = 1;
            }
        }
        if (tie) top_feat = 0;
        
        fprintf(fpcfp, "%s,%zu,%u,%u\n", barcode, n_feats_in_bc, top_feat, total_umi);
    }
    
    fclose(fpcfp);
    return 0;
}

int mex_write_heatmaps(
    const char *output_dir,
    feature_arrays *features,
    pf_counts_result *counts,
    khash_t(strptr) *filtered_barcodes_hash,
    int min_heatmap_counts
) {
    return build_and_write_heatmaps(output_dir, features, counts, 
                                    filtered_barcodes_hash, min_heatmap_counts);
}

/* ============================================================================
 * Internal helper functions
 * ============================================================================ */

static int mex_write_barcodes_matrix(
    const char *output_dir,
    feature_arrays *features,
    data_structures *hashes,
    pf_counts_result *counts,
    khash_t(strptr) *filtered_barcodes_hash,
    int *out_total_barcodes,
    int *out_total_excluded,
    int *out_total_deduped,
    int *out_total_raw
) {
    char barcodes_file[FILENAME_LENGTH];
    char matrix_file[FILENAME_LENGTH];
    
    snprintf(barcodes_file, FILENAME_LENGTH, "%s/barcodes.txt", output_dir);
    snprintf(matrix_file, FILENAME_LENGTH, "%s/matrix.mtx", output_dir);
    
    FILE *barcodesfp = fopen(barcodes_file, "w");
    if (!barcodesfp) {
        fprintf(stderr, "Error opening barcodes file %s\n", barcodes_file);
        return -1;
    }
    
    FILE *matrixfp = fopen(matrix_file, "w");
    if (!matrixfp) {
        fprintf(stderr, "Error opening matrix file %s\n", matrix_file);
        fclose(barcodesfp);
        return -1;
    }
    
    /* First pass: count barcodes and non-zeros for matrix header */
    size_t number_of_barcode_entries = 0;
    size_t number_of_features_seen = 0;
    
    khint_t k;
    for (k = kh_begin(counts->barcode_to_deduped_hash); 
         k != kh_end(counts->barcode_to_deduped_hash); ++k) {
        if (!kh_exist(counts->barcode_to_deduped_hash, k)) continue;
        
        uint32_t bcode_key = kh_key(counts->barcode_to_deduped_hash, k);
        khash_t(u32u32) *deduped_hash = (khash_t(u32u32)*)kh_val(counts->barcode_to_deduped_hash, k);
        
        if (filtered_barcodes_hash) {
            char barcode[barcode_length + 1];
            code2string((unsigned char *)&bcode_key, barcode, barcode_code_length);
            if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
            khint_t kf = kh_get(strptr, filtered_barcodes_hash, barcode);
            if (kf == kh_end(filtered_barcodes_hash)) continue;
        }
        
        number_of_barcode_entries++;
        number_of_features_seen += kh_size(deduped_hash);
    }
    
    /* Write Matrix Market header - exact format match */
    fprintf(matrixfp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(matrixfp, "%%metadata_json: {\"software_version\": \"assignBarcodes-0.1\", \"format_version\": 1}\n");
    fprintf(matrixfp, "%d %ld %ld\n", features->number_of_features, 
            number_of_barcode_entries, number_of_features_seen);
    
    /* Second pass: write barcodes and matrix entries */
    int total_barcodes = 0;
    int total_excluded = 0;
    int total_deduped = 0;
    int total_raw = 0;
    int skipped_barcodes = 0;
    int processed_barcodes = 0;
    
    /* IMPORTANT: Clear counts arrays and histograms before populating.
     * This ensures correct behavior when called multiple times (raw + filtered). */
    memset(counts->total_deduped_counts, 0, (counts->n_features + 1) * sizeof(int));
    memset(counts->total_barcoded_counts, 0, (counts->n_features + 1) * sizeof(int));
    for (int i = 0; i <= counts->n_features; i++) {
        if (counts->feature_hist[i]) {
            vec_u32_destroy(counts->feature_hist[i]);
            counts->feature_hist[i] = NULL;
        }
    }
    
    /* Iterate through filtered_hash to match original iteration order */
    for (k = kh_begin(hashes->filtered_hash); k != kh_end(hashes->filtered_hash); ++k) {
        if (!kh_exist(hashes->filtered_hash, k)) continue;
        
        feature_counts *entry = (feature_counts*)kh_val(hashes->filtered_hash, k);
        
        /* Accumulate raw counts */
        khint_t k0 = kh_get(u32u32, entry->counts, 0);
        uint32_t total = (k0 != kh_end(entry->counts)) ? kh_val(entry->counts, k0) : 0;
        total_raw += total;
        
        /* Populate barcoded_counts (raw counts per feature) */
        khint_t kraw;
        for (kraw = kh_begin(entry->counts); kraw != kh_end(entry->counts); ++kraw) {
            if (!kh_exist(entry->counts, kraw)) continue;
            uint32_t feature_index = kh_key(entry->counts, kraw);
            if (feature_index > 0) {
                counts->total_barcoded_counts[feature_index - 1] += kh_val(entry->counts, kraw);
            }
        }
        
        /* Check if barcode has deduped counts */
        uint32_t bcode_key = *(uint32_t*)entry->sequence_code;
        khint_t kdedup = kh_get(u32ptr, counts->barcode_to_deduped_hash, bcode_key);
        if (kdedup != kh_end(counts->barcode_to_deduped_hash)) {
            khash_t(u32u32) *deduped_hash = (khash_t(u32u32)*)kh_val(counts->barcode_to_deduped_hash, kdedup);
            if (deduped_hash && kh_size(deduped_hash) > 0) {
                char barcode[barcode_length + 1];
                code2string(entry->sequence_code, barcode, barcode_code_length);
                if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
                
                if (filtered_barcodes_hash) {
                    khint_t kf = kh_get(strptr, filtered_barcodes_hash, barcode);
                    if (kf == kh_end(filtered_barcodes_hash)) {
                        total_excluded++;
                        skipped_barcodes++;
                        continue;
                    }
                }
                processed_barcodes++;
                
                fprintf(barcodesfp, "%s\n", barcode);
                
                /* Write matrix entries */
                khint_t kdedup_iter;
                for (kdedup_iter = kh_begin(deduped_hash); 
                     kdedup_iter != kh_end(deduped_hash); ++kdedup_iter) {
                    if (!kh_exist(deduped_hash, kdedup_iter)) continue;
                    uint32_t f_idx = kh_key(deduped_hash, kdedup_iter);
                    uint32_t deduped_count = kh_val(deduped_hash, kdedup_iter);
                    
                    counts->total_deduped_counts[f_idx - 1] += deduped_count;
                    update_feature_hist_internal(counts->feature_hist, counts->n_features, 
                                                 f_idx, deduped_count);
                    
                    fprintf(matrixfp, "%d %d %d\n", f_idx, total_barcodes + 1, deduped_count);
                    total_deduped += deduped_count;
                }
                total_barcodes++;
            }
        }
    }
    
    /* Log barcode counts to match original printFeatureCounts() output */
    fprintf(stderr, "Skipped barcodes: %d\n", skipped_barcodes);
    fprintf(stderr, "Processed barcodes: %d\n", processed_barcodes);
    
    fclose(barcodesfp);
    fclose(matrixfp);
    fprintf(stderr, "closing matrix file\n");
    fprintf(stderr, "writing stats file\n");
    
    *out_total_barcodes = total_barcodes;
    *out_total_excluded = total_excluded;
    *out_total_deduped = total_deduped;
    *out_total_raw = total_raw;
    
    return 0;
}

static int mex_write_features_file(
    const char *output_dir,
    feature_arrays *features
) {
    char features_file[FILENAME_LENGTH];
    snprintf(features_file, FILENAME_LENGTH, "%s/features.txt", output_dir);
    
    FILE *featuresfp = fopen(features_file, "w");
    if (!featuresfp) {
        fprintf(stderr, "Error opening features file %s\n", features_file);
        return -1;
    }
    
    for (int idx = 0; idx < features->number_of_features; idx++) {
        fprintf(featuresfp, "%s\n", features->feature_names[idx]);
    }
    
    fclose(featuresfp);
    return 0;
}

static void update_feature_hist_internal(
    vec_u32_t **feature_hist,
    int n_features,
    uint32_t feature_index,
    uint32_t count
) {
    if (feature_index > (uint32_t)n_features || feature_index == 0) return;
    
    if (!feature_hist[feature_index]) {
        feature_hist[feature_index] = vec_u32_init();
    }
    
    vec_u32_t *h = feature_hist[feature_index];
    if (vec_u32_size(h) <= count) {
        vec_u32_set_size(h, count + 1);
    }
    vec_u32_inc(h, count);
}

static void build_cumulative_feature_hist_internal(
    feature_arrays *features,
    vec_u32_t **feature_hist
) {
    /* Build cumulative histogram at index 0 */
    /* IMPORTANT: Clear any existing cumulative histogram first to ensure
     * correct behavior when called multiple times (raw + filtered passes). */
    if (feature_hist[0]) {
        vec_u32_destroy(feature_hist[0]);
    }
    feature_hist[0] = vec_u32_init();
    
    for (int i = 1; i <= features->number_of_features; i++) {
        if (!feature_hist[i]) continue;
        
        for (size_t j = 0; j < vec_u32_size(feature_hist[i]); j++) {
            uint32_t freq = vec_u32_get(feature_hist[i], j);
            if (freq > 0) {
                if (vec_u32_size(feature_hist[0]) <= j) {
                    vec_u32_set_size(feature_hist[0], j + 1);
                }
                /* Add to cumulative */
                vec_u32_set(feature_hist[0], j, vec_u32_get(feature_hist[0], j) + freq);
            }
        }
    }
}

static int build_and_write_heatmaps(
    const char *output_dir,
    feature_arrays *features,
    pf_counts_result *counts,
    khash_t(strptr) *filtered_barcodes_hash,
    int min_heatmap_counts
) {
    khint_t k;
    
    /* Build co-expression histogram */
    khash_t(u32ptr) *feature_to_barcodes_hash = kh_init(u32ptr);
    int global_max_cooccur = 0;
    
    for (k = kh_begin(counts->barcode_to_deduped_hash); 
         k != kh_end(counts->barcode_to_deduped_hash); ++k) {
        if (!kh_exist(counts->barcode_to_deduped_hash, k)) continue;
        
        uint32_t barcode_key = kh_key(counts->barcode_to_deduped_hash, k);
        khash_t(u32u32) *features_in_barcode = (khash_t(u32u32)*)kh_val(counts->barcode_to_deduped_hash, k);
        
        char barcode[barcode_length + 1];
        code2string((unsigned char *)&barcode_key, barcode, barcode_code_length);
        if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
        
        if (filtered_barcodes_hash) {
            khint_t kf = kh_get(strptr, filtered_barcodes_hash, barcode);
            if (kf == kh_end(filtered_barcodes_hash)) continue;
        }
        
        size_t n_feats_in_bc = kh_size(features_in_barcode);
        if (n_feats_in_bc == 0) continue;
        
        if ((int)n_feats_in_bc > global_max_cooccur) {
            global_max_cooccur = n_feats_in_bc;
        }
        
        khint_t kf;
        for (kf = kh_begin(features_in_barcode); kf != kh_end(features_in_barcode); ++kf) {
            if (!kh_exist(features_in_barcode, kf)) continue;
            uint32_t f_idx = kh_key(features_in_barcode, kf);
            if (f_idx <= 0) continue;
            
            khint_t kfb = kh_get(u32ptr, feature_to_barcodes_hash, f_idx);
            barcode_list_t *barcodes_for_feature;
            if (kfb == kh_end(feature_to_barcodes_hash)) {
                barcodes_for_feature = barcode_list_init();
                int ret;
                khint_t kh = kh_put(u32ptr, feature_to_barcodes_hash, f_idx, &ret);
                kh_val(feature_to_barcodes_hash, kh) = barcodes_for_feature;
            } else {
                barcodes_for_feature = (barcode_list_t*)kh_val(feature_to_barcodes_hash, kfb);
            }
            barcode_list_add(barcodes_for_feature, barcode_key);
        }
    }
    
    /* Allocate co-expression histogram matrix */
    int **coexpression_histograms = malloc((features->number_of_features + 1) * sizeof(int*));
    for (int i = 0; i <= features->number_of_features; ++i) {
        coexpression_histograms[i] = calloc(global_max_cooccur + 1, sizeof(int));
    }
    
    /* Populate co-expression histogram */
    for (k = kh_begin(feature_to_barcodes_hash); k != kh_end(feature_to_barcodes_hash); ++k) {
        if (!kh_exist(feature_to_barcodes_hash, k)) continue;
        uint32_t f_idx = kh_key(feature_to_barcodes_hash, k);
        barcode_list_t *barcodes = (barcode_list_t*)kh_val(feature_to_barcodes_hash, k);
        
        for (size_t i = 0; i < barcodes->n; ++i) {
            uint32_t bcode_key = barcodes->keys[i];
            khint_t kb = kh_get(u32ptr, counts->barcode_to_deduped_hash, bcode_key);
            if (kb == kh_end(counts->barcode_to_deduped_hash)) continue;
            khash_t(u32u32) *deduped_hash = (khash_t(u32u32)*)kh_val(counts->barcode_to_deduped_hash, kb);
            int num_feats = kh_size(deduped_hash);
            if (num_feats > 0) {
                coexpression_histograms[f_idx][num_feats]++;
            }
        }
    }
    
    /* Set column 0 to max co-occurrence for each feature */
    for (int i = 1; i <= features->number_of_features; ++i) {
        int row_max = 0;
        for (int j = global_max_cooccur; j >= 1; --j) {
            if (coexpression_histograms[i][j] > 0) {
                row_max = j;
                break;
            }
        }
        coexpression_histograms[i][0] = row_max;
    }
    
    /* Generate co-expression heatmap */
    generate_heatmap(output_dir, features, coexpression_histograms);
    
    /* Clean up feature_to_barcodes_hash */
    for (k = kh_begin(feature_to_barcodes_hash); k != kh_end(feature_to_barcodes_hash); ++k) {
        if (kh_exist(feature_to_barcodes_hash, k)) {
            barcode_list_destroy((barcode_list_t*)kh_val(feature_to_barcodes_hash, k));
        }
    }
    kh_destroy(u32ptr, feature_to_barcodes_hash);
    
    /* Free co-expression histograms */
    for (int i = 0; i <= features->number_of_features; ++i) {
        free(coexpression_histograms[i]);
    }
    free(coexpression_histograms);
    
    /* Build deduped counts heatmap */
    int max_len = 0;
    for (int i = 1; i <= features->number_of_features; i++) {
        if (counts->feature_hist[i] && (int)vec_u32_size(counts->feature_hist[i]) > max_len) {
            max_len = vec_u32_size(counts->feature_hist[i]);
        }
    }
    
    int max_deduped_count = (max_len > 0) ? max_len - 1 : 0;
    
    int **deduped_histograms = malloc((features->number_of_features + 1) * sizeof(int*));
    for (int i = 0; i <= features->number_of_features; ++i) {
        deduped_histograms[i] = calloc(max_deduped_count + 1, sizeof(int));
    }
    
    for (int i = 1; i <= features->number_of_features; ++i) {
        if (counts->feature_hist[i]) {
            for (size_t j = 0; j < vec_u32_size(counts->feature_hist[i]); j++) {
                deduped_histograms[i][j] = vec_u32_get(counts->feature_hist[i], j);
            }
        }
    }
    
    /* Calculate total deduped counts for heatmap threshold */
    int num_features_for_heatmap = 0;
    int temp_total_deduped[features->number_of_features];
    for (int i = 0; i < features->number_of_features; ++i) {
        long feature_total = 0;
        for (int j = 1; j <= max_deduped_count; j++) {
            feature_total += (long)j * deduped_histograms[i + 1][j];
        }
        temp_total_deduped[i] = feature_total;
        if (feature_total > min_heatmap_counts) {
            num_features_for_heatmap++;
        }
    }
    
    if (num_features_for_heatmap >= 2) {
        generate_deduped_heatmap(output_dir, features, deduped_histograms, 
                                 max_deduped_count, temp_total_deduped, min_heatmap_counts);
        
        /* Build and generate richness heatmap */
        feature_to_barcodes_hash = kh_init(u32ptr);
        int global_max_richness = 0;
        
        for (k = kh_begin(counts->barcode_to_deduped_hash); 
             k != kh_end(counts->barcode_to_deduped_hash); ++k) {
            if (!kh_exist(counts->barcode_to_deduped_hash, k)) continue;
            
            uint32_t barcode_key = kh_key(counts->barcode_to_deduped_hash, k);
            khash_t(u32u32) *features_in_barcode = (khash_t(u32u32)*)kh_val(counts->barcode_to_deduped_hash, k);
            
            char barcode[barcode_length + 1];
            code2string((unsigned char *)&barcode_key, barcode, barcode_code_length);
            if (translate_NXT) translate_nxt_inplace(barcode, barcode_length);
            
            if (filtered_barcodes_hash) {
                khint_t kf = kh_get(strptr, filtered_barcodes_hash, barcode);
                if (kf == kh_end(filtered_barcodes_hash)) continue;
            }
            
            size_t n_feats_in_bc = kh_size(features_in_barcode);
            if (n_feats_in_bc == 0) continue;
            if ((int)n_feats_in_bc > global_max_richness) {
                global_max_richness = n_feats_in_bc;
            }
            
            khint_t kf;
            for (kf = kh_begin(features_in_barcode); kf != kh_end(features_in_barcode); ++kf) {
                if (!kh_exist(features_in_barcode, kf)) continue;
                uint32_t f_idx = kh_key(features_in_barcode, kf);
                if (f_idx <= 0) continue;
                
                khint_t kfb = kh_get(u32ptr, feature_to_barcodes_hash, f_idx);
                barcode_list_t *barcodes_for_feature;
                if (kfb == kh_end(feature_to_barcodes_hash)) {
                    barcodes_for_feature = barcode_list_init();
                    int ret;
                    khint_t kh_iter = kh_put(u32ptr, feature_to_barcodes_hash, f_idx, &ret);
                    kh_val(feature_to_barcodes_hash, kh_iter) = barcodes_for_feature;
                } else {
                    barcodes_for_feature = (barcode_list_t*)kh_val(feature_to_barcodes_hash, kfb);
                }
                barcode_list_add(barcodes_for_feature, barcode_key);
            }
        }
        
        int **richness_histograms = malloc((features->number_of_features + 1) * sizeof(int*));
        for (int i = 0; i <= features->number_of_features; ++i) {
            richness_histograms[i] = calloc(global_max_richness + 1, sizeof(int));
        }
        
        for (k = kh_begin(feature_to_barcodes_hash); k != kh_end(feature_to_barcodes_hash); ++k) {
            if (!kh_exist(feature_to_barcodes_hash, k)) continue;
            uint32_t f_idx = kh_key(feature_to_barcodes_hash, k);
            barcode_list_t *barcodes = (barcode_list_t*)kh_val(feature_to_barcodes_hash, k);
            
            for (size_t i = 0; i < barcodes->n; ++i) {
                uint32_t bcode_key = barcodes->keys[i];
                khint_t kb = kh_get(u32ptr, counts->barcode_to_deduped_hash, bcode_key);
                if (kb == kh_end(counts->barcode_to_deduped_hash)) continue;
                khash_t(u32u32) *deduped_hash = (khash_t(u32u32)*)kh_val(counts->barcode_to_deduped_hash, kb);
                int num_feats = kh_size(deduped_hash);
                if (num_feats > 0) {
                    richness_histograms[f_idx][num_feats]++;
                }
            }
        }
        
        for (int i = 1; i <= features->number_of_features; ++i) {
            int row_max = 0;
            for (int j = global_max_richness; j >= 1; --j) {
                if (richness_histograms[i][j] > 0) {
                    row_max = j;
                    break;
                }
            }
            richness_histograms[i][0] = row_max;
        }
        
        generate_heatmap(output_dir, features, richness_histograms);
        
        for (k = kh_begin(feature_to_barcodes_hash); k != kh_end(feature_to_barcodes_hash); ++k) {
            if (kh_exist(feature_to_barcodes_hash, k)) {
                barcode_list_destroy((barcode_list_t*)kh_val(feature_to_barcodes_hash, k));
            }
        }
        kh_destroy(u32ptr, feature_to_barcodes_hash);
        
        for (int i = 0; i <= features->number_of_features; ++i) {
            free(richness_histograms[i]);
        }
        free(richness_histograms);
    } else {
        fprintf(stdout, "Skipping heatmap generation as there are fewer than 2 features with sufficient counts.\n");
    }
    
    for (int i = 0; i <= features->number_of_features; ++i) {
        free(deduped_histograms[i]);
    }
    free(deduped_histograms);
    
    return 0;
}
