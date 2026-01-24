/**
 * @file pf_api.c
 * @brief Implementation of public API for process_features library
 */

#include "../include/pf_api.h"
#include "../include/common.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/io.h"
#include "../include/utils.h"
#include "../include/memory.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <errno.h>

#define PF_VERSION "1.0.0"
#define PF_ERROR_BUF_SIZE 1024

/* ============================================================================
 * Internal structures
 * ============================================================================ */

struct pf_config {
    int barcode_length;
    int umi_length;
    int max_hamming_distance;
    int stringency;
    int min_counts;
    double min_posterior;
    int feature_offset;
    int barcode_offset;
    int max_barcode_mismatches;
    int max_feature_n;
    int max_barcode_n;
    int max_threads;
    int search_threads;
    int consumer_threads;
    int debug_enabled;
    int reverse_complement_whitelist;
    int limit_search;
    long long max_reads;
    int translate_nxt;
    int use_feature_offset_array;
    int use_feature_anchor_search;
    int require_feature_anchor_match;
    int feature_mode_bootstrap_reads;
    
    /* EmptyDrops control */
    int skip_emptydrops;            /* 1 = skip EmptyDrops entirely */
    int emptydrops_failure_fatal;   /* 1 = treat ED failure as error */
    int expected_cells;             /* 0 = auto-detect */
    int emptydrops_use_fdr;         /* 1 = use FDR gate for tail rescue */
};

struct pf_context {
    pf_config *config;
    feature_arrays *features;
    khash_t(u32ptr) *whitelist_hash;
    unsigned char *whitelist_data;
    khash_t(strptr) *filtered_barcodes_hash;
    int initialized;
    char error_buf[PF_ERROR_BUF_SIZE];
};

/* Global initialization flag */
static int g_pf_global_initialized = 0;

/* ============================================================================
 * Configuration API Implementation
 * ============================================================================ */

pf_config* pf_config_create(void) {
    pf_config *config = calloc(1, sizeof(pf_config));
    if (!config) return NULL;
    
    /* Set defaults matching CLI defaults */
    config->barcode_length = 16;
    config->umi_length = 12;
    config->max_hamming_distance = 1;
    config->stringency = 1;
    config->min_counts = 1;
    config->min_posterior = 0.975;
    config->feature_offset = 0;
    config->barcode_offset = 0;
    config->max_barcode_mismatches = 3;
    config->max_feature_n = 1;
    config->max_barcode_n = 1;
    config->max_threads = 8;
    config->search_threads = 4;
    config->consumer_threads = 1;
    config->debug_enabled = 0;
    config->reverse_complement_whitelist = 0;
    config->limit_search = -1;
    config->max_reads = 0;
    config->translate_nxt = 0;
    config->use_feature_offset_array = 0;
    config->use_feature_anchor_search = 0;
    config->require_feature_anchor_match = 0;
    config->feature_mode_bootstrap_reads = 0;
    
    /* EmptyDrops control defaults */
    config->skip_emptydrops = 0;
    config->emptydrops_failure_fatal = 0;
    config->expected_cells = 0;
    config->emptydrops_use_fdr = 0;
    
    return config;
}

void pf_config_destroy(pf_config *config) {
    if (config) {
        free(config);
    }
}

pf_config* pf_config_clone(const pf_config *config) {
    if (!config) return NULL;
    pf_config *clone = malloc(sizeof(pf_config));
    if (!clone) return NULL;
    memcpy(clone, config, sizeof(pf_config));
    return clone;
}

void pf_config_set_barcode_length(pf_config *config, int length) {
    if (config) config->barcode_length = length;
}

void pf_config_set_umi_length(pf_config *config, int length) {
    if (config) config->umi_length = length;
}

void pf_config_set_max_hamming_distance(pf_config *config, int distance) {
    if (config) config->max_hamming_distance = distance;
}

void pf_config_set_stringency(pf_config *config, int stringency) {
    if (config) config->stringency = stringency;
}

void pf_config_set_min_counts(pf_config *config, int min_counts) {
    if (config) config->min_counts = min_counts;
}

void pf_config_set_min_posterior(pf_config *config, double min_posterior) {
    if (config) config->min_posterior = min_posterior;
}

void pf_config_set_feature_offset(pf_config *config, int offset) {
    if (config) config->feature_offset = offset;
}

void pf_config_set_barcode_offset(pf_config *config, int offset) {
    if (config) config->barcode_offset = offset;
}

void pf_config_set_max_barcode_mismatches(pf_config *config, int mismatches) {
    if (config) config->max_barcode_mismatches = mismatches;
}

void pf_config_set_max_feature_n(pf_config *config, int max_n) {
    if (config) config->max_feature_n = max_n;
}

void pf_config_set_max_barcode_n(pf_config *config, int max_n) {
    if (config) config->max_barcode_n = max_n;
}

void pf_config_set_threads(pf_config *config, int threads) {
    if (config) config->max_threads = threads;
}

void pf_config_set_search_threads(pf_config *config, int threads) {
    if (config) config->search_threads = threads;
}

void pf_config_set_consumer_threads(pf_config *config, int threads) {
    if (config) config->consumer_threads = threads;
}

void pf_config_set_debug(pf_config *config, int enable) {
    if (config) config->debug_enabled = enable;
}

void pf_config_set_reverse_complement_whitelist(pf_config *config, int enable) {
    if (config) config->reverse_complement_whitelist = enable;
}

void pf_config_set_limit_search(pf_config *config, int limit) {
    if (config) config->limit_search = limit;
}

void pf_config_set_max_reads(pf_config *config, long long max_reads) {
    if (config) config->max_reads = max_reads;
}

void pf_config_set_translate_nxt(pf_config *config, int enable) {
    if (config) config->translate_nxt = enable;
}
void pf_config_set_use_feature_offset_array(pf_config *config, int enable) {
    if (config) config->use_feature_offset_array = enable;
}
void pf_config_set_use_feature_anchor_search(pf_config *config, int enable) {
    if (config) config->use_feature_anchor_search = enable;
}

void pf_config_set_require_feature_anchor_match(pf_config *config, int enable) {
    if (config) config->require_feature_anchor_match = enable;
}

void pf_config_set_feature_mode_bootstrap_reads(pf_config *config, int n_reads) {
    if (config) config->feature_mode_bootstrap_reads = n_reads;
}

void pf_config_set_skip_emptydrops(pf_config *config, int enable) {
    if (config) config->skip_emptydrops = enable;
}

void pf_config_set_emptydrops_failure_fatal(pf_config *config, int enable) {
    if (config) config->emptydrops_failure_fatal = enable;
}

void pf_config_set_expected_cells(pf_config *config, int n_cells) {
    if (config) config->expected_cells = n_cells;
}

void pf_config_set_emptydrops_use_fdr(pf_config *config, int enable) {
    if (config) config->emptydrops_use_fdr = enable;
}


/* ============================================================================
 * Global Initialization
 * ============================================================================ */

void pf_global_init(void) {
    if (g_pf_global_initialized) return;
    
    initseq2Code();
    initcode2seq();
    initdiff2hamming(diff2Hamming);
    initialize_complement();
    
    /* Initialize feature code hash if not already done */
    if (!feature_code_hash) {
        feature_code_hash = kh_init(codeu32);
    }
    
    g_pf_global_initialized = 1;
}

/* ============================================================================
 * Context Lifecycle Implementation
 * ============================================================================ */

pf_context* pf_init(const pf_config *config) {
    if (!config) return NULL;
    
    pf_context *ctx = calloc(1, sizeof(pf_context));
    if (!ctx) return NULL;
    
    ctx->config = pf_config_clone(config);
    if (!ctx->config) {
        free(ctx);
        return NULL;
    }
    
    /* Apply config to globals */
    barcode_length = ctx->config->barcode_length;
    barcode_code_length = (barcode_length + 3) / 4;
    umi_length = ctx->config->umi_length;
    umi_code_length = (umi_length + 3) / 4;
    max_barcode_mismatches = ctx->config->max_barcode_mismatches;
    max_feature_n = ctx->config->max_feature_n;
    max_barcode_n = ctx->config->max_barcode_n;
    max_reads = ctx->config->max_reads;
    limit_search = ctx->config->limit_search;
    debug = ctx->config->debug_enabled;
    translate_NXT = ctx->config->translate_nxt;
    use_feature_offset_array = ctx->config->use_feature_offset_array;
    use_feature_anchor_search = ctx->config->use_feature_anchor_search;
    require_feature_anchor_match = ctx->config->require_feature_anchor_match;
    if (require_feature_anchor_match) {
        use_feature_anchor_search = 1;
    }
    feature_mode_bootstrap_reads = ctx->config->feature_mode_bootstrap_reads;
    feature_mode_reads_seen = 0;
    feature_mode_bootstrap_done = 0;
    
    /* Global initialization */
    pf_global_init();
    initialize_unit_sizes();
    
    ctx->initialized = 1;
    ctx->error_buf[0] = '\0';
    
    return ctx;
}

void pf_destroy(pf_context *ctx) {
    if (!ctx) return;
    
    if (ctx->config) {
        pf_config_destroy(ctx->config);
    }
    
    if (ctx->features) {
        if (feature_offsets == ctx->features->feature_offsets) {
            feature_offsets = NULL;
            feature_offsets_count = 0;
        }
        if (feature_anchors == ctx->features->feature_anchors) {
            feature_anchors = NULL;
            feature_anchor_lengths = NULL;
            feature_anchor_count = 0;
        }
        free_feature_arrays(ctx->features);
    }
    if (feature_mode_hist) {
        free(feature_mode_hist);
        feature_mode_hist = NULL;
    }
    if (feature_mode_offsets) {
        free(feature_mode_offsets);
        feature_mode_offsets = NULL;
    }
    
    if (ctx->whitelist_hash) {
        kh_destroy(u32ptr, ctx->whitelist_hash);
    }
    
    if (ctx->whitelist_data) {
        free(ctx->whitelist_data);
    }
    
    if (ctx->filtered_barcodes_hash) {
        free_strptr_hash(ctx->filtered_barcodes_hash);
    }
    
    free(ctx);
}

const char* pf_get_error(pf_context *ctx) {
    if (!ctx || ctx->error_buf[0] == '\0') return NULL;
    return ctx->error_buf;
}

/* ============================================================================
 * Reference Loading Implementation
 * ============================================================================ */

pf_error pf_load_feature_ref(pf_context *ctx, const char *feature_csv) {
    if (!ctx || !feature_csv) return PF_ERR_INVALID_ARG;
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;
    
    if (!file_exists(feature_csv)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, 
                 "Feature reference file not found: %s", feature_csv);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Free existing features if any */
    if (ctx->features) {
        if (feature_offsets == ctx->features->feature_offsets) {
            feature_offsets = NULL;
            feature_offsets_count = 0;
        }
        if (feature_anchors == ctx->features->feature_anchors) {
            feature_anchors = NULL;
            feature_anchor_lengths = NULL;
            feature_anchor_count = 0;
        }
        free_feature_arrays(ctx->features);
        ctx->features = NULL;
    }
    if (feature_mode_hist) {
        free(feature_mode_hist);
        feature_mode_hist = NULL;
    }
    if (feature_mode_offsets) {
        free(feature_mode_offsets);
        feature_mode_offsets = NULL;
    }
    
    ctx->features = read_features_file(feature_csv);
    if (!ctx->features) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Failed to parse feature reference: %s", feature_csv);
        return PF_ERR_PARSE_ERROR;
    }

    feature_offsets = ctx->features->feature_offsets;
    feature_offsets_count = ctx->features->number_of_features;
    feature_anchors = ctx->features->feature_anchors;
    feature_anchor_lengths = ctx->features->feature_anchor_lengths;
    feature_anchor_count = ctx->features->number_of_features;

    if (feature_mode_bootstrap_reads > 0) {
        const int count = ctx->features->number_of_features;
        feature_mode_offsets = malloc(sizeof(int) * count);
        feature_mode_hist = calloc((size_t)count * (size_t)feature_mode_max_offset, sizeof(unsigned int));
        if (!feature_mode_offsets || !feature_mode_hist) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to allocate feature mode arrays");
            return PF_ERR_OUT_OF_MEMORY;
        }
        for (int i = 0; i < count; i++) {
            feature_mode_offsets[i] = -1;
        }
    }
    
    /* Update globals */
    number_of_features = ctx->features->number_of_features;
    maximum_feature_length = ctx->features->max_length;
    feature_code_length = (maximum_feature_length + 3) / 4;
    
    return PF_OK;
}

pf_error pf_load_whitelist(pf_context *ctx, const char *whitelist_path) {
    if (!ctx || !whitelist_path) return PF_ERR_INVALID_ARG;
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;
    
    if (!file_exists(whitelist_path)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Whitelist file not found: %s", whitelist_path);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Free existing whitelist if any */
    if (ctx->whitelist_hash) {
        kh_destroy(u32ptr, ctx->whitelist_hash);
    }
    if (ctx->whitelist_data) {
        free(ctx->whitelist_data);
        ctx->whitelist_data = NULL;
    }
    
    ctx->whitelist_hash = kh_init(u32ptr);
    ctx->whitelist_data = read_whiteList((char*)whitelist_path, ctx->whitelist_hash,
                                          ctx->config->reverse_complement_whitelist);
    
    if (!ctx->whitelist_data) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Failed to load whitelist: %s", whitelist_path);
        kh_destroy(u32ptr, ctx->whitelist_hash);
        ctx->whitelist_hash = NULL;
        return PF_ERR_IO_ERROR;
    }
    
    /* Update global pointer for existing code */
    whitelist_hash = ctx->whitelist_hash;
    whitelist = ctx->whitelist_data;
    
    return PF_OK;
}

pf_error pf_load_filtered_barcodes(pf_context *ctx, const char *filtered_path) {
    if (!ctx || !filtered_path) return PF_ERR_INVALID_ARG;
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;
    
    if (!file_exists(filtered_path)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Filtered barcodes file not found: %s", filtered_path);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Free existing filtered barcodes if any */
    if (ctx->filtered_barcodes_hash) {
        free_strptr_hash(ctx->filtered_barcodes_hash);
    }
    
    ctx->filtered_barcodes_hash = kh_init(strptr);
    read_barcodes_into_hash(filtered_path, ctx->filtered_barcodes_hash);
    
    return PF_OK;
}

/* ============================================================================
 * Processing Implementation
 * ============================================================================ */

pf_error pf_process_fastq_dir(pf_context *ctx,
                               const char *fastq_dir,
                               const char *output_dir,
                               pf_stats *stats_out) {
    if (!ctx || !fastq_dir || !output_dir) return PF_ERR_INVALID_ARG;
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;
    if (!ctx->features) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Features not loaded");
        return PF_ERR_NOT_INITIALIZED;
    }
    if (!ctx->whitelist_hash) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Whitelist not loaded");
        return PF_ERR_NOT_INITIALIZED;
    }
    
    if (!is_directory(fastq_dir)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "FASTQ directory not found: %s", fastq_dir);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Create output directory if needed */
    struct stat st = {0};
    if (stat(output_dir, &st) == -1) {
        if (mkdir(output_dir, 0755) != 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to create output directory: %s", output_dir);
            return PF_ERR_IO_ERROR;
        }
    }
    
    /* Organize FASTQ files */
    fastq_files_collection fastq_files;
    memset(&fastq_files, 0, sizeof(fastq_files));
    
    char *argv_fake[2];
    argv_fake[0] = "pf_api";
    argv_fake[1] = (char*)fastq_dir;
    
    organize_fastq_files_by_directory(1, 2, argv_fake, 1, NULL, NULL, NULL,
                                       &fastq_files, "_R1_", "_R2_", "_R3_");
    
    if (fastq_files.nsamples == 0) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "No FASTQ files found in directory: %s", fastq_dir);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Process each sample */
    double start_time = get_time_in_seconds();
    int any_failed = 0;
    
    for (int i = 0; i < fastq_files.nsamples; i++) {
        char sample_directory[FILENAME_LENGTH];
        snprintf(sample_directory, sizeof(sample_directory), "%s/%s/",
                 output_dir, fastq_files.sample_names[i]);
        
        /* Create sample output directory */
        if (stat(sample_directory, &st) == -1) {
            mkdir(sample_directory, 0755);
        }
        
        /* Per-sample error tracking */
        int sample_error = 0;
        
        /* Set up sample args */
        sample_args args;
        memset(&args, 0, sizeof(args));
        args.sample_index = i;
        args.directory = sample_directory;
        args.filtered_barcodes_name = NULL;
        args.fastq_files = &fastq_files;
        args.features = ctx->features;
        args.maxHammingDistance = ctx->config->max_hamming_distance;
        args.nThreads = ctx->config->search_threads;
        args.stringency = ctx->config->stringency;
        args.min_counts = ctx->config->min_counts;
        args.barcode_constant_offset = ctx->config->barcode_offset;
        args.feature_constant_offset = ctx->config->feature_offset;
        args.read_buffer_lines = READ_BUFFER_LINES;
        args.average_read_length = AVERAGE_READ_LENGTH;
        args.min_posterior = ctx->config->min_posterior;
        args.consumer_threads_per_set = ctx->config->consumer_threads;
        args.filtered_barcodes_hash = ctx->filtered_barcodes_hash;
        args.min_prediction = 1;
        args.min_heatmap = 0;
        args.demux_nsamples = 1;
        args.sample_barcodes = NULL;
        args.sample_max_hamming = 1;
        args.sample_max_N = 0;
        args.sample_constant_offset = -1;
        args.sample_offset_relative = 0;
        
        /* EmptyDrops control from config */
        args.skip_emptydrops = ctx->config->skip_emptydrops;
        args.emptydrops_failure_fatal = ctx->config->emptydrops_failure_fatal;
        args.expected_cells = ctx->config->expected_cells;
        args.emptydrops_use_fdr = ctx->config->emptydrops_use_fdr;
        args.error_out = &sample_error;
        
        /* Process the sample */
        process_files_in_sample(&args);
        
        /* Check for errors */
        if (sample_error) {
            fprintf(stderr, "[pf_api] Sample %s reported an error\n", fastq_files.sample_names[i]);
            any_failed = 1;
        }
    }
    
    double end_time = get_time_in_seconds();
    
    /* Fill in stats if requested */
    if (stats_out) {
        memset(stats_out, 0, sizeof(pf_stats));
        stats_out->processing_time_sec = end_time - start_time;
        stats_out->total_features = ctx->features->number_of_features;
    }
    
    free_fastq_files_collection(&fastq_files);
    
    if (any_failed) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "One or more samples failed during processing");
        return PF_ERR_IO_ERROR;
    }
    
    return PF_OK;
}

pf_error pf_process_fastqs(pf_context *ctx,
                            const char **barcode_fastqs,
                            const char **feature_fastqs,
                            int n_files,
                            const char *output_dir,
                            const char *sample_name,
                            pf_stats *stats_out) {
    if (!ctx || !barcode_fastqs || !feature_fastqs || n_files <= 0 || !output_dir) {
        return PF_ERR_INVALID_ARG;
    }
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;
    if (!ctx->features) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Features not loaded");
        return PF_ERR_NOT_INITIALIZED;
    }
    if (!ctx->whitelist_hash) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Whitelist not loaded");
        return PF_ERR_NOT_INITIALIZED;
    }
    
    /* Create output directory */
    struct stat st = {0};
    if (stat(output_dir, &st) == -1) {
        if (mkdir(output_dir, 0755) != 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to create output directory: %s", output_dir);
            return PF_ERR_IO_ERROR;
        }
    }
    
    /* Create sample subdirectory */
    const char *sname = sample_name ? sample_name : "sample";
    char sample_directory[FILENAME_LENGTH];
    snprintf(sample_directory, sizeof(sample_directory), "%s/%s/", output_dir, sname);
    
    if (stat(sample_directory, &st) == -1) {
        if (mkdir(sample_directory, 0755) != 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to create sample directory: %s", sample_directory);
            return PF_ERR_IO_ERROR;
        }
    }
    
    /* Build fastq_files_collection manually */
    fastq_files_collection fastq_files;
    memset(&fastq_files, 0, sizeof(fastq_files));
    
    fastq_files.nsamples = 1;
    fastq_files.nbarcode_files = n_files;
    fastq_files.nforward_files = n_files;
    fastq_files.nreverse_files = 0;
    
    /* Allocate arrays */
    fastq_files.barcode_fastq = malloc(n_files * sizeof(char*));
    fastq_files.forward_fastq = malloc(n_files * sizeof(char*));
    fastq_files.sample_names = malloc(sizeof(char*));
    fastq_files.sample_sizes = malloc(sizeof(int));
    fastq_files.sample_offsets = malloc(sizeof(int));
    fastq_files.sorted_index = malloc(sizeof(int));
    
    if (!fastq_files.barcode_fastq || !fastq_files.forward_fastq ||
        !fastq_files.sample_names || !fastq_files.sample_sizes ||
        !fastq_files.sample_offsets || !fastq_files.sorted_index) {
        /* Cleanup on failure */
        free(fastq_files.barcode_fastq);
        free(fastq_files.forward_fastq);
        free(fastq_files.sample_names);
        free(fastq_files.sample_sizes);
        free(fastq_files.sample_offsets);
        free(fastq_files.sorted_index);
        return PF_ERR_OUT_OF_MEMORY;
    }
    
    for (int i = 0; i < n_files; i++) {
        fastq_files.barcode_fastq[i] = strdup(barcode_fastqs[i]);
        fastq_files.forward_fastq[i] = strdup(feature_fastqs[i]);
    }
    
    fastq_files.sample_names[0] = strdup(sname);
    fastq_files.sample_sizes[0] = n_files;
    fastq_files.sample_offsets[0] = 0;
    fastq_files.sorted_index[0] = 0;
    fastq_files.max_sample_size = n_files;
    
    /* Process */
    double start_time = get_time_in_seconds();
    int sample_error = 0;
    
    sample_args args;
    memset(&args, 0, sizeof(args));
    args.sample_index = 0;
    args.directory = sample_directory;
    args.filtered_barcodes_name = NULL;
    args.fastq_files = &fastq_files;
    args.features = ctx->features;
    args.maxHammingDistance = ctx->config->max_hamming_distance;
    args.nThreads = ctx->config->search_threads;
    args.stringency = ctx->config->stringency;
    args.min_counts = ctx->config->min_counts;
    args.barcode_constant_offset = ctx->config->barcode_offset;
    args.feature_constant_offset = ctx->config->feature_offset;
    args.read_buffer_lines = READ_BUFFER_LINES;
    args.average_read_length = AVERAGE_READ_LENGTH;
    args.min_posterior = ctx->config->min_posterior;
    args.consumer_threads_per_set = ctx->config->consumer_threads;
    args.filtered_barcodes_hash = ctx->filtered_barcodes_hash;
    args.min_prediction = 1;
    args.min_heatmap = 0;
    args.demux_nsamples = 1;
    args.sample_barcodes = NULL;
    args.sample_max_hamming = 1;
    args.sample_max_N = 0;
    args.sample_constant_offset = -1;
    args.sample_offset_relative = 0;
    
    /* EmptyDrops control from config */
    args.skip_emptydrops = ctx->config->skip_emptydrops;
    args.emptydrops_failure_fatal = ctx->config->emptydrops_failure_fatal;
    args.expected_cells = ctx->config->expected_cells;
    args.emptydrops_use_fdr = ctx->config->emptydrops_use_fdr;
    args.error_out = &sample_error;
    
    process_files_in_sample(&args);
    
    double end_time = get_time_in_seconds();
    
    /* Fill in stats if requested */
    if (stats_out) {
        memset(stats_out, 0, sizeof(pf_stats));
        stats_out->processing_time_sec = end_time - start_time;
        stats_out->total_features = ctx->features->number_of_features;
    }
    
    /* Cleanup */
    for (int i = 0; i < n_files; i++) {
        free(fastq_files.barcode_fastq[i]);
        free(fastq_files.forward_fastq[i]);
    }
    free(fastq_files.barcode_fastq);
    free(fastq_files.forward_fastq);
    free(fastq_files.sample_names[0]);
    free(fastq_files.sample_names);
    free(fastq_files.sample_sizes);
    free(fastq_files.sample_offsets);
    free(fastq_files.sorted_index);
    
    if (sample_error) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Sample processing failed");
        return PF_ERR_IO_ERROR;
    }
    
    return PF_OK;
}

/* ============================================================================
 * Output API Implementation
 * ============================================================================ */

int pf_get_num_features(pf_context *ctx) {
    if (!ctx || !ctx->features) return 0;
    return ctx->features->number_of_features;
}

const char* pf_get_feature_name(pf_context *ctx, int index) {
    if (!ctx || !ctx->features) return NULL;
    if (index < 0 || index >= ctx->features->number_of_features) return NULL;
    return ctx->features->feature_names[index];
}

const char* pf_get_feature_sequence(pf_context *ctx, int index) {
    if (!ctx || !ctx->features) return NULL;
    if (index < 0 || index >= ctx->features->number_of_features) return NULL;
    return ctx->features->feature_sequences[index];
}

/* ============================================================================
 * EmptyDrops Filtering API Implementation (via libscrna)
 * ============================================================================ */

#include "scrna_api.h"
#include <dirent.h>
#include <ctype.h>

/* Helper: read barcodes.tsv into array */
static int read_barcodes_tsv(const char *path, char ***barcodes_out, uint32_t *n_barcodes_out) {
    FILE *f = fopen(path, "r");
    if (!f) {
        /* Try .gz */
        char gz_path[1024];
        snprintf(gz_path, sizeof(gz_path), "%s.gz", path);
        f = fopen(gz_path, "r");
        if (!f) return -1;
    }
    
    char **barcodes = NULL;
    uint32_t n = 0, capacity = 1000;
    barcodes = malloc(capacity * sizeof(char*));
    if (!barcodes) { fclose(f); return -1; }
    
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        /* Strip newline */
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;
        
        if (n >= capacity) {
            capacity *= 2;
            barcodes = realloc(barcodes, capacity * sizeof(char*));
        }
        barcodes[n++] = strdup(line);
    }
    fclose(f);
    
    *barcodes_out = barcodes;
    *n_barcodes_out = n;
    return 0;
}

/* Helper: read features.tsv into array */
static int read_features_tsv(const char *path, char ***features_out, uint32_t *n_features_out) {
    FILE *f = fopen(path, "r");
    if (!f) {
        /* Try .gz */
        char gz_path[1024];
        snprintf(gz_path, sizeof(gz_path), "%s.gz", path);
        f = fopen(gz_path, "r");
        if (!f) return -1;
    }
    
    char **features = NULL;
    uint32_t n = 0, capacity = 100;
    features = malloc(capacity * sizeof(char*));
    if (!features) { fclose(f); return -1; }
    
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        /* Extract first column (tab-separated) */
        char *tab = strchr(line, '\t');
        if (tab) *tab = '\0';
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;
        
        if (n >= capacity) {
            capacity *= 2;
            features = realloc(features, capacity * sizeof(char*));
        }
        features[n++] = strdup(line);
    }
    fclose(f);
    
    *features_out = features;
    *n_features_out = n;
    return 0;
}

/* Helper: read matrix.mtx and build sparse arrays */
static int read_matrix_mtx(const char *path,
                           uint32_t n_features, uint32_t n_barcodes,
                           uint32_t **umi_counts_out,
                           uint32_t **sparse_gene_ids_out,
                           uint32_t **sparse_counts_out,
                           uint32_t **sparse_cell_index_out,
                           uint32_t **n_genes_per_cell_out,
                           size_t *nnz_out) {
    FILE *f = fopen(path, "r");
    if (!f) {
        /* Try .gz - would need zlib, skip for now */
        return -1;
    }
    
    char line[1024];
    
    /* Skip header comments */
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '%') break;
    }
    
    /* Parse dimensions: nrows ncols nnz */
    uint32_t nrows, ncols, nnz;
    if (sscanf(line, "%u %u %u", &nrows, &ncols, &nnz) != 3) {
        fclose(f);
        return -1;
    }
    
    if (nrows != n_features || ncols != n_barcodes) {
        fprintf(stderr, "[pf_api] Matrix dimensions (%u x %u) don't match features (%u) and barcodes (%u)\n",
                nrows, ncols, n_features, n_barcodes);
        fclose(f);
        return -1;
    }
    
    /* Allocate sparse arrays */
    uint32_t *sparse_gene_ids = malloc(nnz * sizeof(uint32_t));
    uint32_t *sparse_counts = malloc(nnz * sizeof(uint32_t));
    uint32_t *cell_starts = calloc(n_barcodes + 1, sizeof(uint32_t));  /* Temp: count per cell */
    
    if (!sparse_gene_ids || !sparse_counts || !cell_starts) {
        free(sparse_gene_ids);
        free(sparse_counts);
        free(cell_starts);
        fclose(f);
        return -1;
    }
    
    /* First pass: count entries per cell */
    long pos = ftell(f);
    while (fgets(line, sizeof(line), f)) {
        uint32_t row, col, val;
        if (sscanf(line, "%u %u %u", &row, &col, &val) == 3) {
            if (col >= 1 && col <= n_barcodes) {
                cell_starts[col]++;  /* 1-indexed to 0-indexed mapping */
            }
        }
    }
    
    /* Compute cumulative starts */
    uint32_t *sparse_cell_index = malloc((n_barcodes + 1) * sizeof(uint32_t));
    uint32_t *n_genes_per_cell = calloc(n_barcodes, sizeof(uint32_t));
    sparse_cell_index[0] = 0;
    for (uint32_t i = 0; i < n_barcodes; i++) {
        n_genes_per_cell[i] = cell_starts[i + 1];  /* cell_starts is 1-indexed counts */
        sparse_cell_index[i + 1] = sparse_cell_index[i] + cell_starts[i + 1];
    }
    
    /* Reset cell_starts for use as current write position */
    memset(cell_starts, 0, (n_barcodes + 1) * sizeof(uint32_t));
    
    /* Second pass: fill sparse arrays */
    fseek(f, pos, SEEK_SET);
    uint32_t *umi_counts = calloc(n_barcodes, sizeof(uint32_t));
    
    while (fgets(line, sizeof(line), f)) {
        uint32_t row, col, val;
        if (sscanf(line, "%u %u %u", &row, &col, &val) == 3) {
            if (col >= 1 && col <= n_barcodes && row >= 1 && row <= n_features) {
                uint32_t cell_idx = col - 1;
                uint32_t gene_idx = row - 1;
                uint32_t write_pos = sparse_cell_index[cell_idx] + cell_starts[col];
                
                sparse_gene_ids[write_pos] = gene_idx;
                sparse_counts[write_pos] = val;
                cell_starts[col]++;
                
                umi_counts[cell_idx] += val;
            }
        }
    }
    
    fclose(f);
    free(cell_starts);
    
    *umi_counts_out = umi_counts;
    *sparse_gene_ids_out = sparse_gene_ids;
    *sparse_counts_out = sparse_counts;
    *sparse_cell_index_out = sparse_cell_index;
    *n_genes_per_cell_out = n_genes_per_cell;
    *nnz_out = nnz;
    
    return 0;
}

/* ============================================================================
 * Pre-MEX EmptyDrops (pipeline integration point)
 * ============================================================================ */

pf_error pf_run_emptydrops_premex(
    const uint32_t *umi_counts,
    const char **barcodes,
    uint32_t n_barcodes,
    const char **features,
    uint32_t n_features,
    const uint32_t *sparse_gene_ids,
    const uint32_t *sparse_counts,
    const uint32_t *sparse_cell_index,
    const uint32_t *n_genes_per_cell,
    const char *output_dir,
    int n_expected_cells,
    int use_fdr_gate,
    char ***filtered_barcodes_out,
    uint32_t *n_filtered_out
) {
    if (!umi_counts || !barcodes || n_barcodes == 0 || !output_dir) {
        return PF_ERR_INVALID_ARG;
    }
    
    fprintf(stderr, "[pf_api] Running EmptyDrops on pre-MEX data (%u barcodes)...\n", n_barcodes);
    
    /* Set up EmptyDrops input */
    scrna_matrix_input input;
    memset(&input, 0, sizeof(input));
    input.umi_counts = (uint32_t*)umi_counts;  /* Cast away const for API compatibility */
    input.barcodes = (char**)barcodes;
    input.features = (char**)features;
    input.n_cells = n_barcodes;
    input.n_features = n_features;
    input.sparse_gene_ids = (uint32_t*)sparse_gene_ids;
    input.sparse_counts = (uint32_t*)sparse_counts;
    input.sparse_cell_index = (uint32_t*)sparse_cell_index;
    input.n_genes_per_cell = (uint32_t*)n_genes_per_cell;
    
    /* Create config */
    scrna_ed_config *config = scrna_ed_config_create();
    if (n_expected_cells > 0) {
        config->n_expected_cells = n_expected_cells;
    }
    config->disable_occupancy_filter = 1;  /* Disabled for compat mode */
    config->use_fdr_gate = use_fdr_gate ? 1 : 0;
    
    /* Run EmptyDrops */
    scrna_ed_result result;
    int rc = scrna_emptydrops_run(&input, config, &result);
    
    if (rc != 0) {
        fprintf(stderr, "[pf_api] EmptyDrops failed: %s\n", 
                result.error_message ? result.error_message : "unknown error");
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return PF_ERR_IO_ERROR;
    }
    
    fprintf(stderr, "[pf_api] EmptyDrops: %zu cells pass\n", result.n_barcodes);
    
    /* Write outputs */
    rc = scrna_emptydrops_write_outputs(&result, output_dir);
    if (rc != 0) {
        fprintf(stderr, "[pf_api] Failed to write EmptyDrops outputs\n");
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        return PF_ERR_IO_ERROR;
    }
    
    /* Return filtered barcodes to caller */
    if (filtered_barcodes_out && n_filtered_out) {
        *n_filtered_out = result.n_barcodes;
        *filtered_barcodes_out = (char**)malloc(result.n_barcodes * sizeof(char*));
        if (*filtered_barcodes_out) {
            for (size_t i = 0; i < result.n_barcodes; i++) {
                (*filtered_barcodes_out)[i] = strdup(result.barcodes[i]);
            }
        }
    }
    
    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
    
    return PF_OK;
}

/* ============================================================================
 * MEX-based EmptyDrops (standalone tool only - NOT for pipeline use)
 * ============================================================================ */

pf_error pf_run_emptydrops_mex(const char *mex_dir, const char *output_dir, int n_expected_cells, int use_fdr_gate) {
    if (!mex_dir || !output_dir) {
        return PF_ERR_INVALID_ARG;
    }
    
    fprintf(stderr, "[pf_api] Running EmptyDrops on MEX directory %s (standalone mode)...\n", mex_dir);
    
    /* Build file paths */
    char barcodes_path[1024], features_path[1024], matrix_path[1024];
    snprintf(barcodes_path, sizeof(barcodes_path), "%s/barcodes.tsv", mex_dir);
    snprintf(features_path, sizeof(features_path), "%s/features.tsv", mex_dir);
    snprintf(matrix_path, sizeof(matrix_path), "%s/matrix.mtx", mex_dir);
    
    /* Read barcodes */
    char **barcodes = NULL;
    uint32_t n_barcodes = 0;
    if (read_barcodes_tsv(barcodes_path, &barcodes, &n_barcodes) != 0) {
        fprintf(stderr, "[pf_api] Failed to read barcodes from %s\n", barcodes_path);
        return PF_ERR_FILE_NOT_FOUND;
    }
    fprintf(stderr, "[pf_api] Loaded %u barcodes\n", n_barcodes);
    
    /* Read features */
    char **features = NULL;
    uint32_t n_features = 0;
    if (read_features_tsv(features_path, &features, &n_features) != 0) {
        fprintf(stderr, "[pf_api] Failed to read features from %s\n", features_path);
        for (uint32_t i = 0; i < n_barcodes; i++) free(barcodes[i]);
        free(barcodes);
        return PF_ERR_FILE_NOT_FOUND;
    }
    fprintf(stderr, "[pf_api] Loaded %u features\n", n_features);
    
    /* Read matrix */
    uint32_t *umi_counts = NULL;
    uint32_t *sparse_gene_ids = NULL;
    uint32_t *sparse_counts = NULL;
    uint32_t *sparse_cell_index = NULL;
    uint32_t *n_genes_per_cell = NULL;
    size_t nnz = 0;
    
    if (read_matrix_mtx(matrix_path, n_features, n_barcodes,
                        &umi_counts, &sparse_gene_ids, &sparse_counts,
                        &sparse_cell_index, &n_genes_per_cell, &nnz) != 0) {
        fprintf(stderr, "[pf_api] Failed to read matrix from %s\n", matrix_path);
        for (uint32_t i = 0; i < n_barcodes; i++) free(barcodes[i]);
        free(barcodes);
        for (uint32_t i = 0; i < n_features; i++) free(features[i]);
        free(features);
        return PF_ERR_FILE_NOT_FOUND;
    }
    fprintf(stderr, "[pf_api] Loaded matrix: %zu non-zero entries\n", nnz);
    
    /* Set up EmptyDrops input */
    scrna_matrix_input input;
    memset(&input, 0, sizeof(input));
    input.umi_counts = umi_counts;
    input.barcodes = barcodes;
    input.features = features;
    input.n_cells = n_barcodes;
    input.n_features = n_features;
    input.sparse_gene_ids = sparse_gene_ids;
    input.sparse_counts = sparse_counts;
    input.sparse_cell_index = sparse_cell_index;
    input.n_genes_per_cell = n_genes_per_cell;
    input.sparse_nnz = nnz;
    
    /* Create config */
    scrna_ed_config *config = scrna_ed_config_create();
    if (n_expected_cells > 0) {
        config->n_expected_cells = n_expected_cells;
    }
    config->disable_occupancy_filter = 1;  /* Disabled for compat mode */
    config->use_fdr_gate = use_fdr_gate ? 1 : 0;
    
    /* Run EmptyDrops */
    scrna_ed_result result;
    int rc = scrna_emptydrops_run(&input, config, &result);
    
    if (rc != 0) {
        fprintf(stderr, "[pf_api] EmptyDrops failed: %s\n", 
                result.error_message ? result.error_message : "unknown error");
        scrna_ed_result_free(&result);
        scrna_ed_config_destroy(config);
        goto cleanup;
    }
    
    fprintf(stderr, "[pf_api] EmptyDrops: %zu cells pass\n", result.n_barcodes);
    
    /* Write outputs */
    rc = scrna_emptydrops_write_outputs(&result, output_dir);
    if (rc != 0) {
        fprintf(stderr, "[pf_api] Failed to write EmptyDrops outputs\n");
    }
    
    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
    
cleanup:
    /* Free memory */
    for (uint32_t i = 0; i < n_barcodes; i++) free(barcodes[i]);
    free(barcodes);
    for (uint32_t i = 0; i < n_features; i++) free(features[i]);
    free(features);
    free(umi_counts);
    free(sparse_gene_ids);
    free(sparse_counts);
    free(sparse_cell_index);
    free(n_genes_per_cell);
    
    return (rc == 0) ? PF_OK : PF_ERR_IO_ERROR;
}

/* ============================================================================
 * Utility Functions Implementation
 * ============================================================================ */

const char* pf_version(void) {
    return PF_VERSION;
}
