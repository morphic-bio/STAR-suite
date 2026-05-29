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
#include <pthread.h>

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
    int legacy_cb_rescue;
    int feature_offset;
    int feature_offset_explicit;  /* 1 if user called pf_config_set_feature_offset() */
    int barcode_offset;
    int max_barcode_mismatches;
    int max_feature_n;
    int max_barcode_n;
    int max_threads;
    int search_threads;
    int consumer_threads;
    pf_permit_acquire_fn permit_acquire_cb;
    pf_permit_release_fn permit_release_cb;
    void *permit_hook_ctx;
    int debug_enabled;
    int reverse_complement_whitelist;
    int limit_search;
    int feature_limited_fallback_mode;
    long long max_reads;
    int translate_nxt;
    int use_feature_offset_array;
    int strict_offset_check;
    int use_feature_anchor_search;
    int require_feature_anchor_match;
    int feature_mode_bootstrap_reads;
    int use_hot_hash;
    int skip_heatmaps;
    
    /* Prehash memory budget (0 = auto-detect) */
    unsigned long long prehash_memory_budget;

    /* EmptyDrops control */
    int skip_emptydrops;            /* 1 = skip EmptyDrops entirely */
    int emptydrops_failure_fatal;   /* 1 = treat ED failure as error */
    int expected_cells;             /* 0 = auto-detect */
    int emptydrops_use_fdr;         /* 1 = use FDR gate for tail rescue */

    /* NXT/TRU auto-detection */
    int autodetect_chemistry;           /* 0=off, 1=on */
    int autodetect_chemistry_reads;     /* N reads to sample (default 10000) */
    int autodetect_chemistry_min_hits;  /* minimum total hits for decision (default 50) */
    int probe_only;                     /* 1=lightweight chemistry probe, no outputs */
    int skip_qc_outputs;                /* 1=skip feature histograms/heatmaps */

    /* Union whitelist support (legacy compat) */
    int allow_union_whitelist;          /* 0=strict, 1=accept mixed NXT+TRU */

    /* Namespace normalization */
    pf_namespace_t source_namespace;    /* namespace of filtered barcode file (PF_NS_UNKNOWN = not set) */
    pf_namespace_t target_namespace;    /* namespace of assignment output (PF_NS_UNKNOWN = not set) */
};

struct pf_context {
    pf_config *config;
    feature_arrays *features;
    khash_t(u32ptr) *whitelist_hash;
    unsigned char *whitelist_data;
    khash_t(strptr) *filtered_barcodes_hash;
    int initialized;
    char error_buf[PF_ERROR_BUF_SIZE];

    /* NXT/TRU auto-detection state */
    unsigned long long chem_detect_raw_hits;
    unsigned long long chem_detect_nxt_hits;
    unsigned long long chem_detect_ticket;
    int chem_detect_done;       /* 0=sampling, 1=decided */
    int detected_match_mode;    /* 0=unknown, 1=RAW_MATCH, 2=TRANSLATED_MATCH, 3=AMBIGUOUS */
};

struct pf_record_stream {
    pf_context *ctx;
    sample_args sample_args;
    statistics *stats;
    data_structures *hashes;
    memory_pool_collection **pools;
    fastq_reader_set *reader_sets[1];
    fastq_processor *processors;
    pthread_t *consumer_threads;
    struct chem_detect_state chem_detect;
    char sample_directory[FILENAME_LENGTH];
    double start_time;
    int nconsumers;
    int nreaders;
    int lines_per_block;
    int queue_initialized;
    int consumers_started;
    int consumers_joined;
    int processors_initialized;
    int worker_state_initialized;
    int chem_detect_enabled;
    int failed;
    int closed;
};

/* Global initialization flag */
static int g_pf_global_initialized = 0;
/* process_features legacy core still uses process-global mutable state.
 * Serialize API entrypoints that mutate/read that state to avoid cross-context
 * namespace/config bleed when multiple contexts are used in one process. */
static pthread_mutex_t g_pf_runtime_mutex = PTHREAD_MUTEX_INITIALIZER;

static void pf_apply_context_globals(pf_context *ctx) {
    if (!ctx || !ctx->config) {
        return;
    }

    barcode_length = ctx->config->barcode_length;
    barcode_code_length = (barcode_length + 3) / 4;
    umi_length = ctx->config->umi_length;
    umi_code_length = (umi_length + 3) / 4;
    max_barcode_mismatches = ctx->config->max_barcode_mismatches;
    max_feature_n = ctx->config->max_feature_n;
    max_barcode_n = ctx->config->max_barcode_n;
    max_reads = ctx->config->max_reads;
    limit_search = ctx->config->limit_search;
    feature_limited_fallback_mode = ctx->config->feature_limited_fallback_mode;
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
    feature_mode_search_offsets_reset();
    use_hot_hash = ctx->config->use_hot_hash;
    skip_heatmaps = ctx->config->skip_heatmaps;

    if (ctx->features) {
        /* Set global hash alias from per-instance hash (non-owning) */
        feature_code_hash = ctx->features->feature_code_hash;
        feature_code_hash_mode = ctx->features->code_hash_mode;
        feature_code_hash_fixed_length = ctx->features->code_hash_fixed_length;
        feature_offsets = ctx->features->feature_offsets;
        feature_offsets_count = ctx->features->number_of_features;
        feature_anchors = ctx->features->feature_anchors;
        feature_anchor_lengths = ctx->features->feature_anchor_lengths;
        feature_suffix_anchors = ctx->features->feature_suffix_anchors;
        feature_suffix_anchor_lengths = ctx->features->feature_suffix_anchor_lengths;
        feature_anchor_count = ctx->features->number_of_features;
        number_of_features = ctx->features->number_of_features;
        maximum_feature_length = ctx->features->max_length;
        feature_code_length = (maximum_feature_length + 3) / 4;
        initialize_unit_sizes();
    } else {
        feature_offsets = NULL;
        feature_offsets_count = 0;
        feature_anchors = NULL;
        feature_anchor_lengths = NULL;
        feature_suffix_anchors = NULL;
        feature_suffix_anchor_lengths = NULL;
        feature_anchor_count = 0;
        number_of_features = 0;
        maximum_feature_length = 0;
        feature_code_length = 0;
    }

    whitelist_hash = ctx->whitelist_hash;
    whitelist = ctx->whitelist_data;
}

/* ============================================================================
 * Namespace API Implementation
 * ============================================================================ */

static inline char pf_complement_base(char b) {
    switch ((unsigned char)b) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 'T'; case 't': return 'A';
        case 'c': return 'G'; case 'g': return 'C';
        default: return b;
    }
}

void pf_translate_barcode_inplace(char *barcode, size_t len) {
    if (!barcode || len < 9) return;
    barcode[7] = pf_complement_base(barcode[7]);
    barcode[8] = pf_complement_base(barcode[8]);
}

const char* pf_namespace_to_string(pf_namespace_t ns) {
    switch (ns) {
        case PF_NS_NXT:     return "NXT";
        case PF_NS_TRU:     return "TRU";
        case PF_NS_UNION:   return "UNION";
        case PF_NS_UNKNOWN: return "UNKNOWN";
        default:             return "UNKNOWN";
    }
}

pf_namespace_t pf_namespace_from_string(const char *s) {
    if (!s) return PF_NS_UNKNOWN;
    if (strcmp(s, "NXT") == 0 || strcmp(s, "nxt") == 0) return PF_NS_NXT;
    if (strcmp(s, "TRU") == 0 || strcmp(s, "tru") == 0) return PF_NS_TRU;
    if (strcmp(s, "UNION") == 0 || strcmp(s, "union") == 0) return PF_NS_UNION;
    return PF_NS_UNKNOWN;
}

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
    config->legacy_cb_rescue = 0;
    config->feature_offset = 0;
    config->feature_offset_explicit = 0;
    config->barcode_offset = 0;
    config->max_barcode_mismatches = 3;
    config->max_feature_n = 1;
    config->max_barcode_n = 1;
    config->max_threads = 8;
    config->search_threads = 4;
    config->consumer_threads = 1;
    config->permit_acquire_cb = NULL;
    config->permit_release_cb = NULL;
    config->permit_hook_ctx = NULL;
    config->debug_enabled = 0;
    config->reverse_complement_whitelist = 0;
    config->limit_search = -1;
    config->feature_limited_fallback_mode = 0;
    config->max_reads = 0;
    config->translate_nxt = 0;
    config->use_feature_offset_array = 0;
    config->strict_offset_check = 0;
    config->use_feature_anchor_search = 0;
    config->require_feature_anchor_match = 0;
    config->feature_mode_bootstrap_reads = 0;
    config->use_hot_hash = 0;
    config->skip_heatmaps = 0;
    
    /* Prehash memory budget default (0 = auto-detect) */
    config->prehash_memory_budget = 0;

    /* EmptyDrops control defaults */
    config->skip_emptydrops = 0;
    config->emptydrops_failure_fatal = 0;
    config->expected_cells = 0;
    config->emptydrops_use_fdr = 0;

    /* NXT/TRU auto-detection defaults */
    config->autodetect_chemistry = 0;
    config->autodetect_chemistry_reads = 10000;
    config->autodetect_chemistry_min_hits = 50;
    config->probe_only = 0;
    config->skip_qc_outputs = 0;

    /* Union whitelist: off by default (strict namespace) */
    config->allow_union_whitelist = 0;

    /* Namespace normalization: both unknown = no normalization */
    config->source_namespace = PF_NS_UNKNOWN;
    config->target_namespace = PF_NS_UNKNOWN;
    
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

void pf_config_set_legacy_cb_rescue(pf_config *config, int enable) {
    if (config) config->legacy_cb_rescue = enable;
}

void pf_config_set_feature_offset(pf_config *config, int offset) {
    if (config) {
        config->feature_offset = offset;
        config->feature_offset_explicit = 1;
    }
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

void pf_config_set_permit_hooks(
    pf_config *config,
    pf_permit_acquire_fn acquire_cb,
    pf_permit_release_fn release_cb,
    void *hook_ctx
) {
    if (!config) {
        return;
    }
    config->permit_acquire_cb = acquire_cb;
    config->permit_release_cb = release_cb;
    config->permit_hook_ctx = hook_ctx;
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

void pf_config_set_feature_limited_mode(pf_config *config, int mode) {
    if (!config) return;
    if (mode < 0 || mode > 1) {
        fprintf(stderr, "ERROR: feature_limited_mode must be 0 or 1 (got %d), clamping to 0\n", mode);
        mode = 0;
    }
    config->feature_limited_fallback_mode = mode;
}

void pf_config_set_feature_limited_fallback(pf_config *config, int mode) {
    pf_config_set_feature_limited_mode(config, mode);
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
void pf_config_set_strict_offset_check(pf_config *config, int enable) {
    if (config) config->strict_offset_check = enable;
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

void pf_config_set_use_hot_hash(pf_config *config, int enable) {
    if (config) config->use_hot_hash = enable;
}

void pf_config_set_skip_heatmaps(pf_config *config, int enable) {
    if (config) config->skip_heatmaps = enable;
}

void pf_config_set_prehash_memory_budget(pf_config *config, unsigned long long budget) {
    if (config) config->prehash_memory_budget = budget;
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

void pf_config_set_autodetect_chemistry(pf_config *config, int enabled) {
    if (config) config->autodetect_chemistry = enabled;
}

void pf_config_set_autodetect_chemistry_reads(pf_config *config, int n_reads) {
    if (!config) return;
    if (n_reads < 1) {
        fprintf(stderr, "ERROR: autodetect_chemistry_reads must be >= 1 (got %d), clamping to 1\n", n_reads);
        n_reads = 1;
    }
    config->autodetect_chemistry_reads = n_reads;
}

void pf_config_set_autodetect_chemistry_min_hits(pf_config *config, int min_hits) {
    if (!config) return;
    if (min_hits < 1) {
        fprintf(stderr, "ERROR: autodetect_chemistry_min_hits must be >= 1 (got %d), clamping to 1\n", min_hits);
        min_hits = 1;
    }
    config->autodetect_chemistry_min_hits = min_hits;
}

void pf_config_set_probe_only(pf_config *config, int enabled) {
    if (config) config->probe_only = enabled;
}

void pf_config_set_skip_qc_outputs(pf_config *config, int enabled) {
    if (config) config->skip_qc_outputs = enabled;
}

void pf_config_set_allow_union_whitelist(pf_config *config, int enable) {
    if (config) config->allow_union_whitelist = (enable != 0);
}

void pf_config_set_source_namespace(pf_config *config, pf_namespace_t ns) {
    if (config) config->source_namespace = ns;
}

void pf_config_set_target_namespace(pf_config *config, pf_namespace_t ns) {
    if (config) config->target_namespace = ns;
}

const char* pf_get_detected_match_mode(pf_context *ctx) {
    if (!ctx) return "UNKNOWN";
    switch (ctx->detected_match_mode) {
        case 1: return "RAW_MATCH";
        case 2: return "TRANSLATED_MATCH";
        case 3: return "AMBIGUOUS";
        default: return "UNKNOWN";
    }
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
    
    /* Feature code hash is now per-instance (owned by feature_arrays).
     * Global alias is set by pf_apply_context_globals() after feature load. */
    
    g_pf_global_initialized = 1;
}

/* ============================================================================
 * Context Lifecycle Implementation
 * ============================================================================ */

pf_context* pf_init(const pf_config *config) {
    if (!config) return NULL;

    pthread_mutex_lock(&g_pf_runtime_mutex);

    pf_context *ctx = calloc(1, sizeof(pf_context));
    if (!ctx) {
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return NULL;
    }
    
    ctx->config = pf_config_clone(config);
    if (!ctx->config) {
        free(ctx);
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return NULL;
    }
    
    /* Apply config to globals for legacy internals. */
    pf_apply_context_globals(ctx);
    
    /* NXT/TRU auto-detection state */
    ctx->chem_detect_raw_hits = 0;
    ctx->chem_detect_nxt_hits = 0;
    ctx->chem_detect_ticket = 0;
    ctx->chem_detect_done = 0;
    ctx->detected_match_mode = 0;

    /* Global initialization */
    pf_global_init();
    
    ctx->initialized = 1;
    ctx->error_buf[0] = '\0';
    
    pthread_mutex_unlock(&g_pf_runtime_mutex);
    return ctx;
}

void pf_destroy(pf_context *ctx) {
    if (!ctx) return;

    pthread_mutex_lock(&g_pf_runtime_mutex);

    if (ctx->config) {
        pf_config_destroy(ctx->config);
    }
    
    if (ctx->features) {
        /* Nullify global alias if it points into this instance */
        if (feature_code_hash_mode == ctx->features->code_hash_mode) {
            if (feature_code_hash_mode == SEQ_KEY_64 &&
                feature_code_hash.h64 == ctx->features->feature_code_hash.h64) {
                feature_code_hash.h64 = NULL;
            } else if (feature_code_hash_mode == SEQ_KEY_128 &&
                       feature_code_hash.h128 == ctx->features->feature_code_hash.h128) {
                feature_code_hash.h128 = NULL;
            }
        }
        if (feature_offsets == ctx->features->feature_offsets) {
            feature_offsets = NULL;
            feature_offsets_count = 0;
        }
        if (feature_anchors == ctx->features->feature_anchors) {
            feature_anchors = NULL;
            feature_anchor_lengths = NULL;
            feature_suffix_anchors = NULL;
            feature_suffix_anchor_lengths = NULL;
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
    feature_mode_search_offsets_reset();
    
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
    pthread_mutex_unlock(&g_pf_runtime_mutex);
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

    pthread_mutex_lock(&g_pf_runtime_mutex);
    pf_apply_context_globals(ctx);

    if (!pf_file_exists(feature_csv)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, 
                 "Feature reference file not found: %s", feature_csv);
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Free existing features if any */
    if (ctx->features) {
        /* Nullify global hash alias before free */
        if (feature_code_hash_mode == ctx->features->code_hash_mode) {
            if (feature_code_hash_mode == SEQ_KEY_64 &&
                feature_code_hash.h64 == ctx->features->feature_code_hash.h64) {
                feature_code_hash.h64 = NULL;
            } else if (feature_code_hash_mode == SEQ_KEY_128 &&
                       feature_code_hash.h128 == ctx->features->feature_code_hash.h128) {
                feature_code_hash.h128 = NULL;
            }
        }
        if (feature_offsets == ctx->features->feature_offsets) {
            feature_offsets = NULL;
            feature_offsets_count = 0;
        }
        if (feature_anchors == ctx->features->feature_anchors) {
            feature_anchors = NULL;
            feature_anchor_lengths = NULL;
            feature_suffix_anchors = NULL;
            feature_suffix_anchor_lengths = NULL;
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
    feature_mode_search_offsets_reset();
    
    /* Wire prehash memory budget from config */
    if (ctx->config->prehash_memory_budget > 0) {
        feature_prehash_memory_budget = ctx->config->prehash_memory_budget;
    } else if (feature_prehash_memory_budget == 0) {
        feature_prehash_memory_budget = prehash_detect_memory_budget();
    }

    int prehash_max_saved = feature_prehash_max_hamming;
    int prehash_max_effective = feature_prehash_max_hamming;
    if (ctx->config && ctx->config->max_hamming_distance >= 0 &&
        prehash_max_effective > ctx->config->max_hamming_distance) {
        prehash_max_effective = ctx->config->max_hamming_distance;
    }
    if (prehash_max_effective < 0) {
        prehash_max_effective = 0;
    }
    if (prehash_max_effective != feature_prehash_max_hamming) {
        fprintf(stderr,
                "Feature prehash clamp: requested max_hamming=%d, assignment max_hamming=%d, effective=%d\n",
                feature_prehash_max_hamming,
                ctx->config ? ctx->config->max_hamming_distance : -1,
                prehash_max_effective);
        feature_prehash_max_hamming = prehash_max_effective;
    }

    ctx->features = read_features_file(feature_csv);
    feature_prehash_max_hamming = prehash_max_saved;
    if (!ctx->features) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Failed to parse feature reference: %s", feature_csv);
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_PARSE_ERROR;
    }

    feature_offsets = ctx->features->feature_offsets;
    feature_offsets_count = ctx->features->number_of_features;
    feature_anchors = ctx->features->feature_anchors;
    feature_anchor_lengths = ctx->features->feature_anchor_lengths;
    feature_suffix_anchors = ctx->features->feature_suffix_anchors;
    feature_suffix_anchor_lengths = ctx->features->feature_suffix_anchor_lengths;
    feature_anchor_count = ctx->features->number_of_features;

    if (feature_mode_bootstrap_reads > 0) {
        const int count = ctx->features->number_of_features;
        feature_mode_offsets = malloc(sizeof(int) * count);
        feature_mode_hist = calloc((size_t)count * (size_t)feature_mode_max_offset, sizeof(unsigned int));
        if (!feature_mode_offsets || !feature_mode_hist) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to allocate feature mode arrays");
            pthread_mutex_unlock(&g_pf_runtime_mutex);
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
    initialize_unit_sizes();

    pthread_mutex_unlock(&g_pf_runtime_mutex);
    return PF_OK;
}

pf_error pf_load_whitelist(pf_context *ctx, const char *whitelist_path) {
    if (!ctx || !whitelist_path) return PF_ERR_INVALID_ARG;
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;

    pthread_mutex_lock(&g_pf_runtime_mutex);
    pf_apply_context_globals(ctx);

    if (!pf_file_exists(whitelist_path)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Whitelist file not found: %s", whitelist_path);
        pthread_mutex_unlock(&g_pf_runtime_mutex);
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
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_IO_ERROR;
    }
    
    /* Update global pointer for existing code */
    whitelist_hash = ctx->whitelist_hash;
    whitelist = ctx->whitelist_data;

    pthread_mutex_unlock(&g_pf_runtime_mutex);
    return PF_OK;
}

pf_error pf_load_filtered_barcodes(pf_context *ctx, const char *filtered_path) {
    if (!ctx || !filtered_path) return PF_ERR_INVALID_ARG;
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;

    pthread_mutex_lock(&g_pf_runtime_mutex);
    pf_apply_context_globals(ctx);

    if (!pf_file_exists(filtered_path)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Filtered barcodes file not found: %s", filtered_path);
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Free existing filtered barcodes if any */
    if (ctx->filtered_barcodes_hash) {
        free_strptr_hash(ctx->filtered_barcodes_hash);
    }
    
    ctx->filtered_barcodes_hash = kh_init(strptr);
    read_barcodes_into_hash(filtered_path, ctx->filtered_barcodes_hash);

    /* Ingress namespace normalization: if source and target are both known
     * single-namespace (NXT vs TRU) and differ, translate every barcode. */
    pf_namespace_t src_ns = ctx->config->source_namespace;
    pf_namespace_t tgt_ns = ctx->config->target_namespace;
    if ((src_ns == PF_NS_NXT || src_ns == PF_NS_TRU) &&
        (tgt_ns == PF_NS_NXT || tgt_ns == PF_NS_TRU) &&
        src_ns != tgt_ns &&
        !ctx->config->allow_union_whitelist) {
        int rc = pf_normalize_hash_namespace(ctx->filtered_barcodes_hash);
        if (rc < 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to normalize filtered barcodes from %s to %s "
                     "(allocation failure)",
                     pf_namespace_to_string(src_ns),
                     pf_namespace_to_string(tgt_ns));
            pthread_mutex_unlock(&g_pf_runtime_mutex);
            return PF_ERR_ALLOC;
        }
        fprintf(stderr,
            "NOTICE: Normalized %d filtered barcodes from %s to %s namespace.\n",
            rc, pf_namespace_to_string(src_ns), pf_namespace_to_string(tgt_ns));
    }

    if (ctx->config->allow_union_whitelist) {
        int added = expand_hash_union_namespace(ctx->filtered_barcodes_hash);
        if (added < 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to expand filtered barcodes for union whitelist "
                     "(allocation failure)");
            pthread_mutex_unlock(&g_pf_runtime_mutex);
            return PF_ERR_ALLOC;
        }
        if (added > 0) {
            fprintf(stderr,
                "WARNING: --allow_union_whitelist is active. Added %d NXT/TRU "
                "translated barcodes to filtered set (total now %u).\n"
                "         This is a legacy compatibility mode for union whitelists\n"
                "         (e.g. raw 3M-february-2018.txt). Migrate to namespace-split\n"
                "         whitelist files for new workflows.\n",
                added, kh_size(ctx->filtered_barcodes_hash));
        } else {
            fprintf(stderr,
                "NOTICE: --allow_union_whitelist active but no additional translations "
                "needed (filtered set already covers both namespaces or is single-namespace).\n");
        }
    }

    /* Hard error when namespace metadata is incomplete or missing and union
     * mode is off.  With exact-only matching, any gap means silent drops. */
    if (!ctx->config->allow_union_whitelist &&
        kh_size(ctx->filtered_barcodes_hash) > 0) {
        int src_known = (src_ns == PF_NS_NXT || src_ns == PF_NS_TRU);
        int tgt_known = (tgt_ns == PF_NS_NXT || tgt_ns == PF_NS_TRU);
        if (!src_known || !tgt_known) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                "Filtered barcodes loaded with incomplete namespace metadata "
                "(source=%s, target=%s) and without --allow_union_whitelist.  "
                "With exact-only matching, a namespace mismatch will silently "
                "drop barcodes.  Set both --source_namespace and "
                "--target_namespace, or use --allow_union_whitelist.",
                pf_namespace_to_string(src_ns),
                pf_namespace_to_string(tgt_ns));
            pthread_mutex_unlock(&g_pf_runtime_mutex);
            return PF_ERR_NAMESPACE;
        }
    }

    pthread_mutex_unlock(&g_pf_runtime_mutex);
    return PF_OK;
}

/* ============================================================================
 * Processing Implementation
 * ============================================================================ */

static pf_error pf_prepare_feature_offset_config(pf_context *ctx) {
    if (ctx->config->use_feature_offset_array && ctx->config->feature_offset_explicit) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Cannot specify both use_feature_offset_array and explicit feature_offset");
        return PF_ERR_OFFSET_CONFLICT;
    }

    if (!ctx->config->use_feature_offset_array && !ctx->config->feature_offset_explicit &&
        ctx->features && feature_offsets && feature_offsets_count > 0) {
        int offset_counts[256] = {0};
        int max_offset_seen = -1;
        int valid_offsets = 0;

        for (int i = 0; i < feature_offsets_count; i++) {
            int off = feature_offsets[i];
            if (off >= 0 && off < 256) {
                offset_counts[off]++;
                valid_offsets++;
                if (off > max_offset_seen) max_offset_seen = off;
            }
        }

        if (valid_offsets > 0) {
            int dominant_offset = 0;
            int dominant_count = 0;
            int second_count = 0;

            for (int i = 0; i <= max_offset_seen; i++) {
                if (offset_counts[i] > dominant_count) {
                    second_count = dominant_count;
                    dominant_count = offset_counts[i];
                    dominant_offset = i;
                } else if (offset_counts[i] > second_count) {
                    second_count = offset_counts[i];
                }
            }

            double heterogeneity_threshold = 0.05;
            if (second_count > 0 &&
                (double)second_count / (double)dominant_count > heterogeneity_threshold) {
                if (ctx->config->strict_offset_check) {
                    snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                             "Multiple feature offsets detected (strict mode). "
                             "Dominant: %d (%d features), second: %d features. "
                             "Use pf_config_set_feature_offset() or pf_config_set_use_feature_offset_array(1).",
                             dominant_offset, dominant_count, second_count);
                    return PF_ERR_MULTI_OFFSET_DETECTED;
                }
                fprintf(stderr, "\nWARNING: Multiple feature offsets detected in pattern column.\n");
                fprintf(stderr, "         Dominant offset: %d (used by %d features)\n",
                        dominant_offset, dominant_count);
                fprintf(stderr, "         Other offsets detected:\n");
                for (int i = 0; i <= max_offset_seen; i++) {
                    if (i != dominant_offset && offset_counts[i] > 0) {
                        double pct = 100.0 * offset_counts[i] / dominant_count;
                        fprintf(stderr, "           offset %d: %d features (%.1f%%)\n",
                                i, offset_counts[i], pct);
                    }
                }
                fprintf(stderr, "\n         Proceeding with dominant offset %d.\n\n",
                        dominant_offset);
            }

            ctx->config->feature_offset = dominant_offset;
            fprintf(stderr, "[offset-detect] Auto-detected global offset: %d (from %d features with pattern)\n",
                    dominant_offset, valid_offsets);
        } else {
            fprintf(stderr, "[offset-detect] No pattern offsets found, using default offset: 0\n");
        }
    }

    return PF_OK;
}

static pf_error pf_create_directory_if_needed(pf_context *ctx, const char *path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        if (mkdir(path, 0755) != 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to create output directory: %.900s", path);
            return PF_ERR_IO_ERROR;
        }
    }
    return PF_OK;
}

static void pf_default_quality(char *dst, size_t sequence_len) {
    if (sequence_len >= LINE_LENGTH - 1) {
        sequence_len = LINE_LENGTH - 2;
    }
    memset(dst, 'I', sequence_len);
    dst[sequence_len] = '\n';
    dst[sequence_len + 1] = '\0';
}

static pf_sequence_view pf_view_from_cstr(const char *src) {
    pf_sequence_view view;
    view.data = src;
    view.length = src ? strlen(src) : 0;
    return view;
}

static int pf_copy_record_view_line(pf_context *ctx,
                                    char *dst,
                                    pf_sequence_view view,
                                    const char *field_name,
                                    int required) {
    if (!view.data) {
        if (required || view.length > 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "In-memory record is missing %s", field_name);
            return 0;
        }
        dst[0] = '\0';
        return 1;
    }
    const int has_newline = view.length > 0 && view.data[view.length - 1] == '\n';
    if ((!has_newline && view.length >= LINE_LENGTH - 1) ||
        (has_newline && view.length >= LINE_LENGTH)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "In-memory record %s length %zu exceeds LINE_LENGTH=%d "
                 "after adding FASTQ-line newline",
                 field_name, view.length, LINE_LENGTH);
        return 0;
    }
    memcpy(dst, view.data, view.length);
    if (has_newline) {
        dst[view.length] = '\0';
    } else {
        dst[view.length] = '\n';
        dst[view.length + 1] = '\0';
    }
    return 1;
}

static int pf_validate_record_view_field(pf_context *ctx,
                                         pf_sequence_view view,
                                         const char *field_name,
                                         int required) {
    if (!view.data) {
        if (required || view.length > 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "In-memory record is missing %s", field_name);
            return 0;
        }
        return 1;
    }
    const int has_newline = view.length > 0 && view.data[view.length - 1] == '\n';
    if ((!has_newline && view.length >= LINE_LENGTH - 1) ||
        (has_newline && view.length >= LINE_LENGTH)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "In-memory record %s length %zu exceeds LINE_LENGTH=%d "
                 "after adding FASTQ-line newline",
                 field_name, view.length, LINE_LENGTH);
        return 0;
    }
    return 1;
}

static fastq_reader *pf_allocate_decoded_reader(int filetype) {
    fastq_reader *reader = (fastq_reader *)calloc(1, sizeof(fastq_reader));
    if (reader) {
        reader->filetype = filetype;
    }
    return reader;
}

static fastq_reader_set *pf_allocate_decoded_reader_set(pf_context *ctx,
                                                        int nreaders,
                                                        size_t read_size,
                                                        size_t read_buffer_lines) {
    if (nreaders != 2 && nreaders != 3) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Decoded process_features input supports 2 or 3 readers, got %d",
                 nreaders);
        return NULL;
    }

    fastq_reader_set *set = (fastq_reader_set *)calloc(1, sizeof(fastq_reader_set));
    if (!set) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Failed to allocate decoded process_features reader set");
        return NULL;
    }

    set->barcode_reader = pf_allocate_decoded_reader(1);
    set->forward_reader = pf_allocate_decoded_reader(2);
    set->reverse_reader = (nreaders == 3) ? pf_allocate_decoded_reader(3) : NULL;
    set->read_buffer_lines = read_buffer_lines;
    set->buffer_storage = (char *)malloc(read_buffer_lines * (read_size + 1));
    set->buffer = (char **)malloc(read_buffer_lines * sizeof(char *));

    if (!set->barcode_reader || !set->forward_reader ||
        (nreaders == 3 && !set->reverse_reader) ||
        !set->buffer_storage || !set->buffer) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Failed to allocate decoded process_features buffers");
        free_fastq_reader(set->barcode_reader);
        free_fastq_reader(set->forward_reader);
        free_fastq_reader(set->reverse_reader);
        free(set->buffer_storage);
        free(set->buffer);
        free(set);
        return NULL;
    }

    for (size_t i = 0; i < read_buffer_lines; ++i) {
        set->buffer[i] = set->buffer_storage + i * (read_size + 1);
    }

    pthread_mutex_init(&set->mutex, NULL);
    pthread_cond_init(&set->can_produce, NULL);
    pthread_cond_init(&set->can_consume, NULL);
    return set;
}

static void pf_stream_signal_done(pf_record_stream *stream) {
    if (!stream || !stream->queue_initialized || !stream->reader_sets[0]) {
        return;
    }
    fastq_reader_set *set = stream->reader_sets[0];
    pthread_mutex_lock(&set->mutex);
    set->done = 1;
    pthread_cond_broadcast(&set->can_consume);
    pthread_cond_broadcast(&set->can_produce);
    pthread_mutex_unlock(&set->mutex);
}

static void pf_stream_join_consumers(pf_record_stream *stream) {
    if (!stream || stream->consumers_joined) {
        return;
    }
    pf_stream_signal_done(stream);
    for (int i = 0; i < stream->consumers_started; ++i) {
        pthread_join(stream->consumer_threads[i], NULL);
    }
    stream->consumers_joined = 1;

    if (stream->chem_detect_enabled) {
        pf_context *ctx = stream->ctx;
        ctx->chem_detect_raw_hits = stream->chem_detect.raw_hits;
        ctx->chem_detect_nxt_hits = stream->chem_detect.nxt_hits;
        ctx->chem_detect_ticket = stream->chem_detect.ticket;
        ctx->chem_detect_done = stream->chem_detect.done;
        ctx->detected_match_mode = stream->chem_detect.match_mode;
    }
}

static void pf_stream_cleanup_runtime(pf_record_stream *stream) {
    if (!stream) {
        return;
    }

    if (stream->queue_initialized && !stream->consumers_joined) {
        pf_stream_join_consumers(stream);
    }

    if (stream->reader_sets[0]) {
        free_fastq_reader_set(stream->reader_sets[0]);
        stream->reader_sets[0] = NULL;
    }

    if (stream->processors) {
        for (int i = 0; i < stream->processors_initialized; ++i) {
            pthread_mutex_destroy(&stream->processors[i].process_mutex);
        }
        free(stream->processors);
        stream->processors = NULL;
    }

    free(stream->consumer_threads);
    stream->consumer_threads = NULL;

    if (stream->hashes) {
        for (int i = 0; i < stream->worker_state_initialized; ++i) {
            destroy_data_structures(&stream->hashes[i]);
        }
        free(stream->hashes);
        stream->hashes = NULL;
    }
    if (stream->pools) {
        for (int i = 0; i < stream->worker_state_initialized; ++i) {
            if (stream->pools[i]) {
                free_memory_pool_collection(stream->pools[i]);
            }
        }
        free(stream->pools);
        stream->pools = NULL;
    }
    free(stream->stats);
    stream->stats = NULL;
}

static int pf_stream_initialize_queue(pf_record_stream *stream, int nreaders) {
    if (stream->queue_initialized) {
        if (stream->nreaders != nreaders) {
            snprintf(stream->ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Mixed decoded record layouts in one process_features stream");
            return 0;
        }
        return 1;
    }

    pf_context *ctx = stream->ctx;
    const int nconsumers = ctx->config->consumer_threads > 0
        ? ctx->config->consumer_threads
        : 1;
    stream->nconsumers = nconsumers;
    stream->nreaders = nreaders;
    stream->lines_per_block = 2 * nreaders;

    stream->reader_sets[0] = pf_allocate_decoded_reader_set(
        ctx, nreaders, LINE_LENGTH, READ_BUFFER_LINES);
    if (!stream->reader_sets[0]) {
        return 0;
    }
    stream->reader_sets[0]->thread_id = 0;

    stream->stats = (statistics *)calloc((size_t)nconsumers, sizeof(statistics));
    stream->hashes = (data_structures *)calloc((size_t)nconsumers, sizeof(data_structures));
    stream->pools = (memory_pool_collection **)calloc((size_t)nconsumers,
                                                      sizeof(memory_pool_collection *));
    stream->processors = (fastq_processor *)calloc((size_t)nconsumers,
                                                   sizeof(fastq_processor));
    stream->consumer_threads = (pthread_t *)calloc((size_t)nconsumers,
                                                   sizeof(pthread_t));
    if (!stream->stats || !stream->hashes || !stream->pools ||
        !stream->processors || !stream->consumer_threads) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Failed to allocate decoded process_features consumer state");
        pf_stream_cleanup_runtime(stream);
        return 0;
    }

    for (int i = 0; i < nconsumers; ++i) {
        initialize_statistics(&stream->stats[i]);
        initialize_data_structures(&stream->hashes[i]);
        stream->pools[i] = initialize_memory_pool_collection();
    }
    stream->worker_state_initialized = nconsumers;

    memset(&stream->sample_args, 0, sizeof(stream->sample_args));
    stream->sample_args.sample_index = 0;
    stream->sample_args.directory = stream->sample_directory;
    stream->sample_args.filtered_barcodes_name = NULL;
    stream->sample_args.fastq_files = NULL;
    stream->sample_args.features = ctx->features;
    stream->sample_args.maxHammingDistance = ctx->config->max_hamming_distance;
    stream->sample_args.nThreads = ctx->config->search_threads;
    stream->sample_args.pools = stream->pools;
    stream->sample_args.stats = stream->stats;
    stream->sample_args.hashes = stream->hashes;
    stream->sample_args.stringency = ctx->config->stringency;
    stream->sample_args.min_counts = ctx->config->min_counts;
    stream->sample_args.read_buffer_lines = READ_BUFFER_LINES;
    stream->sample_args.average_read_length = AVERAGE_READ_LENGTH;
    stream->sample_args.barcode_constant_offset = ctx->config->barcode_offset;
    stream->sample_args.feature_constant_offset = ctx->config->feature_offset;
    stream->sample_args.min_posterior = ctx->config->min_posterior;
    stream->sample_args.legacy_cb_rescue = ctx->config->legacy_cb_rescue;
    stream->sample_args.consumer_threads_per_set = nconsumers;
    stream->sample_args.permit_acquire_hook = ctx->config->permit_acquire_cb;
    stream->sample_args.permit_release_hook = ctx->config->permit_release_cb;
    stream->sample_args.permit_hook_ctx = ctx->config->permit_hook_ctx;
    stream->sample_args.permit_hooks_enabled =
        (ctx->config->permit_acquire_cb != NULL &&
         ctx->config->permit_release_cb != NULL);
    stream->sample_args.filtered_barcodes_hash = ctx->filtered_barcodes_hash;
    stream->sample_args.min_prediction = 1;
    stream->sample_args.min_heatmap = 0;
    stream->sample_args.demux_nsamples = 1;
    stream->sample_args.sample_barcodes = NULL;
    stream->sample_args.sample_max_hamming = 1;
    stream->sample_args.sample_max_N = 0;
    stream->sample_args.sample_constant_offset = -1;
    stream->sample_args.sample_offset_relative = 0;
    stream->sample_args.skip_emptydrops = ctx->config->skip_emptydrops;
    stream->sample_args.emptydrops_failure_fatal = ctx->config->emptydrops_failure_fatal;
    stream->sample_args.expected_cells = ctx->config->expected_cells;
    stream->sample_args.emptydrops_use_fdr = ctx->config->emptydrops_use_fdr;
    stream->sample_args.probe_only = ctx->config->probe_only;
    stream->sample_args.skip_qc_outputs = ctx->config->skip_qc_outputs;

    if (ctx->config->autodetect_chemistry) {
        memset(&stream->chem_detect, 0, sizeof(stream->chem_detect));
        stream->chem_detect.max_reads = ctx->config->autodetect_chemistry_reads;
        stream->chem_detect.min_hits = ctx->config->autodetect_chemistry_min_hits;
        stream->sample_args.chem_detect = &stream->chem_detect;
        stream->reader_sets[0]->chem_detect = &stream->chem_detect;
        stream->chem_detect_enabled = 1;
    }

    stream->reader_sets[0]->probe_only = ctx->config->probe_only;

    for (int i = 0; i < nconsumers; ++i) {
        stream->processors[i].sample_args = &stream->sample_args;
        stream->processors[i].reader_sets = stream->reader_sets;
        stream->processors[i].nsets = 1;
        stream->processors[i].thread_id = i;
        stream->processors[i].nreaders = nreaders;
        pthread_mutex_init(&stream->processors[i].process_mutex, NULL);
        stream->processors_initialized++;
    }

    stream->queue_initialized = 1;

    for (int i = 0; i < nconsumers; ++i) {
        if (pthread_create(&stream->consumer_threads[i], NULL,
                           consume_reads, &stream->processors[i]) != 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to create decoded process_features consumer thread");
            pf_stream_join_consumers(stream);
            pf_stream_cleanup_runtime(stream);
            return 0;
        }
        stream->consumers_started++;
    }

    return 1;
}

static int pf_record_view_reader_count(pf_context *ctx,
                                       const pf_read_record_view *record,
                                       int *nreaders_out) {
    const int has_feature2 =
        (record->feature_sequence2.data || record->feature_sequence2.length > 0);

    if (!pf_validate_record_view_field(ctx, record->barcode_sequence,
                                       "barcode_sequence", 1) ||
        !pf_validate_record_view_field(ctx, record->barcode_quality,
                                       "barcode_quality", 0) ||
        !pf_validate_record_view_field(ctx, record->feature_sequence,
                                       "feature_sequence", 1) ||
        !pf_validate_record_view_field(ctx, record->feature_quality,
                                       "feature_quality", 0)) {
        return 0;
    }

    if (has_feature2) {
        if (!pf_validate_record_view_field(ctx, record->feature_sequence2,
                                           "feature_sequence2", 1) ||
            !pf_validate_record_view_field(ctx, record->feature_quality2,
                                           "feature_quality2", 0)) {
            return 0;
        }
        *nreaders_out = 3;
    } else {
        if (record->feature_quality2.data || record->feature_quality2.length > 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "feature_quality2 was provided without feature_sequence2");
            return 0;
        }
        *nreaders_out = 2;
    }

    return 1;
}

static void pf_copy_quality_or_default(char *dst,
                                       pf_sequence_view quality,
                                       size_t sequence_len) {
    if (quality.data || quality.length > 0) {
        memcpy(dst, quality.data, quality.length);
        if (quality.length > 0 && quality.data[quality.length - 1] == '\n') {
            dst[quality.length] = '\0';
        } else {
            dst[quality.length] = '\n';
            dst[quality.length + 1] = '\0';
        }
    } else {
        pf_default_quality(dst, sequence_len);
    }
}

static int pf_stream_push_record_view(pf_record_stream *stream,
                                      const pf_read_record_view *record) {
    pf_context *ctx = stream->ctx;
    int nreaders = 0;
    if (!pf_record_view_reader_count(ctx, record, &nreaders)) {
        return 0;
    }
    if (!pf_stream_initialize_queue(stream, nreaders)) {
        return 0;
    }

    if (stream->sample_args.probe_only &&
        stream->sample_args.chem_detect &&
        stream->sample_args.chem_detect->done) {
        return 1;
    }

    fastq_reader_set *set = stream->reader_sets[0];
    pthread_mutex_lock(&set->mutex);
    while (set->filled >= set->read_buffer_lines - (size_t)stream->lines_per_block) {
        pthread_cond_wait(&set->can_produce, &set->mutex);
    }

    const size_t p = set->produce_index;
    pf_copy_record_view_line(ctx, set->buffer[p],
                             record->barcode_sequence, "barcode_sequence", 1);
    pf_copy_quality_or_default(set->buffer[(p + 1) % set->read_buffer_lines],
                               record->barcode_quality,
                               record->barcode_sequence.length);

    pf_copy_record_view_line(ctx, set->buffer[(p + 2) % set->read_buffer_lines],
                             record->feature_sequence, "feature_sequence", 1);
    pf_copy_quality_or_default(set->buffer[(p + 3) % set->read_buffer_lines],
                               record->feature_quality,
                               record->feature_sequence.length);

    if (nreaders == 3) {
        pf_copy_record_view_line(ctx,
                                 set->buffer[(p + 4) % set->read_buffer_lines],
                                 record->feature_sequence2,
                                 "feature_sequence2", 1);
        pf_copy_quality_or_default(
            set->buffer[(p + 5) % set->read_buffer_lines],
            record->feature_quality2,
            record->feature_sequence2.length);
    }

    set->produce_index = (p + (size_t)stream->lines_per_block) % set->read_buffer_lines;
    set->filled += (size_t)stream->lines_per_block;
    pthread_cond_signal(&set->can_consume);
    pthread_mutex_unlock(&set->mutex);
    return 1;
}

pf_error pf_process_fastq_dir(pf_context *ctx,
                               const char *fastq_dir,
                               const char *output_dir,
                               pf_stats *stats_out) {
    if (!ctx || !fastq_dir || !output_dir) return PF_ERR_INVALID_ARG;
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;

    pthread_mutex_lock(&g_pf_runtime_mutex);
    pf_apply_context_globals(ctx);

    /* Reset detection state unconditionally to prevent stale leakage */
    ctx->chem_detect_raw_hits = 0;
    ctx->chem_detect_nxt_hits = 0;
    ctx->chem_detect_ticket = 0;
    ctx->chem_detect_done = 0;
    ctx->detected_match_mode = 0;
    if (!ctx->features) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Features not loaded");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_NOT_INITIALIZED;
    }
    if (!ctx->whitelist_hash) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Whitelist not loaded");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_NOT_INITIALIZED;
    }
    
    /* Feature offset preflight detection */
    if (ctx->config->use_feature_offset_array && ctx->config->feature_offset_explicit) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Cannot specify both use_feature_offset_array and explicit feature_offset");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_OFFSET_CONFLICT;
    }
    
    if (!ctx->config->use_feature_offset_array && !ctx->config->feature_offset_explicit &&
        ctx->features && feature_offsets && feature_offsets_count > 0) {
        /* Auto-detect: scan feature_offsets for heterogeneity */
        int offset_counts[256] = {0};
        int max_offset_seen = -1;
        int valid_offsets = 0;
        
        for (int i = 0; i < feature_offsets_count; i++) {
            int off = feature_offsets[i];
            if (off >= 0 && off < 256) {
                offset_counts[off]++;
                valid_offsets++;
                if (off > max_offset_seen) max_offset_seen = off;
            }
        }
        
        if (valid_offsets > 0) {
            /* Find dominant offset */
            int dominant_offset = 0;
            int dominant_count = 0;
            int second_count = 0;
            
            for (int i = 0; i <= max_offset_seen; i++) {
                if (offset_counts[i] > dominant_count) {
                    second_count = dominant_count;
                    dominant_count = offset_counts[i];
                    dominant_offset = i;
                } else if (offset_counts[i] > second_count) {
                    second_count = offset_counts[i];
                }
            }
            
            /* Check for heterogeneity: second offset > 5% of dominant */
            double heterogeneity_threshold = 0.05;
            if (second_count > 0 && (double)second_count / (double)dominant_count > heterogeneity_threshold) {
                if (ctx->config->strict_offset_check) {
                    snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                             "Multiple feature offsets detected (strict mode). "
                             "Dominant: %d (%d features), second: %d features. "
                             "Use pf_config_set_feature_offset() or pf_config_set_use_feature_offset_array(1).",
                             dominant_offset, dominant_count, second_count);
                    pthread_mutex_unlock(&g_pf_runtime_mutex);
                    return PF_ERR_MULTI_OFFSET_DETECTED;
                }
                /* Non-strict: warn and proceed with dominant */
                fprintf(stderr, "\nWARNING: Multiple feature offsets detected in pattern column.\n");
                fprintf(stderr, "         Dominant offset: %d (used by %d features)\n", dominant_offset, dominant_count);
                fprintf(stderr, "         Other offsets detected:\n");
                for (int i = 0; i <= max_offset_seen; i++) {
                    if (i != dominant_offset && offset_counts[i] > 0) {
                        double pct = 100.0 * offset_counts[i] / dominant_count;
                        fprintf(stderr, "           offset %d: %d features (%.1f%%)\n", i, offset_counts[i], pct);
                    }
                }
                fprintf(stderr, "\n         Proceeding with dominant offset %d.\n\n", dominant_offset);
            }
            
            /* Use dominant offset as global */
            ctx->config->feature_offset = dominant_offset;
            fprintf(stderr, "[offset-detect] Auto-detected global offset: %d (from %d features with pattern)\n",
                    dominant_offset, valid_offsets);
        } else {
            /* No valid offsets - default to 0 */
            fprintf(stderr, "[offset-detect] No pattern offsets found, using default offset: 0\n");
        }
    }
    
    if (!pf_is_directory(fastq_dir)) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "FASTQ directory not found: %s", fastq_dir);
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Create output directory if needed */
    struct stat st = {0};
    if (stat(output_dir, &st) == -1) {
        if (mkdir(output_dir, 0755) != 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to create output directory: %s", output_dir);
            pthread_mutex_unlock(&g_pf_runtime_mutex);
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
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_FILE_NOT_FOUND;
    }
    
    /* Process each sample */
    /* NXT/TRU auto-detection shared state */
    struct chem_detect_state chem_detect_buf;
    struct chem_detect_state *chem_detect_ptr = NULL;
    if (ctx->config->autodetect_chemistry) {
        memset(&chem_detect_buf, 0, sizeof(chem_detect_buf));
        chem_detect_buf.max_reads = ctx->config->autodetect_chemistry_reads;
        chem_detect_buf.min_hits = ctx->config->autodetect_chemistry_min_hits;
        chem_detect_ptr = &chem_detect_buf;
        fprintf(stderr, "NOTICE: chemistry auto-detect enabled for directory "
                "(aggregates across all samples; reads=%d, min_hits=%d)\n",
                chem_detect_buf.max_reads, chem_detect_buf.min_hits);
    }

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
        args.legacy_cb_rescue = ctx->config->legacy_cb_rescue;
        args.consumer_threads_per_set = ctx->config->consumer_threads;
        args.permit_acquire_hook = ctx->config->permit_acquire_cb;
        args.permit_release_hook = ctx->config->permit_release_cb;
        args.permit_hook_ctx = ctx->config->permit_hook_ctx;
        args.permit_hooks_enabled = (ctx->config->permit_acquire_cb != NULL &&
                                     ctx->config->permit_release_cb != NULL);
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
        args.chem_detect = chem_detect_ptr;
        args.probe_only = ctx->config->probe_only;
        args.skip_qc_outputs = ctx->config->skip_qc_outputs;
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

    /* Copy auto-detection results back to context */
    if (chem_detect_ptr) {
        ctx->chem_detect_raw_hits = chem_detect_ptr->raw_hits;
        ctx->chem_detect_nxt_hits = chem_detect_ptr->nxt_hits;
        ctx->chem_detect_ticket = chem_detect_ptr->ticket;
        ctx->chem_detect_done = chem_detect_ptr->done;
        ctx->detected_match_mode = chem_detect_ptr->match_mode;
    }
    
    if (any_failed) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "One or more samples failed during processing");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_IO_ERROR;
    }

    pthread_mutex_unlock(&g_pf_runtime_mutex);
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

    pthread_mutex_lock(&g_pf_runtime_mutex);
    pf_apply_context_globals(ctx);

    /* Reset detection state unconditionally to prevent stale leakage */
    ctx->chem_detect_raw_hits = 0;
    ctx->chem_detect_nxt_hits = 0;
    ctx->chem_detect_ticket = 0;
    ctx->chem_detect_done = 0;
    ctx->detected_match_mode = 0;

    if (!ctx->features) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Features not loaded");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_NOT_INITIALIZED;
    }
    if (!ctx->whitelist_hash) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Whitelist not loaded");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_NOT_INITIALIZED;
    }
    
    /* Feature offset preflight detection (same as pf_process_fastq_dir) */
    if (ctx->config->use_feature_offset_array && ctx->config->feature_offset_explicit) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Cannot specify both use_feature_offset_array and explicit feature_offset");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_OFFSET_CONFLICT;
    }
    
    if (!ctx->config->use_feature_offset_array && !ctx->config->feature_offset_explicit &&
        ctx->features && feature_offsets && feature_offsets_count > 0) {
        /* Auto-detect: scan feature_offsets for heterogeneity */
        int offset_counts[256] = {0};
        int max_offset_seen = -1;
        int valid_offsets = 0;
        
        for (int i = 0; i < feature_offsets_count; i++) {
            int off = feature_offsets[i];
            if (off >= 0 && off < 256) {
                offset_counts[off]++;
                valid_offsets++;
                if (off > max_offset_seen) max_offset_seen = off;
            }
        }
        
        if (valid_offsets > 0) {
            /* Find dominant offset */
            int dominant_offset = 0;
            int dominant_count = 0;
            int second_count = 0;
            
            for (int i = 0; i <= max_offset_seen; i++) {
                if (offset_counts[i] > dominant_count) {
                    second_count = dominant_count;
                    dominant_count = offset_counts[i];
                    dominant_offset = i;
                } else if (offset_counts[i] > second_count) {
                    second_count = offset_counts[i];
                }
            }
            
            /* Check for heterogeneity: second offset > 5% of dominant */
            double heterogeneity_threshold = 0.05;
            if (second_count > 0 && (double)second_count / (double)dominant_count > heterogeneity_threshold) {
                if (ctx->config->strict_offset_check) {
                    snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                             "Multiple feature offsets detected (strict mode). "
                             "Dominant: %d (%d features), second: %d features. "
                             "Use pf_config_set_feature_offset() or pf_config_set_use_feature_offset_array(1).",
                             dominant_offset, dominant_count, second_count);
                    pthread_mutex_unlock(&g_pf_runtime_mutex);
                    return PF_ERR_MULTI_OFFSET_DETECTED;
                }
                /* Non-strict: warn and proceed with dominant */
                fprintf(stderr, "\nWARNING: Multiple feature offsets detected in pattern column.\n");
                fprintf(stderr, "         Dominant offset: %d (used by %d features)\n", dominant_offset, dominant_count);
                fprintf(stderr, "         Other offsets detected:\n");
                for (int i = 0; i <= max_offset_seen; i++) {
                    if (i != dominant_offset && offset_counts[i] > 0) {
                        double pct = 100.0 * offset_counts[i] / dominant_count;
                        fprintf(stderr, "           offset %d: %d features (%.1f%%)\n", i, offset_counts[i], pct);
                    }
                }
                fprintf(stderr, "\n         Proceeding with dominant offset %d.\n\n", dominant_offset);
            }
            
            /* Use dominant offset as global */
            ctx->config->feature_offset = dominant_offset;
            fprintf(stderr, "[offset-detect] Auto-detected global offset: %d (from %d features with pattern)\n",
                    dominant_offset, valid_offsets);
        } else {
            /* No valid offsets - default to 0 */
            fprintf(stderr, "[offset-detect] No pattern offsets found, using default offset: 0\n");
        }
    }
    
    /* Create output directory */
    struct stat st = {0};
    if (stat(output_dir, &st) == -1) {
        if (mkdir(output_dir, 0755) != 0) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Failed to create output directory: %s", output_dir);
            pthread_mutex_unlock(&g_pf_runtime_mutex);
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
            pthread_mutex_unlock(&g_pf_runtime_mutex);
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
        pthread_mutex_unlock(&g_pf_runtime_mutex);
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
    args.legacy_cb_rescue = ctx->config->legacy_cb_rescue;
    args.consumer_threads_per_set = ctx->config->consumer_threads;
    args.permit_acquire_hook = ctx->config->permit_acquire_cb;
    args.permit_release_hook = ctx->config->permit_release_cb;
    args.permit_hook_ctx = ctx->config->permit_hook_ctx;
    args.permit_hooks_enabled = (ctx->config->permit_acquire_cb != NULL &&
                                 ctx->config->permit_release_cb != NULL);
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
    struct chem_detect_state chem_detect_buf2;
    struct chem_detect_state *chem_detect_ptr2 = NULL;
    if (ctx->config->autodetect_chemistry) {
        memset(&chem_detect_buf2, 0, sizeof(chem_detect_buf2));
        chem_detect_buf2.max_reads = ctx->config->autodetect_chemistry_reads;
        chem_detect_buf2.min_hits = ctx->config->autodetect_chemistry_min_hits;
        chem_detect_ptr2 = &chem_detect_buf2;
    }
    args.chem_detect = chem_detect_ptr2;
    args.probe_only = ctx->config->probe_only;
    args.skip_qc_outputs = ctx->config->skip_qc_outputs;
    args.error_out = &sample_error;
    
    process_files_in_sample(&args);

    if (chem_detect_ptr2) {
        ctx->chem_detect_raw_hits = chem_detect_ptr2->raw_hits;
        ctx->chem_detect_nxt_hits = chem_detect_ptr2->nxt_hits;
        ctx->chem_detect_ticket = chem_detect_ptr2->ticket;
        ctx->chem_detect_done = chem_detect_ptr2->done;
        ctx->detected_match_mode = chem_detect_ptr2->match_mode;
    }
    
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
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_IO_ERROR;
    }

    pthread_mutex_unlock(&g_pf_runtime_mutex);
    return PF_OK;
}

pf_error pf_process_records_begin(pf_context *ctx,
                                  const char *output_dir,
                                  const char *sample_name,
                                  pf_record_stream **stream_out) {
    if (stream_out) {
        *stream_out = NULL;
    }
    if (!ctx || !output_dir || !stream_out) {
        return PF_ERR_INVALID_ARG;
    }
    if (!ctx->initialized) return PF_ERR_NOT_INITIALIZED;

    pthread_mutex_lock(&g_pf_runtime_mutex);
    pf_apply_context_globals(ctx);

    ctx->chem_detect_raw_hits = 0;
    ctx->chem_detect_nxt_hits = 0;
    ctx->chem_detect_ticket = 0;
    ctx->chem_detect_done = 0;
    ctx->detected_match_mode = 0;

    if (!ctx->features) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Features not loaded");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_NOT_INITIALIZED;
    }
    if (!ctx->whitelist_hash) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE, "Whitelist not loaded");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_NOT_INITIALIZED;
    }

    pf_error prep = pf_prepare_feature_offset_config(ctx);
    if (prep != PF_OK) {
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return prep;
    }

    pf_error mkdir_err = pf_create_directory_if_needed(ctx, output_dir);
    if (mkdir_err != PF_OK) {
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return mkdir_err;
    }

    const char *sname = sample_name ? sample_name : "sample";
    pf_record_stream *stream = (pf_record_stream *)calloc(1, sizeof(pf_record_stream));
    if (!stream) {
        snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                 "Failed to allocate process_features record stream");
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return PF_ERR_OUT_OF_MEMORY;
    }

    stream->ctx = ctx;
    snprintf(stream->sample_directory, sizeof(stream->sample_directory),
             "%s/%s/", output_dir, sname);
    mkdir_err = pf_create_directory_if_needed(ctx, stream->sample_directory);
    if (mkdir_err != PF_OK) {
        free(stream);
        pthread_mutex_unlock(&g_pf_runtime_mutex);
        return mkdir_err;
    }

    stream->start_time = get_time_in_seconds();
    *stream_out = stream;
    return PF_OK;
}

pf_error pf_process_record_views(pf_record_stream *stream,
                                 const pf_read_record_view *records,
                                 size_t n_records) {
    if (!stream || (!records && n_records > 0) || stream->closed) {
        return PF_ERR_INVALID_ARG;
    }
    if (stream->failed) {
        return PF_ERR_INVALID_ARG;
    }

    for (size_t i = 0; i < n_records; ++i) {
        if (!pf_stream_push_record_view(stream, &records[i])) {
            stream->failed = 1;
            return PF_ERR_INVALID_ARG;
        }
    }

    return PF_OK;
}

pf_error pf_process_record_batch(pf_record_stream *stream,
                                 const pf_read_record *records,
                                 size_t n_records) {
    if (!stream || (!records && n_records > 0)) {
        return PF_ERR_INVALID_ARG;
    }
    for (size_t i = 0; i < n_records; ++i) {
        pf_read_record_view view;
        view.barcode_sequence = pf_view_from_cstr(records[i].barcode_sequence);
        view.barcode_quality = pf_view_from_cstr(records[i].barcode_quality);
        view.feature_sequence = pf_view_from_cstr(records[i].feature_sequence);
        view.feature_quality = pf_view_from_cstr(records[i].feature_quality);
        view.feature_sequence2 = pf_view_from_cstr(records[i].feature_sequence2);
        view.feature_quality2 = pf_view_from_cstr(records[i].feature_quality2);
        pf_error err = pf_process_record_views(stream, &view, 1);
        if (err != PF_OK) {
            return err;
        }
    }
    return PF_OK;
}

pf_error pf_process_records_end(pf_record_stream *stream,
                                pf_stats *stats_out) {
    if (!stream || stream->closed) {
        return PF_ERR_INVALID_ARG;
    }
    pf_context *ctx = stream->ctx;
    int sample_error = 0;
    pf_error result = PF_OK;

    stream->closed = 1;

    if (stream->failed) {
        result = PF_ERR_INVALID_ARG;
    } else if (!stream->queue_initialized) {
        if (!pf_stream_initialize_queue(stream, 2)) {
            result = PF_ERR_OUT_OF_MEMORY;
        }
    }

    if (result == PF_OK) {
        pf_stream_join_consumers(stream);
    }

    if (result == PF_OK && stream->nconsumers > 1) {
        for (int i = 1; i < stream->nconsumers; ++i) {
            merge_process_feature_thread_data(&stream->hashes[0],
                                              stream->pools[0],
                                              &stream->stats[0],
                                              &stream->hashes[i],
                                              stream->pools[i],
                                              &stream->stats[i]);
        }
    }

    if (result == PF_OK && !ctx->config->probe_only) {
        finalize_processing(ctx->features, &stream->hashes[0],
                            stream->sample_directory, stream->pools[0],
                            &stream->stats[0], ctx->config->stringency,
                            ctx->config->min_counts, ctx->config->min_posterior,
                            ctx->config->legacy_cb_rescue,
                            ctx->filtered_barcodes_hash,
                            ctx->config->skip_emptydrops,
                            ctx->config->emptydrops_failure_fatal,
                            ctx->config->expected_cells,
                            ctx->config->emptydrops_use_fdr,
                            ctx->config->skip_qc_outputs,
                            &sample_error);
        if (sample_error) {
            snprintf(ctx->error_buf, PF_ERROR_BUF_SIZE,
                     "Sample processing failed");
            result = PF_ERR_IO_ERROR;
        }
    }

    double end_time = get_time_in_seconds();
    if (stats_out) {
        memset(stats_out, 0, sizeof(pf_stats));
        if (stream->stats) {
            stats_out->total_reads = stream->stats[0].number_of_reads;
            stats_out->matched_reads = stream->stats[0].valid;
            stats_out->unmatched_reads = stream->stats[0].total_unmatched_features;
        }
        stats_out->total_features = ctx->features->number_of_features;
        stats_out->processing_time_sec = end_time - stream->start_time;
    }

    pf_stream_cleanup_runtime(stream);
    free(stream);
    pthread_mutex_unlock(&g_pf_runtime_mutex);
    return result;
}

void pf_process_records_abort(pf_record_stream *stream) {
    if (!stream || stream->closed) {
        return;
    }
    stream->closed = 1;
    pf_stream_cleanup_runtime(stream);
    free(stream);
    pthread_mutex_unlock(&g_pf_runtime_mutex);
}

pf_error pf_process_records(pf_context *ctx,
                            const pf_read_record *records,
                            size_t n_records,
                            const char *output_dir,
                            const char *sample_name,
                            pf_stats *stats_out) {
    if (!ctx || !records || n_records == 0 || !output_dir) {
        return PF_ERR_INVALID_ARG;
    }

    pf_record_stream *stream = NULL;
    pf_error err = pf_process_records_begin(ctx, output_dir, sample_name, &stream);
    if (err != PF_OK) {
        return err;
    }

    err = pf_process_record_batch(stream, records, n_records);
    if (err != PF_OK) {
        pf_process_records_abort(stream);
        return err;
    }

    return pf_process_records_end(stream, stats_out);
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

int pf_get_feature_no_ambiguity(pf_context *ctx, int level, int index) {
    if (!ctx || !ctx->features) return -1;
    if (index < 0 || index >= ctx->features->number_of_features) return -1;
    const unsigned char *arr = NULL;
    if (level == 1) arr = ctx->features->feature_no_ambiguity_le1;
    else if (level == 2) arr = ctx->features->feature_no_ambiguity_le2;
    if (!arr) return -1;
    return (int)arr[index];
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
