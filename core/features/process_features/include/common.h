#ifndef COMMON_H
#define COMMON_H

#include <ctype.h>
#include <getopt.h>
#include "khash_wrapper.h"
#include <pthread.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/wait.h>
#include <dirent.h>
#include <stdatomic.h>
#include <sys/mman.h>
#include <errno.h>

// defines that affect structs
#define BARCODE_LENGTH 16
#define UMI_LENGTH 12
#define BARCODE_CODE_LENGTH ((BARCODE_LENGTH + 3) / 4)
#define UMI_CODE_LENGTH ((UMI_LENGTH + 3) / 4)
#define MAX_FEATURES 128
#define MAX_BARCODE_MISMATCHES 3

/* NXT/TRU auto-detection shared state (carried via sample_args) */
struct chem_detect_state {
    unsigned long long raw_hits;
    unsigned long long nxt_hits;
    unsigned long long ticket;
    int done;           /* 0=sampling, 1=decided */
    int match_mode;     /* 0=unknown, 1=RAW_MATCH, 2=TRANSLATED_MATCH, 3=AMBIGUOUS */
    int max_reads;      /* N reads to sample */
    int min_hits;       /* minimum total hits for decision */
};

/* Lookup strategy for tiered feature search (set via PF_LOOKUP_STRATEGY env). */
typedef enum {
    PF_STRATEGY_DEFAULT = 0,      /* d0 exact -> d1/d2 integer prehash */
    PF_STRATEGY_HOT_D0  = 1,      /* hot d0 -> d0 exact -> d1/d2 integer prehash */
    PF_STRATEGY_HOT_D0_BLOOM = 2  /* reserved for future bloom-assisted variants */
} pf_lookup_strategy_t;

// other defines
#define MIN_POSTERIOR 0.975
#define MAX_FEATURE_N 1
#define MAX_BARCODE_N 1
#define MAX_MALLOCS 1000
#define MAX_FEATURE_SEQUENCE_LENGTH 64
#define READ_BUFFER_LINES 8192
#define AVERAGE_READ_LENGTH 300
#define FILENAME_LENGTH 2048
#define LINE_LENGTH 1024
#define NAME_LENGTH 40
#define MAX_FEATURE_CODE_LENGTH 40
#define FEATURE_NAME_LENGTH 32
#define LOOKUP_STRING "AAAAAAACAAAGAAATAACAAACCAACGAACTAAGAAAGCAAGGAAGTAATAAATCAATGAATTACAAACACACAGACATACCAACCCACCGACCTACGAACGCACGGACGTACTAACTCACTGACTTAGAAAGACAGAGAGATAGCAAGCCAGCGAGCTAGGAAGGCAGGGAGGTAGTAAGTCAGTGAGTTATAAATACATAGATATATCAATCCATCGATCTATGAATGCATGGATGTATTAATTCATTGATTTCAAACAACCAAGCAATCACACACCCACGCACTCAGACAGCCAGGCAGTCATACATCCATGCATTCCAACCACCCAGCCATCCCACCCCCCCGCCCTCCGACCGCCCGGCCGTCCTACCTCCCTGCCTTCGAACGACCGAGCGATCGCACGCCCGCGCGCTCGGACGGCCGGGCGGTCGTACGTCCGTGCGTTCTAACTACCTAGCTATCTCACTCCCTCGCTCTCTGACTGCCTGGCTGTCTTACTTCCTTGCTTTGAAAGAACGAAGGAATGACAGACCGACGGACTGAGAGAGCGAGGGAGTGATAGATCGATGGATTGCAAGCACGCAGGCATGCCAGCCCGCCGGCCTGCGAGCGCGCGGGCGTGCTAGCTCGCTGGCTTGGAAGGACGGAGGGATGGCAGGCCGGCGGGCTGGGAGGGCGGGGGGGTGGTAGGTCGGTGGGTTGTAAGTACGTAGGTATGTCAGTCCGTCGGTCTGTGAGTGCGTGGGTGTGTTAGTTCGTTGGTTTTAAATAACTAAGTAATTACATACCTACGTACTTAGATAGCTAGGTAGTTATATATCTATGTATTTCAATCACTCAGTCATTCCATCCCTCCGTCCTTCGATCGCTCGGTCGTTCTATCTCTCTGTCTTTGAATGACTGAGTGATTGCATGCCTGCGTGCTTGGATGGCTGGGTGGTTGTATGTCTGTGTGTTTTAATTACTTAGTTATTTCATTCCTTCGTTCTTTGATTGCTTGGTTGTTTTATTTCTTTGTTTT"
#define VISITED 255
#define CELL_SIZE 10
#define BAR_WIDTH 20
#define BASE_PADDING 10
#define BAR_GRAPH_HEIGHT 100

// Debug print function
#define DEBUG_PRINT(fmt, ...) \
    do { \
        const char *env_debug = getenv("DEBUG"); \
        if ((env_debug && atoi(env_debug) != 0) || debug) { \
            fprintf(stderr, fmt, ##__VA_ARGS__); \
        } \
    } while (0)



// Struct
//forward declare Queue - defined in queue.h
typedef struct _queue Queue;

// Structs for memory management
typedef struct _storage_block {
    struct _storage_block *next;
    unsigned char *storage;
} storage_block;

/* NBSignalCut struct removed - EM functionality no longer needed */

typedef struct feature_search_tables{
    struct feature_arrays *features;
    void *feature_code;
    unsigned char feature_code_length;
    seq_hash_t feature_code_hash;
    seq_key_mode_t code_hash_mode;
} feature_search_tables;

typedef struct feature_arrays {
    int number_of_features;
    int max_length;
    int common_length;
    char **feature_names;
    char *feature_names_storage;
    unsigned int *feature_lengths;
    unsigned char *feature_code_lengths;
    char **feature_sequences;
    char *feature_sequences_storage;
    unsigned char **feature_codes;
    unsigned char *feature_codes_storage;
    char **feature_anchors;
    char *feature_anchors_storage;
    unsigned int *feature_anchor_lengths;
    char **feature_suffix_anchors;
    char *feature_suffix_anchors_storage;
    unsigned int *feature_suffix_anchor_lengths;
    int *feature_offsets; /* 0-based array; entry i corresponds to feature index (i+1). -1 = unknown */
    int number_of_mismatched_features;
    int *mismatched_feature_indices;
    /* Per-library cumulative prehash tables (integer-keyed via seq_hash_t).
     * Keys are 2-bit encoded DNA sequences; no strdup, no per-key free. */
    seq_hash_t feature_hamming_le1_hash; /* <=1 */
    seq_hash_t feature_hamming_le2_hash; /* <=2 */
    unsigned char *feature_no_ambiguity_le1;   /* per feature (1=true) */
    unsigned char *feature_no_ambiguity_le2;   /* per feature (1=true) */
    int feature_hamming_le1_enabled;
    int feature_hamming_le2_enabled;
    /* Per-instance exact code hash (owned by this feature_arrays) */
    seq_hash_t feature_code_hash;
    seq_key_mode_t code_hash_mode;
    int code_hash_fixed_length;
} feature_arrays;

typedef struct feature_counts {
    unsigned char sequence_code[4];
    khash_t(u32u32) *counts;
} feature_counts;

typedef struct feature_sequences {
    uint32_t counts;
    char hamming_distance;
    uint32_t feature_index;
    uint16_t match_position;
    char sequence[];
} feature_sequences;

typedef struct feature_umi_counts {
    unsigned char sequence_umi_code[8];
    khash_t(u32u32) *counts;
} feature_umi_counts;

typedef struct unmatched_barcodes_features_block {
    struct unmatched_barcodes_features_block *next;
    unsigned char storage[];
} unmatched_barcodes_features_block;

typedef struct unmatched_barcodes_features {
    struct unmatched_barcodes_features_block *next;
    uint32_t feature_index;
    unsigned char *barcode;
    unsigned char *umi;
    unsigned char number_of_closest_barcodes;
    unsigned char *closest_barcodes;
    unsigned char *Qscores;
    uint16_t match_position;
} unmatched_barcodes_features;

typedef struct unmatched_barcodes_features_block_list {
    unmatched_barcodes_features_block *first_entry;
    unmatched_barcodes_features_block *last_entry;
} unmatched_barcodes_features_block_list;

typedef struct unit_sizes {
    size_t feature_counts;
    size_t feature_umi_counts;
    size_t feature_sequences;
    size_t unmatched_barcodes_features_block;
    size_t cb_probe_counts; /* per-CB probe-count array */
} unit_sizes;

typedef struct statistics {
    double start_time;
    size_t nMismatches;
    size_t recovered;
    size_t pending;
    size_t valid;
    size_t pending_recovered;
    size_t total_unmatched_features;
    size_t number_of_reads;
    size_t limited_exact_checks;
    size_t limited_exact_hits;
    size_t limited_simple_fallback_calls;
    size_t limited_simple_fallback_hits;
    size_t limited_full_fallback_calls;
    size_t limited_full_fallback_hits;
    size_t resolve_calls_total;
    size_t resolve_calls_multi_alt;
    size_t resolve_too_many_n;
    size_t resolve_no_hit;
    size_t resolve_ambiguous_tie;
    size_t resolve_ambiguous_no_feature;
    size_t resolve_resolved;
    size_t resolve_multi_alt_too_many_n;
    size_t resolve_multi_alt_no_hit;
    size_t resolve_multi_alt_ambiguous_tie;
    size_t resolve_multi_alt_ambiguous_no_feature;
    size_t resolve_multi_alt_resolved;
    unmatched_barcodes_features_block_list unmatched_list;
} statistics;


typedef struct data_structures {
    khash_t(u32ptr) *filtered_hash;
    khash_t(u64ptr) *sequence_umi_hash;
    khash_t(strptr) *unique_features_match;
    Queue *neighbors_queue;
} data_structures;

typedef struct fastq_files_collection {
    char **barcode_fastq;
    char **forward_fastq;
    char **reverse_fastq;
    char *concatenated_files;
    int nbarcode_files;
    int nforward_files;
    int nreverse_files;
    char **sample_names;
    char *concatenated_sample_names;
    int *sample_sizes;
    int *sample_offsets;
    int nsamples;
    int max_sample_size;
    int *sorted_index;
} fastq_files_collection;

//need to have this here for sample_args
typedef struct _memory_pool {
    size_t block_size;
    size_t blocks_per_pool;
    size_t free_blocks;
    void *next_free;
    storage_block *first_block;
    storage_block *current_block;
} memory_pool;

typedef struct _memory_pool_collection {
    memory_pool *feature_counts_pool;
    memory_pool *feature_umi_counts_pool;
    memory_pool *feature_sequences_pool;
    memory_pool *unmatched_barcodes_features_block_pool;
    memory_pool *cb_counts_pool; /* new pool for CB×probe count arrays */
} memory_pool_collection;

typedef struct sample_args {
    int sample_index;
    char *directory;
    char *filtered_barcodes_name;
    fastq_files_collection *fastq_files;
    feature_arrays *features;
    int maxHammingDistance;
    int nThreads;
    memory_pool_collection **pools;
    statistics *stats;
    data_structures *hashes;
    uint16_t stringency;
    uint16_t min_counts;
    int read_buffer_lines;
    int average_read_length;
    int barcode_constant_offset;
    int feature_constant_offset;
    int parallel_by_file;
    double min_posterior;
    int legacy_cb_rescue;            /* 1 = use legacy order-dependent pending rescue */
    int consumer_threads_per_set;
    uint64_t (*permit_acquire_hook)(void *hook_ctx);
    void (*permit_release_hook)(void *hook_ctx, uint64_t wait_ns, uint64_t work_units, uint64_t work_bytes, uint64_t work_ns);
    void *permit_hook_ctx;
    int permit_hooks_enabled;
    khash_t(strptr) *filtered_barcodes_hash;
    int heatmap_minimum_counts;
    int min_prediction;
    int min_heatmap;
    int demux_nsamples; // number of demultiplexed samples (1 = legacy)
    data_structures **sample_hashes; // [demux_nsamples][threads]
    statistics **sample_stats;       // per sample, per thread
    memory_pool_collection **sample_pools; // per sample, per thread
    /* Sample barcode demux config (Phase 1-4) */
    int sample_max_hamming;          /* default 1 */
    int sample_max_N;                /* default 0 */
    int sample_constant_offset;      /* >=0 absolute offset; -1 unused */
    int sample_offset_relative;      /* negative/positive offset from feature end; 0 unused */
    feature_arrays *sample_barcodes; /* loaded sample barcode list */
    
    /* EmptyDrops control (copied from pf_config) */
    int skip_emptydrops;             /* 1 = skip EmptyDrops entirely */
    int emptydrops_failure_fatal;    /* 1 = treat ED failure as error */
    int expected_cells;              /* 0 = auto-detect */
    int emptydrops_use_fdr;          /* 1 = use FDR gate for tail rescue */
    
    /* NXT/TRU auto-detection state (shared with consumer threads) */
    struct chem_detect_state *chem_detect;
    int probe_only;                  /* 1 = chemistry probe pass, suppress outputs */
    int skip_qc_outputs;             /* 1 = skip feature histograms/heatmaps */

    /* Error propagation */
    int *error_out;                  /* Set to non-zero if fatal error occurred */
} sample_args;

typedef struct fastq_reader {
    char *concatenated_filenames;
    char **filenames;
    int nfiles;
    gzFile gz_pointer;
    int filetype;
} fastq_reader;

typedef struct fastq_reader_set {
    int thread_id;
    struct fastq_reader *barcode_reader;
    struct fastq_reader *forward_reader;
    struct fastq_reader *reverse_reader;
    char   **buffer;
    char    *buffer_storage;
    size_t   read_buffer_lines;
    size_t   produce_index;
    size_t   consume_index;
    size_t   filled;
    pthread_mutex_t mutex;
    pthread_cond_t  can_produce;
    pthread_cond_t  can_consume;
    int done;
    struct chem_detect_state *chem_detect;
    int probe_only;
} fastq_reader_set;

typedef struct fastq_processor {
    int thread_id; // Add thread_id to identify the consumer thread
    int nsets;
    int nreaders;
    struct fastq_reader_set **reader_sets;
    struct sample_args *sample_args;
    pthread_mutex_t process_mutex;
    pthread_cond_t can_process;
} fastq_processor;





#endif // COMMON_H
