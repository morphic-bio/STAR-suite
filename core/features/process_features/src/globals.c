#include "../include/globals.h"

// Global variables definitions
unsigned char seq2code[256];
char code2seq[256][4];
unsigned char diff2Hamming[256];
unsigned char match[256];
unit_sizes dynamic_struct_sizes;

int debug;

unsigned char *whitelist;
khash_t(u32ptr) *whitelist_hash; 
seq_hash_t feature_code_hash = { .h64 = NULL };
seq_key_mode_t feature_code_hash_mode = SEQ_KEY_64;
int feature_code_hash_fixed_length = 0;
int feature_prehash_max_hamming = 2;
unsigned long long feature_prehash_max_entries = 50000000ULL;
unsigned long long feature_prehash_memory_budget = 0;

unsigned long long prehash_detect_memory_budget(void) {
    unsigned long long mem_bytes = 0;

    /* Try cgroup v2 memory.max */
    FILE *f = fopen("/sys/fs/cgroup/memory.max", "r");
    if (f) {
        char buf[64];
        if (fgets(buf, sizeof(buf), f)) {
            if (strncmp(buf, "max", 3) != 0) {
                mem_bytes = strtoull(buf, NULL, 10);
            }
        }
        fclose(f);
    }

    /* Try cgroup v1 memory.limit_in_bytes */
    if (mem_bytes == 0) {
        f = fopen("/sys/fs/cgroup/memory/memory.limit_in_bytes", "r");
        if (f) {
            char buf[64];
            if (fgets(buf, sizeof(buf), f)) {
                unsigned long long limit = strtoull(buf, NULL, 10);
                /* cgroup v1 returns a very large number when unlimited */
                if (limit < (1ULL << 62)) {
                    mem_bytes = limit;
                }
            }
            fclose(f);
        }
    }

    /* Fallback to /proc/meminfo MemAvailable */
    if (mem_bytes == 0) {
        f = fopen("/proc/meminfo", "r");
        if (f) {
            char line[256];
            while (fgets(line, sizeof(line), f)) {
                unsigned long long val;
                if (sscanf(line, "MemAvailable: %llu kB", &val) == 1) {
                    mem_bytes = val * 1024ULL;
                    break;
                }
            }
            fclose(f);
        }
    }

    /* 75% safety factor */
    if (mem_bytes > 0) {
        return mem_bytes * 3ULL / 4ULL;
    }

    /* If all detection fails, return 0 (no budget enforcement) */
    return 0;
}


int barcode_length=BARCODE_LENGTH;
int barcode_code_length=BARCODE_CODE_LENGTH;
int number_of_features;
int maximum_feature_length;
int feature_code_length;

int max_feature_n = MAX_FEATURE_N;
int max_barcode_n = MAX_BARCODE_N;
int max_barcode_mismatches = MAX_BARCODE_MISMATCHES;
int umi_length = UMI_LENGTH;
int umi_code_length = UMI_CODE_LENGTH;
long long max_reads = 0;
int limit_search = -1;
int feature_limited_fallback_mode = 0;
int min_heatmap = -1;
int min_em_counts = 100;
int translate_NXT = 0;

int use_feature_offset_array = 0;
int use_feature_anchor_search = 0;
int require_feature_anchor_match = 0;
int *feature_offsets = NULL;
int feature_offsets_count = 0;
char **feature_anchors = NULL;
unsigned int *feature_anchor_lengths = NULL;
char **feature_suffix_anchors = NULL;
unsigned int *feature_suffix_anchor_lengths = NULL;
int feature_anchor_count = 0;
int feature_mode_bootstrap_reads = 0;
unsigned long long feature_mode_reads_seen = 0;
int feature_mode_bootstrap_done = 0;
int feature_mode_max_offset = LINE_LENGTH;
int *feature_mode_offsets = NULL;
unsigned int *feature_mode_hist = NULL;

void clear_feature_lookup_hashes(void) {
    /* Global hash is now a non-owning alias; just null the pointer.
     * The actual hash is owned by feature_arrays and freed there. */
    feature_code_hash.h64 = NULL;
    feature_code_hash.h128 = NULL;
}
