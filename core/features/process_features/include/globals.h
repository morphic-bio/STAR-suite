#ifndef GLOBALS_H
#define GLOBALS_H

#include "common.h"

// Use 'extern' to declare global variables, making them accessible across files.
extern unsigned char seq2code[256];
extern char code2seq[256][4];
extern unsigned char diff2Hamming[256];
extern unsigned char match[256];
extern unit_sizes dynamic_struct_sizes;

extern int debug;

// Valid barcode list and hash
extern unsigned char *whitelist;
extern khash_t(u32ptr) *whitelist_hash; 

extern seq_hash_t feature_code_hash;
extern seq_key_mode_t feature_code_hash_mode;
extern int feature_code_hash_fixed_length;
extern int feature_prehash_max_hamming;
extern unsigned long long feature_prehash_max_entries;
extern unsigned long long feature_prehash_memory_budget;

unsigned long long prehash_detect_memory_budget(void);

// Size globals to replace constants
extern int barcode_length;
extern int barcode_code_length;
extern int number_of_features;
extern int maximum_feature_length;
extern int feature_code_length;

// Default values for the program
extern int max_feature_n;
extern int max_barcode_n;
extern int max_barcode_mismatches;
extern int umi_length;
extern int umi_code_length;
extern long long max_reads;
extern int limit_search;
/* 0=full fallback (checkAndCorrectFeature), 1=simple fallback (simpleCorrectFeature at limited offsets) */
extern int feature_limited_fallback_mode;
extern int min_heatmap;

// Translation flag for final-stage barcode output/filters
extern int translate_NXT;

// Feature offset control
extern int use_feature_offset_array;
extern int use_feature_anchor_search;
extern int require_feature_anchor_match;
extern int *feature_offsets;
extern int feature_offsets_count;
extern char **feature_anchors;
extern unsigned int *feature_anchor_lengths;
extern char **feature_suffix_anchors;
extern unsigned int *feature_suffix_anchor_lengths;
extern int feature_anchor_count;
extern int feature_mode_bootstrap_reads;
extern unsigned long long feature_mode_reads_seen;
extern int feature_mode_bootstrap_done;
extern int feature_mode_max_offset;
extern int *feature_mode_offsets;
extern unsigned int *feature_mode_hist;
extern int use_hot_hash;
extern int skip_heatmaps;

/* ADT/protein 10x MEX output mode */
extern int adt_mex_output;
extern char adt_feature_ref_path[4096];
extern char adt_command_line[8192];

/* Hash / HTO / CMO demux settings (adt_mex extension) */
extern int hash_demux_mode;
extern char hash_feature_selector[256];
extern char hash_demux_method[64];
extern char library_feature_type_cli[128];
extern int hash_min_total;
extern int hash_min_top;
extern double hash_min_ratio;

/* Hash lifecycle helpers for feature matching tables */
void clear_feature_lookup_hashes(void);


#endif // GLOBALS_H
