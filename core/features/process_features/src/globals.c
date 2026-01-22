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
khash_t(codeu32) *feature_code_hash;


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
int min_heatmap = -1;
int min_em_counts = 100;
int translate_NXT = 0;