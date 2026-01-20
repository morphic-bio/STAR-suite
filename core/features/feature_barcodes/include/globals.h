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
extern GHashTable *whitelist_hash; 

extern GHashTable *feature_code_hash;

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
extern int min_heatmap;

// Translation flag for final-stage barcode output/filters
extern int translate_NXT;


#endif // GLOBALS_H