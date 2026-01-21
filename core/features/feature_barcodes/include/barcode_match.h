#ifndef BARCODE_MATCH_H
#define BARCODE_MATCH_H

#include "common.h"
#include "prototypes.h"
#include "globals.h"

/* Initializer to match main.c order */
static inline void barcode_match_init(void) {
    initseq2Code();
    initcode2seq();
    initdiff2hamming(diff2Hamming);
}

/* Shared helpers (moved from assignBarcodes.c) */
int split_line(char *line, char *fields[], const char *split_string);
char check_sequence(char *sequence, int sequence_length);
int string2code(char *string, int sequence_length, unsigned char *code);
void code2string(unsigned char *code, char *string, int length);
void initseq2Code(void);
void initcode2seq(void);
void initdiff2hamming(unsigned char *difference);
void free_feature_arrays(feature_arrays *features);

int feature_lookup_kmer(const char *seq, int len, const struct feature_arrays *fa, int direct_search);
int feature_lookup_code(const unsigned char *code, int code_len);

#endif /* BARCODE_MATCH_H */


