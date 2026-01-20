#include "../include/common.h"
#include "../include/prototypes.h"
#include "../include/globals.h"

/*
  Canonical implementations moved from assignBarcodes.c so they can be
  shared by both the main pipeline and demux_fastq.
*/

int split_line(char *line, char *fields[], const char *split_string) {
    int count = 0;
    char *token;

    token = strtok(line, split_string);
    while (token != NULL) {
        fields[count++] = token;
        token = strtok(NULL, split_string);
    }
    return count;
}

char check_sequence(char *sequence, int sequence_length){
    for (int i=0; i<sequence_length; i++){
        if (sequence[i] == 'A' || sequence[i] == 'C' || sequence[i] == 'G' || sequence[i] == 'T'){
            continue;
        }
        return 0;
    }
    return 1;
}

int string2code(char *string, int sequence_length, unsigned char *code){
    const int last_element=sequence_length/4;
    for (int i=0; i<last_element; i++){
        code[i]=seq2code[(unsigned char)string[4*i]]<<6 | seq2code[(unsigned char)string[4*i+1]]<<4 | seq2code[(unsigned char)string[4*i+2]]<<2 | seq2code[(unsigned char)string[4*i+3]];
    }
    char leftover=sequence_length%4;
    if (!leftover){
        return last_element;
    }
    switch (leftover){
        case 1:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6;
            break;
        case 2:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6 | seq2code[(unsigned char)string[4*last_element+1]]<<4;
            break;
        case 3:
            code[last_element]=seq2code[(unsigned char)string[4*last_element]]<<6 | seq2code[(unsigned char)string[4*last_element+1]]<<4 | seq2code[(unsigned char)string[4*last_element+2]]<<2;
            break;
    }
    return last_element+1;
}

void code2string(unsigned char *code, char *string, int length){
    for (int i=0; i<length; i++){
        string[4*i]=code2seq[code[i]][0];
        string[4*i+1]=code2seq[code[i]][1];
        string[4*i+2]=code2seq[code[i]][2];
        string[4*i+3]=code2seq[code[i]][3];
    }
    string[4*length]='\0';
}

void initseq2Code(){
    for (int i=0; i<256; i++){
        seq2code[i]=0;
    }
    seq2code['A']=0;
    seq2code['C']=1;
    seq2code['G']=2;
    seq2code['T']=3;
}

void initcode2seq(){
    char temp[1025];
    strcpy(temp, LOOKUP_STRING);
    for (int i=0; i<256; i++){
        code2seq[i][0]=temp[4*i];
        code2seq[i][1]=temp[4*i+1];
        code2seq[i][2]=temp[4*i+2];
        code2seq[i][3]=temp[4*i+3];
    }
}

void initdiff2hamming(unsigned char *difference){
    memset(difference,0,256);
    for (int i=0; i<256; i++){
        difference[i]=0;
        unsigned char mask=3;
        for (int j=0; j<4; j++){
            if (i & mask){
                difference[i]++;
            }
            mask=mask << 2;
        }
    }
}

void free_feature_arrays(feature_arrays *features) {
    free(features->feature_names_storage);
    free(features->feature_lengths);
    free(features->feature_code_lengths);
    free(features->feature_sequences_storage);
    free(features->feature_codes_storage);
    free(features->feature_names);
    free(features->feature_sequences);
    free(features->feature_codes);
    free(features->mismatched_feature_indices);
    free(features);
}

void initialize_complement(){
    match['A']='T';
    match['T']='A';
    match['C']='G';
    match['G']='C';
    match['N']='N';
}

int feature_lookup_code(const unsigned char *code, int code_len) {
    // Hash look-up path (current behavior)
    GBytes *key = g_bytes_new_static(code, code_len);
    gpointer val = g_hash_table_lookup(feature_code_hash, key);
    g_bytes_unref(key);
    if (val) return GPOINTER_TO_UINT(val);
    return 0;
}

// Precomputed 8-mer table for direct search
static const struct feature_arrays *fa_cached_8 = NULL;
static uint64_t *fa_u64 = NULL;
static int fa_u64_n = 0;

static void ensure_feature_u64(const struct feature_arrays *fa) {
    if (!fa || fa->common_length != 8) return;
    if (fa_u64 && fa_cached_8 == fa && fa_u64_n == fa->number_of_features) return;

    free(fa_u64);
    fa_u64_n = fa->number_of_features;
    fa_u64 = (uint64_t*)malloc((size_t)fa_u64_n * sizeof(uint64_t));
    if (!fa_u64) { perror("malloc fa_u64"); exit(EXIT_FAILURE); }

    for (int i = 0; i < fa_u64_n; ++i) {
        uint64_t v = 0;
        if (fa->feature_lengths[i] == 8) memcpy(&v, fa->feature_sequences[i], 8);
        fa_u64[i] = v; // 0 for non-8 lengths; they won’t match
    }
    fa_cached_8 = fa;
}

int feature_lookup_kmer(const char *seq, int len, const struct feature_arrays *fa, int direct_search) {
    int use_direct = (direct_search == 1) ||
                     (direct_search < 0 && fa && fa->number_of_features <= 128);

    if (use_direct && fa) {
        if (len == 8 && fa->common_length == 8) {
            ensure_feature_u64(fa);
            uint64_t needle; memcpy(&needle, seq, 8); // safe load
            for (int i = 0; i < fa_u64_n; ++i) {
                if (fa->feature_lengths[i] != 8) continue;
                if (fa_u64[i] == needle) return i + 1;
            }
            return 0;
        } else {
            for (int i = 0; i < fa->number_of_features; ++i) {
                if (fa->feature_lengths[i] != (unsigned)len) continue;
                if (memcmp(fa->feature_sequences[i], seq, (size_t)len) == 0) return i + 1;
            }
            return 0;
        }
    }

    unsigned char code[(MAX_FEATURE_CODE_LENGTH > 0 ? MAX_FEATURE_CODE_LENGTH : 40)];
    int clen = string2code((char*)seq, len, code);
    return feature_lookup_code(code, clen);
}


