#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include "common.h"
#include <zlib.h>
#include <stdatomic.h>
#include "queue.h"

// Forward declarations

// Function prototypes
void initialize_complement();
void initialize_statistics(statistics *stats);
double get_time_in_seconds();
void lowerCaseDifferences(char *ref, char *query, int length);
int compare_feature_sequences(const void *a, const void *b);
void initseq2Code();
void initcode2seq();
void initdiff2hamming(unsigned char *difference);
void free_feature_arrays(feature_arrays *features);
void initialize_unit_sizes();
int is_directory(const char *path);
void read_unmatched_barcodes_features_block(unmatched_barcodes_features_block *entry_block, unmatched_barcodes_features *entry);
int insert_feature_sequence(char *sequence, uint32_t feature_index, unsigned char hamming_distance, uint16_t match_position, data_structures *hashes, memory_pool_collection *pools);
int print_feature_sequences(feature_arrays *features, int *total_counts, char *directory, data_structures *hashes);
unmatched_barcodes_features_block* add_unmatched_barcode_store_feature(unsigned char *barcodes, unsigned char* corrected_barcodes, char *umi, unsigned char *qscores, uint32_t feature_index, int number_of_variants, uint16_t match_position, memory_pool_collection *pools, statistics *stats);
unsigned char* read_whiteList(char *whitelist_filename,GHashTable *hash, int reverse_complement_flag);
int split_line(char *line, char *fields[], const char *split_string);
char check_sequence(char *sequence, int sequence_length);
int string2code(char *string, int sequence_length, unsigned char *code);
void string2all_codes(char *string, unsigned char codes[][LINE_LENGTH/2+1], int *lengths);
int check_barcode(char *barcode, unsigned char *code);
int find_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t feature_code_length, size_t feature_length);
int fuzzy_matching_codes(unsigned char *sequence_code, unsigned char *feature_code, size_t feature_code_length, size_t feature_length, int maxHammingDistance);
int find_matches_in_sub_arrays(unsigned char *sequence_code, unsigned char *feature_code, size_t sequence_code_length, size_t feature_code_length, size_t feature_length, int maxHammingDistance, int *best_offset);
int checkSequenceAndCorrectForN(char *line, char *corrected_lines[], char *buffer,int sequence_length, int maxN);
char* printToBuffer(char *string, int sequence_length, char *buffer);
int generate_sequences(char* string, int string_length, int *indices, char **output, int index, int number_of_indices, char *output_indices[], int number_of_outputs);
int simple_code_search(feature_arrays *features, char *sequence);
int simple_hash_search(feature_arrays *features, char *sequence);
int simple_search(feature_arrays *features, char *line);
int simple_hamming_search(feature_arrays *features, char *line, int maxHammingDistance, int *hamming_distance);
int find_feature_match_single(feature_arrays *features, char *lineR2, int maxHammingDistance,int *bestScore, char **matching_sequence, uint16_t *match_position);
int find_feature_match_parallel(feature_arrays *features, char *lineR2, int maxHammingDistance, int nThreads, int *bestScore, char **matching_sequence, uint16_t *match_position);
void update_feature_counts(char *barcodeString, char *umi, uint32_t feature_index, data_structures *hashes,  memory_pool_collection *pools);
void update_feature_counts_from_code(unsigned char *code, char *umi, uint32_t feature_index, data_structures *hashes, memory_pool_collection *pools);
void update_umi_counts(unsigned char *code, char *umi,  uint32_t feature_index,data_structures *hashes, memory_pool_collection *pools);
char check_neighbor(uint64_t code64,uint32_t *counts, data_structures *hashes);
int find_neighbors(uint64_t key64, uint64_t *neighbors, uint32_t *counts, data_structures *hashes);
void find_deduped_counts(data_structures *hashes, GHashTable* barcode_to_deduped_counts, uint16_t stringency, uint16_t min_counts);
void find_connected_component(gpointer start_key, uint32_t *counts, data_structures *hashes);
void add_deduped_count(GHashTable* temp_deduped_hash, uint32_t *clique_counts, uint16_t stringency, uint16_t min_counts);
void code2string(unsigned char *code, char *string, int length);
void printFeatureCounts(feature_arrays *features, int *deduped_counts, int *barcoded_counts, GArray **feature_hist, char *directory, data_structures *hashes, statistics *stats, GHashTable *barcode_to_deduped_hash, GHashTable *filtered_barcodes_hash);
int find_closest_barcodes(unsigned char* code,unsigned char *corrected_codes, unsigned char *indices);
int find_variant_match(unsigned char *code, int sequence_index, unsigned char *corrected_bases );
void process_pending_barcodes( data_structures *hashes, memory_pool_collection *pools, statistics *stats, double min_posterior);
unsigned char* find_best_posterior_match (unmatched_barcodes_features_block *entry_block, int number_of_features, double min_posterior, statistics *stats, data_structures *hashes);
int simpleCorrectFeature(char *line, feature_arrays *features, int maxN, int maxHammingDistance, int *hamming_distance);
int checkAndCorrectFeature(char *line, feature_arrays *features,int maxHammingDistance, int nThreads, int *hamming_distance, char *matching_sequence, int maxN,char *ambiguous, uint16_t *match_position);
size_t barcode_code2number(unsigned char *code);
int checkAndCorrectBarcode(char **lines, int maxN, uint32_t feature_index, uint16_t match_position, data_structures *hashes, memory_pool_collection *pools, statistics *stats, int barcode_constant_offset);
void finalize_processing(feature_arrays *features, data_structures *hashes,  char *directory, memory_pool_collection *pools, statistics *stats, uint16_t stringency, uint16_t min_counts, double min_posterior, GHashTable *filtered_barcodes_hash);
void open_fastq_files(const char *barcode_fastq, const char *forward_fastq, const char *reverse_fastq, gzFile *barcode_fastqgz, gzFile *forward_fastqgz, gzFile *reverse_fastqgz);
fastq_reader* allocate_fastq_reader( char **filenames, int nfiles, int filetype, size_t read_size, size_t read_buffer_lines);
fastq_reader_set *  allocate_fastq_reader_set( char **barcode_filenames, char **forward_filenames, char **reverse_filenames, int nfiles, size_t read_size, size_t read_buffer_lines);
void *read_fastqs_by_set(void *arg);
void *consume_reads(void *arg);
void free_fastq_reader(fastq_reader *reader);
void free_fastq_reader_set(fastq_reader_set *reader_set);
void free_fastq_processor(fastq_processor *processor_args);
void merge_stats(statistics *merged_stats, statistics *thread_stats);
void merge_feature_counts(gpointer key, gpointer value, gpointer user_data);
void merge_feature_umi_counts(gpointer key, gpointer value, gpointer user_data);
void merge_feature_sequences(gpointer key, gpointer value, gpointer user_data);
void merge_queues(Queue *dest_q, Queue *src_q);
void merge_unmatched_barcodes(unmatched_barcodes_features_block_list *merged_list, unmatched_barcodes_features_block_list *thread_list, memory_pool_collection *merged_pool);
void process_files_in_sample(sample_args *args);
void initialize_data_structures(data_structures *hashes);
void destroy_data_structures(data_structures *hashes);
int existing_output_skip(char keep_existing, char *directory);
void process_feature_sequence(char *sequence, feature_arrays *features, int maxHammingDistance, int nThreads, int feature_constant_offset, int max_feature_n, uint32_t *feature_index, int *hamming_distance, char *matching_sequence, uint16_t *match_position);
void process_multiple_feature_sequences(int nsequences, char **sequences, int *orientations, feature_arrays *features, int maxHammingDistance, int nThreads, int feature_constant_offset, int max_feature_n, uint32_t *feature_index, int *hamming_distance, char *matching_sequence, uint16_t *match_position);
void reverse_in_place(char *str);
char complement(char base);
void reverse_complement_in_place(char *seq);
void reverse_complement_sequence(char *sequence,  char *reverse, int length);
void translate_nxt_inplace(char *barcode, int len);
void cleanup_sample(memory_pool_collection *pools, data_structures *hashes);
void sort_samples_by_size(fastq_files_collection *fastq_files, int *sample_order);

// EM functions removed - no longer needed
void generate_heatmap(const char *directory, feature_arrays *features, int **coexpression_histograms);
void generate_deduped_heatmap(const char *directory, feature_arrays *features, int **deduped_histograms, int max_deduped_count, int *total_deduped_counts, int histogram_minimum_counts);
char* grab_sample_barcode(char *read, uint16_t feature_end, const sample_args *args);
int find_sample_index(char *seq, const sample_args *args, int *hamming_distance);


#endif // PROTOTYPES_H