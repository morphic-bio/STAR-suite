#ifndef IO_H
#define IO_H

#include "common.h"
#include "prototypes.h"
#include "utils.h"
#include "memory.h"
#include "globals.h"
// Function prototypes for file I/O operations
feature_arrays* read_features_file(const char* filename);
int get_feature_line_sizes(char *line, int nameIndex, int seqIndex, int *name_size, int *seq_size, int *code_size, int *maxFeatureLength);
void process_feature_line(char *line, int nameIndex, int seqIndex, feature_arrays *myfeatures, int count);
feature_arrays* allocate_feature_arrays(int name_size, int seq_size, int code_size, int count, int maxFeatureLength);
void find_name_and_sequence_fields(char *line, int *nameIndex, int *seqIndex);
int put_fastq_files_string_into_collection(char *fastqFilesString, char **fastq_files, int *nFiles, char *concatenated_fastq);
void check_filecounts(fastq_files_collection *fastq_files);
int count_character(char *string, char character);
int compare_filenames(const void *a, const void *b);
int count_files_with_pattern(const char *directory_path, const char *pattern);
void organize_fastq_files_by_directory(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern);
void organize_fastq_files_by_type(int positional_arg_count, int argc, char *argv[], int optind, char *barcodeFastqFilesString, char *forwardFastqFilesString, char *reverseFastqFilesString, fastq_files_collection *fastq_files, char *barcode_pattern, char *forward_pattern, char *reverse_pattern, int sample_flag);
void sort_samples_by_size(fastq_files_collection *fastq_files, int *sample_order);
size_t get_file_size(char *filepath);
int file_exists(const char *filename);
void read_barcodes_into_hash(char *filename, GHashTable *hash);
const char* get_basename(const char *path);

#endif // IO_H