#ifndef UTILS_H
#define UTILS_H

#include "common.h"

// These functions are no longer needed - khash provides built-in hash functions
void free_fastq_files_collection(fastq_files_collection *fastq_files);
int mkdir_p(const char *path);
int pf_copy_file(const char *src_path, const char *dst_path);
int pf_file_fingerprint(const char *path, char *out_hex, size_t out_hex_len);
#endif // UTILS_H