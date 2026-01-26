#ifndef UTILS_H
#define UTILS_H

#include "common.h"

// These functions are no longer needed - khash provides built-in hash functions
void free_fastq_files_collection(fastq_files_collection *fastq_files);
int mkdir_p(const char *path);
#endif // UTILS_H