#ifndef UTILS_H
#define UTILS_H

#include "common.h"

guint hash_int64(gconstpointer v);
gboolean equal_int64(gconstpointer v1, gconstpointer v2);
guint hash_int32(gconstpointer v);
gboolean equal_int32(gconstpointer v1, gconstpointer v2);
void free_fastq_files_collection(fastq_files_collection *fastq_files);
int mkdir_p(const char *path);
#endif // UTILS_H