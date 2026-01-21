#include "../include/utils.h"

guint hash_int64(gconstpointer v) {
    guint64 k = *(const guint64*)v;
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    return k;
}

gboolean equal_int64(gconstpointer v1, gconstpointer v2) {
    return *(const long long*)v1 == *(const long long*)v2;
}

// Helper functions for GHashTable, moved from assignBarcodes.c
guint hash_int32(gconstpointer v) {
    guint32 k = *(const guint32*) v;
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    return k;
}
gboolean equal_int32(gconstpointer v1, gconstpointer v2) {
     return *((const gint*)v1) == *((const gint*)v2);
}
void free_fastq_files_collection(fastq_files_collection *fastq_files){
    //TODO make memory management consistent between organize by directory and type
    //for now just freeing what is obviously allocated
    free(fastq_files->concatenated_files);
    free(fastq_files->concatenated_sample_names);
    if(fastq_files->barcode_fastq)
        free(fastq_files->barcode_fastq);
    if(fastq_files->forward_fastq)
        free(fastq_files->forward_fastq);
    if(fastq_files->reverse_fastq)
        free(fastq_files->reverse_fastq);
    if(fastq_files->sample_names)
        free(fastq_files->sample_names);
    if(fastq_files->sample_sizes)
        free(fastq_files->sample_sizes);
    if(fastq_files->sample_offsets)
        free(fastq_files->sample_offsets);
    if(fastq_files->sorted_index)
        free(fastq_files->sorted_index);
}
int mkdir_p(const char *path) {
    char temp[1024];
    char *p = NULL;
    size_t len;

    // Copy path and ensure it ends with '/'
    snprintf(temp, sizeof(temp), "%s", path);
    len = strlen(temp);
    if (temp[len - 1] == '/')
        temp[len - 1] = 0;

    // Iterate through each directory in the path
    for (p = temp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;

            // Create directory if it doesn't exist
            if (mkdir(temp, S_IRWXU) != 0 && errno != EEXIST) {
                perror("mkdir");
                return -1;
            }
            *p = '/';
        }
    }
    // Create the final directory
    if (mkdir(temp, S_IRWXU) != 0 && errno != EEXIST) {
        perror("mkdir");
        return -1;
    }
    return 0;
}