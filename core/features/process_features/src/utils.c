#include "../include/utils.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

// Hash functions removed - khash provides built-in hash functions
void free_fastq_files_collection(fastq_files_collection *fastq_files){
    //TODO make memory management consistent between organize by directory and type
    //for now just freeing what is obviously allocated
    free(fastq_files->concatenated_files);
    free(fastq_files->concatenated_sample_names);
    // Free individual file path strings before freeing the arrays
    if(fastq_files->barcode_fastq) {
        for (int i = 0; i < fastq_files->nbarcode_files; i++) {
            free(fastq_files->barcode_fastq[i]);
        }
        free(fastq_files->barcode_fastq);
    }
    if(fastq_files->forward_fastq) {
        for (int i = 0; i < fastq_files->nforward_files; i++) {
            free(fastq_files->forward_fastq[i]);
        }
        free(fastq_files->forward_fastq);
    }
    if(fastq_files->reverse_fastq) {
        for (int i = 0; i < fastq_files->nreverse_files; i++) {
            free(fastq_files->reverse_fastq[i]);
        }
        free(fastq_files->reverse_fastq);
    }
    // Free individual sample name strings before freeing the array
    if(fastq_files->sample_names) {
        for (int i = 0; i < fastq_files->nsamples; i++) {
            free(fastq_files->sample_names[i]);
        }
        free(fastq_files->sample_names);
    }
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

int pf_copy_file(const char *src_path, const char *dst_path) {
    if (!src_path || !dst_path) return -1;
    FILE *src = fopen(src_path, "rb");
    if (!src) return -1;
    FILE *dst = fopen(dst_path, "wb");
    if (!dst) {
        fclose(src);
        return -1;
    }
    char buf[8192];
    size_t nread;
    while ((nread = fread(buf, 1, sizeof(buf), src)) > 0) {
        if (fwrite(buf, 1, nread, dst) != nread) {
            fclose(src);
            fclose(dst);
            return -1;
        }
    }
    fclose(src);
    fclose(dst);
    return 0;
}

int pf_file_fingerprint(const char *path, char *out_hex, size_t out_hex_len) {
    if (!path || !out_hex || out_hex_len < 17) return -1;
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;
    unsigned long long hash = 1469598103934665603ULL;
    char buf[8192];
    size_t nread;
    while ((nread = fread(buf, 1, sizeof(buf), fp)) > 0) {
        for (size_t i = 0; i < nread; ++i) {
            hash ^= (unsigned char)buf[i];
            hash *= 1099511628211ULL;
        }
    }
    fclose(fp);
    snprintf(out_hex, out_hex_len, "%016llx", hash);
    return 0;
}