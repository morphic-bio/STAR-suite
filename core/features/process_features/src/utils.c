#include "../include/utils.h"

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