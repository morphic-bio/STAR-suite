/**
 * @file call_features.c
 * @brief Implementation of dominant feature calling
 */

#include "../include/call_features.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include <inttypes.h>
#include <zlib.h>

#define MAX_LINE_LENGTH 4096

typedef struct {
    int is_gz;
    FILE *fp;
    gzFile gz;
} cf_file;

static int file_exists(const char *path) {
    struct stat st;
    return (stat(path, &st) == 0);
}

static int file_exists_with_gz(const char *path) {
    if (file_exists(path)) return 1;
    char gz_path[MAX_LINE_LENGTH];
    snprintf(gz_path, sizeof(gz_path), "%s.gz", path);
    return file_exists(gz_path);
}

static int open_text_file(const char *path, cf_file *out) {
    if (!out) return -1;
    memset(out, 0, sizeof(*out));

    FILE *fp = fopen(path, "r");
    if (fp) {
        out->fp = fp;
        out->is_gz = 0;
        return 0;
    }

    gzFile gz = gzopen(path, "rb");
    if (gz) {
        out->gz = gz;
        out->is_gz = 1;
        return 0;
    }

    if (!strstr(path, ".gz")) {
        char gz_path[MAX_LINE_LENGTH];
        snprintf(gz_path, sizeof(gz_path), "%s.gz", path);
        gz = gzopen(gz_path, "rb");
        if (gz) {
            out->gz = gz;
            out->is_gz = 1;
            return 0;
        }
    }

    return -1;
}

static char *cf_gets(cf_file *fh, char *buf, int size) {
    if (fh->is_gz) {
        return gzgets(fh->gz, buf, size);
    }
    return fgets(buf, size, fh->fp);
}

static void close_text_file(cf_file *fh) {
    if (!fh) return;
    if (fh->is_gz) {
        if (fh->gz) gzclose(fh->gz);
    } else {
        if (fh->fp) fclose(fh->fp);
    }
}

/* ============================================================================
 * Configuration Implementation
 * ============================================================================ */

cf_config* cf_config_create(void) {
    cf_config *config = calloc(1, sizeof(cf_config));
    if (!config) return NULL;
    
    /* Default values */
    config->min_deduped_counts = 2;
    config->dominance_fraction = 0.8;
    config->dominance_margin = 1;
    config->include_ambiguous = 1;
    
    return config;
}

void cf_config_destroy(cf_config *config) {
    if (config) free(config);
}

/* ============================================================================
 * MEX Matrix Reader Implementation
 * ============================================================================ */

/* Helper: count lines in a file */
static int count_lines(const char *path) {
    cf_file fh;
    if (open_text_file(path, &fh) != 0) return -1;
    
    int count = 0;
    char line[MAX_LINE_LENGTH];
    while (cf_gets(&fh, line, sizeof(line))) {
        count++;
    }
    close_text_file(&fh);
    return count;
}

/* Helper: load names from a file (one per line) */
static int load_names(const char *path, char ***names_out, char **storage_out, int *count_out) {
    cf_file fh;
    if (open_text_file(path, &fh) != 0) return -1;
    
    /* Count lines first */
    int count = 0;
    char line[MAX_LINE_LENGTH];
    size_t total_len = 0;
    while (cf_gets(&fh, line, sizeof(line))) {
        size_t len = strlen(line);
        if (len > 0 && line[len-1] == '\n') line[--len] = '\0';
        if (len > 0 && line[len-1] == '\r') line[--len] = '\0';
        total_len += len + 1;
        count++;
    }
    
    if (count == 0) {
        close_text_file(&fh);
        *names_out = NULL;
        *storage_out = NULL;
        *count_out = 0;
        return 0;
    }
    
    /* Allocate storage */
    char *storage = malloc(total_len);
    char **names = malloc(count * sizeof(char*));
    if (!storage || !names) {
        free(storage);
        free(names);
        close_text_file(&fh);
        return -1;
    }
    
    /* Read names */
    close_text_file(&fh);
    if (open_text_file(path, &fh) != 0) {
        free(storage);
        free(names);
        return -1;
    }
    int i = 0;
    char *ptr = storage;
    while (cf_gets(&fh, line, sizeof(line)) && i < count) {
        size_t len = strlen(line);
        if (len > 0 && line[len-1] == '\n') line[--len] = '\0';
        if (len > 0 && line[len-1] == '\r') line[--len] = '\0';
        
        /* Handle tab-separated (features.tsv format) - take first column */
        char *tab = strchr(line, '\t');
        if (tab) *tab = '\0';
        len = strlen(line);
        
        /* Note: Keep barcode suffixes like -1 intact for CR9 compatibility */
        
        strcpy(ptr, line);
        names[i] = ptr;
        ptr += len + 1;
        i++;
    }
    
    close_text_file(&fh);
    *names_out = names;
    *storage_out = storage;
    *count_out = i;
    return 0;
}

cf_sparse_matrix* cf_load_mex(const char *mex_dir) {
    if (!mex_dir) return NULL;
    
    /* Build file paths */
    char matrix_path[MAX_LINE_LENGTH];
    char barcodes_path[MAX_LINE_LENGTH];
    char features_path[MAX_LINE_LENGTH];
    
    snprintf(matrix_path, sizeof(matrix_path), "%s/matrix.mtx", mex_dir);
    snprintf(barcodes_path, sizeof(barcodes_path), "%s/barcodes.txt", mex_dir);
    snprintf(features_path, sizeof(features_path), "%s/features.txt", mex_dir);
    
    /* Check if files exist (try .tsv variants) */
    if (!file_exists_with_gz(barcodes_path)) {
        snprintf(barcodes_path, sizeof(barcodes_path), "%s/barcodes.tsv", mex_dir);
    }
    if (!file_exists_with_gz(features_path)) {
        snprintf(features_path, sizeof(features_path), "%s/features.tsv", mex_dir);
    }
    
    /* Allocate matrix */
    cf_sparse_matrix *matrix = calloc(1, sizeof(cf_sparse_matrix));
    if (!matrix) return NULL;
    
    /* Load barcodes */
    if (load_names(barcodes_path, &matrix->col_names, &matrix->col_names_storage, &matrix->num_cols) != 0) {
        fprintf(stderr, "Failed to load barcodes from %s\n", barcodes_path);
        cf_free_matrix(matrix);
        return NULL;
    }
    
    /* Load features */
    if (load_names(features_path, &matrix->row_names, &matrix->row_names_storage, &matrix->num_rows) != 0) {
        fprintf(stderr, "Failed to load features from %s\n", features_path);
        cf_free_matrix(matrix);
        return NULL;
    }
    
    /* Open matrix file */
    cf_file fh;
    if (open_text_file(matrix_path, &fh) != 0) {
        fprintf(stderr, "Failed to open matrix file %s (or .gz)\n", matrix_path);
        cf_free_matrix(matrix);
        return NULL;
    }
    
    /* Skip header comments */
    char line[MAX_LINE_LENGTH];
    while (cf_gets(&fh, line, sizeof(line))) {
        if (line[0] != '%') break;
    }
    
    /* Parse dimensions line: rows cols nnz */
    int rows, cols, nnz;
    if (sscanf(line, "%d %d %d", &rows, &cols, &nnz) != 3) {
        fprintf(stderr, "Invalid matrix header: %s\n", line);
        close_text_file(&fh);
        cf_free_matrix(matrix);
        return NULL;
    }
    
    /* Verify dimensions match */
    if (rows != matrix->num_rows) {
        fprintf(stderr, "Warning: matrix rows (%d) != features count (%d)\n", rows, matrix->num_rows);
    }
    if (cols != matrix->num_cols) {
        fprintf(stderr, "Warning: matrix cols (%d) != barcodes count (%d)\n", cols, matrix->num_cols);
    }
    
    /* Allocate entries */
    matrix->entries = malloc(nnz * sizeof(cf_matrix_entry));
    if (!matrix->entries) {
        close_text_file(&fh);
        cf_free_matrix(matrix);
        return NULL;
    }
    matrix->num_entries = 0;
    
    /* Read entries (1-based indices in file, convert to 0-based) */
    int row, col, val;
    while (cf_gets(&fh, line, sizeof(line))) {
        if (sscanf(line, "%d %d %d", &row, &col, &val) == 3) {
            if (matrix->num_entries < nnz) {
                matrix->entries[matrix->num_entries].row = row - 1;  /* Convert to 0-based */
                matrix->entries[matrix->num_entries].col = col - 1;
                matrix->entries[matrix->num_entries].value = val;
                matrix->num_entries++;
            }
        }
    }
    
    close_text_file(&fh);
    return matrix;
}

void cf_free_matrix(cf_sparse_matrix *matrix) {
    if (!matrix) return;
    free(matrix->entries);
    free(matrix->row_names);
    free(matrix->col_names);
    free(matrix->row_names_storage);
    free(matrix->col_names_storage);
    free(matrix);
}

/* ============================================================================
 * Feature Calling Implementation
 * ============================================================================ */

/* Helper: comparison for sorting matrix entries by column (barcode) */
static int compare_by_col(const void *a, const void *b) {
    const cf_matrix_entry *ea = (const cf_matrix_entry *)a;
    const cf_matrix_entry *eb = (const cf_matrix_entry *)b;
    if (ea->col != eb->col) return (int)ea->col - (int)eb->col;
    return (int)ea->row - (int)eb->row;
}

cf_call_results* cf_call_features(const cf_sparse_matrix *matrix, const cf_config *config) {
    if (!matrix) return NULL;
    
    /* Use default config if none provided */
    cf_config default_config;
    if (!config) {
        default_config.min_deduped_counts = 2;
        default_config.dominance_fraction = 0.8;
        default_config.dominance_margin = 1;
        default_config.include_ambiguous = 1;
        config = &default_config;
    }
    
    /* Allocate results */
    cf_call_results *results = calloc(1, sizeof(cf_call_results));
    if (!results) return NULL;
    
    results->feature_names = matrix->row_names;
    results->num_features = matrix->num_rows;
    
    /* Sort entries by column (barcode) for efficient processing */
    cf_matrix_entry *sorted_entries = malloc(matrix->num_entries * sizeof(cf_matrix_entry));
    if (!sorted_entries) {
        free(results);
        return NULL;
    }
    memcpy(sorted_entries, matrix->entries, matrix->num_entries * sizeof(cf_matrix_entry));
    qsort(sorted_entries, matrix->num_entries, sizeof(cf_matrix_entry), compare_by_col);
    
    /* Allocate per-barcode feature counts */
    int *feature_counts = calloc(matrix->num_rows, sizeof(int));
    if (!feature_counts) {
        free(sorted_entries);
        free(results);
        return NULL;
    }
    
    /* Allocate calls array (one per barcode) */
    results->calls = calloc(matrix->num_cols, sizeof(cf_cell_call));
    if (!results->calls) {
        free(feature_counts);
        free(sorted_entries);
        free(results);
        return NULL;
    }
    
    /* Process each barcode */
    int current_col = -1;
    int entry_start = 0;
    
    for (int i = 0; i <= matrix->num_entries; i++) {
        int col = (i < matrix->num_entries) ? (int)sorted_entries[i].col : -1;
        
        /* When we move to a new barcode (or finish), process the previous one */
        if (col != current_col && current_col >= 0) {
            /* Sum up counts for this barcode */
            memset(feature_counts, 0, matrix->num_rows * sizeof(int));
            for (int j = entry_start; j < i; j++) {
                int row = sorted_entries[j].row;
                if (row >= 0 && row < matrix->num_rows) {
                    feature_counts[row] += sorted_entries[j].value;
                }
            }
            
            /* Find top features */
            int max_idx = -1, max_count = 0;
            int second_idx = -1, second_count = 0;
            int total_count = 0;
            int num_passing = 0;
            
            for (int f = 0; f < matrix->num_rows; f++) {
                if (feature_counts[f] > 0) {
                    total_count += feature_counts[f];
                    
                    if (feature_counts[f] >= config->min_deduped_counts) {
                        num_passing++;
                        
                        if (feature_counts[f] > max_count) {
                            second_idx = max_idx;
                            second_count = max_count;
                            max_idx = f;
                            max_count = feature_counts[f];
                        } else if (feature_counts[f] > second_count) {
                            second_idx = f;
                            second_count = feature_counts[f];
                        }
                    }
                }
            }
            
            /* Make call decision */
            cf_cell_call *call = &results->calls[results->num_calls];
            call->barcode = matrix->col_names[current_col];
            call->total_umi_count = total_count;
            call->num_features = num_passing;
            call->second_feature_index = second_idx;
            call->second_umi_count = second_count;
            
            if (num_passing == 0) {
                /* No features pass threshold */
                call->call_type = CF_CALL_NONE;
                call->feature_index = -1;
                call->feature_name = NULL;
                call->umi_count = 0;
                results->num_none++;
            } else if (num_passing == 1) {
                /* Single feature - assign it */
                call->call_type = CF_CALL_ASSIGNED;
                call->feature_index = max_idx;
                call->feature_name = matrix->row_names[max_idx];
                call->umi_count = max_count;
                results->num_assigned++;
            } else {
                /* Multiple features - check dominance criteria */
                double fraction = (total_count > 0) ? (double)max_count / total_count : 0.0;
                int margin = max_count - second_count;
                
                if (fraction >= config->dominance_fraction && margin >= config->dominance_margin) {
                    call->call_type = CF_CALL_ASSIGNED;
                    call->feature_index = max_idx;
                    call->feature_name = matrix->row_names[max_idx];
                    call->umi_count = max_count;
                    results->num_assigned++;
                } else {
                    call->call_type = CF_CALL_AMBIGUOUS;
                    call->feature_index = -1;
                    call->feature_name = NULL;
                    call->umi_count = max_count;
                    results->num_ambiguous++;
                }
            }
            
            results->num_calls++;
        }
        
        if (i < matrix->num_entries) {
            if (col != current_col) {
                current_col = col;
                entry_start = i;
            }
        }
    }
    
    free(feature_counts);
    free(sorted_entries);
    
    return results;
}

void cf_free_results(cf_call_results *results) {
    if (!results) return;
    free(results->calls);
    free(results);
}

/* ============================================================================
 * Output Writers
 * ============================================================================ */

int cf_write_calls_csv(const cf_call_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) {
        fprintf(stderr, "Failed to open output file: %s\n", output_path);
        return -1;
    }
    
    /* Write header */
    fprintf(fp, "barcode,feature_call,num_features,num_umis\n");
    
    /* Write calls */
    for (int i = 0; i < results->num_calls; i++) {
        const cf_cell_call *call = &results->calls[i];
        
        const char *call_str;
        switch (call->call_type) {
            case CF_CALL_ASSIGNED:
                call_str = call->feature_name ? call->feature_name : "Unknown";
                break;
            case CF_CALL_AMBIGUOUS:
                call_str = "Multiplet";
                break;
            case CF_CALL_NONE:
            default:
                call_str = "Unassigned";
                break;
        }
        
        fprintf(fp, "%s,%s,%d,%d\n",
                call->barcode,
                call_str,
                call->num_features,
                call->umi_count);
    }
    
    fclose(fp);
    return 0;
}

int cf_write_summary(const cf_call_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) {
        fprintf(stderr, "Failed to open summary file: %s\n", output_path);
        return -1;
    }
    
    fprintf(fp, "Feature Calling Summary\n");
    fprintf(fp, "=======================\n\n");
    fprintf(fp, "Total cells:      %d\n", results->num_calls);
    fprintf(fp, "Assigned:         %d (%.1f%%)\n", 
            results->num_assigned,
            results->num_calls > 0 ? 100.0 * results->num_assigned / results->num_calls : 0.0);
    fprintf(fp, "Ambiguous:        %d (%.1f%%)\n",
            results->num_ambiguous,
            results->num_calls > 0 ? 100.0 * results->num_ambiguous / results->num_calls : 0.0);
    fprintf(fp, "Unassigned:       %d (%.1f%%)\n",
            results->num_none,
            results->num_calls > 0 ? 100.0 * results->num_none / results->num_calls : 0.0);
    fprintf(fp, "\nTotal features:   %d\n", results->num_features);
    
    /* Count cells per feature */
    if (results->num_features > 0 && results->calls) {
        int *feature_counts = calloc(results->num_features, sizeof(int));
        if (feature_counts) {
            for (int i = 0; i < results->num_calls; i++) {
                if (results->calls[i].call_type == CF_CALL_ASSIGNED &&
                    results->calls[i].feature_index >= 0 &&
                    results->calls[i].feature_index < results->num_features) {
                    feature_counts[results->calls[i].feature_index]++;
                }
            }
            
            fprintf(fp, "\nCells per feature:\n");
            for (int f = 0; f < results->num_features; f++) {
                if (results->feature_names && results->feature_names[f]) {
                    fprintf(fp, "  %-30s %d\n", results->feature_names[f], feature_counts[f]);
                } else {
                    fprintf(fp, "  Feature_%d                       %d\n", f, feature_counts[f]);
                }
            }
            free(feature_counts);
        }
    }
    
    fclose(fp);
    return 0;
}

/* ============================================================================
 * Convenience Function
 * ============================================================================ */

int cf_process_mex_dir(const char *mex_dir, const char *output_dir, const cf_config *config) {
    if (!mex_dir || !output_dir) return -1;
    
    /* Create output directory if needed */
    struct stat st = {0};
    if (stat(output_dir, &st) == -1) {
        if (mkdir(output_dir, 0755) != 0) {
            fprintf(stderr, "Failed to create output directory: %s\n", output_dir);
            return -1;
        }
    }
    
    /* Load matrix */
    printf("Loading MEX from: %s\n", mex_dir);
    cf_sparse_matrix *matrix = cf_load_mex(mex_dir);
    if (!matrix) {
        fprintf(stderr, "Failed to load MEX directory\n");
        return -1;
    }
    printf("  Features: %d, Barcodes: %d, Non-zero entries: %d\n",
           matrix->num_rows, matrix->num_cols, matrix->num_entries);
    
    /* Call features */
    printf("Calling features...\n");
    cf_call_results *results = cf_call_features(matrix, config);
    if (!results) {
        fprintf(stderr, "Failed to call features\n");
        cf_free_matrix(matrix);
        return -1;
    }
    
    /* Write outputs */
    char calls_path[MAX_LINE_LENGTH];
    char summary_path[MAX_LINE_LENGTH];
    snprintf(calls_path, sizeof(calls_path), "%s/feature_calls.csv", output_dir);
    snprintf(summary_path, sizeof(summary_path), "%s/feature_calls_summary.txt", output_dir);
    
    printf("Writing calls to: %s\n", calls_path);
    if (cf_write_calls_csv(results, calls_path) != 0) {
        cf_free_results(results);
        cf_free_matrix(matrix);
        return -1;
    }
    
    printf("Writing summary to: %s\n", summary_path);
    if (cf_write_summary(results, summary_path) != 0) {
        cf_free_results(results);
        cf_free_matrix(matrix);
        return -1;
    }
    
    printf("Done. Assigned: %d, Ambiguous: %d, Unassigned: %d\n",
           results->num_assigned, results->num_ambiguous, results->num_none);
    
    cf_free_results(results);
    cf_free_matrix(matrix);
    return 0;
}

/* ============================================================================
 * CR9-Compatible GMM Calling Implementation
 * ============================================================================ */

#include "../include/gmm.h"
#include <limits.h>

cf_gmm_config* cf_gmm_config_create(void) {
    cf_gmm_config *config = calloc(1, sizeof(cf_gmm_config));
    if (!config) return NULL;
    
    config->min_umi_threshold = 3;
    config->n_init = 10;
    
    return config;
}

void cf_gmm_config_destroy(cf_gmm_config *config) {
    if (config) free(config);
}

/* Build per-feature count vectors from sparse matrix */
static int** build_per_feature_counts(const cf_sparse_matrix *matrix, int *num_features, int *num_cells) {
    *num_features = matrix->num_rows;
    *num_cells = matrix->num_cols;
    
    /* Allocate count arrays (one per feature) */
    int **counts = calloc(*num_features, sizeof(int*));
    if (!counts) return NULL;
    
    for (int f = 0; f < *num_features; f++) {
        counts[f] = calloc(*num_cells, sizeof(int));
        if (!counts[f]) {
            for (int i = 0; i < f; i++) free(counts[i]);
            free(counts);
            return NULL;
        }
    }
    
    /* Fill in counts from sparse matrix */
    for (int i = 0; i < matrix->num_entries; i++) {
        int row = matrix->entries[i].row;
        int col = matrix->entries[i].col;
        int val = matrix->entries[i].value;
        if (row >= 0 && row < *num_features && col >= 0 && col < *num_cells) {
            counts[row][col] = val;
        }
    }
    
    return counts;
}

static void free_per_feature_counts(int **counts, int num_features) {
    if (counts) {
        for (int f = 0; f < num_features; f++) {
            free(counts[f]);
        }
        free(counts);
    }
}

cf_gmm_results* cf_call_features_gmm(const cf_sparse_matrix *matrix, const cf_gmm_config *config) {
    if (!matrix) return NULL;
    
    /* Use default config if none provided */
    cf_gmm_config default_config;
    if (!config) {
        default_config.min_umi_threshold = 3;
        default_config.n_init = 10;
        config = &default_config;
    }
    
    int num_features, num_cells;
    int **feature_counts = build_per_feature_counts(matrix, &num_features, &num_cells);
    if (!feature_counts) return NULL;
    
    /* Allocate results */
    cf_gmm_results *results = calloc(1, sizeof(cf_gmm_results));
    if (!results) {
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    results->num_cells = num_cells;
    results->num_features_total = num_features;
    results->feature_names = matrix->row_names;
    results->cell_barcodes = matrix->col_names;
    
    /* Allocate per-cell arrays */
    results->feature_calls = calloc(num_cells, sizeof(char*));
    results->num_features = calloc(num_cells, sizeof(int));
    results->num_umis = calloc(num_cells, sizeof(int));
    
    if (!results->feature_calls || !results->num_features || !results->num_umis) {
        cf_free_gmm_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    /* Allocate per-feature results */
    results->feature_results = calloc(num_features, sizeof(cf_feature_gmm_result));
    if (!results->feature_results) {
        cf_free_gmm_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    /* Allocate per-cell positive feature tracking */
    int **cell_positive_features = calloc(num_cells, sizeof(int*));
    int *cell_num_positive = calloc(num_cells, sizeof(int));
    int *cell_positive_cap = calloc(num_cells, sizeof(int));
    
    if (!cell_positive_features || !cell_num_positive || !cell_positive_cap) {
        free(cell_positive_features);
        free(cell_num_positive);
        free(cell_positive_cap);
        cf_free_gmm_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    /* Run GMM for each feature */
    int *positive_calls = malloc(num_cells * sizeof(int));
    if (!positive_calls) {
        free(cell_positive_features);
        free(cell_num_positive);
        free(cell_positive_cap);
        cf_free_gmm_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    for (int f = 0; f < num_features; f++) {
        int umi_threshold = 0;
        
        /* Run GMM calling for this feature */
        gmm_call_feature(feature_counts[f], num_cells, config->min_umi_threshold,
                         positive_calls, &umi_threshold);
        
        /* Store feature result */
        results->feature_results[f].feature_index = f;
        results->feature_results[f].feature_name = matrix->row_names ? matrix->row_names[f] : NULL;
        results->feature_results[f].umi_threshold = umi_threshold;
        
        /* Calculate total UMIs for this feature */
        int total_umis = 0;
        int num_positive = 0;
        for (int c = 0; c < num_cells; c++) {
            total_umis += feature_counts[f][c];
            if (positive_calls[c]) {
                num_positive++;
                
                /* Track positive feature for this cell */
                if (cell_num_positive[c] >= cell_positive_cap[c]) {
                    int new_cap = cell_positive_cap[c] ? cell_positive_cap[c] * 2 : 4;
                    int *new_arr = realloc(cell_positive_features[c], new_cap * sizeof(int));
                    if (new_arr) {
                        cell_positive_features[c] = new_arr;
                        cell_positive_cap[c] = new_cap;
                    }
                }
                if (cell_num_positive[c] < cell_positive_cap[c]) {
                    cell_positive_features[c][cell_num_positive[c]++] = f;
                }
            }
        }
        results->feature_results[f].total_umis = total_umis;
        results->feature_results[f].num_positive = num_positive;
    }
    
    free(positive_calls);
    
    /* Build per-cell feature call strings */
    for (int c = 0; c < num_cells; c++) {
        int n_pos = cell_num_positive[c];
        results->num_features[c] = n_pos;
        
        /* Calculate total UMIs for called features */
        int total_umi = 0;
        for (int i = 0; i < n_pos; i++) {
            total_umi += feature_counts[cell_positive_features[c][i]][c];
        }
        results->num_umis[c] = total_umi;
        
        /* Check if cell has any UMIs at all */
        int has_any_umi = 0;
        for (int f = 0; f < num_features && !has_any_umi; f++) {
            if (feature_counts[f][c] > 0) has_any_umi = 1;
        }
        
        if (n_pos == 0) {
            if (!has_any_umi) {
                results->cells_no_molecules++;
            } else {
                results->cells_no_call++;
            }
            results->feature_calls[c] = strdup("None");
        } else if (n_pos == 1) {
            results->cells_1_feature++;
            int f = cell_positive_features[c][0];
            results->feature_calls[c] = strdup(matrix->row_names ? matrix->row_names[f] : "Unknown");
        } else {
            results->cells_multi_feature++;
            /* Build pipe-joined string of feature names */
            size_t total_len = 0;
            for (int i = 0; i < n_pos; i++) {
                int f = cell_positive_features[c][i];
                if (matrix->row_names && matrix->row_names[f]) {
                    total_len += strlen(matrix->row_names[f]) + 1;
                }
            }
            char *call_str = malloc(total_len + 1);
            if (call_str) {
                call_str[0] = '\0';
                for (int i = 0; i < n_pos; i++) {
                    int f = cell_positive_features[c][i];
                    if (i > 0) strcat(call_str, "|");
                    if (matrix->row_names && matrix->row_names[f]) {
                        strcat(call_str, matrix->row_names[f]);
                    }
                }
                results->feature_calls[c] = call_str;
            } else {
                results->feature_calls[c] = strdup("Multiple");
            }
        }
    }
    
    /* Cleanup */
    for (int c = 0; c < num_cells; c++) {
        free(cell_positive_features[c]);
    }
    free(cell_positive_features);
    free(cell_num_positive);
    free(cell_positive_cap);
    free_per_feature_counts(feature_counts, num_features);
    
    return results;
}

void cf_free_gmm_results(cf_gmm_results *results) {
    if (!results) return;
    
    if (results->feature_calls) {
        for (int i = 0; i < results->num_cells; i++) {
            free(results->feature_calls[i]);
        }
        free(results->feature_calls);
    }
    free(results->num_features);
    free(results->num_umis);
    free(results->feature_results);
    free(results);
}

int cf_write_protospacer_calls_per_cell(const cf_gmm_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;
    
    /* CR9 format: cell_barcode,num_features,feature_call,num_umis */
    fprintf(fp, "cell_barcode,num_features,feature_call,num_umis\n");
    
    for (int i = 0; i < results->num_cells; i++) {
        fprintf(fp, "%s,%d,%s,%d\n",
                results->cell_barcodes[i],
                results->num_features[i],
                results->feature_calls[i],
                results->num_umis[i]);
    }
    
    fclose(fp);
    return 0;
}

int cf_write_protospacer_calls_summary(const cf_gmm_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;
    
    /* CR9 format header */
    fprintf(fp, "Category,Metric,Value\n");
    
    /* Overall summary */
    fprintf(fp, "All,Cells with 0 molecules,%" PRId64 "\n", (int64_t)results->cells_no_molecules);
    fprintf(fp, "All,Cells with no confident call,%" PRId64 "\n", (int64_t)results->cells_no_call);
    fprintf(fp, "All,Cells with 1 feature,%" PRId64 "\n", (int64_t)results->cells_1_feature);
    fprintf(fp, "All,Cells with >1 features,%" PRId64 "\n", (int64_t)results->cells_multi_feature);
    
    /* Per-feature breakdown */
    for (int f = 0; f < results->num_features_total; f++) {
        const char *name = results->feature_names ? results->feature_names[f] : "Unknown";
        fprintf(fp, "%s,Cells,%" PRId64 "\n", name, (int64_t)results->feature_results[f].num_positive);
    }
    
    fclose(fp);
    return 0;
}

int cf_write_protospacer_umi_thresholds(const cf_gmm_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;
    
    /* CR9 format: only output features with positive cells (threshold > 0) */
    fprintf(fp, "Protospacer,UMI threshold\n");
    
    for (int f = 0; f < results->num_features_total; f++) {
        /* Skip features with no positive cells */
        if (results->feature_results[f].umi_threshold <= 0) continue;
        
        const char *name = results->feature_names ? results->feature_names[f] : "Unknown";
        fprintf(fp, "%s,%d\n",
                name,
                results->feature_results[f].umi_threshold);
    }
    
    fclose(fp);
    return 0;
}

int cf_write_protospacer_umi_thresholds_json(const cf_gmm_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;
    
    /* CR9 format: simple key-value pairs, only features with positive cells */
    fprintf(fp, "{\n");
    int first = 1;
    for (int f = 0; f < results->num_features_total; f++) {
        /* Skip features with no positive cells */
        if (results->feature_results[f].umi_threshold <= 0) continue;
        
        const char *name = results->feature_names ? results->feature_names[f] : "Unknown";
        if (!first) fprintf(fp, ",\n");
        fprintf(fp, "    \"%s\": %d", name, results->feature_results[f].umi_threshold);
        first = 0;
    }
    if (!first) fprintf(fp, "\n");
    fprintf(fp, "}\n");
    
    fclose(fp);
    return 0;
}

int cf_process_mex_dir_gmm(const char *mex_dir, const char *output_dir, const cf_gmm_config *config) {
    if (!mex_dir || !output_dir) return -1;
    
    /* Create output directory if needed */
    struct stat st = {0};
    if (stat(output_dir, &st) == -1) {
        if (mkdir(output_dir, 0755) != 0) {
            fprintf(stderr, "Failed to create output directory: %s\n", output_dir);
            return -1;
        }
    }
    
    /* Load matrix */
    printf("Loading MEX from: %s\n", mex_dir);
    cf_sparse_matrix *matrix = cf_load_mex(mex_dir);
    if (!matrix) {
        fprintf(stderr, "Failed to load MEX directory\n");
        return -1;
    }
    printf("  Features: %d, Barcodes: %d, Non-zero entries: %d\n",
           matrix->num_rows, matrix->num_cols, matrix->num_entries);
    
    /* Call features with GMM */
    printf("Calling features with CR9-style GMM...\n");
    cf_gmm_results *results = cf_call_features_gmm(matrix, config);
    if (!results) {
        fprintf(stderr, "Failed to call features\n");
        cf_free_matrix(matrix);
        return -1;
    }
    
    /* Write CR9-compatible outputs */
    char path[MAX_LINE_LENGTH];
    
    snprintf(path, sizeof(path), "%s/protospacer_calls_per_cell.csv", output_dir);
    printf("Writing: %s\n", path);
    cf_write_protospacer_calls_per_cell(results, path);
    
    snprintf(path, sizeof(path), "%s/protospacer_calls_summary.csv", output_dir);
    printf("Writing: %s\n", path);
    cf_write_protospacer_calls_summary(results, path);
    
    snprintf(path, sizeof(path), "%s/protospacer_umi_thresholds.csv", output_dir);
    printf("Writing: %s\n", path);
    cf_write_protospacer_umi_thresholds(results, path);
    
    snprintf(path, sizeof(path), "%s/protospacer_umi_thresholds.json", output_dir);
    printf("Writing: %s\n", path);
    cf_write_protospacer_umi_thresholds_json(results, path);
    
    printf("\n=== GMM Calling Summary ===\n");
    printf("Total cells:              %d\n", results->num_cells);
    printf("Cells with 0 molecules:   %d\n", results->cells_no_molecules);
    printf("Cells with no call:       %d\n", results->cells_no_call);
    printf("Cells with 1 feature:     %d\n", results->cells_1_feature);
    printf("Cells with >1 features:   %d\n", results->cells_multi_feature);
    
    cf_free_gmm_results(results);
    cf_free_matrix(matrix);
    return 0;
}
