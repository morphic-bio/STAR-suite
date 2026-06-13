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
#include <math.h>
#include <ctype.h>
#include <htslib/kfunc.h>

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

static void cf_gmm_rstrip_line(char *line) {
    if (!line) return;
    size_t len = strlen(line);
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
        line[--len] = '\0';
    }
}

static int cf_copy_name_subset(char ***names_out, char **storage_out,
                               char **source_names, const int *selected_rows,
                               int n_selected) {
    size_t total_len = 0;
    for (int i = 0; i < n_selected; i++) {
        const char *name = source_names[selected_rows ? selected_rows[i] : i];
        total_len += strlen(name ? name : "") + 1;
    }

    char **names = calloc(n_selected > 0 ? n_selected : 1, sizeof(char *));
    char *storage = calloc(total_len > 0 ? total_len : 1, 1);
    if (!names || !storage) {
        free(names);
        free(storage);
        return -1;
    }

    char *ptr = storage;
    for (int i = 0; i < n_selected; i++) {
        const char *name = source_names[selected_rows ? selected_rows[i] : i];
        size_t len = strlen(name ? name : "");
        memcpy(ptr, name ? name : "", len + 1);
        names[i] = ptr;
        ptr += len + 1;
    }

    *names_out = names;
    *storage_out = storage;
    return 0;
}

static cf_sparse_matrix *cf_filter_matrix_rows(const cf_sparse_matrix *matrix,
                                               const int *row_map,
                                               const int *selected_rows,
                                               int n_selected) {
    if (!matrix || !row_map || !selected_rows || n_selected <= 0) return NULL;

    cf_sparse_matrix *filtered = calloc(1, sizeof(cf_sparse_matrix));
    if (!filtered) return NULL;

    filtered->num_rows = n_selected;
    filtered->num_cols = matrix->num_cols;

    if (cf_copy_name_subset(&filtered->row_names, &filtered->row_names_storage,
                            matrix->row_names, selected_rows, n_selected) != 0 ||
        cf_copy_name_subset(&filtered->col_names, &filtered->col_names_storage,
                            matrix->col_names, NULL, matrix->num_cols) != 0) {
        cf_free_matrix(filtered);
        return NULL;
    }

    int kept_entries = 0;
    for (int i = 0; i < matrix->num_entries; i++) {
        uint32_t row = matrix->entries[i].row;
        if (row < (uint32_t)matrix->num_rows && row_map[row] >= 0) {
            kept_entries++;
        }
    }

    filtered->entries = calloc(kept_entries > 0 ? kept_entries : 1, sizeof(cf_matrix_entry));
    if (!filtered->entries) {
        cf_free_matrix(filtered);
        return NULL;
    }

    for (int i = 0; i < matrix->num_entries; i++) {
        uint32_t row = matrix->entries[i].row;
        if (row >= (uint32_t)matrix->num_rows || row_map[row] < 0) continue;
        filtered->entries[filtered->num_entries].row = (uint32_t)row_map[row];
        filtered->entries[filtered->num_entries].col = matrix->entries[i].col;
        filtered->entries[filtered->num_entries].value = matrix->entries[i].value;
        filtered->num_entries++;
    }

    return filtered;
}

static cf_sparse_matrix *cf_filter_matrix_to_crispr_guides(const char *mex_dir,
                                                          cf_sparse_matrix *matrix,
                                                          int *filtered_out) {
    if (filtered_out) *filtered_out = 0;
    if (!mex_dir || !matrix) return NULL;

    char features_path[MAX_LINE_LENGTH];
    snprintf(features_path, sizeof(features_path), "%s/features.tsv", mex_dir);
    if (!file_exists_with_gz(features_path)) {
        snprintf(features_path, sizeof(features_path), "%s/features.txt", mex_dir);
    }
    if (!file_exists_with_gz(features_path)) {
        return matrix;
    }

    cf_file fh;
    if (open_text_file(features_path, &fh) != 0) {
        return matrix;
    }

    int *row_map = malloc(matrix->num_rows * sizeof(int));
    int *selected_rows = malloc(matrix->num_rows * sizeof(int));
    if (!row_map || !selected_rows) {
        free(row_map);
        free(selected_rows);
        close_text_file(&fh);
        return NULL;
    }
    for (int i = 0; i < matrix->num_rows; i++) row_map[i] = -1;

    char line[MAX_LINE_LENGTH];
    int row = 0;
    int n_guides = 0;
    while (row < matrix->num_rows && cf_gets(&fh, line, sizeof(line))) {
        cf_gmm_rstrip_line(line);
        char *name = NULL;
        char *type = NULL;
        char *tab1 = strchr(line, '\t');
        if (tab1) {
            *tab1 = '\0';
            name = tab1 + 1;
            char *tab2 = strchr(name, '\t');
            if (tab2) {
                *tab2 = '\0';
                type = tab2 + 1;
            }
        }
        if (type && strcmp(type, "CRISPR Guide Capture") == 0) {
            row_map[row] = n_guides;
            selected_rows[n_guides] = row;
            n_guides++;
        }
        row++;
    }
    close_text_file(&fh);

    if (n_guides == 0 || n_guides == matrix->num_rows) {
        free(row_map);
        free(selected_rows);
        return matrix;
    }

    cf_sparse_matrix *filtered = cf_filter_matrix_rows(matrix, row_map, selected_rows, n_guides);
    free(row_map);
    free(selected_rows);
    if (!filtered) return NULL;
    if (filtered_out) *filtered_out = 1;
    return filtered;
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

    int filtered_to_guides = 0;
    cf_sparse_matrix *guide_matrix = cf_filter_matrix_to_crispr_guides(mex_dir, matrix, &filtered_to_guides);
    if (!guide_matrix) {
        fprintf(stderr, "Failed to filter MEX to CRISPR Guide Capture rows\n");
        cf_free_matrix(matrix);
        return -1;
    }
    if (guide_matrix != matrix) {
        cf_free_matrix(matrix);
        matrix = guide_matrix;
    }
    if (filtered_to_guides) {
        printf("  CRISPR Guide Capture rows: %d, guide non-zero entries: %d\n",
               matrix->num_rows, matrix->num_entries);
    }
    
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

/* ============================================================================
 * Ambient-FDR Guide Calling Implementation
 * ============================================================================ */

typedef struct cf_feature_meta {
    char *id;
    char *name;
    char *type;
} cf_feature_meta;

typedef struct cf_feature_table {
    cf_feature_meta *features;
    int count;
} cf_feature_table;

typedef struct cf_id_map_entry {
    const char *id;
    int guide_index;
} cf_id_map_entry;

typedef struct cf_ambient_entry {
    int cell;
    int guide;
    uint64_t count;
    double pvalue;
    double qvalue;
} cf_ambient_entry;

typedef struct cf_pvalue_rank {
    int entry_index;
    double pvalue;
} cf_pvalue_rank;

static int cf_mkdir_p(const char *path) {
    if (!path || !*path) return -1;

    char tmp[MAX_LINE_LENGTH];
    snprintf(tmp, sizeof(tmp), "%s", path);
    size_t len = strlen(tmp);
    if (len == 0) return -1;
    if (tmp[len - 1] == '/') tmp[len - 1] = '\0';

    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
            *p = '/';
        }
    }
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
    return 0;
}

static char *cf_strdup_safe(const char *s) {
    return strdup(s ? s : "");
}

static void cf_rstrip_line(char *line) {
    if (!line) return;
    size_t len = strlen(line);
    while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
        line[--len] = '\0';
    }
}

static void cf_json_write_string(FILE *fp, const char *s) {
    fputc('"', fp);
    if (s) {
        for (const unsigned char *p = (const unsigned char *)s; *p; ++p) {
            switch (*p) {
                case '\\': fputs("\\\\", fp); break;
                case '"': fputs("\\\"", fp); break;
                case '\b': fputs("\\b", fp); break;
                case '\f': fputs("\\f", fp); break;
                case '\n': fputs("\\n", fp); break;
                case '\r': fputs("\\r", fp); break;
                case '\t': fputs("\\t", fp); break;
                default:
                    if (*p < 0x20) {
                        fprintf(fp, "\\u%04x", *p);
                    } else {
                        fputc(*p, fp);
                    }
                    break;
            }
        }
    }
    fputc('"', fp);
}

static int cf_resolve_features_path(const char *mex_dir, char *path, size_t path_size) {
    snprintf(path, path_size, "%s/features.txt", mex_dir);
    if (file_exists_with_gz(path)) return 0;
    snprintf(path, path_size, "%s/features.tsv", mex_dir);
    if (file_exists_with_gz(path)) return 0;
    return -1;
}

static int cf_load_feature_table(const char *mex_dir, cf_feature_table *table) {
    if (!mex_dir || !table) return -1;
    memset(table, 0, sizeof(*table));

    char path[MAX_LINE_LENGTH];
    if (cf_resolve_features_path(mex_dir, path, sizeof(path)) != 0) {
        fprintf(stderr, "Failed to find features.tsv/features.txt in %s\n", mex_dir);
        return -1;
    }

    cf_file fh;
    if (open_text_file(path, &fh) != 0) {
        fprintf(stderr, "Failed to open features file %s\n", path);
        return -1;
    }

    int capacity = 128;
    cf_feature_meta *features = calloc(capacity, sizeof(cf_feature_meta));
    if (!features) {
        close_text_file(&fh);
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    int count = 0;
    while (cf_gets(&fh, line, sizeof(line))) {
        cf_rstrip_line(line);
        if (line[0] == '\0') continue;

        char *id = line;
        char *name = NULL;
        char *type = NULL;

        char *tab1 = strchr(line, '\t');
        if (tab1) {
            *tab1 = '\0';
            name = tab1 + 1;
            char *tab2 = strchr(name, '\t');
            if (tab2) {
                *tab2 = '\0';
                type = tab2 + 1;
            }
        }
        if (!name || !*name) name = id;
        if (!type) type = "";

        if (count >= capacity) {
            capacity *= 2;
            cf_feature_meta *grown = realloc(features, capacity * sizeof(cf_feature_meta));
            if (!grown) {
                close_text_file(&fh);
                for (int i = 0; i < count; i++) {
                    free(features[i].id);
                    free(features[i].name);
                    free(features[i].type);
                }
                free(features);
                return -1;
            }
            memset(grown + count, 0, (capacity - count) * sizeof(cf_feature_meta));
            features = grown;
        }

        features[count].id = cf_strdup_safe(id);
        features[count].name = cf_strdup_safe(name);
        features[count].type = cf_strdup_safe(type);
        if (!features[count].id || !features[count].name || !features[count].type) {
            close_text_file(&fh);
            for (int i = 0; i <= count; i++) {
                free(features[i].id);
                free(features[i].name);
                free(features[i].type);
            }
            free(features);
            return -1;
        }
        count++;
    }

    close_text_file(&fh);
    table->features = features;
    table->count = count;
    return 0;
}

static void cf_free_feature_table(cf_feature_table *table) {
    if (!table || !table->features) return;
    for (int i = 0; i < table->count; i++) {
        free(table->features[i].id);
        free(table->features[i].name);
        free(table->features[i].type);
    }
    free(table->features);
    table->features = NULL;
    table->count = 0;
}

static int cf_has_crispr_feature_type(const cf_feature_table *table) {
    if (!table) return 0;
    for (int i = 0; i < table->count; i++) {
        if (table->features[i].type && strcmp(table->features[i].type, "CRISPR Guide Capture") == 0) {
            return 1;
        }
    }
    return 0;
}

static int cf_feature_is_guide(const cf_feature_table *table, int row, int require_crispr_type) {
    if (!table || row < 0 || row >= table->count) return 0;
    if (!require_crispr_type) return 1;
    return table->features[row].type && strcmp(table->features[row].type, "CRISPR Guide Capture") == 0;
}

static int cf_compare_id_map(const void *a, const void *b) {
    const cf_id_map_entry *ea = (const cf_id_map_entry *)a;
    const cf_id_map_entry *eb = (const cf_id_map_entry *)b;
    return strcmp(ea->id, eb->id);
}

static int cf_compare_id_key(const void *key, const void *elem) {
    const char *id = *(const char * const *)key;
    const cf_id_map_entry *entry = (const cf_id_map_entry *)elem;
    return strcmp(id, entry->id);
}

static int cf_lookup_guide_id(const cf_id_map_entry *map, int n, const char *id) {
    if (!map || n <= 0 || !id) return -1;
    const char *key = id;
    cf_id_map_entry *found = bsearch(&key, map, n, sizeof(cf_id_map_entry), cf_compare_id_key);
    return found ? found->guide_index : -1;
}

static int cf_compare_string_ptr(const void *a, const void *b) {
    const char *sa = *(const char * const *)a;
    const char *sb = *(const char * const *)b;
    return strcmp(sa, sb);
}

static int cf_string_set_contains(char **sorted, int n, const char *value) {
    if (!sorted || n <= 0 || !value) return 0;
    char *key = (char *)value;
    return bsearch(&key, sorted, n, sizeof(char *), cf_compare_string_ptr) != NULL;
}

static int cf_compare_ambient_entry_cell_guide(const void *a, const void *b) {
    const cf_ambient_entry *ea = (const cf_ambient_entry *)a;
    const cf_ambient_entry *eb = (const cf_ambient_entry *)b;
    if (ea->cell != eb->cell) return ea->cell - eb->cell;
    return ea->guide - eb->guide;
}

static int cf_compare_pvalue_rank(const void *a, const void *b) {
    const cf_pvalue_rank *ra = (const cf_pvalue_rank *)a;
    const cf_pvalue_rank *rb = (const cf_pvalue_rank *)b;
    if (ra->pvalue < rb->pvalue) return -1;
    if (ra->pvalue > rb->pvalue) return 1;
    return ra->entry_index - rb->entry_index;
}

static int cf_append_ambient_entry(cf_ambient_entry **entries, size_t *n, size_t *cap,
                                   int cell, int guide, uint64_t count) {
    if (*n >= *cap) {
        size_t new_cap = *cap ? (*cap * 2) : 1024;
        cf_ambient_entry *grown = realloc(*entries, new_cap * sizeof(cf_ambient_entry));
        if (!grown) return -1;
        *entries = grown;
        *cap = new_cap;
    }
    (*entries)[*n].cell = cell;
    (*entries)[*n].guide = guide;
    (*entries)[*n].count = count;
    (*entries)[*n].pvalue = 1.0;
    (*entries)[*n].qvalue = 1.0;
    (*n)++;
    return 0;
}

static size_t cf_collapse_ambient_entries(cf_ambient_entry *entries, size_t n) {
    if (!entries || n == 0) return 0;
    qsort(entries, n, sizeof(cf_ambient_entry), cf_compare_ambient_entry_cell_guide);

    size_t out = 0;
    for (size_t i = 0; i < n; i++) {
        if (out > 0 &&
            entries[out - 1].cell == entries[i].cell &&
            entries[out - 1].guide == entries[i].guide) {
            entries[out - 1].count += entries[i].count;
        } else {
            entries[out++] = entries[i];
        }
    }
    return out;
}

static double cf_poisson_upper_tail(uint64_t observed, double lambda) {
    if (observed == 0) return 1.0;
    if (lambda <= 0.0) return 0.0;

    double p = kf_gammap((double)observed, lambda);
    if (isnan(p)) return 1.0;
    if (p < 0.0) return 0.0;
    if (p > 1.0) return 1.0;
    return p;
}

static void cf_apply_bh(cf_ambient_entry *entries, size_t n_entries, uint64_t total_tests) {
    if (!entries || n_entries == 0 || total_tests == 0) return;

    cf_pvalue_rank *ranks = malloc(n_entries * sizeof(cf_pvalue_rank));
    if (!ranks) {
        for (size_t i = 0; i < n_entries; i++) entries[i].qvalue = 1.0;
        return;
    }

    for (size_t i = 0; i < n_entries; i++) {
        ranks[i].entry_index = (int)i;
        ranks[i].pvalue = entries[i].pvalue;
    }
    qsort(ranks, n_entries, sizeof(cf_pvalue_rank), cf_compare_pvalue_rank);

    double running_min = 1.0;
    double m = (double)total_tests;
    for (size_t rev = n_entries; rev > 0; rev--) {
        size_t i = rev - 1;
        double rank = (double)(i + 1);
        double q = ranks[i].pvalue * m / rank;
        if (q > 1.0) q = 1.0;
        if (q > running_min) q = running_min;
        running_min = q;
        entries[ranks[i].entry_index].qvalue = q;
    }

    free(ranks);
}

static int cf_append_call_name(char **call, const char *name) {
    const char *safe = (name && *name) ? name : "Unknown";
    size_t old_len = (*call) ? strlen(*call) : 0;
    size_t add_len = strlen(safe);
    size_t sep = old_len > 0 ? 1 : 0;
    char *grown = realloc(*call, old_len + sep + add_len + 1);
    if (!grown) return -1;
    if (sep) grown[old_len++] = '|';
    memcpy(grown + old_len, safe, add_len + 1);
    *call = grown;
    return 0;
}

static int cf_build_cell_calls(const cf_ambient_entry *entries, size_t n_entries,
                               const char **guide_names, int n_cells,
                               double threshold, int min_umi,
                               char ***cell_calls_out, int **cell_num_features_out,
                               uint64_t **cell_num_umis_out, double **cell_min_q_out,
                               uint64_t **cell_min_called_umi_out,
                               uint64_t **cell_max_called_umi_out,
                               int *cells_no_call_out, int *cells_1_feature_out,
                               int *cells_multi_feature_out) {
    char **cell_calls = calloc(n_cells, sizeof(char *));
    int *cell_num_features = calloc(n_cells, sizeof(int));
    uint64_t *cell_num_umis = calloc(n_cells, sizeof(uint64_t));
    double *cell_min_q = malloc(n_cells * sizeof(double));
    uint64_t *cell_min_called_umi = malloc(n_cells * sizeof(uint64_t));
    uint64_t *cell_max_called_umi = calloc(n_cells, sizeof(uint64_t));
    if (!cell_calls || !cell_num_features || !cell_num_umis || !cell_min_q ||
        !cell_min_called_umi || !cell_max_called_umi) {
        free(cell_calls);
        free(cell_num_features);
        free(cell_num_umis);
        free(cell_min_q);
        free(cell_min_called_umi);
        free(cell_max_called_umi);
        return -1;
    }

    for (int c = 0; c < n_cells; c++) {
        cell_min_q[c] = 1.0;
        cell_min_called_umi[c] = UINT64_MAX;
    }

    for (size_t i = 0; i < n_entries; i++) {
        int cell = entries[i].cell;
        int guide = entries[i].guide;
        if (cell < 0 || cell >= n_cells) continue;
        if (entries[i].qvalue < cell_min_q[cell]) cell_min_q[cell] = entries[i].qvalue;

        if (entries[i].count >= (uint64_t)min_umi && entries[i].qvalue <= threshold) {
            if (cf_append_call_name(&cell_calls[cell], guide_names ? guide_names[guide] : "Unknown") != 0) {
                for (int c = 0; c < n_cells; c++) free(cell_calls[c]);
                free(cell_calls);
                free(cell_num_features);
                free(cell_num_umis);
                free(cell_min_q);
                return -1;
            }
            cell_num_features[cell]++;
            cell_num_umis[cell] += entries[i].count;
            if (entries[i].count < cell_min_called_umi[cell]) {
                cell_min_called_umi[cell] = entries[i].count;
            }
            if (entries[i].count > cell_max_called_umi[cell]) {
                cell_max_called_umi[cell] = entries[i].count;
            }
        }
    }

    int no_call = 0, one = 0, multi = 0;
    for (int c = 0; c < n_cells; c++) {
        if (cell_num_features[c] == 0) {
            cell_min_called_umi[c] = 0;
            no_call++;
        } else if (cell_num_features[c] == 1) {
            one++;
        } else {
            multi++;
        }
    }

    *cell_calls_out = cell_calls;
    *cell_num_features_out = cell_num_features;
    *cell_num_umis_out = cell_num_umis;
    *cell_min_q_out = cell_min_q;
    *cell_min_called_umi_out = cell_min_called_umi;
    *cell_max_called_umi_out = cell_max_called_umi;
    if (cells_no_call_out) *cells_no_call_out = no_call;
    if (cells_1_feature_out) *cells_1_feature_out = one;
    if (cells_multi_feature_out) *cells_multi_feature_out = multi;
    return 0;
}

static void cf_free_cell_calls(char **cell_calls, int n_cells,
                               int *cell_num_features, uint64_t *cell_num_umis,
                               double *cell_min_q,
                               uint64_t *cell_min_called_umi,
                               uint64_t *cell_max_called_umi) {
    if (cell_calls) {
        for (int c = 0; c < n_cells; c++) free(cell_calls[c]);
    }
    free(cell_calls);
    free(cell_num_features);
    free(cell_num_umis);
    free(cell_min_q);
    free(cell_min_called_umi);
    free(cell_max_called_umi);
}

static void cf_count_calls_at_threshold(const cf_ambient_entry *entries, size_t n_entries,
                                        int n_cells, double threshold, int min_umi,
                                        int *no_call, int *one, int *multi) {
    int *counts = calloc(n_cells, sizeof(int));
    if (!counts) {
        *no_call = n_cells;
        *one = 0;
        *multi = 0;
        return;
    }

    for (size_t i = 0; i < n_entries; i++) {
        if (entries[i].count >= (uint64_t)min_umi && entries[i].qvalue <= threshold &&
            entries[i].cell >= 0 && entries[i].cell < n_cells) {
            counts[entries[i].cell]++;
        }
    }

    *no_call = *one = *multi = 0;
    for (int c = 0; c < n_cells; c++) {
        if (counts[c] == 0) (*no_call)++;
        else if (counts[c] == 1) (*one)++;
        else (*multi)++;
    }
    free(counts);
}

static int cf_write_ambient_skip_summary(const char *output_dir, const char *status,
                                         const char *message) {
    if (cf_mkdir_p(output_dir) != 0) return -1;

    char path[MAX_LINE_LENGTH];
    snprintf(path, sizeof(path), "%s/guide_fdr_summary.json", output_dir);
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;
    fprintf(fp, "{\n");
    fprintf(fp, "  \"caller\": \"ambient-fdr\",\n");
    fprintf(fp, "  \"status\": ");
    cf_json_write_string(fp, status);
    fprintf(fp, ",\n  \"message\": ");
    cf_json_write_string(fp, message);
    fprintf(fp, "\n}\n");
    fclose(fp);
    return 0;
}

static int cf_write_ambient_rates(const char *output_dir,
                                  const char **guide_ids,
                                  const char **guide_names,
                                  const uint64_t *ambient_counts,
                                  const double *ambient_rates,
                                  int n_guides) {
    char path[MAX_LINE_LENGTH];
    snprintf(path, sizeof(path), "%s/guide_ambient_rates.tsv", output_dir);
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;

    fprintf(fp, "feature_id\tfeature_name\tambient_umis\tambient_rate\n");
    for (int g = 0; g < n_guides; g++) {
        fprintf(fp, "%s\t%s\t%" PRIu64 "\t%.17g\n",
                guide_ids[g] ? guide_ids[g] : "",
                guide_names[g] ? guide_names[g] : "",
                ambient_counts[g],
                ambient_rates[g]);
    }
    fclose(fp);
    return 0;
}

static int cf_write_ambient_qvalues(const char *output_dir,
                                    const cf_ambient_entry *entries,
                                    size_t n_entries,
                                    const cf_sparse_matrix *filtered_matrix,
                                    const char **guide_ids,
                                    const char **guide_names,
                                    int n_guides) {
    char path[MAX_LINE_LENGTH];

    snprintf(path, sizeof(path), "%s/guide_qvalues.mtx", output_dir);
    FILE *mtx = fopen(path, "w");
    if (!mtx) return -1;
    fprintf(mtx, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(mtx, "%% rows are filtered cell barcodes; columns are CRISPR guides; missing entries imply qvalue=1\n");
    fprintf(mtx, "%d %d %zu\n", filtered_matrix->num_cols, n_guides, n_entries);
    for (size_t i = 0; i < n_entries; i++) {
        fprintf(mtx, "%d %d %.17g\n",
                entries[i].cell + 1,
                entries[i].guide + 1,
                entries[i].qvalue);
    }
    fclose(mtx);

    snprintf(path, sizeof(path), "%s/guide_qvalues_barcodes.tsv", output_dir);
    FILE *bc = fopen(path, "w");
    if (!bc) return -1;
    for (int c = 0; c < filtered_matrix->num_cols; c++) {
        fprintf(bc, "%s\n", filtered_matrix->col_names[c]);
    }
    fclose(bc);

    snprintf(path, sizeof(path), "%s/guide_qvalues_features.tsv", output_dir);
    FILE *feat = fopen(path, "w");
    if (!feat) return -1;
    for (int g = 0; g < n_guides; g++) {
        fprintf(feat, "%s\t%s\tCRISPR Guide Capture\n",
                guide_ids[g] ? guide_ids[g] : "",
                guide_names[g] ? guide_names[g] : "");
    }
    fclose(feat);

    return 0;
}

static int cf_write_ambient_calls(const char *output_dir,
                                  const cf_sparse_matrix *filtered_matrix,
                                  char **cell_calls,
                                  const int *cell_num_features,
                                  const uint64_t *cell_num_umis,
                                  const uint64_t *cell_min_called_umi,
                                  const uint64_t *cell_max_called_umi,
                                  const double *cell_min_q,
                                  double fdr_threshold) {
    char path[MAX_LINE_LENGTH];
    snprintf(path, sizeof(path), "%s/guide_fdr_calls_per_cell.csv", output_dir);
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;

    fprintf(fp, "cell_barcode,num_features,feature_call,num_umis,min_called_umi,max_called_umi,min_qvalue,num_features_at_default_fdr,call_status,default_fdr,caller\n");
    for (int c = 0; c < filtered_matrix->num_cols; c++) {
        const char *call = cell_calls[c] ? cell_calls[c] : "None";
        const char *status = "none";
        if (cell_num_features[c] == 1) status = "singlet";
        else if (cell_num_features[c] > 1) status = "multiplet";

        fprintf(fp, "%s,%d,%s,%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%.17g,%d,%s,%.17g,ambient-fdr\n",
                filtered_matrix->col_names[c],
                cell_num_features[c],
                call,
                cell_num_umis[c],
                cell_min_called_umi[c],
                cell_max_called_umi[c],
                cell_min_q[c],
                cell_num_features[c],
                status,
                fdr_threshold);
    }

    fclose(fp);
    return 0;
}

static int cf_write_ambient_summary(const char *output_dir,
                                    const cf_ambient_entry *entries,
                                    size_t n_entries,
                                    int n_cells,
                                    int n_guides,
                                    uint64_t total_tests,
                                    int raw_barcodes,
                                    int ambient_barcodes,
                                    int filtered_cells_absent_from_raw,
                                    uint64_t ambient_total,
                                    int min_umi,
                                    double fdr_threshold,
                                    int emit_sparse_qvalues,
                                    int cells_no_call,
                                    int cells_1_feature,
                                    int cells_multi_feature) {
    char path[MAX_LINE_LENGTH];
    snprintf(path, sizeof(path), "%s/guide_fdr_summary.json", output_dir);
    FILE *fp = fopen(path, "w");
    if (!fp) return -1;

    int assigned = cells_1_feature + cells_multi_feature;
    fprintf(fp, "{\n");
    fprintf(fp, "  \"caller\": \"ambient-fdr\",\n");
    fprintf(fp, "  \"status\": \"ok\",\n");
    fprintf(fp, "  \"default_fdr\": %.17g,\n", fdr_threshold);
    fprintf(fp, "  \"min_umi\": %d,\n", min_umi);
    fprintf(fp, "  \"emit_qvalues\": \"%s\",\n", emit_sparse_qvalues ? "sparse" : "none");
    fprintf(fp, "  \"num_filtered_cells\": %d,\n", n_cells);
    fprintf(fp, "  \"num_guides\": %d,\n", n_guides);
    fprintf(fp, "  \"total_tests\": %" PRIu64 ",\n", total_tests);
    fprintf(fp, "  \"observed_filtered_entries\": %zu,\n", n_entries);
    fprintf(fp, "  \"qvalue_entries\": %zu,\n", emit_sparse_qvalues ? n_entries : 0);
    fprintf(fp, "  \"raw_barcodes\": %d,\n", raw_barcodes);
    fprintf(fp, "  \"ambient_barcodes\": %d,\n", ambient_barcodes);
    fprintf(fp, "  \"filtered_cells_absent_from_raw\": %d,\n", filtered_cells_absent_from_raw);
    fprintf(fp, "  \"ambient_total_umis\": %" PRIu64 ",\n", ambient_total);
    fprintf(fp, "  \"cells_no_call\": %d,\n", cells_no_call);
    fprintf(fp, "  \"cells_1_feature\": %d,\n", cells_1_feature);
    fprintf(fp, "  \"cells_multi_feature\": %d,\n", cells_multi_feature);
    fprintf(fp, "  \"assigned_cells\": %d,\n", assigned);
    fprintf(fp, "  \"assignment_rate\": %.17g,\n", n_cells > 0 ? (double)assigned / (double)n_cells : 0.0);
    fprintf(fp, "  \"fdr_grid\": [\n");

    const double grid[] = {0.001, 0.005, 0.01, 0.05, 0.1};
    const int n_grid = (int)(sizeof(grid) / sizeof(grid[0]));
    for (int i = 0; i < n_grid; i++) {
        int no_call = 0, one = 0, multi = 0;
        cf_count_calls_at_threshold(entries, n_entries, n_cells, grid[i], min_umi,
                                    &no_call, &one, &multi);
        int grid_assigned = one + multi;
        fprintf(fp, "    {\"fdr\": %.17g, \"cells_no_call\": %d, \"cells_1_feature\": %d, \"cells_multi_feature\": %d, \"assigned_cells\": %d, \"assignment_rate\": %.17g}%s\n",
                grid[i], no_call, one, multi, grid_assigned,
                n_cells > 0 ? (double)grid_assigned / (double)n_cells : 0.0,
                (i + 1 < n_grid) ? "," : "");
    }
    fprintf(fp, "  ]\n");
    fprintf(fp, "}\n");
    fclose(fp);

    snprintf(path, sizeof(path), "%s/guide_fdr_qc.json", output_dir);
    FILE *qc = fopen(path, "w");
    if (qc) {
        fprintf(qc, "{\n");
        fprintf(qc, "  \"caller\": \"ambient-fdr\",\n");
        fprintf(qc, "  \"status\": \"ok\",\n");
        fprintf(qc, "  \"default_fdr\": %.17g,\n", fdr_threshold);
        fprintf(qc, "  \"assigned_cells\": %d,\n", assigned);
        fprintf(qc, "  \"assignment_rate\": %.17g,\n", n_cells > 0 ? (double)assigned / (double)n_cells : 0.0);
        fprintf(qc, "  \"cells_no_call\": %d,\n", cells_no_call);
        fprintf(qc, "  \"cells_1_feature\": %d,\n", cells_1_feature);
        fprintf(qc, "  \"cells_multi_feature\": %d,\n", cells_multi_feature);
        fprintf(qc, "  \"ambient_barcodes\": %d,\n", ambient_barcodes);
        fprintf(qc, "  \"ambient_total_umis\": %" PRIu64 "\n", ambient_total);
        fprintf(qc, "}\n");
        fclose(qc);
    }
    return 0;
}

cf_ambient_fdr_config* cf_ambient_fdr_config_create(void) {
    cf_ambient_fdr_config *config = calloc(1, sizeof(cf_ambient_fdr_config));
    if (!config) return NULL;
    config->fdr_threshold = 0.01;
    config->min_umi = 1;
    config->emit_sparse_qvalues = 1;
    return config;
}

void cf_ambient_fdr_config_destroy(cf_ambient_fdr_config *config) {
    if (config) free(config);
}

int cf_process_mex_dir_ambient_fdr(const char *raw_mex_dir,
                                   const char *filtered_mex_dir,
                                   const char *output_dir,
                                   const cf_ambient_fdr_config *config) {
    if (!raw_mex_dir || !filtered_mex_dir || !output_dir) return -1;

    cf_ambient_fdr_config default_config;
    if (!config) {
        default_config.fdr_threshold = 0.01;
        default_config.min_umi = 1;
        default_config.emit_sparse_qvalues = 1;
        config = &default_config;
    }
    double fdr_threshold = config->fdr_threshold;
    if (!(fdr_threshold > 0.0 && fdr_threshold <= 1.0)) {
        fprintf(stderr, "Invalid Ambient-FDR threshold %.17g; expected 0 < FDR <= 1\n", fdr_threshold);
        return -1;
    }
    int min_umi = config->min_umi > 0 ? config->min_umi : 1;

    if (cf_mkdir_p(output_dir) != 0) {
        fprintf(stderr, "Failed to create output directory: %s\n", output_dir);
        return -1;
    }

    printf("Loading raw MEX from: %s\n", raw_mex_dir);
    cf_sparse_matrix *raw_matrix = cf_load_mex(raw_mex_dir);
    if (!raw_matrix) {
        fprintf(stderr, "Failed to load raw MEX directory\n");
        return -1;
    }
    printf("Loading filtered MEX from: %s\n", filtered_mex_dir);
    cf_sparse_matrix *filtered_matrix = cf_load_mex(filtered_mex_dir);
    if (!filtered_matrix) {
        cf_free_matrix(raw_matrix);
        fprintf(stderr, "Failed to load filtered MEX directory\n");
        return -1;
    }

    cf_feature_table raw_features = {0};
    cf_feature_table filtered_features = {0};
    if (cf_load_feature_table(raw_mex_dir, &raw_features) != 0) {
        cf_free_matrix(raw_matrix);
        cf_free_matrix(filtered_matrix);
        return -1;
    }
    if (cf_load_feature_table(filtered_mex_dir, &filtered_features) != 0) {
        cf_free_feature_table(&raw_features);
        cf_free_matrix(raw_matrix);
        cf_free_matrix(filtered_matrix);
        return -1;
    }

    if (raw_features.count != raw_matrix->num_rows) {
        fprintf(stderr, "Warning: raw features count (%d) != matrix rows (%d)\n",
                raw_features.count, raw_matrix->num_rows);
    }
    if (filtered_features.count != filtered_matrix->num_rows) {
        fprintf(stderr, "Warning: filtered features count (%d) != matrix rows (%d)\n",
                filtered_features.count, filtered_matrix->num_rows);
    }

    int filtered_has_crispr_type = cf_has_crispr_feature_type(&filtered_features);
    int raw_has_crispr_type = cf_has_crispr_feature_type(&raw_features);

    int n_guides = 0;
    for (int r = 0; r < filtered_features.count && r < filtered_matrix->num_rows; r++) {
        if (cf_feature_is_guide(&filtered_features, r, filtered_has_crispr_type)) n_guides++;
    }
    if (n_guides == 0 || filtered_matrix->num_cols == 0) {
        cf_free_feature_table(&raw_features);
        cf_free_feature_table(&filtered_features);
        cf_free_matrix(raw_matrix);
        cf_free_matrix(filtered_matrix);
        return cf_write_ambient_skip_summary(output_dir, "skipped_no_guides",
                                             "No CRISPR Guide Capture features or filtered cells were found.");
    }

    const char **guide_ids = calloc(n_guides, sizeof(char *));
    const char **guide_names = calloc(n_guides, sizeof(char *));
    int *filtered_row_to_guide = malloc(filtered_matrix->num_rows * sizeof(int));
    int *raw_row_to_guide = malloc(raw_matrix->num_rows * sizeof(int));
    cf_id_map_entry *id_map = calloc(n_guides, sizeof(cf_id_map_entry));
    uint64_t *ambient_counts = calloc(n_guides, sizeof(uint64_t));
    double *ambient_rates = calloc(n_guides, sizeof(double));
    uint64_t *cell_depth = calloc(filtered_matrix->num_cols, sizeof(uint64_t));

    if (!guide_ids || !guide_names || !filtered_row_to_guide || !raw_row_to_guide ||
        !id_map || !ambient_counts || !ambient_rates || !cell_depth) {
        free(guide_ids); free(guide_names); free(filtered_row_to_guide); free(raw_row_to_guide);
        free(id_map); free(ambient_counts); free(ambient_rates); free(cell_depth);
        cf_free_feature_table(&raw_features);
        cf_free_feature_table(&filtered_features);
        cf_free_matrix(raw_matrix);
        cf_free_matrix(filtered_matrix);
        return -1;
    }

    for (int r = 0; r < filtered_matrix->num_rows; r++) filtered_row_to_guide[r] = -1;
    for (int r = 0; r < raw_matrix->num_rows; r++) raw_row_to_guide[r] = -1;

    int guide_index = 0;
    for (int r = 0; r < filtered_features.count && r < filtered_matrix->num_rows; r++) {
        if (!cf_feature_is_guide(&filtered_features, r, filtered_has_crispr_type)) continue;
        filtered_row_to_guide[r] = guide_index;
        guide_ids[guide_index] = filtered_features.features[r].id;
        guide_names[guide_index] = filtered_features.features[r].name;
        id_map[guide_index].id = filtered_features.features[r].id;
        id_map[guide_index].guide_index = guide_index;
        guide_index++;
    }
    qsort(id_map, n_guides, sizeof(cf_id_map_entry), cf_compare_id_map);

    for (int r = 0; r < raw_features.count && r < raw_matrix->num_rows; r++) {
        if (!cf_feature_is_guide(&raw_features, r, raw_has_crispr_type)) continue;
        raw_row_to_guide[r] = cf_lookup_guide_id(id_map, n_guides, raw_features.features[r].id);
    }

    char **filtered_barcodes_sorted = malloc(filtered_matrix->num_cols * sizeof(char *));
    char **raw_barcodes_sorted = malloc(raw_matrix->num_cols * sizeof(char *));
    unsigned char *raw_is_filtered_cell = calloc(raw_matrix->num_cols, sizeof(unsigned char));
    if (!filtered_barcodes_sorted || !raw_barcodes_sorted || !raw_is_filtered_cell) {
        free(filtered_barcodes_sorted); free(raw_barcodes_sorted); free(raw_is_filtered_cell);
        free(guide_ids); free(guide_names); free(filtered_row_to_guide); free(raw_row_to_guide);
        free(id_map); free(ambient_counts); free(ambient_rates); free(cell_depth);
        cf_free_feature_table(&raw_features);
        cf_free_feature_table(&filtered_features);
        cf_free_matrix(raw_matrix);
        cf_free_matrix(filtered_matrix);
        return -1;
    }
    for (int c = 0; c < filtered_matrix->num_cols; c++) filtered_barcodes_sorted[c] = filtered_matrix->col_names[c];
    for (int c = 0; c < raw_matrix->num_cols; c++) raw_barcodes_sorted[c] = raw_matrix->col_names[c];
    qsort(filtered_barcodes_sorted, filtered_matrix->num_cols, sizeof(char *), cf_compare_string_ptr);
    qsort(raw_barcodes_sorted, raw_matrix->num_cols, sizeof(char *), cf_compare_string_ptr);

    int ambient_barcodes = 0;
    for (int c = 0; c < raw_matrix->num_cols; c++) {
        if (cf_string_set_contains(filtered_barcodes_sorted, filtered_matrix->num_cols, raw_matrix->col_names[c])) {
            raw_is_filtered_cell[c] = 1;
        } else {
            ambient_barcodes++;
        }
    }

    int filtered_cells_absent_from_raw = 0;
    for (int c = 0; c < filtered_matrix->num_cols; c++) {
        if (!cf_string_set_contains(raw_barcodes_sorted, raw_matrix->num_cols, filtered_matrix->col_names[c])) {
            filtered_cells_absent_from_raw++;
        }
    }

    uint64_t ambient_total = 0;
    for (int i = 0; i < raw_matrix->num_entries; i++) {
        int row = (int)raw_matrix->entries[i].row;
        int col = (int)raw_matrix->entries[i].col;
        if (row < 0 || row >= raw_matrix->num_rows || col < 0 || col >= raw_matrix->num_cols) continue;
        int guide = raw_row_to_guide[row];
        if (guide < 0) continue;
        if (raw_is_filtered_cell[col]) continue;
        ambient_counts[guide] += raw_matrix->entries[i].value;
        ambient_total += raw_matrix->entries[i].value;
    }

    if (ambient_barcodes == 0 || ambient_total == 0) {
        cf_write_ambient_skip_summary(output_dir, "skipped_no_ambient",
                                      "No non-cell raw guide UMIs were available to estimate ambient rates.");
        free(filtered_barcodes_sorted); free(raw_barcodes_sorted); free(raw_is_filtered_cell);
        free(guide_ids); free(guide_names); free(filtered_row_to_guide); free(raw_row_to_guide);
        free(id_map); free(ambient_counts); free(ambient_rates); free(cell_depth);
        cf_free_feature_table(&raw_features);
        cf_free_feature_table(&filtered_features);
        cf_free_matrix(raw_matrix);
        cf_free_matrix(filtered_matrix);
        return 0;
    }

    for (int g = 0; g < n_guides; g++) {
        ambient_rates[g] = (double)ambient_counts[g] / (double)ambient_total;
    }

    cf_ambient_entry *entries = NULL;
    size_t n_entries = 0;
    size_t entries_cap = 0;
    for (int i = 0; i < filtered_matrix->num_entries; i++) {
        int row = (int)filtered_matrix->entries[i].row;
        int col = (int)filtered_matrix->entries[i].col;
        if (row < 0 || row >= filtered_matrix->num_rows || col < 0 || col >= filtered_matrix->num_cols) continue;
        int guide = filtered_row_to_guide[row];
        if (guide < 0 || filtered_matrix->entries[i].value == 0) continue;
        cell_depth[col] += filtered_matrix->entries[i].value;
        if (cf_append_ambient_entry(&entries, &n_entries, &entries_cap,
                                    col, guide, filtered_matrix->entries[i].value) != 0) {
            free(entries);
            free(filtered_barcodes_sorted); free(raw_barcodes_sorted); free(raw_is_filtered_cell);
            free(guide_ids); free(guide_names); free(filtered_row_to_guide); free(raw_row_to_guide);
            free(id_map); free(ambient_counts); free(ambient_rates); free(cell_depth);
            cf_free_feature_table(&raw_features);
            cf_free_feature_table(&filtered_features);
            cf_free_matrix(raw_matrix);
            cf_free_matrix(filtered_matrix);
            return -1;
        }
    }
    n_entries = cf_collapse_ambient_entries(entries, n_entries);

    for (size_t i = 0; i < n_entries; i++) {
        double lambda = ambient_rates[entries[i].guide] * (double)cell_depth[entries[i].cell];
        entries[i].pvalue = cf_poisson_upper_tail(entries[i].count, lambda);
    }

    uint64_t total_tests = (uint64_t)filtered_matrix->num_cols * (uint64_t)n_guides;
    cf_apply_bh(entries, n_entries, total_tests);

    char **cell_calls = NULL;
    int *cell_num_features = NULL;
    uint64_t *cell_num_umis = NULL;
    uint64_t *cell_min_called_umi = NULL;
    uint64_t *cell_max_called_umi = NULL;
    double *cell_min_q = NULL;
    int cells_no_call = 0, cells_1_feature = 0, cells_multi_feature = 0;
    if (cf_build_cell_calls(entries, n_entries, guide_names, filtered_matrix->num_cols,
                            fdr_threshold, min_umi,
                            &cell_calls, &cell_num_features, &cell_num_umis, &cell_min_q,
                            &cell_min_called_umi, &cell_max_called_umi,
                            &cells_no_call, &cells_1_feature, &cells_multi_feature) != 0) {
        free(entries);
        free(filtered_barcodes_sorted); free(raw_barcodes_sorted); free(raw_is_filtered_cell);
        free(guide_ids); free(guide_names); free(filtered_row_to_guide); free(raw_row_to_guide);
        free(id_map); free(ambient_counts); free(ambient_rates); free(cell_depth);
        cf_free_feature_table(&raw_features);
        cf_free_feature_table(&filtered_features);
        cf_free_matrix(raw_matrix);
        cf_free_matrix(filtered_matrix);
        return -1;
    }

    int ret = 0;
    if (cf_write_ambient_rates(output_dir, guide_ids, guide_names, ambient_counts, ambient_rates, n_guides) != 0) ret = -1;
    if (cf_write_ambient_calls(output_dir, filtered_matrix, cell_calls, cell_num_features,
                               cell_num_umis, cell_min_called_umi, cell_max_called_umi,
                               cell_min_q, fdr_threshold) != 0) ret = -1;
    if (config->emit_sparse_qvalues) {
        if (cf_write_ambient_qvalues(output_dir, entries, n_entries, filtered_matrix,
                                     guide_ids, guide_names, n_guides) != 0) ret = -1;
    }
    if (cf_write_ambient_summary(output_dir, entries, n_entries,
                                 filtered_matrix->num_cols, n_guides, total_tests,
                                 raw_matrix->num_cols, ambient_barcodes, filtered_cells_absent_from_raw,
                                 ambient_total, min_umi, fdr_threshold,
                                 config->emit_sparse_qvalues,
                                 cells_no_call, cells_1_feature, cells_multi_feature) != 0) ret = -1;

    printf("\n=== Ambient-FDR Guide Calling Summary ===\n");
    printf("Filtered cells:          %d\n", filtered_matrix->num_cols);
    printf("Guides:                  %d\n", n_guides);
    printf("Raw barcodes:            %d\n", raw_matrix->num_cols);
    printf("Ambient barcodes:        %d\n", ambient_barcodes);
    printf("Ambient guide UMIs:      %" PRIu64 "\n", ambient_total);
    printf("Observed guide entries:  %zu\n", n_entries);
    printf("Total BH tests:          %" PRIu64 "\n", total_tests);
    printf("Calls at FDR %.4g:       %d singlet, %d multiplet, %d none\n",
           fdr_threshold, cells_1_feature, cells_multi_feature, cells_no_call);

    cf_free_cell_calls(cell_calls, filtered_matrix->num_cols, cell_num_features,
                       cell_num_umis, cell_min_q,
                       cell_min_called_umi, cell_max_called_umi);
    free(entries);
    free(filtered_barcodes_sorted);
    free(raw_barcodes_sorted);
    free(raw_is_filtered_cell);
    free(guide_ids);
    free(guide_names);
    free(filtered_row_to_guide);
    free(raw_row_to_guide);
    free(id_map);
    free(ambient_counts);
    free(ambient_rates);
    free(cell_depth);
    cf_free_feature_table(&raw_features);
    cf_free_feature_table(&filtered_features);
    cf_free_matrix(raw_matrix);
    cf_free_matrix(filtered_matrix);
    return ret;
}

/* ============================================================================
 * NB-EM Feature Calling Implementation (SCEPTRE-style)
 * ============================================================================ */

#include "../include/nbem.h"

/* Forward declaration for internal function */
static cf_nbem_results* cf_call_features_nbem_with_output(const cf_sparse_matrix *matrix, 
                                                           const cf_nbem_config *config,
                                                           const char *posteriors_temp_path,
                                                           double posteriors_threshold);

cf_nbem_config* cf_nbem_config_create(void) {
    cf_nbem_config *config = calloc(1, sizeof(cf_nbem_config));
    if (!config) return NULL;
    
    config->moi_mode = CF_MOI_AUTO;
    config->moi_min_umi = 1;
    config->moi_pmulti_threshold = 0.05;
    config->prob_threshold = 0.8;
    config->backup_threshold = 5;
    config->max_iter = 100;
    config->tol = 1e-4;
    config->use_poisson = 0;
    config->sceptre_parity = 0;
    config->global_phi = 0.0;  /* Auto */
    config->phi_min = 0.01;
    config->phi_max = 1000.0;
    config->covariate_tsv = NULL;
    config->debug = 0;
    config->debug_max_features = 0;
    config->debug_csv = NULL;
    config->debug_feature = NULL;
    config->debug_iter_csv = NULL;
    
    return config;
}

void cf_nbem_config_destroy(cf_nbem_config *config) {
    if (config) free(config);
}

cf_nbem_results* cf_call_features_nbem(const cf_sparse_matrix *matrix, const cf_nbem_config *config) {
    return cf_call_features_nbem_with_output(matrix, config, NULL, 0.01);
}

/* Internal implementation that supports streaming posteriors to file */
static cf_nbem_results* cf_call_features_nbem_with_output(const cf_sparse_matrix *matrix, 
                                                    const cf_nbem_config *config,
                                                    const char *posteriors_temp_path,
                                                    double posteriors_threshold) {
    if (!matrix) return NULL;
    
    /* Use default config if none provided */
    cf_nbem_config default_config;
    if (!config) {
        memset(&default_config, 0, sizeof(default_config));
        default_config.moi_mode = CF_MOI_AUTO;
        default_config.moi_min_umi = 1;
        default_config.moi_pmulti_threshold = 0.05;
        default_config.prob_threshold = 0.8;
        default_config.backup_threshold = 5;
        default_config.max_iter = 100;
        default_config.tol = 1e-4;
        default_config.use_poisson = 0;
        default_config.sceptre_parity = 0;
        default_config.global_phi = 0.0;
        default_config.phi_min = 0.01;
        default_config.phi_max = 1000.0;
        config = &default_config;
    }
    
    if (posteriors_threshold <= 0) posteriors_threshold = 0.01;
    
    int num_features, num_cells;
    int **feature_counts = build_per_feature_counts(matrix, &num_features, &num_cells);
    if (!feature_counts) return NULL;
    
    /* Allocate results */
    cf_nbem_results *results = calloc(1, sizeof(cf_nbem_results));
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
    
    /* Allocate per-cell best/second posterior tracking for LOW-MOI disambiguation */
    results->best_post = calloc(num_cells, sizeof(double));
    results->second_post = calloc(num_cells, sizeof(double));
    results->best_feature = malloc(num_cells * sizeof(int));
    results->second_feature = malloc(num_cells * sizeof(int));
    
    if (!results->feature_calls || !results->num_features || !results->num_umis ||
        !results->best_post || !results->second_post || !results->best_feature || !results->second_feature) {
        cf_free_nbem_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    /* Initialize best/second tracking */
    for (int c = 0; c < num_cells; c++) {
        results->best_feature[c] = -1;
        results->second_feature[c] = -1;
        results->best_post[c] = -1.0;
        results->second_post[c] = -1.0;
    }
    
    /* Allocate per-feature results */
    results->feature_results = calloc(num_features, sizeof(cf_feature_nbem_result));
    if (!results->feature_results) {
        cf_free_nbem_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }

    /* Step 1: Estimate MOI */
    nbem_moi_result moi_raw = nbem_estimate_moi(feature_counts, num_features, num_cells,
                                                  config->moi_min_umi, config->moi_pmulti_threshold);
    results->moi.p0 = moi_raw.p0;
    results->moi.lambda = moi_raw.lambda;
    results->moi.p_multi_exp = moi_raw.p_multi_exp;
    results->moi.p_multi_obs = moi_raw.p_multi_obs;
    results->moi.n_cells_0 = moi_raw.n_cells_0;
    results->moi.n_cells_1 = moi_raw.n_cells_1;
    results->moi.n_cells_multi = moi_raw.n_cells_multi;
    
    cf_moi_mode effective_moi = config->moi_mode;
    if (effective_moi == CF_MOI_AUTO) {
        effective_moi = (moi_raw.classification == NBEM_MOI_HIGH) ? CF_MOI_HIGH : CF_MOI_LOW;
    }
    if (config->sceptre_parity) {
        effective_moi = CF_MOI_HIGH;
    }
    results->moi.classification = effective_moi;
    
    printf("MOI estimation: p0=%.3f, lambda=%.3f, p_multi_obs=%.3f, p_multi_exp=%.3f\n",
           moi_raw.p0, moi_raw.lambda, moi_raw.p_multi_obs, moi_raw.p_multi_exp);
    printf("MOI classification: %s\n", (effective_moi == CF_MOI_HIGH) ? "HIGH" : "LOW");
    
    /* Step 2: Compute covariates */
    nbem_covariates *covariates = nbem_compute_grna_covariates(feature_counts, num_features, num_cells);
    if (!covariates) {
        cf_free_nbem_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    /* Load GEX covariates if provided */
    if (config->covariate_tsv) {
        nbem_load_gex_covariates(covariates, config->covariate_tsv, matrix->col_names);
    }
    
    /* Step 3: Build design matrix (shared across all features) */
    int n_covars = covariates->has_gex ? 5 : 3;
    double *X = malloc(num_cells * n_covars * sizeof(double));
    double *beta = malloc(n_covars * sizeof(double));
    double *mu0 = malloc(num_cells * sizeof(double));  /* Per-feature background, reused */
    
    if (!X || !beta || !mu0) {
        free(X); free(beta); free(mu0);
        nbem_free_covariates(covariates);
        cf_free_nbem_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    nbem_build_design_matrix(covariates, X);
    
    /* Step 4: Estimate global phi using streaming accumulator (no O(features*cells) alloc) */
    double global_phi = config->global_phi;
    if (global_phi <= 0) {
        nbem_phi_accum phi_acc;
        nbem_phi_accum_init(&phi_acc);
        for (int f = 0; f < num_features; f++) {
            nbem_phi_accum_add(&phi_acc, feature_counts[f], num_cells);
        }
        global_phi = nbem_phi_accum_finish(&phi_acc, config->phi_min, config->phi_max);
    }
    results->global_phi = global_phi;
    printf("Global dispersion (phi): %.4f\n", global_phi);
    
    /* Step 5: Run NB-EM for each feature with PER-FEATURE background model */
    printf("Running NB-EM for %d features (per-feature background)...\n", num_features);
    
    /* Scratch arrays (reused per feature) */
    double *posteriors_scratch = malloc(num_cells * sizeof(double));
    int *positive_calls = malloc(num_cells * sizeof(int));
    
    if (!posteriors_scratch || !positive_calls) {
        free(posteriors_scratch);
        free(positive_calls);
        free(X); free(beta); free(mu0);
        nbem_free_covariates(covariates);
        cf_free_nbem_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    /* Track per-cell positive features for final call assembly */
    int **cell_positive_features = calloc(num_cells, sizeof(int*));
    int *cell_num_positive = calloc(num_cells, sizeof(int));
    int *cell_positive_cap = calloc(num_cells, sizeof(int));
    
    if (!cell_positive_features || !cell_num_positive || !cell_positive_cap) {
        free(posteriors_scratch);
        free(positive_calls);
        free(X); free(beta); free(mu0);
        free(cell_positive_features);
        free(cell_num_positive);
        free(cell_positive_cap);
        nbem_free_covariates(covariates);
        cf_free_nbem_results(results);
        free_per_feature_counts(feature_counts, num_features);
        return NULL;
    }
    
    /* Open temp file for streaming posteriors */
    FILE *posteriors_fp = NULL;
    int posteriors_nnz = 0;
    if (posteriors_temp_path) {
        posteriors_fp = fopen(posteriors_temp_path, "w");
        if (!posteriors_fp) {
            fprintf(stderr, "Warning: failed to open posteriors temp file: %s\n", posteriors_temp_path);
        }
        /* Note: we write raw entries first, header added later */
    }
    
    int num_em_failures = 0;
    FILE *debug_fp = NULL;
    int debug_written = 0;
    if (config->debug) {
        if (config->debug_csv) {
            debug_fp = fopen(config->debug_csv, "w");
            if (!debug_fp) {
                fprintf(stderr, "Warning: failed to open NB-EM debug CSV: %s\n", config->debug_csv);
            }
        }
        if (debug_fp) {
            fprintf(debug_fp, "feature,pi,delta,phi,log_likelihood,converged,n_iter,used_fallback,failure_reason,num_positive,total_umis\n");
        }
    }
    
    for (int f = 0; f < num_features; f++) {
        /* === Per-feature background model === */
        /* Fit Poisson GLM for THIS feature's counts (not total counts) */
        int irls_ok = poisson_irls(X, feature_counts[f], num_cells, n_covars, beta, 25, 1e-6);
        
        if (irls_ok != 0) {
            /* Fallback: use feature mean as constant background */
            double feat_mean = 0.0;
            for (int c = 0; c < num_cells; c++) feat_mean += feature_counts[f][c];
            feat_mean /= num_cells;
            if (feat_mean < 1e-10) feat_mean = 1e-10;
            for (int c = 0; c < num_cells; c++) mu0[c] = feat_mean;
        } else {
            /* Predict per-feature background: mu0_f[i] = exp(X_i * beta_f) */
            poisson_predict(X, beta, num_cells, n_covars, mu0);
        }
        
        /* Fit NB-EM for this feature (with optional per-iteration debug) */
        nbem_fit_result fit;
        const char *feat_name = matrix->row_names ? matrix->row_names[f] : NULL;
        double pi_seed = 0.1;
        double delta_seed = 1.0;
        if (config->debug && config->debug_feature && feat_name &&
            strcmp(config->debug_feature, feat_name) == 0) {
            FILE *iter_fp = NULL;
            if (config->debug_iter_csv) {
                iter_fp = fopen(config->debug_iter_csv, "w");
                if (!iter_fp) {
                    fprintf(stderr, "Warning: failed to open NB-EM iter debug CSV: %s\n", config->debug_iter_csv);
                }
            }
            fit = nbem_fit_feature_debug(feature_counts[f], mu0, num_cells,
                                         global_phi, config->max_iter, config->tol,
                                         pi_seed, delta_seed, config->use_poisson, feat_name, iter_fp);
            if (iter_fp) {
                fclose(iter_fp);
            }
        } else {
            fit = nbem_fit_feature(feature_counts[f], mu0, num_cells,
                                   global_phi, config->use_poisson, config->max_iter, config->tol,
                                   pi_seed, delta_seed);
        }
        
        /* Track EM failures */
        if (fit.used_fallback) {
            num_em_failures++;
        }
        
        /* Compute posteriors into scratch buffer */
        nbem_posteriors(feature_counts[f], mu0, num_cells,
                        fit.pi, fit.delta, global_phi, config->use_poisson, posteriors_scratch);
        
        /* Make calls (uses fallback threshold if EM failed) */
        int use_fallback = fit.used_fallback;
        if (config->sceptre_parity) {
            use_fallback = 0;
        }
        int n_positive = nbem_make_calls(posteriors_scratch, feature_counts[f], num_cells,
                                          config->prob_threshold, config->backup_threshold,
                                          use_fallback, positive_calls);
        
        /* Update per-cell best/second posterior tracking (only for positive calls) */
        for (int c = 0; c < num_cells; c++) {
            if (!positive_calls[c]) continue;
            double post = posteriors_scratch[c];
            if (post > results->best_post[c]) {
                /* New best: demote current best to second */
                results->second_post[c] = results->best_post[c];
                results->second_feature[c] = results->best_feature[c];
                results->best_post[c] = post;
                results->best_feature[c] = f;
            } else if (post > results->second_post[c]) {
                /* New second best */
                results->second_post[c] = post;
                results->second_feature[c] = f;
            }
        }
        
        /* Stream posteriors to temp file */
        if (posteriors_fp) {
            for (int c = 0; c < num_cells; c++) {
                if (posteriors_scratch[c] >= posteriors_threshold) {
                    fprintf(posteriors_fp, "%d %d %.6f\n", f + 1, c + 1, posteriors_scratch[c]);
                    posteriors_nnz++;
                }
            }
        }
        
        /* Store feature result */
        results->feature_results[f].feature_index = f;
        results->feature_results[f].feature_name = matrix->row_names ? matrix->row_names[f] : NULL;
        results->feature_results[f].pi = fit.pi;
        results->feature_results[f].delta = fit.delta;
        results->feature_results[f].phi = global_phi;
        results->feature_results[f].log_likelihood = fit.log_likelihood;
        results->feature_results[f].converged = fit.converged;
        results->feature_results[f].n_iter = fit.n_iter;
        results->feature_results[f].used_fallback = fit.used_fallback;
        results->feature_results[f].num_positive = n_positive;
        
        /* Calculate total UMIs */
        int total_umis = 0;
        for (int c = 0; c < num_cells; c++) {
            total_umis += feature_counts[f][c];
        }
        results->feature_results[f].total_umis = total_umis;

        if (debug_fp && (config->debug_max_features <= 0 || debug_written < config->debug_max_features)) {
            const char *name = matrix->row_names ? matrix->row_names[f] : "Unknown";
            fprintf(debug_fp, "%s,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%d\n",
                    name,
                    fit.pi,
                    fit.delta,
                    global_phi,
                    fit.log_likelihood,
                    fit.converged,
                    fit.n_iter,
                    fit.used_fallback,
                    fit.failure_reason,
                    n_positive,
                    total_umis);
            debug_written++;
        }
        
        /* Track positive cells for this feature */
        for (int c = 0; c < num_cells; c++) {
            if (positive_calls[c]) {
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
    }
    
    if (posteriors_fp) {
        fclose(posteriors_fp);
    }
    if (debug_fp) {
        fclose(debug_fp);
    }
    
    results->num_em_failures = num_em_failures;
    results->posteriors_nnz = posteriors_nnz;
    
    free(posteriors_scratch);
    free(positive_calls);
    free(X);
    free(beta);
    free(mu0);
    nbem_free_covariates(covariates);
    
    /* Step 6: Build per-cell call strings (handle MOI mode) */
    /* Now uses best_post/second_post arrays instead of dense posteriors matrix */
    for (int c = 0; c < num_cells; c++) {
        int n_pos = cell_num_positive[c];
        
        /* Check if cell has any UMIs */
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
            results->num_features[c] = 0;
            results->num_umis[c] = 0;
        } else if (n_pos == 1 || effective_moi == CF_MOI_LOW) {
            /* Single assignment (or LOW MOI forces single) */
            int best_f = results->best_feature[c];
            
            if (effective_moi == CF_MOI_LOW && n_pos > 1) {
                /* Use pre-tracked best/second posteriors for disambiguation */
                double best_post_val = results->best_post[c];
                double second_post_val = results->second_post[c];
                
                /* If second best is close, mark ambiguous */
                if (second_post_val >= 0 && (best_post_val - second_post_val) < 0.05) {
                    results->feature_calls[c] = strdup("Multiplet");
                    results->num_features[c] = n_pos;
                    results->num_umis[c] = 0;
                    for (int i = 0; i < n_pos; i++) {
                        results->num_umis[c] += feature_counts[cell_positive_features[c][i]][c];
                    }
                    results->cells_multi_feature++;
                    continue;
                }
            } else if (n_pos == 1) {
                best_f = cell_positive_features[c][0];
            }
            
            if (best_f < 0) best_f = cell_positive_features[c][0];  /* Safety fallback */
            
            results->cells_1_feature++;
            results->feature_calls[c] = strdup(matrix->row_names ? matrix->row_names[best_f] : "Unknown");
            results->num_features[c] = 1;
            results->num_umis[c] = feature_counts[best_f][c];
        } else {
            /* Multiple assignments (HIGH MOI) */
            results->cells_multi_feature++;
            
            /* Build pipe-joined string */
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
            
            results->num_features[c] = n_pos;
            results->num_umis[c] = 0;
            for (int i = 0; i < n_pos; i++) {
                results->num_umis[c] += feature_counts[cell_positive_features[c][i]][c];
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

void cf_free_nbem_results(cf_nbem_results *results) {
    if (!results) return;
    
    if (results->feature_calls) {
        for (int i = 0; i < results->num_cells; i++) {
            free(results->feature_calls[i]);
        }
        free(results->feature_calls);
    }
    free(results->num_features);
    free(results->num_umis);
    
    /* Free per-cell tracking arrays */
    free(results->best_post);
    free(results->second_post);
    free(results->best_feature);
    free(results->second_feature);
    
    free(results->feature_results);
    free(results);
}

int cf_write_nbem_calls_per_cell(const cf_nbem_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;
    
    /* Same format as GMM: cell_barcode,num_features,feature_call,num_umis */
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

int cf_write_nbem_calls_summary(const cf_nbem_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;
    
    /* Same format as GMM */
    fprintf(fp, "Category,Metric,Value\n");
    
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

int cf_write_nbem_feature_params(const cf_nbem_results *results, const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;
    
    fprintf(fp, "feature,pi,delta,phi,log_likelihood,converged,n_iter,used_fallback,num_positive,total_umis\n");
    
    for (int f = 0; f < results->num_features_total; f++) {
        const cf_feature_nbem_result *fr = &results->feature_results[f];
        const char *name = results->feature_names ? results->feature_names[f] : "Unknown";
        
        fprintf(fp, "%s,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d\n",
                name,
                fr->pi,
                fr->delta,
                fr->phi,
                fr->log_likelihood,
                fr->converged,
                fr->n_iter,
                fr->used_fallback,
                fr->num_positive,
                fr->total_umis);
    }
    
    fclose(fp);
    return 0;
}

int cf_finalize_nbem_posteriors_mtx(const char *temp_path, const char *output_path,
                                     int num_features, int num_cells, int nnz) {
    if (!temp_path || !output_path) return -1;
    
    /* Read temp file contents */
    FILE *temp_fp = fopen(temp_path, "r");
    if (!temp_fp) return -1;
    
    FILE *out_fp = fopen(output_path, "w");
    if (!out_fp) {
        fclose(temp_fp);
        return -1;
    }
    
    /* Write MatrixMarket header */
    fprintf(out_fp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(out_fp, "%%\n");
    fprintf(out_fp, "%d %d %d\n", num_features, num_cells, nnz);
    
    /* Copy temp file contents */
    char line[256];
    while (fgets(line, sizeof(line), temp_fp)) {
        fputs(line, out_fp);
    }
    
    fclose(temp_fp);
    fclose(out_fp);
    
    /* Remove temp file */
    remove(temp_path);
    
    return 0;
}

int cf_write_nbem_summary_json(const cf_nbem_results *results, const cf_nbem_config *config,
                                const char *output_path) {
    if (!results || !output_path) return -1;
    
    FILE *fp = fopen(output_path, "w");
    if (!fp) return -1;
    
    fprintf(fp, "{\n");
    
    /* MOI section */
    fprintf(fp, "  \"moi\": {\n");
    fprintf(fp, "    \"classification\": \"%s\",\n", 
            (results->moi.classification == CF_MOI_HIGH) ? "high" : "low");
    fprintf(fp, "    \"p0\": %.6f,\n", results->moi.p0);
    fprintf(fp, "    \"lambda\": %.6f,\n", results->moi.lambda);
    fprintf(fp, "    \"p_multi_obs\": %.6f,\n", results->moi.p_multi_obs);
    fprintf(fp, "    \"p_multi_exp\": %.6f,\n", results->moi.p_multi_exp);
    fprintf(fp, "    \"n_cells_0\": %d,\n", results->moi.n_cells_0);
    fprintf(fp, "    \"n_cells_1\": %d,\n", results->moi.n_cells_1);
    fprintf(fp, "    \"n_cells_multi\": %d\n", results->moi.n_cells_multi);
    fprintf(fp, "  },\n");
    
    /* Parameters section */
    fprintf(fp, "  \"parameters\": {\n");
    fprintf(fp, "    \"global_phi\": %.6f,\n", results->global_phi);
    if (config) {
        fprintf(fp, "    \"prob_threshold\": %.6f,\n", config->prob_threshold);
        fprintf(fp, "    \"backup_threshold\": %d,\n", config->backup_threshold);
        fprintf(fp, "    \"max_iter\": %d,\n", config->max_iter);
        fprintf(fp, "    \"tol\": %.6e,\n", config->tol);
        fprintf(fp, "    \"use_poisson\": %d,\n", config->use_poisson);
        fprintf(fp, "    \"sceptre_parity\": %d\n", config->sceptre_parity);
    } else {
        fprintf(fp, "    \"prob_threshold\": 0.8,\n");
        fprintf(fp, "    \"backup_threshold\": 5,\n");
        fprintf(fp, "    \"max_iter\": 100,\n");
        fprintf(fp, "    \"tol\": 1e-4,\n");
        fprintf(fp, "    \"use_poisson\": 0,\n");
        fprintf(fp, "    \"sceptre_parity\": 0\n");
    }
    fprintf(fp, "  },\n");
    
    /* Summary section */
    fprintf(fp, "  \"summary\": {\n");
    fprintf(fp, "    \"num_cells\": %d,\n", results->num_cells);
    fprintf(fp, "    \"num_features\": %d,\n", results->num_features_total);
    fprintf(fp, "    \"cells_no_molecules\": %d,\n", results->cells_no_molecules);
    fprintf(fp, "    \"cells_no_call\": %d,\n", results->cells_no_call);
    fprintf(fp, "    \"cells_1_feature\": %d,\n", results->cells_1_feature);
    fprintf(fp, "    \"cells_multi_feature\": %d,\n", results->cells_multi_feature);
    fprintf(fp, "    \"num_em_failures\": %d\n", results->num_em_failures);
    fprintf(fp, "  }\n");
    
    fprintf(fp, "}\n");
    
    fclose(fp);
    return 0;
}

int cf_process_mex_dir_nbem(const char *mex_dir, const char *output_dir, const cf_nbem_config *config) {
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
    
    /* Set up temp file for streaming posteriors */
    char temp_path[MAX_LINE_LENGTH];
    char output_path[MAX_LINE_LENGTH];
    snprintf(temp_path, sizeof(temp_path), "%s/.nbem_posteriors_temp.txt", output_dir);
    snprintf(output_path, sizeof(output_path), "%s/nbem_cell_posteriors.mtx", output_dir);
    
    double posteriors_threshold = 0.01;
    
    /* Call features with NB-EM (with streaming posteriors) */
    printf("Calling features with NB-EM (SCEPTRE-style)...\n");
    cf_nbem_config local_cfg;
    const cf_nbem_config *cfg_ptr = config;
    if (config && config->debug && !config->debug_csv) {
        local_cfg = *config;
        static char debug_path[MAX_LINE_LENGTH];
        snprintf(debug_path, sizeof(debug_path), "%s/nbem_debug.csv", output_dir);
        local_cfg.debug_csv = debug_path;
        cfg_ptr = &local_cfg;
    }
    cf_nbem_results *results = cf_call_features_nbem_with_output(matrix, cfg_ptr,
                                                                  temp_path, posteriors_threshold);
    if (!results) {
        fprintf(stderr, "Failed to call features\n");
        cf_free_matrix(matrix);
        return -1;
    }
    
    /* Write outputs */
    char path[MAX_LINE_LENGTH];
    
    /* Standard outputs (same as GMM) */
    snprintf(path, sizeof(path), "%s/protospacer_calls_per_cell.csv", output_dir);
    printf("Writing: %s\n", path);
    cf_write_nbem_calls_per_cell(results, path);
    
    snprintf(path, sizeof(path), "%s/protospacer_calls_summary.csv", output_dir);
    printf("Writing: %s\n", path);
    cf_write_nbem_calls_summary(results, path);
    
    /* NB-EM specific outputs */
    snprintf(path, sizeof(path), "%s/nbem_feature_params.csv", output_dir);
    printf("Writing: %s\n", path);
    cf_write_nbem_feature_params(results, path);
    
    /* Finalize posteriors MTX (add header to temp file) */
    if (file_exists(temp_path)) {
        printf("Writing: %s\n", output_path);
        if (cf_finalize_nbem_posteriors_mtx(temp_path, output_path,
                                             results->num_features_total, results->num_cells,
                                             results->posteriors_nnz) != 0) {
            fprintf(stderr, "Warning: failed to finalize posteriors MTX: %s\n", output_path);
        }
    } else {
        fprintf(stderr, "Warning: missing posteriors temp file; skipping MTX output\n");
    }
    
    snprintf(path, sizeof(path), "%s/nbem_summary.json", output_dir);
    printf("Writing: %s\n", path);
    cf_write_nbem_summary_json(results, config, path);
    
    printf("\n=== NB-EM Calling Summary ===\n");
    printf("Total cells:              %d\n", results->num_cells);
    printf("Cells with 0 molecules:   %d\n", results->cells_no_molecules);
    printf("Cells with no call:       %d\n", results->cells_no_call);
    printf("Cells with 1 feature:     %d\n", results->cells_1_feature);
    printf("Cells with >1 features:   %d\n", results->cells_multi_feature);
    printf("Global phi:               %.4f\n", results->global_phi);
    printf("MOI classification:       %s\n", 
           (results->moi.classification == CF_MOI_HIGH) ? "HIGH" : "LOW");
    printf("EM failures (fallback):   %d/%d features\n", 
           results->num_em_failures, results->num_features_total);
    
    cf_free_nbem_results(results);
    cf_free_matrix(matrix);
    return 0;
}
