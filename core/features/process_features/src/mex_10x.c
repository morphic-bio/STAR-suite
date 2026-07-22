#include "../include/mex_10x.h"
#include "../include/common.h"
#include "../include/globals.h"
#include "../include/prototypes.h"
#include "../include/pf_counts.h"
#include "../include/io.h"
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>

typedef struct {
    int gzip;
    FILE *fp;
    gzFile gz;
} out_handle;

static int ensure_dir(const char *path) {
    struct stat st;
    if (stat(path, &st) == 0 && S_ISDIR(st.st_mode)) {
        return 0;
    }
    if (mkdir(path, 0755) == 0) {
        return 0;
    }
    return -1;
}

static int out_open(out_handle *out, const char *path, int gzip) {
    memset(out, 0, sizeof(*out));
    out->gzip = gzip;
    if (gzip) {
        out->gz = gzopen(path, "wb");
        return out->gz ? 0 : -1;
    }
    out->fp = fopen(path, "w");
    return out->fp ? 0 : -1;
}

static void out_close(out_handle *out) {
    if (!out) return;
    if (out->gzip) {
        if (out->gz) gzclose(out->gz);
    } else {
        if (out->fp) fclose(out->fp);
    }
}

static int out_vprintf(out_handle *out, const char *fmt, va_list args) {
    char buf[4096];
    int n = vsnprintf(buf, sizeof(buf), fmt, args);
    if (n < 0) return -1;
    if (out->gzip) {
        return gzputs(out->gz, buf) >= 0 ? 0 : -1;
    }
    return fputs(buf, out->fp) >= 0 ? 0 : -1;
}

static int out_printf(out_handle *out, const char *fmt, ...) {
    int rc;
    va_list args;
    va_start(args, fmt);
    rc = out_vprintf(out, fmt, args);
    va_end(args);
    return rc;
}

static int copy_lines(FILE *in, out_handle *out) {
    char line[4096];
    while (fgets(line, sizeof(line), in)) {
        if (out_printf(out, "%s", line) != 0) return -1;
    }
    return 0;
}

static int write_barcodes(const char *input_dir, const char *output_dir, int gzip) {
    char in_path[4096];
    char out_path[4096];
    snprintf(in_path, sizeof(in_path), "%s/barcodes.txt", input_dir);
    snprintf(out_path, sizeof(out_path), "%s/barcodes.tsv%s", output_dir, gzip ? ".gz" : "");

    FILE *in = fopen(in_path, "r");
    if (!in) return -1;

    out_handle out;
    if (out_open(&out, out_path, gzip) != 0) {
        fclose(in);
        return -1;
    }

    int rc = copy_lines(in, &out);
    fclose(in);
    out_close(&out);
    return rc;
}

static int write_features_from_ref(const char *output_dir, feature_arrays *features, const char *default_feature_type, int gzip) {
    char out_path[4096];
    snprintf(out_path, sizeof(out_path), "%s/features.tsv%s", output_dir, gzip ? ".gz" : "");

    out_handle out;
    if (out_open(&out, out_path, gzip) != 0) {
        return -1;
    }

    for (int i = 0; i < features->number_of_features; ++i) {
        const char *id = features->feature_ids ? features->feature_ids[i] : features->feature_names[i];
        const char *name = features->feature_names[i];
        const char *ftype = default_feature_type;
        if (!id || !id[0]) id = name;
        if (!name || !name[0]) name = id;
        if (out_printf(&out, "%s\t%s\t%s\n", id, name, ftype) != 0) {
            out_close(&out);
            return -1;
        }
    }

    out_close(&out);
    return 0;
}

static int write_features(const char *input_dir, const char *output_dir, const char *feature_type, int gzip) {
    char in_path[4096];
    char out_path[4096];
    snprintf(in_path, sizeof(in_path), "%s/features.txt", input_dir);
    snprintf(out_path, sizeof(out_path), "%s/features.tsv%s", output_dir, gzip ? ".gz" : "");

    FILE *in = fopen(in_path, "r");
    if (!in) return -1;

    out_handle out;
    if (out_open(&out, out_path, gzip) != 0) {
        fclose(in);
        return -1;
    }

    char line[4096];
    while (fgets(line, sizeof(line), in)) {
        char *name = line;
        size_t len = strlen(name);
        while (len > 0 && (name[len - 1] == '\n' || name[len - 1] == '\r')) {
            name[--len] = '\0';
        }
        if (len == 0) continue;
        if (out_printf(&out, "%s\t%s\t%s\n", name, name, feature_type) != 0) {
            fclose(in);
            out_close(&out);
            return -1;
        }
    }

    fclose(in);
    out_close(&out);
    return 0;
}

static int write_matrix(const char *input_dir, const char *output_dir, int gzip) {
    char in_path[4096];
    char out_path[4096];
    snprintf(in_path, sizeof(in_path), "%s/matrix.mtx", input_dir);
    snprintf(out_path, sizeof(out_path), "%s/matrix.mtx%s", output_dir, gzip ? ".gz" : "");

    FILE *in = fopen(in_path, "r");
    if (!in) return -1;

    out_handle out;
    if (out_open(&out, out_path, gzip) != 0) {
        fclose(in);
        return -1;
    }

    char line[4096];
    int header_written = 0;
    while (fgets(line, sizeof(line), in)) {
        if (!header_written) {
            if (out_printf(&out, "%%%%MatrixMarket matrix coordinate integer general\n") != 0) {
                fclose(in);
                out_close(&out);
                return -1;
            }
            header_written = 1;
            continue;
        }
        if (line[0] == '%') {
            if (out_printf(&out, "%s", line) != 0) {
                fclose(in);
                out_close(&out);
                return -1;
            }
            continue;
        }
        break;
    }

    if (line[0] != '\0') {
        if (out_printf(&out, "%s", line) != 0) {
            fclose(in);
            out_close(&out);
            return -1;
        }
    }

    while (fgets(line, sizeof(line), in)) {
        long row = 0, col = 0;
        double val = 0.0;
        if (sscanf(line, "%ld %ld %lf", &row, &col, &val) != 3) continue;
        if (out_printf(&out, "%ld %ld %ld\n", row, col, (long)(val + 0.5)) != 0) {
            fclose(in);
            out_close(&out);
            return -1;
        }
    }

    fclose(in);
    out_close(&out);
    return 0;
}

int pf_write_mex_10x_from_features(const char *input_dir,
                                   const char *output_dir,
                                   feature_arrays *features,
                                   const char *default_feature_type,
                                   int gzip_output) {
    if (!input_dir || !output_dir || !features || !default_feature_type) return -1;
    if (ensure_dir(output_dir) != 0) return -1;

    int same_dir = (strcmp(input_dir, output_dir) == 0);

    if (write_barcodes(input_dir, output_dir, 0) != 0) return -1;
    if (write_features_from_ref(output_dir, features, default_feature_type, 0) != 0) return -1;
    if (!same_dir && write_matrix(input_dir, output_dir, 0) != 0) return -1;

    if (gzip_output) {
        if (write_barcodes(input_dir, output_dir, 1) != 0) return -1;
        if (write_features_from_ref(output_dir, features, default_feature_type, 1) != 0) return -1;
        if (write_matrix(input_dir, output_dir, 1) != 0) return -1;
    }

    return 0;
}

static uint32_t pf_barcode_string_to_key(const char *barcode) {
    unsigned char code[8];
    memset(code, 0, sizeof(code));
    string2code((char *)barcode, barcode_length, code);
    return *(uint32_t *)code;
}

int pf_write_mex_10x_from_counts_subset(const char *assign_output_dir,
                                          const char *output_dir,
                                          feature_arrays *features,
                                          pf_counts_result *counts,
                                          const unsigned char *feature_include_mask,
                                          const char *default_feature_type,
                                          int gzip_output) {
    if (!assign_output_dir || !output_dir || !features || !counts ||
        !feature_include_mask || !default_feature_type) {
        return -1;
    }
    if (ensure_dir(output_dir) != 0) return -1;

    const int n_features = features->number_of_features;
    int n_selected = 0;
    int *old_to_new = calloc((size_t)n_features + 1, sizeof(int));
    if (!old_to_new) return -1;
    for (int i = 0; i < n_features; ++i) {
        if (feature_include_mask[i]) {
            n_selected++;
            old_to_new[i + 1] = n_selected;
        }
    }
    if (n_selected == 0) {
        free(old_to_new);
        return -1;
    }

    char barcodes_in[4096];
    snprintf(barcodes_in, sizeof(barcodes_in), "%s/barcodes.txt", assign_output_dir);
    FILE *bc_in = fopen(barcodes_in, "r");
    if (!bc_in) {
        free(old_to_new);
        return -1;
    }

    char **barcodes = NULL;
    size_t n_barcodes = 0;
    size_t barcodes_cap = 0;
    char line[4096];
    while (fgets(line, sizeof(line), bc_in)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[--len] = '\0';
        }
        if (len == 0) continue;
        if (n_barcodes >= barcodes_cap) {
            barcodes_cap = barcodes_cap ? barcodes_cap * 2 : 64;
            char **tmp = realloc(barcodes, barcodes_cap * sizeof(char *));
            if (!tmp) {
                fclose(bc_in);
                for (size_t i = 0; i < n_barcodes; ++i) free(barcodes[i]);
                free(barcodes);
                free(old_to_new);
                return -1;
            }
            barcodes = tmp;
        }
        barcodes[n_barcodes] = strdup(line);
        if (!barcodes[n_barcodes]) {
            fclose(bc_in);
            for (size_t i = 0; i < n_barcodes; ++i) free(barcodes[i]);
            free(barcodes);
            free(old_to_new);
            return -1;
        }
        n_barcodes++;
    }
    fclose(bc_in);

    char barcodes_out[4096];
    char features_out[4096];
    char matrix_out[4096];
    snprintf(barcodes_out, sizeof(barcodes_out), "%s/barcodes.tsv", output_dir);
    snprintf(features_out, sizeof(features_out), "%s/features.tsv", output_dir);
    snprintf(matrix_out, sizeof(matrix_out), "%s/matrix.mtx", output_dir);

    FILE *bc_fp = fopen(barcodes_out, "w");
    FILE *feat_fp = fopen(features_out, "w");
    if (!bc_fp || !feat_fp) {
        if (bc_fp) fclose(bc_fp);
        if (feat_fp) fclose(feat_fp);
        for (size_t i = 0; i < n_barcodes; ++i) free(barcodes[i]);
        free(barcodes);
        free(old_to_new);
        return -1;
    }

    for (int i = 0; i < n_features; ++i) {
        if (!feature_include_mask[i]) continue;
        const char *id = features->feature_ids ? features->feature_ids[i] : features->feature_names[i];
        const char *name = features->feature_names[i];
        if (!id || !id[0]) id = name;
        if (!name || !name[0]) name = id;
        fprintf(feat_fp, "%s\t%s\t%s\n", id, name, default_feature_type);
    }

    for (size_t col = 0; col < n_barcodes; ++col) {
        fprintf(bc_fp, "%s\n", barcodes[col]);
    }
    fclose(bc_fp);
    fclose(feat_fp);

    char matrix_tmp[4096];
    snprintf(matrix_tmp, sizeof(matrix_tmp), "%s/matrix_body.tmp", output_dir);
    FILE *body_fp = fopen(matrix_tmp, "w");
    if (!body_fp) {
        for (size_t i = 0; i < n_barcodes; ++i) free(barcodes[i]);
        free(barcodes);
        free(old_to_new);
        return -1;
    }

    long nnz = 0;
    for (size_t col = 0; col < n_barcodes; ++col) {
        uint32_t bkey = pf_barcode_string_to_key(barcodes[col]);
        khint_t kb = kh_get(u32ptr, counts->barcode_to_deduped_hash, bkey);
        if (kb == kh_end(counts->barcode_to_deduped_hash)) continue;
        khash_t(u32u32) *feat_hash = (khash_t(u32u32) *)kh_val(counts->barcode_to_deduped_hash, kb);
        for (khint_t kf = kh_begin(feat_hash); kf != kh_end(feat_hash); ++kf) {
            if (!kh_exist(feat_hash, kf)) continue;
            uint32_t f_idx = kh_key(feat_hash, kf);
            if (f_idx == 0 || f_idx > (uint32_t)n_features) continue;
            if (!feature_include_mask[f_idx - 1]) continue;
            int new_row = old_to_new[f_idx];
            if (new_row <= 0) continue;
            uint32_t val = kh_val(feat_hash, kf);
            if (val == 0) continue;
            fprintf(body_fp, "%d %zu %u\n", new_row, col + 1, val);
            nnz++;
        }
    }
    fclose(body_fp);

    FILE *mtx_fp = fopen(matrix_out, "w");
    if (!mtx_fp) {
        remove(matrix_tmp);
        for (size_t i = 0; i < n_barcodes; ++i) free(barcodes[i]);
        free(barcodes);
        free(old_to_new);
        return -1;
    }
    fprintf(mtx_fp, "%%%%MatrixMarket matrix coordinate integer general\n");
    fprintf(mtx_fp, "%%metadata_json: {\"software_version\": \"assignBarcodes-0.1\", \"format_version\": 1}\n");
    fprintf(mtx_fp, "%d %zu %ld\n", n_selected, n_barcodes, nnz);
    FILE *body_in = fopen(matrix_tmp, "r");
    if (!body_in) {
        fclose(mtx_fp);
        remove(matrix_tmp);
        for (size_t i = 0; i < n_barcodes; ++i) free(barcodes[i]);
        free(barcodes);
        free(old_to_new);
        return -1;
    }
    while (fgets(line, sizeof(line), body_in)) {
        fputs(line, mtx_fp);
    }
    fclose(body_in);
    fclose(mtx_fp);
    remove(matrix_tmp);

    if (gzip_output) {
        char barcodes_gz[4096];
        char features_gz[4096];
        char matrix_gz[4096];
        snprintf(barcodes_gz, sizeof(barcodes_gz), "%s/barcodes.tsv.gz", output_dir);
        snprintf(features_gz, sizeof(features_gz), "%s/features.tsv.gz", output_dir);
        snprintf(matrix_gz, sizeof(matrix_gz), "%s/matrix.mtx.gz", output_dir);

        FILE *plain = fopen(barcodes_out, "r");
        out_handle gz_out;
        if (plain && out_open(&gz_out, barcodes_gz, 1) == 0) {
            copy_lines(plain, &gz_out);
            fclose(plain);
            out_close(&gz_out);
        } else if (plain) {
            fclose(plain);
        }

        plain = fopen(features_out, "r");
        if (plain && out_open(&gz_out, features_gz, 1) == 0) {
            copy_lines(plain, &gz_out);
            fclose(plain);
            out_close(&gz_out);
        } else if (plain) {
            fclose(plain);
        }

        plain = fopen(matrix_out, "r");
        if (plain && out_open(&gz_out, matrix_gz, 1) == 0) {
            copy_lines(plain, &gz_out);
            fclose(plain);
            out_close(&gz_out);
        } else if (plain) {
            fclose(plain);
        }
    }

    for (size_t i = 0; i < n_barcodes; ++i) free(barcodes[i]);
    free(barcodes);
    free(old_to_new);
    return 0;
}

int pf_write_mex_10x(const char *input_dir,
                     const char *output_dir,
                     const char *feature_type,
                     int gzip_output) {
    if (!input_dir || !output_dir || !feature_type) return -1;
    if (ensure_dir(output_dir) != 0) return -1;

    /* Always write plain TSV/MTX for local consumers; optionally add .gz */
    if (write_barcodes(input_dir, output_dir, 0) != 0) return -1;
    if (write_features(input_dir, output_dir, feature_type, 0) != 0) return -1;
    if (write_matrix(input_dir, output_dir, 0) != 0) return -1;

    if (gzip_output) {
        if (write_barcodes(input_dir, output_dir, 1) != 0) return -1;
        if (write_features(input_dir, output_dir, feature_type, 1) != 0) return -1;
        if (write_matrix(input_dir, output_dir, 1) != 0) return -1;
    }

    return 0;
}
