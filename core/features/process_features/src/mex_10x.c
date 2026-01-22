#include "../include/mex_10x.h"
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
