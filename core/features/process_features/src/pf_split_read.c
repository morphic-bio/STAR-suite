#include "../include/pf_split_read.h"
#include "../include/common.h"

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void pf_split_set_error(char *error_buf, size_t error_cap, const char *fmt, ...) {
    if (!error_buf || error_cap == 0 || !fmt) {
        return;
    }
    va_list args;
    va_start(args, fmt);
    vsnprintf(error_buf, error_cap, fmt, args);
    va_end(args);
}

static char pf_complement_base(char base) {
    switch (toupper((unsigned char)base)) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return (char)toupper((unsigned char)base);
    }
}

void pf_revcomp_seq(char *seq, size_t len) {
    for (size_t i = 0; i < len / 2; ++i) {
        char a = pf_complement_base(seq[i]);
        char b = pf_complement_base(seq[len - 1 - i]);
        seq[i] = b;
        seq[len - 1 - i] = a;
    }
    if (len % 2) {
        seq[len / 2] = pf_complement_base(seq[len / 2]);
    }
}

void pf_revcomp_qual(char *qual, size_t len) {
    for (size_t i = 0; i < len / 2; ++i) {
        char tmp = qual[i];
        qual[i] = qual[len - 1 - i];
        qual[len - 1 - i] = tmp;
    }
}

static int pf_hamming_at_most(const char *haystack, size_t hay_len,
                              const char *needle, size_t needle_len,
                              int max_hamming) {
    if (needle_len == 0 || hay_len < needle_len) {
        return 0;
    }
    for (size_t offset = 0; offset + needle_len <= hay_len; ++offset) {
        int dist = 0;
        for (size_t i = 0; i < needle_len; ++i) {
            if (toupper((unsigned char)haystack[offset + i]) !=
                toupper((unsigned char)needle[i])) {
                if (++dist > max_hamming) {
                    break;
                }
            }
        }
        if (dist <= max_hamming) {
            return 1;
        }
    }
    return 0;
}

static int pf_substring_search(const char *haystack, size_t hay_len,
                               const char *needle, size_t needle_len) {
    if (needle_len == 0 || hay_len < needle_len) {
        return 0;
    }
    for (size_t offset = 0; offset + needle_len <= hay_len; ++offset) {
        if (memcmp(haystack + offset, needle, needle_len) == 0) {
            return 1;
        }
    }
    return 0;
}

static int pf_capture_match(const char *seq, size_t seq_len,
                            const pf_capture_pattern_view *patterns, int n_patterns,
                            int max_hamming) {
    for (int i = 0; i < n_patterns; ++i) {
        if (max_hamming <= 0) {
            if (pf_substring_search(seq, seq_len, patterns[i].seq, patterns[i].len)) {
                return 1;
            }
        } else if (pf_hamming_at_most(seq, seq_len, patterns[i].seq,
                                      patterns[i].len, max_hamming)) {
            return 1;
        }
    }
    return 0;
}

static int pf_tokenize_capture_patterns(char *buf,
                                        pf_capture_pattern_view *patterns,
                                        int max_patterns) {
    if (!buf || !buf[0]) {
        return 0;
    }
    int count = 0;
    char *saveptr = NULL;
    for (char *token = strtok_r(buf, "|", &saveptr);
         token != NULL;
         token = strtok_r(NULL, "|", &saveptr)) {
        while (*token && isspace((unsigned char)*token)) {
            ++token;
        }
        char *end = token + strlen(token);
        while (end > token && isspace((unsigned char)end[-1])) {
            --end;
        }
        *end = '\0';
        if (*token == '\0') {
            continue;
        }
        if (count >= max_patterns) {
            return -1;
        }
        patterns[count].seq = token;
        patterns[count].len = strlen(token);
        ++count;
    }
    return count;
}

static int pf_split_read_valid_read_idx(int idx) {
    return idx >= 0 && idx <= 2;
}

int pf_split_read_layout_prepare_ex(pf_split_read_layout *layout,
                                    char *error_buf,
                                    size_t error_cap) {
    if (!layout) {
        pf_split_set_error(error_buf, error_cap, "Split-read layout is null");
        return -1;
    }
    layout->n_capture_patterns = 0;
    layout->capture_prepared = 0;
    memset(layout->capture_patterns, 0, sizeof(layout->capture_patterns));

    int bc_start = 0;
    int bc_end_inclusive = 0;
    int bc_revcomp = 0;
    if (!pf_parse_bc_format(layout->barcode_format, &bc_start, &bc_end_inclusive, &bc_revcomp)) {
        pf_split_set_error(error_buf, error_cap,
                           "Invalid split-read barcode format: '%s'",
                           layout->barcode_format);
        return -1;
    }
    (void)bc_revcomp;

    if (!pf_split_read_valid_read_idx(layout->barcode_read_idx)) {
        pf_split_set_error(error_buf, error_cap,
                           "Invalid split-read barcode read index: %d",
                           layout->barcode_read_idx);
        return -1;
    }
    if (!pf_split_read_valid_read_idx(layout->umi_read_idx)) {
        pf_split_set_error(error_buf, error_cap,
                           "Invalid split-read UMI read index: %d",
                           layout->umi_read_idx);
        return -1;
    }
    if (!pf_split_read_valid_read_idx(layout->feature_read_idx)) {
        pf_split_set_error(error_buf, error_cap,
                           "Invalid split-read feature read index: %d",
                           layout->feature_read_idx);
        return -1;
    }
    if (layout->umi_start < 0 || layout->umi_length <= 0) {
        pf_split_set_error(error_buf, error_cap,
                           "Invalid split-read UMI window: start=%d length=%d",
                           layout->umi_start, layout->umi_length);
        return -1;
    }
    const int bc_len = bc_end_inclusive - bc_start + 1;
    if ((size_t)bc_len + (size_t)layout->umi_length + 1 > LINE_LENGTH) {
        pf_split_set_error(error_buf, error_cap,
                           "Split-read barcode+UMI length exceeds buffer: barcode=%d UMI=%d",
                           bc_len, layout->umi_length);
        return -1;
    }

    if (!layout->capture_sequences[0]) {
        layout->capture_prepared = 1;
        return 0;
    }
    if (!pf_split_read_valid_read_idx(layout->capture_read_idx)) {
        pf_split_set_error(error_buf, error_cap,
                           "Invalid split-read capture read index: %d",
                           layout->capture_read_idx);
        return -1;
    }
    if (layout->capture_max_hamming < 0) {
        pf_split_set_error(error_buf, error_cap,
                           "Invalid split-read capture hamming distance: %d",
                           layout->capture_max_hamming);
        return -1;
    }
    strncpy(layout->capture_parse_buf, layout->capture_sequences,
            sizeof(layout->capture_parse_buf) - 1);
    layout->capture_parse_buf[sizeof(layout->capture_parse_buf) - 1] = '\0';
    layout->n_capture_patterns = pf_tokenize_capture_patterns(
        layout->capture_parse_buf,
        layout->capture_patterns,
        PF_SPLIT_MAX_CAPTURES);
    if (layout->n_capture_patterns < 0) {
        pf_split_set_error(error_buf, error_cap,
                           "Too many split-read capture sequences; maximum is %d",
                           PF_SPLIT_MAX_CAPTURES);
        layout->n_capture_patterns = 0;
        return -1;
    }
    if (layout->n_capture_patterns == 0) {
        pf_split_set_error(error_buf, error_cap,
                           "No non-empty split-read capture sequences were parsed");
        return -1;
    }
    layout->capture_prepared = 1;
    return 0;
}

int pf_split_read_layout_prepare(pf_split_read_layout *layout) {
    return pf_split_read_layout_prepare_ex(layout, NULL, 0);
}

int pf_parse_bc_format(const char *format, int *start_out, int *end_inclusive_out, int *revcomp_out) {
    if (!format || !start_out || !end_inclusive_out || !revcomp_out) {
        return 0;
    }
    int start = -1;
    int end_inclusive = -1;
    char strand = '+';
    if (sscanf(format, "bc:%d:%d:%c", &start, &end_inclusive, &strand) != 3) {
        return 0;
    }
    if (start < 0 || end_inclusive < start) {
        return 0;
    }
    *start_out = start;
    *end_inclusive_out = end_inclusive;
    *revcomp_out = (strand == '-' || strand == 'r' || strand == 'R');
    return 1;
}

static int pf_extract_substring(const char *seq, const char *qual, int start, int length,
                                char *out_seq, size_t out_seq_cap,
                                char *out_qual, size_t out_qual_cap) {
    if (!seq || !out_seq || out_seq_cap == 0 || length <= 0 ||
        (size_t)length + 1 > out_seq_cap ||
        (out_qual && out_qual_cap > 0 && (size_t)length + 1 > out_qual_cap)) {
        return 0;
    }
    size_t seq_len = strlen(seq);
    size_t qual_len = qual ? strlen(qual) : 0;
    if (start < 0 || (size_t)start + (size_t)length > seq_len) {
        return 0;
    }
    memcpy(out_seq, seq + start, (size_t)length);
    out_seq[length] = '\0';
    if (out_qual && out_qual_cap > 0) {
        if ((size_t)start + (size_t)length <= qual_len) {
            memcpy(out_qual, qual + start, (size_t)length);
        } else {
            memset(out_qual, 'I', (size_t)length);
        }
        out_qual[length] = '\0';
    }
    return 1;
}

int pf_split_read_synthesize_record(const pf_split_read_layout *layout,
                                    const char *seqs[3],
                                    const char *quals[3],
                                    char *synth_seq_out,
                                    size_t synth_seq_cap,
                                    char *synth_qual_out,
                                    size_t synth_qual_cap,
                                    const char **feature_seq_out,
                                    const char **feature_qual_out,
                                    pf_split_read_metrics *metrics) {
    if (!layout || !seqs || !quals || !synth_seq_out || synth_seq_cap == 0 ||
        !synth_qual_out || synth_qual_cap == 0 ||
        !feature_seq_out || !feature_qual_out) {
        return -1;
    }
    if (!layout->capture_prepared) {
        return -1;
    }

    if (metrics) {
        metrics->total_reads++;
        metrics->capture_total++;
    }

    const int capture_enabled = layout->capture_sequences[0] != '\0';
    int capture_hit = 1;
    if (capture_enabled) {
        if (layout->capture_read_idx < 0 || layout->capture_read_idx > 2 ||
            !seqs[layout->capture_read_idx]) {
            return -1;
        }
        const char *capture_seq = seqs[layout->capture_read_idx];
        const size_t capture_len = strlen(capture_seq);
        if (layout->n_capture_patterns == 0) {
            capture_hit = 1;
        } else {
            capture_hit = pf_capture_match(capture_seq, capture_len,
                                           layout->capture_patterns,
                                           layout->n_capture_patterns,
                                           layout->capture_max_hamming);
            if (metrics && layout->n_capture_patterns >= 1 &&
                pf_capture_match(capture_seq, capture_len,
                                 &layout->capture_patterns[0], 1,
                                 layout->capture_max_hamming)) {
                metrics->capture_cs1_hits++;
            }
            if (metrics && layout->n_capture_patterns >= 2 &&
                pf_capture_match(capture_seq, capture_len,
                                 &layout->capture_patterns[1], 1,
                                 layout->capture_max_hamming)) {
                metrics->capture_cs2_hits++;
            }
        }
        if (!capture_hit) {
            if (metrics) {
                metrics->capture_filtered_out++;
            }
            return 0;
        }
    }
    if (metrics) {
        metrics->capture_either_hits++;
    }

    int bc_start = 0;
    int bc_end_inclusive = 0;
    int bc_revcomp = 0;
    if (!pf_parse_bc_format(layout->barcode_format, &bc_start, &bc_end_inclusive, &bc_revcomp)) {
        return -1;
    }
    const int bc_len = bc_end_inclusive - bc_start + 1;
    if (layout->barcode_read_idx < 0 || layout->barcode_read_idx > 2 ||
        layout->umi_read_idx < 0 || layout->umi_read_idx > 2 ||
        layout->feature_read_idx < 0 || layout->feature_read_idx > 2 ||
        !seqs[layout->barcode_read_idx] || !seqs[layout->umi_read_idx] ||
        !seqs[layout->feature_read_idx]) {
        return -1;
    }
    if (bc_len <= 0 || layout->umi_length <= 0 ||
        (size_t)bc_len + (size_t)layout->umi_length + 1 > synth_seq_cap ||
        (size_t)bc_len + (size_t)layout->umi_length + 1 > synth_qual_cap) {
        return -1;
    }

    char barcode_part[LINE_LENGTH];
    char barcode_qual_part[LINE_LENGTH];
    char umi_part[LINE_LENGTH];
    char umi_qual_part[LINE_LENGTH];
    if (!pf_extract_substring(seqs[layout->barcode_read_idx],
                              quals[layout->barcode_read_idx],
                              bc_start, bc_len,
                              barcode_part, sizeof(barcode_part),
                              barcode_qual_part, sizeof(barcode_qual_part))) {
        return 0;
    }
    if (!pf_extract_substring(seqs[layout->umi_read_idx],
                              quals[layout->umi_read_idx],
                              layout->umi_start, layout->umi_length,
                              umi_part, sizeof(umi_part),
                              umi_qual_part, sizeof(umi_qual_part))) {
        return 0;
    }
    if (bc_revcomp) {
        pf_revcomp_seq(barcode_part, (size_t)bc_len);
        pf_revcomp_qual(barcode_qual_part, (size_t)bc_len);
    }

    int wrote = snprintf(synth_seq_out, synth_seq_cap, "%s%s", barcode_part, umi_part);
    int wrote_qual = snprintf(synth_qual_out, synth_qual_cap, "%s%s",
                              barcode_qual_part, umi_qual_part);
    if (wrote < 0 || wrote_qual < 0 ||
        (size_t)wrote >= synth_seq_cap || (size_t)wrote_qual >= synth_qual_cap) {
        return -1;
    }
    *feature_seq_out = seqs[layout->feature_read_idx];
    *feature_qual_out = quals[layout->feature_read_idx];

    if (metrics) {
        metrics->umi_valid++;
        metrics->barcode_synth_ok++;
    }
    return 1;
}

void pf_write_split_read_metrics(const char *output_dir,
                                 const pf_split_read_metrics *metrics) {
    if (!output_dir || !metrics) {
        return;
    }
    char path[1024];
    snprintf(path, sizeof(path), "%s/split_read_metrics.tsv", output_dir);
    FILE *fp = fopen(path, "w");
    if (!fp) {
        return;
    }
    fprintf(fp, "metric\tvalue\n");
    fprintf(fp, "total_reads\t%llu\n", metrics->total_reads);
    fprintf(fp, "capture_total\t%llu\n", metrics->capture_total);
    fprintf(fp, "capture_cs1_hits\t%llu\n", metrics->capture_cs1_hits);
    fprintf(fp, "capture_cs2_hits\t%llu\n", metrics->capture_cs2_hits);
    fprintf(fp, "capture_either_hits\t%llu\n", metrics->capture_either_hits);
    fprintf(fp, "capture_filtered_out\t%llu\n", metrics->capture_filtered_out);
    fprintf(fp, "umi_valid\t%llu\n", metrics->umi_valid);
    fprintf(fp, "barcode_synth_ok\t%llu\n", metrics->barcode_synth_ok);
    fclose(fp);
}
