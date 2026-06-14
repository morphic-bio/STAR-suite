#ifndef PF_SPLIT_READ_H
#define PF_SPLIT_READ_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PF_SPLIT_MAX_CAPTURES 8

typedef struct {
    const char *seq;
    size_t len;
} pf_capture_pattern_view;

typedef struct {
    int enabled;
    char barcode_format[64];
    int barcode_read_idx;
    int umi_read_idx;
    int umi_start;
    int umi_length;
    int feature_read_idx;
    int capture_read_idx;
    char capture_sequences[1024];
    int capture_max_hamming;
    char capture_parse_buf[1024];
    int n_capture_patterns;
    pf_capture_pattern_view capture_patterns[PF_SPLIT_MAX_CAPTURES];
    int capture_prepared;
    char fastq_r1_pattern[64];
    char fastq_r2_pattern[64];
    char fastq_r3_pattern[64];
} pf_split_read_layout;

typedef struct {
    unsigned long long total_reads;
    unsigned long long capture_total;
    unsigned long long capture_cs1_hits;
    unsigned long long capture_cs2_hits;
    unsigned long long capture_either_hits;
    unsigned long long capture_filtered_out;
    unsigned long long umi_valid;
    unsigned long long barcode_synth_ok;
} pf_split_read_metrics;

int pf_parse_bc_format(const char *format, int *start_out, int *end_inclusive_out, int *revcomp_out);

void pf_revcomp_seq(char *seq, size_t len);
void pf_revcomp_qual(char *qual, size_t len);

int pf_split_read_layout_prepare(pf_split_read_layout *layout);
int pf_split_read_layout_prepare_ex(pf_split_read_layout *layout,
                                    char *error_buf,
                                    size_t error_cap);

int pf_split_read_synthesize_record(const pf_split_read_layout *layout,
                                    const char *seqs[3],
                                    const char *quals[3],
                                    char *synth_seq_out,
                                    size_t synth_seq_cap,
                                    char *synth_qual_out,
                                    size_t synth_qual_cap,
                                    const char **feature_seq_out,
                                    const char **feature_qual_out,
                                    pf_split_read_metrics *metrics);

void pf_write_split_read_metrics(const char *output_dir,
                                 const pf_split_read_metrics *metrics);

#ifdef __cplusplus
}
#endif

#endif
