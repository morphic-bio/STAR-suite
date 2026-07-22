#ifndef ADT_MEX_H
#define ADT_MEX_H

#include "common.h"
#include "pf_counts.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PF_ADT_FEATURE_TYPE_DEFAULT "Antibody Capture"

typedef struct {
    const char *assign_output_dir;
    const char *mex_output_dir;
    feature_arrays *features;
    pf_counts_result *counts;
    statistics *stats;
    int barcode_length;
    int umi_length;
    int barcode_offset;
    int feature_offset;
    int max_hamming_distance;
    int stringency;
    const char *command_line;
} pf_adt_mex_config;

int pf_write_adt_protein_outputs(const pf_adt_mex_config *config);

/* Legacy protein-only entry; prefer pf_write_adt_mex_outputs from hash_demux.h. */
int pf_write_adt_protein_outputs_masked(const pf_adt_mex_config *config,
                                        const unsigned char *protein_mask,
                                        const char *mex_output_dir);

#ifdef __cplusplus
}
#endif

#endif
