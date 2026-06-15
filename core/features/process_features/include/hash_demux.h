#ifndef HASH_DEMUX_H
#define HASH_DEMUX_H

#include "adt_mex.h"
#include "pf_counts.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PF_HASH_FEATURE_TYPE_DEFAULT "Multiplexing Capture"
#define PF_HASH_DEMUX_METHOD_RATIO "ratio"

#define PF_HASH_DEMUX_AUTO (-1)
#define PF_HASH_DEMUX_NO 0
#define PF_HASH_DEMUX_YES 1

typedef struct {
    pf_adt_mex_config base;
    int hash_demux_mode;           /* PF_HASH_DEMUX_AUTO|NO|YES */
    const char *hash_feature_selector;
    const char *hash_demux_method;
    const char *library_feature_type; /* pf-multi feature_types hint */
    int hash_min_total;
    int hash_min_top;
    double hash_min_ratio;
} pf_hash_mex_config;

typedef struct {
    int n_hash_features;
    int n_protein_features;
    int n_singlet;
    int n_doublet;
    int n_negative;
} pf_hash_demux_stats;

int pf_is_hash_like_feature_type(const char *feature_type);
int pf_build_hash_feature_mask(feature_arrays *features,
                               const char *selector,
                               const char *library_feature_type,
                               unsigned char *out_mask,
                               int *n_hash_out,
                               int *n_protein_out);

int pf_write_adt_mex_outputs(const pf_hash_mex_config *config);

#ifdef __cplusplus
}
#endif

#endif
