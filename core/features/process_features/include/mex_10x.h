#ifndef MEX_10X_H
#define MEX_10X_H

#include "common.h"
#include "pf_counts.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Write 10x-style MEX files from assignBarcodes outputs.
 * Input directory should contain:
 *   - matrix.mtx
 *   - barcodes.txt
 *   - features.txt
 * Output directory will contain:
 *   - matrix.mtx(.gz) with integer header/values
 *   - barcodes.tsv(.gz)
 *   - features.tsv(.gz) with 3 columns
 */
int pf_write_mex_10x(const char *input_dir,
                     const char *output_dir,
                     const char *feature_type,
                     int gzip_output);

int pf_write_mex_10x_from_features(const char *input_dir,
                                   const char *output_dir,
                                   feature_arrays *features,
                                   const char *default_feature_type,
                                   int gzip_output);

/* Write MEX from deduped counts, including only features where mask[i] != 0. */
int pf_write_mex_10x_from_counts_subset(const char *assign_output_dir,
                                          const char *output_dir,
                                          feature_arrays *features,
                                          pf_counts_result *counts,
                                          const unsigned char *feature_include_mask,
                                          const char *default_feature_type,
                                          int gzip_output);

#ifdef __cplusplus
}
#endif

#endif
