#ifndef MEX_10X_H
#define MEX_10X_H

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

#ifdef __cplusplus
}
#endif

#endif
