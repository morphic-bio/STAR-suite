/**
 * @file mex_writer.h
 * @brief MEX format output writer for feature barcode processing
 * 
 * This module handles writing all output files produced by assignBarcodes:
 * - barcodes.txt
 * - features.txt  
 * - matrix.mtx (Matrix Market format)
 * - stats.txt
 * - feature_per_cell.csv
 * - Heatmaps (via generate_heatmap)
 * - Histogram HTML files
 * 
 * The module is designed to produce byte-identical output to the original
 * printFeatureCounts() function while providing a cleaner interface.
 * 
 * WIRING STATUS (Stage 2 Complete):
 * This module is now wired into finalize_processing() via mex_write_all().
 * The old printFeatureCounts() function is no longer called from the main pipeline.
 */

#ifndef MEX_WRITER_H
#define MEX_WRITER_H

#include "common.h"
#include "pf_counts.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Configuration for MEX writer output
 */
typedef struct mex_writer_config {
    /** Base output directory (filtered outputs go to output_dir/filtered/) */
    const char *output_dir;
    
    /** Feature arrays containing feature names and metadata */
    feature_arrays *features;
    
    /** Data structures containing raw counts and whitelist info */
    data_structures *hashes;
    
    /** Processing statistics (for stats.txt output) */
    statistics *stats;
    
    /** Optional filter hash - if non-NULL, only write barcodes in this set */
    khash_t(strptr) *filtered_barcodes_hash;
    
    /** Minimum counts for heatmap display (feature threshold) */
    int min_heatmap_counts;
    
} mex_writer_config;

/**
 * @brief Write all MEX outputs for the given counts result
 * 
 * This is the main entry point that produces all output files.
 * When filtered_barcodes_hash is NULL, outputs go to output_dir.
 * When filtered_barcodes_hash is non-NULL, outputs go to output_dir/filtered/.
 * 
 * Output files:
 * - barcodes.txt: One barcode per line
 * - features.txt: One feature name per line
 * - matrix.mtx: Matrix Market coordinate format (feature x barcode)
 * - stats.txt: Summary statistics
 * - feature_per_cell.csv: Per-barcode feature summary
 * - feature_richness_histogram.html: Histogram of features per barcode
 * - feature_multiplicity_histogram.html: Histogram of UMI counts
 * - Heatmaps via generate_heatmap() and generate_deduped_heatmap()
 * 
 * @param config Writer configuration
 * @param counts Deduped counts result from pf_build_deduped_counts()
 * @return 0 on success, -1 on error
 */
int mex_write_all(
    const mex_writer_config *config,
    pf_counts_result *counts
);

/**
 * @brief Write only the core MEX files (barcodes.txt, features.txt, matrix.mtx)
 * 
 * A lighter-weight version that skips heatmaps and histograms.
 * Useful for testing or when only core outputs are needed.
 * 
 * @param config Writer configuration  
 * @param counts Deduped counts result
 * @return 0 on success, -1 on error
 */
int mex_write_core(
    const mex_writer_config *config,
    pf_counts_result *counts
);

/**
 * @brief Write stats.txt file
 * 
 * @param output_dir Directory to write to
 * @param stats Processing statistics
 * @param hashes Data structures for additional stats
 * @param total_raw_counts Total raw feature counts
 * @param total_deduped_counts Total deduped counts
 * @param total_barcodes Number of barcodes written
 * @param total_excluded_barcodes Number of barcodes excluded by filter
 * @return 0 on success, -1 on error
 */
int mex_write_stats(
    const char *output_dir,
    const statistics *stats,
    const data_structures *hashes,
    int total_raw_counts,
    int total_deduped_counts,
    int total_barcodes,
    int total_excluded_barcodes
);

/**
 * @brief Write feature_per_cell.csv file
 * 
 * @param output_dir Directory to write to
 * @param counts Deduped counts result
 * @param filtered_barcodes_hash Optional filter (or NULL for all barcodes)
 * @return 0 on success, -1 on error
 */
int mex_write_feature_per_cell(
    const char *output_dir,
    pf_counts_result *counts,
    khash_t(strptr) *filtered_barcodes_hash
);

/**
 * @brief Generate all heatmaps for the given counts
 * 
 * Generates:
 * - Co-expression heatmap (features vs co-occurrence)
 * - Deduped counts heatmap (features vs UMI count distribution)
 * - Richness heatmap (features vs feature richness per barcode)
 * 
 * @param output_dir Directory to write heatmaps to
 * @param features Feature arrays
 * @param counts Deduped counts result
 * @param filtered_barcodes_hash Optional filter
 * @param min_heatmap_counts Minimum counts threshold for display
 * @return 0 on success, -1 on error
 */
int mex_write_heatmaps(
    const char *output_dir,
    feature_arrays *features,
    pf_counts_result *counts,
    khash_t(strptr) *filtered_barcodes_hash,
    int min_heatmap_counts
);

#ifdef __cplusplus
}
#endif

#endif /* MEX_WRITER_H */
