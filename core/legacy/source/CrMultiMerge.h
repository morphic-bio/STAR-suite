#ifndef CR_MULTI_MERGE_H
#define CR_MULTI_MERGE_H

#include "IncludeDefine.h"
#include "MexWriter.h"
#include <vector>
#include <fstream>

/**
 * @file CrMultiMerge.h
 * @brief Merge multiple MEX files into a combined MEX
 * 
 * Reads GEX MEX and feature MEX files, merges by barcode,
 * and writes combined output.
 */

namespace CrMultiMerge {

/**
 * @struct MexData
 * @brief In-memory representation of a MEX file
 */
struct MexData {
    vector<string> features;    // Feature IDs (from features.tsv)
    vector<string> featureNames; // Feature names
    vector<string> featureTypes; // Feature types
    vector<string> barcodes;    // Barcode list (from barcodes.tsv)
    vector<MexWriter::Triplet> triplets; // Sparse matrix entries (row, col, count)
    // Note: triplets use 0-based indices
};

/**
 * @brief Read MEX files from directory
 * @param mexDir Directory containing matrix.mtx, features.tsv, barcodes.tsv
 * @return MexData structure
 * @throws runtime_error on read errors
 */
MexData readMex(const string& mexDir);

/**
 * @brief Filter MEX data to specific feature type
 * @param data Input MEX data
 * @param featureType Feature type to keep (e.g., "Gene Expression")
 * @return Filtered MexData
 */
MexData filterByFeatureType(const MexData& data, const string& featureType);

/**
 * @brief Merge multiple MEX files
 * @param gexData GEX MEX data (determines barcode order)
 * @param featureDataVec Vector of feature MEX data to merge
 * @return Combined MexData
 */
MexData mergeMex(const MexData& gexData, const vector<MexData>& featureDataVec);

/**
 * @brief Compute observed GEX barcodes (barcodes with counts > 0 in GEX triplets)
 * @param gexData GEX MEX data
 * @return Vector of observed barcode strings
 */
vector<string> computeObservedGexBarcodes(const MexData& gexData);

/**
 * @brief Write combined MEX to directory
 * @param outputDir Output directory (will be created)
 * @param data Combined MEX data
 * @param gemWell GEM well suffix (e.g., "1", "2") to append to barcodes
 * @param logStream Log stream for output messages
 * @param gexBarcodes Optional GEX barcode whitelist (if non-empty, filter to GEX-only barcodes)
 * @return 0 on success, -1 on error
 */
int writeCombinedMex(const string& outputDir, const MexData& data, const string& gemWell, ofstream& logStream, const vector<string>& gexBarcodes = vector<string>());

} // namespace CrMultiMerge

#endif // CR_MULTI_MERGE_H
