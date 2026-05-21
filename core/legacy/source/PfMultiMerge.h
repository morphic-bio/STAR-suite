#ifndef PF_MULTI_MERGE_H
#define PF_MULTI_MERGE_H

#include "IncludeDefine.h"
#include "MexWriter.h"
#include "Parameters.h"
#include <vector>
#include <fstream>
#include <chrono>
#include <iostream>

/**
 * @file PfMultiMerge.h
 * @brief Merge multiple MEX files into a combined MEX
 * 
 * Reads GEX MEX and feature MEX files, merges by barcode,
 * and writes combined output.
 */

namespace PfMultiMerge {

struct TimerScope {
    const char* label;
    std::ostream& log;
    std::chrono::steady_clock::time_point t0;
    TimerScope(const char* l, std::ostream& lg)
        : label(l), log(lg), t0(std::chrono::steady_clock::now()) {}
    ~TimerScope() {
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - t0).count();
        log << "pf-timing: " << label << " " << ms << " ms" << std::endl;
    }
};

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
};

/** Feature/barcode axes only (no matrix triplets). */
struct MexAxes {
    vector<string> features;
    vector<string> featureNames;
    vector<string> featureTypes;
    vector<string> barcodes;
};

/** Cell Ranger–compatible sorted/suffixed barcode column layout for MEX output. */
struct CrBarcodeLayout {
    vector<string> sortedBarcodes;
    /** sourceColToSorted[source_col] = output column (0-based), or UINT32_MAX if dropped. */
    vector<uint32_t> sourceColToSorted;
};

/**
 * Apply GEM-well suffixing, optional NXT/TRU namespace translation, duplicate checks,
 * and lexicographic sort (same rules as writeCombinedMex).
 */
CrBarcodeLayout buildCrBarcodeLayoutForColumns(const vector<string>& sourceBarcodes,
                                               const vector<uint32_t>& colIndices,
                                               const string& gemWell,
                                               const string& inputChemistry,
                                               const string& outputChemistry,
                                               ostream& logStream);

string resolveMexFile(const string& mexDir, const string& basename);
vector<string> readLines(const string& path);

MexData readMex(const string& mexDir);

/** Load features + barcodes without reading matrix.mtx into memory. */
MexAxes readMexAxes(const string& mexDir);

/**
 * @brief Build old_col -> new_col remap; vector length = source column count.
 * Only indices listed in colIndices receive dense 0..N-1 targets.
 */
vector<uint32_t> buildColumnRemap(const vector<uint32_t>& colIndices, size_t sourceColumnCount);

/**
 * @brief Stream column-subset of matrix.mtx(.gz) to output matrix.mtx.gz (two-pass for nnz).
 */
uint64_t streamMatrixColumnSubset(const string& inputMatrixPath,
                                  const string& outputMatrixGzPath,
                                  const vector<uint32_t>& oldToNew,
                                  size_t nRows,
                                  size_t nColsOut);

/**
 * @brief Write features/barcodes.gz + streamed matrix subset with CR-compat barcodes.
 */
int writeColumnSubsetMexGz(const string& inputMexDir,
                           const string& outputDir,
                           const MexAxes& sourceAxes,
                           const CrBarcodeLayout& layout,
                           Parameters& P,
                           ostream& logStream);

/**
 * @brief Stream full pool MEX with full raw barcode axis + CR-compat barcodes (multi/count/raw).
 */
int writeStreamedPoolMexGzCrCompat(const string& inputMexDir,
                                   const string& outputDir,
                                   const MexAxes& sourceAxes,
                                   Parameters& P,
                                   ostream& logStream,
                                   const string& gemWell = "1",
                                   const vector<string>& gexBarcodeFilter = vector<string>(),
                                   const string& inputChemistry = "TRU",
                                   const string& outputChemistry = "TRU");

/** Copy pre-formatted gzip MEX files (downstream mirrors; barcodes already CR-compat). */
int copyMexGzDir(const string& inputMexDir, const string& outputDir, Parameters& P);

MexData filterByFeatureType(const MexData& data, const string& featureType);

MexData mergeMex(const MexData& gexData, const vector<MexData>& featureDataVec);

vector<string> computeObservedGexBarcodes(const MexData& gexData);

void pruneZeroCountFeatures(MexData& data, ofstream& logStream);

/** @brief Silent overload for pre-merge pruning. Returns count of pruned features. */
size_t pruneZeroCountFeatures(MexData& data);

/**
 * @brief Write combined MEX to directory with streaming gzip.
 *
 * Uses vector-based O(1) barcode remap and writes directly to gzFile
 * with batched snprintf, eliminating the write-plain-then-compress path.
 */
int writeCombinedMex(const string& outputDir,
                     const MexData& data,
                     const string& gemWell,
                     ofstream& logStream,
                     const vector<string>& gexBarcodes = vector<string>(),
                     const string& inputChemistry = "TRU",
                     const string& outputChemistry = "TRU");

} // namespace PfMultiMerge

#endif // PF_MULTI_MERGE_H
