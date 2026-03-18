#ifndef PF_MULTI_MERGE_H
#define PF_MULTI_MERGE_H

#include "IncludeDefine.h"
#include "MexWriter.h"
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

string resolveMexFile(const string& mexDir, const string& basename);
vector<string> readLines(const string& path);

MexData readMex(const string& mexDir);

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
