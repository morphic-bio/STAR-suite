#ifndef H_FlexFilterIO
#define H_FlexFilterIO

/**
 * @file FlexFilterIO.h
 * @brief Shared I/O helpers for FlexFilter (used by both STAR inline and standalone tool)
 * 
 * Single source of truth for:
 * - Sample whitelist loading (TAB-delimited or tag-only)
 * - Composite MEX loading (barcodes, features, matrix)
 * - Filtered MEX writing (per-sample outputs)
 * - Summary table emission (stdout + TSV)
 */

#include "FlexFilter.h"
#include "../MexWriter.h"
#include <string>
#include <vector>
#include <cstdint>
#include <sys/stat.h>

namespace FlexFilterIO {

//-----------------------------------------------------------------------------
// Data Structures
//-----------------------------------------------------------------------------

/**
 * @struct SampleWhitelist
 * @brief Parsed sample whitelist (labels + tags, sorted by label)
 */
struct SampleWhitelist {
    std::vector<std::string> sampleLabels;  // Sample names (sorted lexicographically)
    std::vector<std::string> sampleTags;    // TAG8 sequences (parallel to labels)
    uint32_t labeledCount;                  // Count of samples with explicit labels
    uint32_t tagOnlyCount;                  // Count of tag-only lines (no label)
};

/**
 * @struct CompositeMexData
 * @brief Loaded composite MEX data (barcodes, features, sparse matrix)
 */
struct CompositeMexData {
    std::vector<std::string> barcodes;                              // Cell barcodes (CB16+TAG8)
    std::vector<MexWriter::Feature> features;                       // Features (id, name, type)
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> cellGeneMap;  // Per-cell (geneIdx, count)
    uint32_t nCells;
    uint32_t nGenes;
    
    CompositeMexData() : nCells(0), nGenes(0) {}
};

/**
 * @struct MexWriteResult
 * @brief Result of writing filtered MEX for a tag
 */
struct MexWriteResult {
    bool success;
    uint32_t cellsWritten;
    uint32_t entriesWritten;
    uint32_t missingBarcodes;  // Barcodes in passingBarcodes but not in MEX
    std::string errorMessage;
    
    MexWriteResult() : success(false), cellsWritten(0), entriesWritten(0), missingBarcodes(0) {}
};

//-----------------------------------------------------------------------------
// Directory Creation Helper
//-----------------------------------------------------------------------------

/**
 * @brief Function pointer type for directory creation
 * Allows STAR to pass its createDirectory wrapper for permission handling
 * @param path Directory path to create
 * @param perm Permission mode (e.g., 0755)
 * @return true on success, false on failure
 */
typedef bool (*CreateDirectoryFunc)(const std::string& path, mode_t perm);

/**
 * @brief Default directory creation using mkdir -p equivalent
 */
bool defaultCreateDirectory(const std::string& path, mode_t perm);

//-----------------------------------------------------------------------------
// Whitelist Loading
//-----------------------------------------------------------------------------

/**
 * @brief Load sample whitelist from TSV file
 * 
 * Format support:
 * - TAB-delimited: sample_name<TAB>tag_sequence
 * - Single-token: tag_sequence only (8bp, auto-labeled as "TAG_<seq>")
 * 
 * @param path Path to whitelist file
 * @param[out] whitelist Parsed whitelist (sorted by label)
 * @param[out] errorMsg Error message on failure
 * @return true on success, false on failure
 * 
 * Logs: "Loaded N samples (M labeled, K tag-only)"
 */
bool loadSampleWhitelist(
    const std::string& path,
    SampleWhitelist& whitelist,
    std::string& errorMsg);

//-----------------------------------------------------------------------------
// MEX Loading
//-----------------------------------------------------------------------------

/**
 * @brief Load composite MEX from directory
 * 
 * Supports:
 * - Standard filenames: matrix.mtx, barcodes.tsv, features.tsv
 * - InlineHashDedup_ prefix fallback
 * 
 * @param mexDir Path to MEX directory
 * @param[out] data Loaded MEX data
 * @param[out] errorMsg Error message on failure
 * @return true on success, false on failure
 * 
 * Warns on malformed feature lines (missing columns)
 */
bool loadCompositeMex(
    const std::string& mexDir,
    CompositeMexData& data,
    std::string& errorMsg);

//-----------------------------------------------------------------------------
// MEX Writing
//-----------------------------------------------------------------------------

/**
 * @brief Write filtered MEX for a single tag/sample
 * 
 * @param outputPrefix Base output prefix (e.g., "/path/to/output/")
 * @param sampleLabel Sample name (used for directory: <prefix>/<sample>/Gene/filtered/)
 * @param passingBarcodes Barcodes that pass filtering
 * @param mexData Preloaded composite MEX data
 * @param createDir Directory creation function (nullptr = use default)
 * @param runDirPerm Permission mode for created directories
 * @param[out] result Write result with statistics
 * @param cb_len Output barcode length. If > 0, barcodes are truncated to this length.
 *               Default: 16 (strip sample tag). Use -1 to keep full CB+TAG barcodes.
 * 
 * Prints: "Writing MEX to: <path>" to stdout
 * Warns: "WARNING: Passing barcode not found in MEX: <bc>"
 * Skips: "Skipping <sample> (no passing barcodes)"
 */
void writeFilteredMexForTag(
    const std::string& outputPrefix,
    const std::string& sampleLabel,
    const std::vector<std::string>& passingBarcodes,
    const CompositeMexData& mexData,
    CreateDirectoryFunc createDir,
    mode_t runDirPerm,
    MexWriteResult& result,
    int cb_len = 16);

//-----------------------------------------------------------------------------
// Summary Emission
//-----------------------------------------------------------------------------

/**
 * @brief Emit summary table to stdout and TSV file
 * 
 * Columns: Sample, Expected, Retain, Simple_ED, Tail_Tested, ED_Pass, Occ_Rem, Final, Total_UMIs
 * Tail_Rescue = ED passers - OrdMag after occupancy
 * 
 * @param outputPrefix Output prefix (writes <prefix>/flexfilter_summary.tsv)
 * @param outputs FlexFilter outputs with per-tag results
 * @param config FlexFilter config (uses emptydropsParams.useFDRGate)
 * 
 * Prints formatted table to stdout with header and TOTAL row
 */
void emitSummary(
    const std::string& outputPrefix,
    const FlexFilter::Outputs& outputs,
    const FlexFilter::Config& config);

//-----------------------------------------------------------------------------
// Utility Functions
//-----------------------------------------------------------------------------

/**
 * @brief Sum UMIs from a matrix.mtx file
 * @param matrixPath Path to matrix.mtx
 * @return Total UMI count, 0 on error
 */
uint64_t sumUMIFromMatrix(const std::string& matrixPath);

/**
 * @brief Trim whitespace from string (in-place)
 */
void trimString(std::string& s);

} // namespace FlexFilterIO

#endif // H_FlexFilterIO

