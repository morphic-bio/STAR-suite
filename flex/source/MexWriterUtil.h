#ifndef CODE_MexWriterUtil
#define CODE_MexWriterUtil

#include "MexWriter.h"
#include <string>
#include <vector>

/**
 * @file MexWriterUtil.h
 * @brief Standalone utility wrapper for calling MexWriter from dedup pipeline
 * 
 * This provides a simple entry point to emit MEX format files from your
 * dedup pipeline, decoupled from Solo and replayer internals.
 * 
 * Usage from dedup pipeline:
 * 
 *   // Prepare your data
 *   std::vector<std::string> composite_barcodes = {"AAACCCAAGAAACACTBC001234", ...};
 *   std::vector<std::string> gene_ids = {"ENSG00000000003", "ENSG00000000005", ...};
 *   std::vector<MexWriter::Triplet> triplets = {
 *       {0, 0, 5},   // cell 0, gene 0, count 5 (0-based)
 *       {1, 2, 12},  // cell 1, gene 2, count 12
 *       ...
 *   };
 *   
 *   // Write MEX
 *   int result = writeMexFromDedup(
 *       "/output/prefix",
 *       composite_barcodes,
 *       gene_ids,
 *       triplets
 *   );
 */

namespace MexWriterUtil {

/**
 * @brief Write MEX format from dedup pipeline data
 * 
 * @param outputPrefix Path prefix for output files (e.g., "/path/to/output")
 *                     Will create: output_matrix.mtx, output_barcodes.tsv, output_features.tsv
 * 
 * @param barcodes Composite barcodes from your dedup pipeline
 *                 Order determines column indices in matrix.mtx (0-based internally, 1-based in output)
 *                 Format: vector of strings (e.g., "AAACCCAAGAAACACTBC001234" for CB16+TAG8)
 *                 You provide the strings and ordering - no synthesis inside STAR
 * 
 * @param geneIds Gene IDs from your dedup pipeline
 *                Order determines row indices in matrix.mtx (0-based internally, 1-based in output)
 *                Format: vector of strings (e.g., "ENSG00000000003")
 * 
 * @param triplets Sparse count matrix entries from your dedup output
 *                 Format: vector of (cell_idx, gene_idx, count) where indices are 0-based
 *                 Will be converted to 1-based for MTX output automatically
 * 
 * @return 0 on success, -1 on error
 * 
 * Notes:
 * - Barcodes and gene IDs must be in the exact order you want in the output
 * - Triplet indices are 0-based (cell_idx < barcodes.size(), gene_idx < geneIds.size())
 * - Output MTX uses 1-based indices (Matrix Market standard)
 * - No dependency on Solo, keys.bin, or replayer structures
 */
int writeMexFromDedup(
    const std::string& outputPrefix,
    const std::vector<std::string>& barcodes,
    const std::vector<std::string>& geneIds,
    const std::vector<MexWriter::Triplet>& triplets
);

/**
 * @brief Overload with explicit gene names and feature types
 * 
 * Use this if you want custom gene names (different from IDs) or feature types.
 */
int writeMexFromDedup(
    const std::string& outputPrefix,
    const std::vector<std::string>& barcodes,
    const std::vector<MexWriter::Feature>& features,
    const std::vector<MexWriter::Triplet>& triplets
);

} // namespace MexWriterUtil

#endif
