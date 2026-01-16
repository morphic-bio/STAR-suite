#ifndef CODE_MexWriter
#define CODE_MexWriter

#include <string>
#include <vector>
#include <cstdint>

/**
 * @file MexWriter.h
 * @brief Minimal MEX (Matrix Exchange Format) writer for single-cell count matrices
 * 
 * This module provides a lightweight utility to write sparse count matrices in the
 * Matrix Market COO format, along with barcodes and features TSV files.
 * 
 * Source: Extracted from SoloKeyReplay::writeMatrixMarket() for reuse across
 * multiple deduplication pipelines (inline-hash, keys-replayer).
 * 
 * Format specification:
 * - matrix.mtx: Matrix Market coordinate format (COO), 1-based indices
 * - barcodes.tsv: One barcode per line, ordered by column index
 * - features.tsv: Tab-separated (ID, name, type), ordered by row index
 * 
 * Assumptions:
 * - CB ordering: Provided barcode list determines column order
 * - Gene ordering: Provided feature list determines row order
 * - Triplets: (cb_idx, gene_idx, count) are 0-based internally, converted to 1-based for MTX
 */

namespace MexWriter {

/**
 * @struct Triplet
 * @brief Sparse matrix entry: (cell_idx, gene_idx, count)
 * Indices are 0-based (internal representation)
 */
struct Triplet {
    uint32_t cell_idx;  // 0-based cell/barcode index
    uint32_t gene_idx;  // 0-based gene/feature index
    uint32_t count;     // UMI count or read count
};

/**
 * @struct Feature
 * @brief Gene/feature metadata for features.tsv
 */
struct Feature {
    std::string id;         // Gene ID or feature ID
    std::string name;       // Gene name (can be same as ID)
    std::string featureType; // e.g., "Gene Expression", "Multiplexing Capture"
    
    Feature(const std::string& id_, const std::string& name_, const std::string& type_ = "Gene Expression")
        : id(id_), name(name_), featureType(type_) {}
};

/**
 * @brief Write MEX format files (matrix.mtx, barcodes.tsv, features.tsv)
 * 
 * @param outputPrefix Path prefix for output files (e.g., "/path/to/output" → output_matrix.mtx)
 * @param barcodes Ordered list of cell barcodes (determines column order, 0-based)
 * @param features Ordered list of features (determines row order, 0-based)
 * @param triplets Sparse matrix entries (cell_idx, gene_idx, count), 0-based
 * @param cb_len Output barcode length. If > 0, barcodes are truncated to this length.
 *               If -1 (default), no truncation is applied.
 *               For per-sample MEX outputs, use cb_len=16 to strip sample tags.
 * @return 0 on success, -1 on error (including duplicate barcodes after truncation)
 * 
 * Example usage:
 * ```cpp
 * std::vector<std::string> barcodes = {"AAACCCAAGAAACACT-1", "AAACCCAAGAAACCAT-1"};
 * std::vector<MexWriter::Feature> features = {
 *     MexWriter::Feature("ENSG00000000003", "TSPAN6"),
 *     MexWriter::Feature("ENSG00000000005", "TNMD")
 * };
 * std::vector<MexWriter::Triplet> triplets = {
 *     {0, 0, 5},  // cell 0, gene 0, count 5
 *     {1, 1, 3}   // cell 1, gene 1, count 3
 * };
 * // Per-sample output with 16bp barcodes (tag stripped)
 * MexWriter::writeMex("/path/to/output/", barcodes, features, triplets, 16);
 * ```
 */
int writeMex(const std::string& outputPrefix,
             const std::vector<std::string>& barcodes,
             const std::vector<Feature>& features,
             const std::vector<Triplet>& triplets,
             int cb_len = -1);

/**
 * @brief Convenience overload: features as simple ID strings (name=ID, type="Gene Expression")
 * @param cb_len Output barcode length. If > 0, barcodes are truncated to this length.
 */
int writeMex(const std::string& outputPrefix,
             const std::vector<std::string>& barcodes,
             const std::vector<std::string>& featureIds,
             const std::vector<Triplet>& triplets,
             int cb_len = -1);

} // namespace MexWriter

#endif
