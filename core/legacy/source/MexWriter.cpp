#include "MexWriter.h"
#include <cstdio>
#include <algorithm>
#include <unordered_set>

namespace MexWriter {

int writeMex(const std::string& outputPrefix,
             const std::vector<std::string>& barcodes,
             const std::vector<Feature>& features,
             const std::vector<Triplet>& triplets,
             int cb_len)
{
    if (barcodes.empty() || features.empty()) {
        std::fprintf(stderr, "[MexWriter] ERROR: empty input (barcodes=%zu, features=%zu)\n",
                     barcodes.size(), features.size());
        return -1; // No data to write
    }
    
    // Normalize barcodes if cb_len > 0 (truncate to cb_len characters)
    std::vector<std::string> outputBarcodes;
    if (cb_len > 0) {
        outputBarcodes.reserve(barcodes.size());
        std::unordered_set<std::string> seen;
        seen.reserve(barcodes.size());
        
        for (size_t i = 0; i < barcodes.size(); i++) {
            const std::string& bc = barcodes[i];
            std::string truncated;
            
            if (static_cast<int>(bc.size()) <= cb_len) {
                // Barcode shorter than or equal to cb_len: use as-is
                truncated = bc;
            } else {
                // Truncate to cb_len
                truncated = bc.substr(0, static_cast<size_t>(cb_len));
            }
            
            // Check for duplicates after truncation
            if (seen.find(truncated) != seen.end()) {
                std::fprintf(stderr, "[MexWriter] ERROR: duplicate barcode after truncation to %d bp\n", cb_len);
                std::fprintf(stderr, "  barcode[%zu] = '%s' -> '%s' (already exists)\n", 
                            i, bc.c_str(), truncated.c_str());
                std::fprintf(stderr, "  This should not happen in per-sample MEX output.\n");
                std::fprintf(stderr, "  SOLUTION: Check for barcode collisions or disable truncation:\n");
                std::fprintf(stderr, "    - External tool: --keep-cb-tag\n");
                std::fprintf(stderr, "    - STAR inline:   --soloFlexKeepCBTag yes\n");
                return -1;
            }
            
            seen.insert(truncated);
            outputBarcodes.push_back(truncated);
        }
    } else {
        // No truncation: use original barcodes
        outputBarcodes = barcodes;
    }
    
    // Create output file paths
    bool dirMode = (!outputPrefix.empty() && outputPrefix.back() == '/');
    std::string base = outputPrefix;
    if (dirMode && base.size() > 1 && base.back() == '/') {
        // leave trailing slash so we emit files inside the directory
    }
    std::string mtx_path = dirMode ? (base + "matrix.mtx") : (outputPrefix + "_matrix.mtx");
    std::string barcodes_path = dirMode ? (base + "barcodes.tsv") : (outputPrefix + "_barcodes.tsv");
    std::string features_path = dirMode ? (base + "features.tsv") : (outputPrefix + "_features.tsv");
    
    FILE *mtx_fp = fopen(mtx_path.c_str(), "w");
    FILE *barcodes_fp = fopen(barcodes_path.c_str(), "w");
    FILE *features_fp = fopen(features_path.c_str(), "w");
    
    if (!mtx_fp || !barcodes_fp || !features_fp) {
        std::fprintf(stderr, "[MexWriter] ERROR: failed to open output files\n");
        if (!mtx_fp) std::fprintf(stderr, "  failed: %s\n", mtx_path.c_str());
        if (!barcodes_fp) std::fprintf(stderr, "  failed: %s\n", barcodes_path.c_str());
        if (!features_fp) std::fprintf(stderr, "  failed: %s\n", features_path.c_str());
        if (mtx_fp) fclose(mtx_fp);
        if (barcodes_fp) fclose(barcodes_fp);
        if (features_fp) fclose(features_fp);
        return -1;
    }
    
    uint32_t num_features = static_cast<uint32_t>(features.size());
    uint32_t num_barcodes = static_cast<uint32_t>(outputBarcodes.size());
    uint64_t num_nonzero = static_cast<uint64_t>(triplets.size());
    
    // Write Matrix Market header
    fprintf(mtx_fp, "%%%%MatrixMarket matrix coordinate integer general\n");
    fprintf(mtx_fp, "%%\n");
    fprintf(mtx_fp, "%u %u %lu\n", num_features, num_barcodes, num_nonzero);
    
    // Write matrix entries (row, col, count) - Matrix Market uses 1-based indexing
    for (const auto& t : triplets) {
        if (t.cell_idx < outputBarcodes.size() && t.gene_idx < features.size()) {
            // Convert 0-based to 1-based
            fprintf(mtx_fp, "%u %u %u\n", 
                   t.gene_idx + 1,  // row (feature)
                   t.cell_idx + 1,  // col (barcode)
                   t.count);
        }
    }
    
    // Write barcodes (ordered by column) - using normalized barcodes
    for (const auto& bc : outputBarcodes) {
        fprintf(barcodes_fp, "%s\n", bc.c_str());
    }
    
    // Write features (ordered by row)
    for (const auto& feat : features) {
        fprintf(features_fp, "%s\t%s\t%s\n", 
               feat.id.c_str(), 
               feat.name.c_str(),
               feat.featureType.c_str());
    }
    
    fclose(mtx_fp);
    fclose(barcodes_fp);
    fclose(features_fp);
    
    return 0;
}

int writeMex(const std::string& outputPrefix,
             const std::vector<std::string>& barcodes,
             const std::vector<std::string>& featureIds,
             const std::vector<Triplet>& triplets,
             int cb_len)
{
    // Convert simple feature IDs to Feature structs
    std::vector<Feature> features;
    features.reserve(featureIds.size());
    for (const auto& id : featureIds) {
        features.emplace_back(id, id, "Gene Expression");
    }
    
    return writeMex(outputPrefix, barcodes, features, triplets, cb_len);
}

} // namespace MexWriter
