#ifndef STAR_SUITE_SIMPLE_ED_CALLER_H
#define STAR_SUITE_SIMPLE_ED_CALLER_H

#include "scrna_api.h"
#include "OrdMagStage.h"
#include "SparseCountView.h"
#include <vector>
#include <string>

// Shared orchestration used by the standalone CLI and in-process Flex caller.
// All sparse arrays and barcode IDs describe the same unmodified matrix.
struct SimpleEDOptions {
    bool legacyRankAmbient = true;
    bool guardedRankAmbient = false;
    uint32_t ambientFallbackMinAbs = 0;
    double ambientFallbackMinFrac = 0.0;
    uint32_t maxExpectedCells = 0;
    uint32_t ordmagRetainCount = 0;
    uint64_t ambientUmiTarget = 0;
    uint32_t bootstrapThreads = 0;
    uint32_t bootstrapWorkers = 0; // Execution cap; bootstrapThreads retains RNG partitioning
    bool invariantChecks = false;
    std::string diagnosticsDir;
};

struct SimpleEDRunInfo {
    uint32_t retainCount = 0;
    uint32_t ordmagCount = 0;
    uint32_t ambientCells = 0;
    OrdMagBootstrapTrace bootstrap;
};

int runSimpleEDWithAmbient(const std::vector<std::string>& barcodes,
                           const std::vector<uint32_t>& umi_counts,
                           const std::vector<uint32_t>& sparse_gene_ids,
                           const std::vector<uint32_t>& sparse_counts,
                           const std::vector<uint32_t>& sparse_cell_index,
                           const std::vector<uint32_t>& n_genes_per_cell,
                           uint32_t n_features,
                           const scrna_ed_config* config,
                           const SimpleEDOptions& options,
                           const std::vector<uint8_t>& mitochondrial_features,
                           scrna_ed_result* result,
                           SimpleEDRunInfo* info = nullptr);

// Borrow the sparse counts for the synchronous call; sample offsets may refer
// to disjoint rows in a shared matrix. Statistics and row identities are unchanged.
int runSimpleEDWithAmbientView(const std::vector<std::string>& barcodes,
                           const std::vector<uint32_t>& umi_counts,
                           const SparseCountView& matrix,
                           uint32_t n_features,
                           const scrna_ed_config* config,
                           const SimpleEDOptions& options,
                           const std::vector<uint8_t>& mitochondrial_features,
                           scrna_ed_result* result,
                           SimpleEDRunInfo* info = nullptr);


#endif
