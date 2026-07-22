#ifndef STAR_MULTI_GENE_UMI_CR_H
#define STAR_MULTI_GENE_UMI_CR_H

#include <cstdint>
#include <string>
#include <vector>

namespace multi_gene_umi_cr {

struct GeneSupport {
    std::uint32_t gene = 0;
    std::uint64_t correctedCount = 0;
    std::uint64_t originalAtCorrectedCount = 0;
    GeneSupport() = default;
    GeneSupport(std::uint32_t geneIn, std::uint64_t correctedCountIn,
                std::uint64_t originalAtCorrectedCountIn)
        : gene(geneIn), correctedCount(correctedCountIn),
          originalAtCorrectedCount(originalAtCorrectedCountIn) {}
};

struct Result {
    bool accepted = false;
    std::uint32_t gene = UINT32_MAX;
    std::string reason;
};

// Resolve one candidate/barcode + corrected-UMI group. The input must contain
// one row per provisional gene after per-gene UMI correction.
Result resolve(const std::vector<GeneSupport> &supports);

} // namespace multi_gene_umi_cr

#endif
