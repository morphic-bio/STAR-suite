#include "MultiGeneUmiCr.h"

#include <set>
#include <stdexcept>

namespace multi_gene_umi_cr {

Result resolve(const std::vector<GeneSupport> &supports)
{
    if (supports.empty()) throw std::invalid_argument("MultiGeneUMI_CR support is empty");
    std::set<std::uint32_t> genes;
    std::uint64_t maximum = 0;
    std::uint32_t winner = UINT32_MAX;
    bool tied = false;
    for (const GeneSupport &support : supports) {
        if (!genes.insert(support.gene).second) {
            throw std::invalid_argument("duplicate gene in MultiGeneUMI_CR support");
        }
        if (support.correctedCount == 0
            || support.originalAtCorrectedCount > support.correctedCount) {
            throw std::invalid_argument("invalid MultiGeneUMI_CR count tuple");
        }
        if (support.correctedCount > maximum) {
            maximum = support.correctedCount;
            winner = support.gene;
            tied = false;
        } else if (support.correctedCount == maximum) {
            tied = true;
        }
    }
    Result result;
    if (tied || winner == UINT32_MAX) {
        result.reason = "corrected_count_tie";
        return result;
    }
    std::uint64_t winnerOriginal = 0;
    for (const GeneSupport &support : supports) {
        if (support.gene == winner) winnerOriginal = support.originalAtCorrectedCount;
    }
    for (const GeneSupport &support : supports) {
        if (support.originalAtCorrectedCount > winnerOriginal) {
            result.reason = "original_umi_dominance";
            return result;
        }
    }
    result.accepted = true;
    result.gene = winner;
    result.reason = "unique_corrected_count_winner";
    return result;
}

} // namespace multi_gene_umi_cr
