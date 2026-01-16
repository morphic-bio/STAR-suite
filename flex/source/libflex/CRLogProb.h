#ifndef H_CRLogProb
#define H_CRLogProb

#include <vector>
#include <cmath>

// DropletUtils-style dense multinomial log-probability:
// logLik = logFac(total) + sum(count_i * log(p_i)) - sum(logFac(count_i))
// Inputs:
//   counts: pointer to counts for this cell (sparse triplets not assumed here)
//   geneIds: pointer to gene IDs for this cell (same length as counts)
//   nGenes: number of genes with nonzero counts for this cell
//   ambientLogP: dense log probabilities for all genes (length = nFeatures)
// Returns: log-likelihood as double
inline double computeDenseMultinomialLogProb(
    const std::vector<uint32_t>& counts,
    const std::vector<uint32_t>& geneIds,
    const std::vector<double>& ambientLogP
) {
    uint32_t total = 0;
    double sumLogFac = 0.0;
    double sumCountLogP = 0.0;
    for (size_t i = 0; i < counts.size(); i++) {
        uint32_t c = counts[i];
        total += c;
        sumLogFac += lgamma(static_cast<double>(c) + 1.0);
        uint32_t gid = geneIds[i];
        if (gid < ambientLogP.size()) {
            sumCountLogP += static_cast<double>(c) * ambientLogP[gid];
        }
    }
    double logFacTotal = lgamma(static_cast<double>(total) + 1.0);
    return logFacTotal - sumLogFac + sumCountLogP;
}

#endif
