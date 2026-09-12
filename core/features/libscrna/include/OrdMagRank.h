#ifndef STAR_SUITE_ORDMAG_RANK_H
#define STAR_SUITE_ORDMAG_RANK_H

#include "scrna_types.h"
#include <algorithm>
#include <cmath>
#include <limits>

struct OrdMagCellQuality {
    uint32 detectedGenes = 0;
    uint64_t nonMitoUMIs = 0;
};

// Count distinct features with positive counts, optionally excluding MT rows.
// The original matrix is never modified. Reuse seen across cells and
// use a different token for each cell, avoiding allocation/sorting per cell.
// Initialize seen to UINT32_MAX; cell tokens must be smaller than UINT32_MAX.
inline OrdMagCellQuality ordMagCellQuality(const uint32* geneIds, const uint32* counts,
                                          std::size_t nEntries, vector<uint32>& seen,
                                          uint32 cellToken, const uint8_t* mitochondrial = nullptr,
                                          std::size_t stride = 1)
{
    OrdMagCellQuality quality;
    for (std::size_t i = 0; i < nEntries; ++i) {
        const std::size_t pos = i * stride;
        const uint32 gene = geneIds[pos];
        if (counts[pos] == 0 || gene >= seen.size() || (mitochondrial && mitochondrial[gene])) continue;
        quality.nonMitoUMIs += counts[pos];
        if (seen[gene] != cellToken) {
            seen[gene] = cellToken;
            ++quality.detectedGenes;
        }
    }
    return quality;
}

// The exact target changes one position at a time. Expanding or dropping an
// entire tied group would create a discontinuity proportional to its size.
inline uint32 ordMagRetainCount(uint32 nNonzero, double meanCount)
{
    if (meanCount <= 0.0) return 0;
    if (meanCount >= nNonzero) return nNonzero;
    return static_cast<uint32>(std::round(meanCount));
}

// Quality ranks only break UMI ties. Barcode identity is the last stable key;
// the original index is used only when identities are absent or duplicated.
inline bool ordMagRankBefore(uint32 a, uint32 b,
                             const vector<uint32>& umiCounts,
                             const vector<uint32>* detectedGenes,
                             const vector<string>* barcodes,
                             const vector<uint64_t>* nonMitoUMIs = nullptr)
{
    if (umiCounts[a] != umiCounts[b]) return umiCounts[a] > umiCounts[b];
    if (nonMitoUMIs && (*nonMitoUMIs)[a] != (*nonMitoUMIs)[b]) {
        return (*nonMitoUMIs)[a] > (*nonMitoUMIs)[b];
    }
    if (detectedGenes && (*detectedGenes)[a] != (*detectedGenes)[b]) {
        return (*detectedGenes)[a] > (*detectedGenes)[b];
    }
    if (barcodes && (*barcodes)[a] != (*barcodes)[b]) {
        return (*barcodes)[a] < (*barcodes)[b];
    }
    return a < b;
}

#endif
