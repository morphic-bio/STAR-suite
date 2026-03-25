#include "SoloFeature.h"
#include "hash_shims_cpp_compat.h"
#include "SoloCommon.h"
#include "SoloReadFeature.h"
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <unordered_set>

void SoloFeature::materializeRGUFromHash() {
    // Note: non-Flex STAR_SOLO_NONFLEX_HASH_BRIDGE path no longer calls this from
    // countCBgeneUMI (direct collapse: core/SoloFeature_collapseUMI_fromBridgeHash.cpp).
    if (!readFeatSum || !readFeatSum->inlineHash_) {
        return;
    }
    
    // Estimate size: use hash size * rguStride
    size_t hashSize = kh_size(readFeatSum->inlineHash_);
    if (hashSize == 0) {
        return;
    }
    
    // Build temporary vector of entries for sorting
    struct HashEntry {
        uint32_t cbIdx;
        uint32_t geneIdx;
        uint32_t umi24;
        uint32_t count;
    };
    std::vector<HashEntry> entries;
    entries.reserve(hashSize);
    
    // Extract all entries from hash
    bool nonFlexHashBridge = pSolo.inlineHashMode
        && !pSolo.flexMode
        && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
    for (khiter_t iter = kh_begin(readFeatSum->inlineHash_); iter != kh_end(readFeatSum->inlineHash_); ++iter) {
        if (!kh_exist(readFeatSum->inlineHash_, iter)) continue;
        
        uint64_t key = kh_key(readFeatSum->inlineHash_, iter);
        uint32_t count = kh_val(readFeatSum->inlineHash_, iter);
        
        HashEntry entry;
        uint16_t compactGeneIdx = 0;
        if (nonFlexHashBridge) {
            unpackBridgeCgAggKey(key, &entry.cbIdx, &entry.umi24, &compactGeneIdx);
            entry.cbIdx = readFeatSum->bridgeCompactToWl(entry.cbIdx);
            if (entry.cbIdx == static_cast<uint32_t>(-1)) {
                continue;
            }
            entry.geneIdx = readFeatSum->bridgeCompactToGene(compactGeneIdx);
            if (entry.geneIdx == static_cast<uint32_t>(-1)) {
                continue;
            }
        } else {
            unpackCgAggKey(key, &entry.cbIdx, &entry.umi24, &compactGeneIdx, nullptr);
            entry.geneIdx = compactGeneIdx;
        }
        entry.count = count;
        
        // CRITICAL: Expand by count (each count represents one read)
        // WARNING: This can cause memory explosion for high-count entries (e.g., count=1000 → 1000 entries)
        // TODO: Replace with direct hash consumption in collapseUMIall() to avoid expansion
        // For now, this maintains compatibility with existing collapse logic that expects per-read entries
        for (uint32_t i = 0; i < count; i++) {
            entries.push_back(entry);
        }
    }
    
    // Sort by (CB, gene, UMI)
    std::sort(entries.begin(), entries.end(), 
        [](const HashEntry &a, const HashEntry &b) {
            if (a.cbIdx != b.cbIdx) return a.cbIdx < b.cbIdx;
            if (a.geneIdx != b.geneIdx) return a.geneIdx < b.geneIdx;
            return a.umi24 < b.umi24;
        });
    
    // Count unique CBs and build indCB/indCBwl
    std::unordered_set<uint32_t> uniqueCBs;
    for (const auto &entry : entries) {
        uniqueCBs.insert(entry.cbIdx);
    }
    
    nCB = uniqueCBs.size();
    indCB.resize(nCB);
    indCBwl.resize(pSolo.cbWLsize, (uint32)-1);
    
    // Build indCB (list of detected CB indices) and indCBwl (reverse mapping)
    uint32_t iCB = 0;
    std::vector<uint32_t> sortedCBs(uniqueCBs.begin(), uniqueCBs.end());
    std::sort(sortedCBs.begin(), sortedCBs.end());
    for (uint32_t cbIdx : sortedCBs) {
        indCB[iCB] = cbIdx;
        indCBwl[cbIdx] = iCB;
        iCB++;
    }
    
    // Allocate rGeneUMI array
    size_t totalEntries = entries.size();
    rGeneUMI = new uint32[totalEntries * rguStride];
    
    // Allocate rCBn and rCBp
    rCBn = new uint32[nCB];
    rCBp = new uint32*[nCB];
    
    // Build rGeneUMI array and rCBp/rCBn
    uint32_t rguOffset = 0;
    uint32_t currentCB = UINT32_MAX;
    uint32_t currentICB = 0;
    uint32_t cbStartOffset = 0;
    
    uint32_t syntheticReadId = 0;
    for (size_t i = 0; i < entries.size(); i++) {
        const auto &entry = entries[i];
        
        // Check if CB changed
        if (entry.cbIdx != currentCB) {
            if (currentCB != UINT32_MAX) {
                // Set rCBn for previous CB (number of entries)
                rCBn[currentICB] = (rguOffset - cbStartOffset) / rguStride;
            }
            currentCB = entry.cbIdx;
            currentICB = indCBwl[entry.cbIdx];
            cbStartOffset = rguOffset;
            rCBp[currentICB] = rGeneUMI + rguOffset;
        }
        
        // Write entry to rGeneUMI
        rGeneUMI[rguOffset + rguG] = entry.geneIdx;
        rGeneUMI[rguOffset + rguU] = entry.umi24;
        if (rguStride == 3) {
            // The bridge does not preserve original read ids, but legacy collapse
            // only requires stable unique ids within this materialized record set.
            rGeneUMI[rguOffset + rguR] = syntheticReadId++;
        }
        rguOffset += rguStride;
    }
    
    // Set rCBn for last CB
    if (currentCB != UINT32_MAX) {
        rCBn[currentICB] = (rguOffset - cbStartOffset) / rguStride;
    }
    
    // Set nReadsMapped
    nReadsMapped = totalEntries;
}
