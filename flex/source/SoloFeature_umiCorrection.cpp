#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "UmiCodec.h"
#include "UMICorrector.h"
#include "ErrorWarning.h"
#include "hash_shims_cpp_compat.h"
#include "SoloReadFeature.h"
#include <algorithm>
#include <atomic>
#include <cstdio>

// Debug flag - controlled by STAR_DEBUG_TAG environment variable
static const bool g_debugURHist = (std::getenv("STAR_DEBUG_TAG") != nullptr);

void SoloFeature::collectURHistogram(uint32_t /*readId*/, uint32_t /*cbIdx*/, uint32_t /*geneIdx*/) {
    // NOTE: This function uses the SoloTagLedger which is not used in the inline flex path.
    // The flex path uses buildHistogramsFromHash() instead, which extracts UMI directly from the hash.
    // This function is kept as a stub for non-flex paths that may still use the ledger.
    return;
}

void SoloFeature::buildHistogramsFromHash() {
    if (!readFeatSum || !readFeatSum->inlineHash_) {
        return;
    }
    
    // Clear existing histograms (if any)
    umiGroupHistograms.clear();
    
    // Iterate over inline hash to build histograms
    for (khiter_t iter = kh_begin(readFeatSum->inlineHash_); iter != kh_end(readFeatSum->inlineHash_); ++iter) {
        if (!kh_exist(readFeatSum->inlineHash_, iter)) continue;
        
        uint64_t key = kh_key(readFeatSum->inlineHash_, iter);
        uint32_t count = kh_val(readFeatSum->inlineHash_, iter);
        
        // Unpack key to extract (CB, UMI, gene, tag)
        uint32_t cbIdx;
        uint32_t umi24;
        uint16_t geneIdx;
        uint8_t tagIdx;
        unpackCgAggKey(key, &cbIdx, &umi24, &geneIdx, &tagIdx);
        
        // Decode UMI24 to UR string
        // Note: decodeUMI12() is static and doesn't require ledger state
        // It simply decodes a packed 24-bit UMI to a 12bp string
        std::string ur = decodeUMI12(umi24);
        if (ur.empty() || ur.length() != 12) continue;
        
        // Get CB16 string for assignment lookup
        // Fix: cbIdx==0 is valid (first whitelist entry), only check bounds
        if (cbIdx >= pSolo.cbWLstr.size()) continue;
        std::string cb16 = pSolo.cbWLstr[cbIdx].substr(0, 16);
        
        // Check cells allow-list (same logic as collectURHistogram); if empty, allow all
        bool allowed = cellsAllowSet.empty();
        if (!allowed) {
            if (cellsAllowSet.find(cb16) != cellsAllowSet.end()) {
                allowed = true;
            } else if (pSolo.cbWLstr[cbIdx].length() >= 24) {
                std::string cb24 = pSolo.cbWLstr[cbIdx].substr(0, 24);
                if (cellsAllowSet.find(cb24) != cellsAllowSet.end()) {
                    allowed = true;
                }
            }
        }
        
        if (!allowed) continue; // Skip if not in allow-list
        
        // Sample dimension disabled
        uint8_t sampleIdx = 0;
        
        // Build group key: (cbIdx << 24) | (sampleIdx << 16) | geneIdx
        uint64_t groupKey = (static_cast<uint64_t>(cbIdx) << 24) | 
                            (static_cast<uint64_t>(sampleIdx) << 16) | 
                            static_cast<uint64_t>(geneIdx);
        
        // Increment histogram (multiply by count from hash)
        umiGroupHistograms[groupKey].urCounts[ur] += count;
        umiGroupHistograms[groupKey].totalReads += count;
        readsURGrouped += count;
    }
}

void SoloFeature::applyCliqueCorrectionsToHash() {
    if (!readFeatSum || !readFeatSum->inlineHash_ || umiCorrections.empty()) {
        return;
    }
    
    // Build a list of updates (can't modify hash while iterating)
    struct HashUpdate {
        uint64_t oldKey;
        uint64_t newKey;
        uint32_t count;
    };
    std::vector<HashUpdate> updates;
    
    // Iterate over hash to find entries that need correction
    for (khiter_t iter = kh_begin(readFeatSum->inlineHash_); iter != kh_end(readFeatSum->inlineHash_); ++iter) {
        if (!kh_exist(readFeatSum->inlineHash_, iter)) continue;
        
        uint64_t key = kh_key(readFeatSum->inlineHash_, iter);
        uint32_t count = kh_val(readFeatSum->inlineHash_, iter);
        
        // Unpack key
        uint32_t cbIdx;
        uint32_t umi24;
        uint16_t geneIdx;
        uint8_t tagIdx;
        unpackCgAggKey(key, &cbIdx, &umi24, &geneIdx, &tagIdx);
        
        // Decode UMI24 to UR string
        std::string ur = decodeUMI12(umi24);
        if (ur.empty() || ur.length() != 12) continue;
        
        // Sample dimension disabled
        uint8_t sampleIdx = 0;
        
        // Build group key
        uint64_t groupKey = (static_cast<uint64_t>(cbIdx) << 24) | 
                            (static_cast<uint64_t>(sampleIdx) << 16) | 
                            static_cast<uint64_t>(geneIdx);
        
        // Check if this group has corrections
        auto corrIt = umiCorrections.find(groupKey);
        if (corrIt == umiCorrections.end()) continue;
        
        const auto &corrections = corrIt->second;
        auto urCorrIt = corrections.find(ur);
        if (urCorrIt == corrections.end()) continue;
        
        // Found correction: UR -> corrected UB
        std::string correctedUb = urCorrIt->second;
        uint32_t correctedUb24 = encodeUMI12(correctedUb);
        if (correctedUb24 == UINT32_MAX) continue; // Invalid corrected UB
        
        // Pack new key with corrected UB
        uint64_t newKey = packCgAggKey(cbIdx, correctedUb24, geneIdx, tagIdx);
        
        // Store update
        updates.push_back({key, newKey, count});
    }
    
    // Apply updates: delete old keys, insert/update new keys
    for (const auto &update : updates) {
        // Delete old key
        khiter_t old_iter = kh_get(cg_agg, readFeatSum->inlineHash_, update.oldKey);
        if (old_iter != kh_end(readFeatSum->inlineHash_)) {
            kh_del(cg_agg, readFeatSum->inlineHash_, old_iter);
        }
        
        // Insert/update new key
        int absent;
        khiter_t new_iter = kh_put(cg_agg, readFeatSum->inlineHash_, update.newKey, &absent);
        if (absent) {
            kh_val(readFeatSum->inlineHash_, new_iter) = update.count;
        } else {
            kh_val(readFeatSum->inlineHash_, new_iter) += update.count;
        }
    }
}

void SoloFeature::runCliqueCorrection() {
    // Only run if UMI correction is enabled
    if (pSolo.umiCorrectionMode == 0) return;
    
    // Allow clique correction when writing keys.bin (even with skipProcessing)
    // This enables the keys-to-MEX replay workflow where we need corrected keys
    // but don't need full Solo counting/matrix construction
    
    // Debug: log skip conditions
    if (g_debugURHist) {
        fprintf(stderr, "[CLIQUE] Starting correction: umiCorrectionMode=%d writeKeysBin=%d skipProcessing=%d histograms=%zu\n",
                pSolo.umiCorrectionMode, pSolo.writeKeysBin ? 1 : 0, pSolo.skipProcessing ? 1 : 0, umiGroupHistograms.size());
    }
    
    P.inOut->logMain << "Running clique-based UMI correction..." << endl;
    
    // Build histograms from inline hash if enabled
    if (pSolo.inlineHashMode && readFeatSum && readFeatSum->inlineHash_) {
        buildHistogramsFromHash();
    }
    
    // Clear previous corrections
    umiCorrections.clear();
    
    // Correction parameters
    UMIParams params(pSolo.umiMinCount, pSolo.umiRatioThresh, pSolo.maxComponentSize);
    
    // Process each group
    for (auto& groupEntry : umiGroupHistograms) {
        uint64_t groupKey = groupEntry.first;
        UMIHistogram& hist = groupEntry.second;
        
        if (hist.urCounts.empty()) continue;
        
        // Convert histogram to UMICount vector
        std::vector<UMICount> counts;
        counts.reserve(hist.urCounts.size());
        uint64_t groupUMIsBefore = 0;
        uint64_t groupReadsBefore = 0;
        for (const auto& urEntry : hist.urCounts) {
            counts.push_back(UMICount(urEntry.first, urEntry.second));
            groupUMIsBefore++;
            groupReadsBefore += urEntry.second;
        }
        
        // Accumulate before totals
        umisBeforeTotal += groupUMIsBefore;
        readsBeforeTotal += groupReadsBefore;
        
        // Run clique correction
        UMICorrectionResult result = UMICorrector::correctClique(counts, params);
        
        // Accumulate after totals (count unique corrected UBs)
        uint64_t groupUMIsAfter = groupUMIsBefore - result.merges;
        uint64_t groupReadsAfter = groupReadsBefore; // Reads don't change, only UMIs merge
        umisAfterTotal += groupUMIsAfter;
        readsAfterTotal += groupReadsAfter;
        
        // Accumulate metrics
        mergesTotal += result.merges;
        componentsTotal += result.components;
        componentsCappedTotal += result.componentsCapped;
        componentsBelowThresholdTotal += result.componentsBelowThreshold;
        
        // Store corrections for this group
        umiCorrections[groupKey] = result.urToUb;
        
        // Track component size histogram properly
        for (uint32_t compSize : result.componentSizes) {
            if (compSize == 1) componentSizeHist[0]++;
            else if (compSize == 2) componentSizeHist[1]++;
            else if (compSize == 3) componentSizeHist[2]++;
            else if (compSize == 4) componentSizeHist[3]++;
            else if (compSize > 4) componentSizeHist[4]++;
            
            if (compSize > maxComponentSeen) {
                maxComponentSeen = compSize;
            }
        }
        
        // Clear histogram for this group to free memory
        hist.urCounts.clear();
    }
    
    // Apply corrections to reads
    // This will be done in collapseUMIperCB() by looking up corrections
    
    P.inOut->logMain << "Clique correction complete: merges=" << mergesTotal
                     << " components=" << componentsTotal
                     << " umis_before=" << umisBeforeTotal
                     << " umis_after=" << umisAfterTotal
                     << " reads_before=" << readsBeforeTotal
                     << " reads_after=" << readsAfterTotal
                     << endl;
    
    // Clear histograms after processing (corrections stored in umiCorrections)
    umiGroupHistograms.clear();
    
    // Apply corrections to inline hash if enabled
    if (pSolo.inlineHashMode && readFeatSum && readFeatSum->inlineHash_ && !umiCorrections.empty()) {
        applyCliqueCorrectionsToHash();
    }
}
