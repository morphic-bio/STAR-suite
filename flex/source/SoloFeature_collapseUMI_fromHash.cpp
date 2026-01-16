#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "hash_shims_cpp_compat.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"
#include <unordered_map>
#include <unordered_set>
#include "SoloCommon.h"
#include "UmiCodec.h"
#include "SoloReadFeature.h"
#include "UMICorrector.h"
#include "ErrorWarning.h"
#include "solo/CbBayesianResolver.h"
#include "solo/CbCorrector.h"
#include <cstdio>
#include <vector>
#include <algorithm>

// Debug counters from SoloReadFeature_record.cpp
extern "C" uint64_t solo_probe_align_count();
extern "C" uint64_t solo_genomic_align_count();
extern "C" uint64_t solo_probe_resolved_count();
extern "C" uint64_t solo_genomic_resolved_count();
extern "C" uint64_t solo_resolver_dropped_count();
extern "C" uint64_t solo_probe_missing_idx_count();
extern "C" uint64_t solo_resolver_drop_probe_disagree_count();
extern "C" uint64_t solo_resolver_drop_genomic_disagree_count();
extern "C" uint64_t solo_resolver_drop_mixed_count();
extern "C" uint64_t solo_resolver_drop_no_candidates_count();
extern "C" uint64_t solo_resolver_keep_probe_count();
extern "C" uint64_t solo_resolver_keep_genomic_count();
extern "C" uint64_t solo_genomic_align_with_probe_genes_count();
extern "C" uint64_t solo_genomic_only_reads_with_probe_genes_count();
extern "C" uint64_t solo_genomic_only_probe_gene_count();
extern "C" uint64_t solo_genomic_dropped_mapq();
extern "C" uint64_t solo_genomic_dropped_nm();
extern "C" uint64_t solo_genomic_dropped_mmrate();

/**
 * @brief Direct hash consumption for collapse/dedup
 * 
 * This function replaces materializeRGUFromHash + collapseUMIall.
 * It iterates directly over the khash_t(cg_agg) and:
 * - Groups entries by (cbIdx, geneIdx)
 * - Deduplicates UMIs within each group
 * - Applies UMI correction methods (CR, Directional, etc.)
 * - Populates count matrices
 * 
 * Benefits:
 * - No memory expansion (counts stored as-is, not per-read)
 * - No readId dependency
 * - Direct consumption, no materialization step
 */
void SoloFeature::collapseUMIall_fromHash()
{
    if (!readFeatSum || !readFeatSum->inlineHash_) {
        P.inOut->logMain << "ERROR: collapseUMIall_fromHash called but inlineHash_ is null" << endl;
        return;
    }
    
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Starting direct hash collapse (no materialization)" << endl;
    static bool warnedMultiMap = false;
    if (pSolo.multiMap.yes.multi && !warnedMultiMap) {
        P.inOut->logMain << "WARNING: inline-hash direct collapse does not emit multimapper weights; countMatMult will remain empty" << endl;
        warnedMultiMap = true;
    }
    
    // Hash structure: key = [CB20][UMI24][GENE15][TAG5], value = count
    khash_t(cg_agg) *hash = readFeatSum->inlineHash_;
    size_t hashSize = kh_size(hash);
    
    if (hashSize == 0) {
        P.inOut->logMain << "WARNING: inlineHash_ is empty, no reads to collapse" << endl;
        return;
    }

    // ========== Phase 2: Resolve accumulated ambiguous CBs ==========
    if (!readFeatSum->pendingAmbiguous_.empty() && pSolo.cbCorrector) {
        const std::vector<std::string> &whitelistSeqs = pSolo.cbCorrector->whitelist();
        CbBayesianResolver resolver(whitelistSeqs.size(), &whitelistSeqs);
        
        uint64_t resolved = 0, stillAmbiguous = 0, addedToHash = 0;
        
        for (auto &kv : readFeatSum->pendingAmbiguous_) {
            SoloReadFeature::ExtendedAmbiguousEntry &entry = kv.second;
            
            if (entry.candidateIdx.empty() || entry.umiCounts.empty()) {
                stillAmbiguous++;
                continue;
            }
            
            // Build context and candidates for resolver
            CBContext context(entry.cbSeq, entry.cbQual);
            
            std::vector<Candidate> candidates;
            candidates.reserve(entry.candidateIdx.size());
            for (uint32_t idx : entry.candidateIdx) {
                if (idx > 0 && idx <= whitelistSeqs.size()) {
                    candidates.emplace_back(idx, whitelistSeqs[idx - 1], 0.0);
                }
            }
            
            if (candidates.empty()) {
                stillAmbiguous++;
                continue;
            }
            
            // Run Bayesian resolution
            BayesianResult result = resolver.resolve(context, candidates, entry.umiCounts);
            
            if (result.status == BayesianResult::Resolved && result.bestIdx > 0) {
                resolved++;
                
                // Add resolved observations to hash with resolved CB index
                uint32_t resolvedCbIdx = result.bestIdx - 1; // Convert to 0-based
                for (const auto &obs : entry.observations) {
                    uint64_t newKey = packCgAggKey(resolvedCbIdx, obs.umi24, obs.geneIdx, obs.tagIdx);
                    int absent;
                    khiter_t iter = kh_put(cg_agg, hash, newKey, &absent);
                    if (absent) {
                        kh_val(hash, iter) = obs.count;
                    } else {
                        kh_val(hash, iter) += obs.count;
                    }
                    addedToHash++;
                }
            } else {
                stillAmbiguous++;
            }
        }
        
        P.inOut->logMain << "[AMBIG-CB-RESOLVE] pending=" << readFeatSum->pendingAmbiguous_.size()
                         << " resolved=" << resolved
                         << " still_ambiguous=" << stillAmbiguous
                         << " added_to_hash=" << addedToHash << endl;
        
        // Update hash size after resolution
        hashSize = kh_size(hash);
    } else if (!readFeatSum->pendingAmbiguous_.empty()) {
        P.inOut->logMain << "[AMBIG-CB-RESOLVE] " << readFeatSum->pendingAmbiguous_.size() 
                         << " pending ambiguous CBs but CbCorrector not available, skipping" << endl;
    }

    // Instrumentation: total entries and counts pre-dedup
    uint64_t totalCountsPre = 0;
    for (khiter_t iter = kh_begin(hash); iter != kh_end(hash); ++iter) {
        if (!kh_exist(hash, iter)) continue;
        totalCountsPre += kh_val(hash, iter);
    }
    P.inOut->logMain << "[INLINE-HASH] pre_dedup entries=" << hashSize
                     << " total_counts=" << totalCountsPre << endl;

    P.inOut->logMain << "[INLINE-STATS] probe_align=" << solo_probe_align_count()
                     << " genomic_align=" << solo_genomic_align_count()
                     << " probe_resolved=" << solo_probe_resolved_count()
                     << " genomic_resolved=" << solo_genomic_resolved_count()
                     << " resolver_dropped=" << solo_resolver_dropped_count()
                     << " probe_missing_idx=" << solo_probe_missing_idx_count()
                     << " resolver_keep_probe=" << solo_resolver_keep_probe_count()
                     << " resolver_keep_genomic=" << solo_resolver_keep_genomic_count()
                     << " resolver_drop_probe_disagree=" << solo_resolver_drop_probe_disagree_count()
                     << " resolver_drop_genomic_disagree=" << solo_resolver_drop_genomic_disagree_count()
                     << " resolver_drop_mixed=" << solo_resolver_drop_mixed_count()
                     << " resolver_drop_no_candidates=" << solo_resolver_drop_no_candidates_count()
                     << " genomic_align_with_probe_genes=" << solo_genomic_align_with_probe_genes_count()
                     << " genomic_only_reads_with_probe_genes=" << solo_genomic_only_reads_with_probe_genes_count()
                     << " genomic_only_probe_gene_count=" << solo_genomic_only_probe_gene_count()
                     << " genomic_dropped_mapq=" << solo_genomic_dropped_mapq()
                     << " genomic_dropped_nm=" << solo_genomic_dropped_nm()
                     << " genomic_dropped_mmrate=" << solo_genomic_dropped_mmrate()
                     << endl;
    
    // Deduped UMI counts per (CB,TAG,gene) using khash: key = packCgAggKey(cb, 0, gene, tag)
    khash_t(cg_agg)* cbTagGeneCounts = kh_init(cg_agg);

    // Track observed CBs and per-(CB,gene) read totals (for Exact)
    std::unordered_set<uint32_t> uniqueCBs;
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> cbGeneReadCounts;
    
    P.inOut->logMain << "Grouping hash entries by CB/gene/tag..." << endl;
    
    for (khiter_t iter = kh_begin(hash); iter != kh_end(hash); ++iter) {
        if (!kh_exist(hash, iter)) continue;
        
        uint64_t key = kh_key(hash, iter);
        uint32_t count = kh_val(hash, iter);

        // Fast-path: extract tag via mask, skip tagIdx==0, and zero UMI bits for dedup aggregation
        uint8_t tagIdxFast = static_cast<uint8_t>(key & 0x1Fu);
        if (tagIdxFast > 0) {
            static const uint64_t UMI_CLEAR_MASK = ~(0xFFFFFFULL << 20); // zero bits 43..20
            uint64_t tripletKey = key & UMI_CLEAR_MASK; // keeps CB/gene/tag, zeros UMI
            int absent;
            khiter_t titer = kh_put(cg_agg, cbTagGeneCounts, tripletKey, &absent);
            if (absent) {
                kh_val(cbTagGeneCounts, titer) = 1u;
            } else {
                kh_val(cbTagGeneCounts, titer) += 1u;
            }
        }

        uint32_t cbIdx, umi24;
        uint16_t geneIdx;
        uint8_t tagIdx;
        unpackCgAggKey(key, &cbIdx, &umi24, &geneIdx, &tagIdx);
        
        uniqueCBs.insert(cbIdx);
        cbGeneReadCounts[cbIdx][geneIdx] += count;
    }
    
    nCB = uniqueCBs.size();
    P.inOut->logMain << "Found " << nCB << " unique CBs and " << kh_size(cbTagGeneCounts) << " (CB, gene, tag) groups" << endl;
    
    // Step 2: Build indCB (dense list of observed whitelist CB indices) and indCBwl (reverse map WL idx -> dense idx)
    indCB.resize(nCB);                      // dense list of CB whitelist indices observed in this run
    indCBwl.resize(pSolo.cbWLsize, (uint32)-1); // reverse map: whitelist idx -> dense column idx (or -1 if absent)
    
    std::vector<uint32_t> sortedCBs(uniqueCBs.begin(), uniqueCBs.end());
    std::sort(sortedCBs.begin(), sortedCBs.end());
    
    for (uint32_t iCB = 0; iCB < nCB; iCB++) {
        uint32_t cbIdx = sortedCBs[iCB];
        indCB[iCB] = cbIdx;
        indCBwl[cbIdx] = iCB;
    }
    
    // For inline-MEX path we do not build Solo dense matrices; keep minimal bookkeeping only
    nReadPerCB.assign(nCB, 0);
    nReadPerCBunique.assign(nCB, 0);
    nReadPerCBtotal.assign(nCB, 0);
    nUMIperCB.assign(nCB, 0);
    nGenePerCB.assign(nCB, 0);
    countMatStride = pSolo.umiDedup.yes.N + 1;
    countCellGeneUMI.clear();
    countCellGeneUMIindex.assign(nCB + 1, 0);
    if (pSolo.multiMap.yes.multi) {
        countMatMult.s = 1 + pSolo.multiMap.yes.N * pSolo.umiDedup.yes.N;
        countMatMult.m.clear();
        countMatMult.i.assign(nCB + 1, 0);
    }

    // Build per-(CB,gene) dedup counts and per-(CB,TAG) gene counts for MEX
    std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>> cbGeneCounts;
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>> cbTagGeneCountsVec;
    size_t totalCbGeneEntries = 0;
    for (khiter_t iter = kh_begin(cbTagGeneCounts); iter != kh_end(cbTagGeneCounts); ++iter) {
        if (!kh_exist(cbTagGeneCounts, iter)) continue;
        uint64_t tripletKey = kh_key(cbTagGeneCounts, iter);
        uint32_t cbIdx, umiDummy;
        uint16_t geneIdx;
        uint8_t tagIdx;
        unpackCgAggKey(tripletKey, &cbIdx, &umiDummy, &geneIdx, &tagIdx);
        if (tagIdx == 0) continue; // should already be filtered
        uint32_t val = kh_val(cbTagGeneCounts, iter);
        cbGeneCounts[cbIdx][geneIdx] += val;
        uint64_t cbTagKey = (static_cast<uint64_t>(cbIdx) << 8) | tagIdx;
        cbTagGeneCountsVec[cbTagKey].emplace_back(geneIdx, val);
    }
    kh_destroy(cg_agg, cbTagGeneCounts);

    // We skip building Solo dense matrices entirely in this path
    countCellGeneUMI.clear();
    countCellGeneUMIindex.assign(nCB + 1, 0);
    if (pSolo.multiMap.yes.multi) {
        countMatMult.i.assign(nCB + 1, 0);
    }
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished direct hash collapse (MEX-only)" << endl;
    
    P.inOut->logMain << "Found " << cbTagGeneCountsVec.size() << " unique (CB, TAG) combinations" << endl;
    
    // Write MEX format directly from dedup data to the standard Solo raw path
    std::string mexDir = P.outFileNamePrefix + pSolo.outFileNames[0] + SoloFeatureTypes::Names[featureType] + "/raw/";
    createDirectory(mexDir, P.runDirPerm, "Solo raw MEX directory", P);
    // Pass trailing-slash prefix so writer emits matrix.mtx/barcodes.tsv/features.tsv
    InlineMatrixBundle inlineMatrix = buildInlineMatrixFromHash(cbTagGeneCountsVec);
    writeMexFromInlineHashDedup(mexDir, inlineMatrix);
    
    // Run flexfilter inline if enabled.
    // Execute once per STAR run: prefer the Gene feature if present, otherwise the first feature encountered.
    // Note: flexfilter reads from on-disk MEX files (FlexFilter::runFromFiles()), not from the hash.
    static bool flexfilterRan = false;
    const bool geneRequested = pSolo.featureYes[SoloFeatureTypes::Gene];
    const bool shouldRunHere = !geneRequested || featureType == SoloFeatureTypes::Gene;

    if (pSolo.runFlexFilter) {
        
        if (flexfilterRan) {
            P.inOut->logMain << "  Skipping flexfilter for " << SoloFeatureTypes::Names[featureType]
                             << " (already ran for another feature)" << endl;
        } else if (!shouldRunHere) {
            P.inOut->logMain << "  Skipping flexfilter for " << SoloFeatureTypes::Names[featureType]
                             << " (will run once for Gene feature)" << endl;
        } else {
            std::string flexOutputPrefix = pSolo.flexFilterOutputPrefix;
            if (flexOutputPrefix.back() != '/') {
                flexOutputPrefix += '/';
            }
            createDirectory(flexOutputPrefix, P.runDirPerm, "FlexFilter output directory", P);
            flexfilterRan = true;
            runFlexFilterInline(inlineMatrix, flexOutputPrefix);
        }
    }
    
    // Destroy merged hash after MEX write and flexfilter completes
    // 
    // Hash usage timeline:
    // 1. Hash is populated during read capture (SoloReadFeature_record.cpp)
    // 2. Hash is merged from per-thread hashes in sumThreads() (SoloFeature_sumThreads.cpp)
    // 3. Hash is used during collapse (lines 68-206): iterate, dedup UMIs, build cbTagGeneCountsVec
    // 4. MEX is written to disk from cbTagGeneCountsVec (line 206: writeMexFromInlineHashDedup)
    // 5. Flexfilter reads from on-disk MEX files via FlexFilter::runFromFiles() (SoloFeature_flexfilter.cpp:187),
    //    NOT from the hash. Flexfilter uses inputs.matrixDir pointing to the on-disk MEX directory.
    // 6. Hash is no longer needed after MEX write - safe to destroy immediately
    //
    // Note: Each feature has its own hash (readFeatSum is per-feature). If flexfilter will run later
    // for Gene feature (shouldRunHere=false), that's Gene feature's hash, not this feature's hash.
    // Since flexfilter uses on-disk MEX files and not the hash, we can safely destroy this feature's
    // hash regardless of flexfilter timing.
    if (pSolo.soloFlexMinimalMemory && pSolo.inlineHashMode) {
        if (readFeatSum->inlineHash_) {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
        }
    }
}
