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
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
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
    resolvePendingAmbiguousToHash(false);
    hashSize = kh_size(hash);

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
    
    // --- Flat-sort collapse: extract, sort, count runs ---
    // Each inline hash entry is a unique (CB, UMI, gene, tag) combination.
    // UMI dedup = count entries sharing the same (CB, gene, tag) after zeroing UMI bits.
    auto collapseStart = std::chrono::steady_clock::now();

    static const uint64_t UMI_CLEAR_MASK = ~(0xFFFFFFULL << 20); // zero bits 43..20

    P.inOut->logMain << "Grouping hash entries by CB/gene/tag (flat sort)..." << endl;

    // Phase A: extract triplet keys (UMI-zeroed) and unique CBs in one pass
    std::vector<uint64_t> tripletKeys;
    tripletKeys.reserve(hashSize);
    std::unordered_set<uint32_t> uniqueCBs;

    for (khiter_t iter = kh_begin(hash); iter != kh_end(hash); ++iter) {
        if (!kh_exist(hash, iter)) continue;
        uint64_t key = kh_key(hash, iter);

        uint32_t cbIdx = static_cast<uint32_t>((key >> 44) & 0xFFFFF);
        uniqueCBs.insert(cbIdx);

        uint8_t tagIdx = static_cast<uint8_t>(key & 0x1F);
        if (tagIdx > 0) {
            tripletKeys.push_back(key & UMI_CLEAR_MASK);
        }
    }

    auto extractEnd = std::chrono::steady_clock::now();

    // Phase B: sort triplet keys
    std::sort(tripletKeys.begin(), tripletKeys.end());

    auto sortEnd = std::chrono::steady_clock::now();

    // Phase C: count consecutive runs → build cbTagGeneCountsVec directly
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>> cbTagGeneCountsVec;
    size_t nTripletGroups = 0;

    if (!tripletKeys.empty()) {
        size_t runStart = 0;
        for (size_t i = 1; i <= tripletKeys.size(); ++i) {
            if (i == tripletKeys.size() || tripletKeys[i] != tripletKeys[i-1]) {
                uint64_t tripletKey = tripletKeys[runStart];
                uint32_t dedupCount = static_cast<uint32_t>(i - runStart);

                uint32_t cbIdx;
                uint16_t geneIdx;
                uint8_t tagIdx;
                unpackCgAggKey(tripletKey, &cbIdx, nullptr, &geneIdx, &tagIdx);

                uint64_t cbTagKey = (static_cast<uint64_t>(cbIdx) << 8) | tagIdx;
                cbTagGeneCountsVec[cbTagKey].emplace_back(geneIdx, dedupCount);
                ++nTripletGroups;

                runStart = i;
            }
        }
    }

    { std::vector<uint64_t>().swap(tripletKeys); }

    auto groupEnd = std::chrono::steady_clock::now();
    double extractSec = std::chrono::duration<double>(extractEnd - collapseStart).count();
    double sortSec = std::chrono::duration<double>(sortEnd - extractEnd).count();
    double groupSec = std::chrono::duration<double>(groupEnd - sortEnd).count();
    double totalCollapseSec = std::chrono::duration<double>(groupEnd - collapseStart).count();
    P.inOut->logMain << "Collapse timing: extract=" << std::fixed << std::setprecision(1)
                     << extractSec << "s sort=" << sortSec << "s group=" << groupSec
                     << "s total=" << totalCollapseSec << "s" << endl;

    nCB = uniqueCBs.size();
    P.inOut->logMain << "Found " << nCB << " unique CBs and " << nTripletGroups
                     << " (CB, gene, tag) groups" << endl;

    // Build indCB (dense list of observed whitelist CB indices) and indCBwl (reverse map)
    indCB.resize(nCB);
    indCBwl.resize(pSolo.cbWLsize, (uint32)-1);

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
    // 3. Hash is used during collapse: flat-sort dedup → build cbTagGeneCountsVec
    // 4. MEX is written to disk from cbTagGeneCountsVec (writeMexFromInlineHashDedup)
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
