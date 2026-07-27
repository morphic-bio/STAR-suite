#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "hash_shims_cpp_compat.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "serviceFuns.cpp"
#include "SoloCommon.h"
#include "UmiCodec.h"
#include "UMICorrector.h"
#include "ErrorWarning.h"
#include "MexWriter.h"
#include "SampleDetector.h"
#include "FlexGdna.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <omp.h>

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

    const char *flexSnapIn = std::getenv("STAR_SOLO_FLEX_HASH_SNAPSHOT_IN");
    const bool flexSnapReplay = (flexSnapIn != nullptr && flexSnapIn[0] != '\0');

    // ========== Phase 2: Resolve accumulated ambiguous CBs ==========
    // Must run BEFORE the empty-hash early return: in no-align mode all reads
    // may still be in pendingAmbiguous_ and only land in inlineHash_ after
    // resolution.
    if (!flexSnapReplay) {
        resolvePendingAmbiguousToHash(false);
    }
    size_t hashSize = kh_size(hash);

    if (!flexSnapReplay) {
        if (const char *snapPath = std::getenv("STAR_SOLO_FLEX_HASH_SNAPSHOT_OUT")) {
            if (snapPath[0] != '\0')
                flexHashSnapshotWrite(snapPath);
        }
    }

    if (hashSize == 0) {
        P.inOut->logMain << "WARNING: inlineHash_ is empty after ambiguous resolution, no reads to collapse" << endl;
        return;
    }

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
    
    // --- Collapse: CB-first CSR layout + parallel per-CB local dedup ---
    // Solo-style architecture: flat CB counter → CSR offsets → scatter into per-CB slices →
    // parallel per-CB sort + serial run-length scan. No second hash table, no global sort.
    // Each CB slice (~323 entries avg) fits in L1; sort + scan is microseconds per CB.
    auto collapseStart = std::chrono::steady_clock::now();

    P.inOut->logMain << "Grouping hash entries by CB (CSR scatter + per-CB local dedup)..." << endl;

    // Phase 1: Count entries per CB + total counts (fused — single hash scan for both)
    const uint32_t wlSize = pSolo.cbWLsize;
    std::vector<size_t> cbCounts(wlSize, 0);
    uint64_t totalCountsPre = 0;
    for (khiter_t iter = kh_begin(hash); iter != kh_end(hash); ++iter) {
        if (!kh_exist(hash, iter)) continue;
        uint32_t cbIdx = static_cast<uint32_t>((kh_key(hash, iter) >> 44) & 0xFFFFF);
        if (cbIdx < wlSize)
            ++cbCounts[cbIdx];
        totalCountsPre += flexGdnaValueCount(kh_val(hash, iter));
    }
    P.inOut->logMain << "[INLINE-HASH] pre_dedup entries=" << hashSize
                     << " total_counts=" << totalCountsPre << endl;

    auto cbCountEnd = std::chrono::steady_clock::now();

    // Build sorted CB list (sequential scan → already sorted) + indCB / indCBwl
    std::vector<uint32_t> sortedCBs;
    sortedCBs.reserve(wlSize / 8);
    for (uint32_t i = 0; i < wlSize; ++i) {
        if (cbCounts[i] > 0)
            sortedCBs.push_back(i);
    }

    nCB = static_cast<uint32_t>(sortedCBs.size());
    indCB.resize(nCB);
    indCBwl.assign(wlSize, static_cast<uint32_t>(-1));
    for (uint32_t i = 0; i < nCB; ++i) {
        indCB[i] = sortedCBs[i];
        indCBwl[sortedCBs[i]] = i;
    }

    // CSR offsets (prefix sum over per-CB counts)
    std::vector<size_t> cbOffsets(nCB + 1, 0);
    for (uint32_t iCB = 0; iCB < nCB; ++iCB)
        cbOffsets[iCB + 1] = cbOffsets[iCB] + cbCounts[indCB[iCB]];
    { std::vector<size_t>().swap(cbCounts); }

    const size_t totalEntries = cbOffsets[nCB];

    // Phase 2: Scatter hash entries into flat CSR array
    // Each slot packs [tag(5)][gene(15)][umi(24)] into uint64 for per-CB sorting.
    std::vector<uint64_t> flat(totalEntries);
    {
        std::vector<size_t> cbWrite(nCB);
        for (uint32_t iCB = 0; iCB < nCB; ++iCB)
            cbWrite[iCB] = cbOffsets[iCB];

        for (khiter_t iter = kh_begin(hash); iter != kh_end(hash); ++iter) {
            if (!kh_exist(hash, iter)) continue;
            uint64_t key = kh_key(hash, iter);
            uint32_t cbIdx = static_cast<uint32_t>((key >> 44) & 0xFFFFF);
            if (cbIdx >= wlSize) continue;
            uint32_t iCB = indCBwl[cbIdx];
            uint8_t  tag   = static_cast<uint8_t>(key & 0x1F);
            uint16_t gene  = static_cast<uint16_t>((key >> 5) & 0x7FFF);
            uint32_t umi24 = static_cast<uint32_t>((key >> 20) & 0xFFFFFF);
            flat[cbWrite[iCB]++] = ((uint64_t)tag << 39) | ((uint64_t)gene << 24) | umi24;
        }
    }

    auto scatterEnd = std::chrono::steady_clock::now();

    // Phase 3: Parallel per-CB sort (avg ~323 entries/CB → L1-resident, microseconds each)
#pragma omp parallel for schedule(dynamic, 64)
    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        size_t beg = cbOffsets[iCB];
        size_t end = cbOffsets[iCB + 1];
        if (end > beg + 1)
            std::sort(flat.data() + beg, flat.data() + end);
    }

    auto sortEnd = std::chrono::steady_clock::now();

    // Phase 4: Serial scan — emit triplets and composite barcodes directly.
    // The CSR scan visits (CB, tag, gene) groups in sorted order. Composite barcodes
    // and triplets are built incrementally — no intermediate unordered_map, no re-sort.
    // Load probe/feature list first (needed for gene validation and MEX output).
    std::vector<std::string> geneIds;
    uint32_t maxGeneIdx = 0;
    {
        bool usedProbeList = false;
        if (!P.pSolo.probeListPath.empty() && P.pSolo.probeListPath != "-") {
            std::ifstream probeFile(P.pSolo.probeListPath);
            if (probeFile.is_open()) {
                std::string line;
                while (std::getline(probeFile, line)) {
                    if (!line.empty() && line[0] == '#') continue;
                    size_t b = line.find_first_not_of(" \t\r\n");
                    size_t e = line.find_last_not_of(" \t\r\n");
                    if (b == std::string::npos) continue;
                    geneIds.push_back(line.substr(b, e - b + 1));
                }
                probeFile.close();
                usedProbeList = true;
            }
        }
        if (!usedProbeList) {
            std::ostringstream err;
            err << "ERROR: Probe list unavailable (path: " << P.pSolo.probeListPath << ")";
            P.inOut->logMain << err.str() << endl;
            throw std::runtime_error(err.str());
        }
    }

    size_t nTripletGroups = 0;
    InlineMatrixBundle inlineMatrix;
    std::vector<MexWriter::Triplet> &triplets = inlineMatrix.triplets;
    std::vector<std::string> &compositeBarcodes = inlineMatrix.matrixData.barcodes;
    std::vector<uint64_t> &compositeCbTagKeys = inlineMatrix.cbTagKeys;
    triplets.reserve(totalEntries / 2);
    compositeBarcodes.reserve(nCB * 4);
    compositeCbTagKeys.reserve(nCB * 4);

    uint64_t prevCbTag = UINT64_MAX;
    uint32_t cellIdx = 0;
    std::vector<uint32_t> cellUMIs;
    std::vector<uint32_t> cellGenes;
    cellUMIs.reserve(nCB * 4);
    cellGenes.reserve(nCB * 4);

    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        const size_t beg = cbOffsets[iCB];
        const size_t end = cbOffsets[iCB + 1];
        const uint32_t cbIdx = indCB[iCB];

        size_t j = beg;
        while (j < end) {
            uint64_t packed = flat[j];
            uint8_t  tag  = static_cast<uint8_t>((packed >> 39) & 0x1F);
            uint16_t gene = static_cast<uint16_t>((packed >> 24) & 0x7FFF);

            if (tag == 0) {
                do { ++j; } while (j < end && (flat[j] >> 39) == 0);
                continue;
            }

            uint32_t dedupCount = 0;
            while (j < end) {
                uint64_t p = flat[j];
                if (static_cast<uint8_t>((p >> 39) & 0x1F) != tag ||
                    static_cast<uint16_t>((p >> 24) & 0x7FFF) != gene)
                    break;
                ++dedupCount;
                ++j;
            }

            uint64_t cbTagKey = (static_cast<uint64_t>(cbIdx) << 8) | tag;
            if (cbTagKey != prevCbTag) {
                if (cbIdx < pSolo.cbWLstr.size() &&
                    tag < gCanonicalTags.size() && !gCanonicalTags[tag].empty()) {
                    compositeBarcodes.push_back(pSolo.cbWLstr[cbIdx] + gCanonicalTags[tag]);
                    compositeCbTagKeys.push_back(cbTagKey);
                    cellIdx = static_cast<uint32_t>(compositeBarcodes.size() - 1);
                    cellUMIs.push_back(0);
                    cellGenes.push_back(0);
                } else {
                    prevCbTag = cbTagKey;
                    continue;
                }
                prevCbTag = cbTagKey;
            }

            if (gene > 0 && gene <= geneIds.size()) {
                triplets.push_back({cellIdx, gene - 1, dedupCount});
                cellUMIs.back() += dedupCount;
                cellGenes.back() += 1;
                if (gene > maxGeneIdx) maxGeneIdx = gene;
            }
            ++nTripletGroups;
        }
    }

    { std::vector<uint64_t>().swap(flat); }
    { std::vector<size_t>().swap(cbOffsets); }

    // Pad gene list if hash references beyond probe list
    if (maxGeneIdx > geneIds.size()) {
        P.inOut->logMain << "WARNING: probe list has " << geneIds.size()
                         << " entries but hash references probe index " << maxGeneIdx
                         << " — padding with placeholders" << endl;
        for (uint32_t p = static_cast<uint32_t>(geneIds.size()); p < maxGeneIdx; ++p)
            geneIds.push_back("UNKNOWN_PROBE_" + std::to_string(p + 1));
    }

    auto groupEnd = std::chrono::steady_clock::now();

    double cbCountSec  = std::chrono::duration<double>(cbCountEnd - collapseStart).count();
    double scatterSec  = std::chrono::duration<double>(scatterEnd - cbCountEnd).count();
    double sortSec     = std::chrono::duration<double>(sortEnd - scatterEnd).count();
    double scanSec     = std::chrono::duration<double>(groupEnd - sortEnd).count();
    double totalCollapseSec = std::chrono::duration<double>(groupEnd - collapseStart).count();
    P.inOut->logMain << "Collapse timing: cb_count=" << std::fixed << std::setprecision(1)
                     << cbCountSec << "s scatter=" << scatterSec << "s sort=" << sortSec
                     << "s scan=" << scanSec << "s total=" << totalCollapseSec << "s" << endl;

    P.inOut->logMain << "Found " << nCB << " unique CBs and " << nTripletGroups
                     << " (CB, gene, tag) groups" << endl;

    // Populate remaining InlineMatrixBundle / SampleMatrixData fields
    const uint32_t nCells = static_cast<uint32_t>(compositeBarcodes.size());
    SampleMatrixData &matrix = inlineMatrix.matrixData;
    matrix.nCells = nCells;
    matrix.nGenes = static_cast<uint32_t>(geneIds.size());
    matrix.countMatStride = 3;
    matrix.features = std::move(geneIds);
    matrix.nUMIperCB = std::move(cellUMIs);
    matrix.nGenePerCB = std::move(cellGenes);
    matrix.countCellGeneUMIindex.assign(nCells + 1, 0);
    {
        uint32_t offset = 0;
        size_t ti = 0;
        for (uint32_t cell = 0; cell < nCells; ++cell) {
            matrix.countCellGeneUMIindex[cell] = offset;
            while (ti < triplets.size() && triplets[ti].cell_idx == cell) {
                matrix.countCellGeneUMI.push_back(triplets[ti].gene_idx);
                matrix.countCellGeneUMI.push_back(triplets[ti].count);
                matrix.countCellGeneUMI.push_back(0);
                offset += matrix.countMatStride;
                ++ti;
            }
        }
        matrix.countCellGeneUMIindex[nCells] = offset;
    }

    P.inOut->logMain << "Found " << nCells << " unique (CB, TAG) combinations" << endl;
    P.inOut->logMain << "  Genes: " << matrix.nGenes << ", Entries: " << triplets.size() << endl;

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
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished direct hash collapse (MEX-only)" << endl;
    
    // Write MEX directly — triplets are already in (cell_idx, gene_idx) order from the sorted scan
    std::string mexDir = P.outFileNamePrefix + pSolo.outFileNames[0] + SoloFeatureTypes::Names[featureType] + "/raw/";
    createDirectory(mexDir, P.runDirPerm, "Solo raw MEX directory", P);
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
