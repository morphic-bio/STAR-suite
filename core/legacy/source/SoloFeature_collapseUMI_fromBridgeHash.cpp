// Non-Flex Solo inline-hash bridge: collapse UMIs directly from aggregated bridge hash
// (tagless key [bridgeCB24][UMI24][GENE16]) without materializeRGUFromHash() / rGeneUMI.
//
// Minimal Unique (+ optional CR multimap rescue at alignment) path (v8+ / v9 tidy):
// - Per-thread: tuple khash -> stable slot id; 8-byte packed slots (no insert-time CB→slot list).
// - Collapse: merge threads into one global tuple hash + compact map wl→slot ids for observed CBs only.
// - Per observed whitelist CB (global sorted order): unpack slots to scratch, sort by (gene, umi), collapseOneBarcodeRows.
// - countMatMult / legacy multimapper redistribution not populated
//
// MultiGeneUMI_CR (v10): per-CB flat aggregates only — no nested umi→gene maps and no
// `vector<unordered_map>` umiCorrected. Corrected counts live in
// `unordered_map<uint64_t,uint32_t> corrUmiGeneCount` keyed by 64-bit pack `(umi<<32)|geneIdx`
// (sorts by UMI then gene); first-seen keys appended to `corrPackedKeys`, sorted, then scanned
// in runs by corrected UMI. Original per-(origUmi,geneSlice) counts for the legacy tie-break
// are merged into `origByOrigUmi[origUmi]` as sparse `(geneIdx,count)` vectors (same information
// as legacy `umiGeneMapCount0`). `umiArrayCorrect_CR(..., readInfoRec=false, ..., dummyMap)` —
// correction path is unchanged; only the correction map recording is skipped on this path.

#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "SoloCommon.h"
#include "ParametersSolo.h"
#include "hash_shims_cpp_compat.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "IncludeDefine.h"

#include <algorithm>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

// Per-CB scratch row (single-barcode slice unpacked from shard khash).
struct ShardRow {
    uint32_t wl;
    uint32_t gene;
    uint32_t umi24;
    uint32_t count;
};

// Packed shard-aggregate key: [wlCb:22][geneFull:18][umi24:24] (64 bits).

static inline uint64_t packShardCollapseAggKey(uint32_t wlCb, uint32_t geneFull, uint32_t umi24)
{
    return (uint64_t(wlCb & ((1u << 22) - 1)) << 42) | (uint64_t(geneFull & ((1u << 18) - 1)) << 24)
           | uint64_t(umi24 & 0xFFFFFFu);
}

static inline void unpackShardCollapseAggKey(uint64_t key, uint32_t *wlCb, uint32_t *geneFull, uint32_t *umi24)
{
    if (wlCb)
        *wlCb = static_cast<uint32_t>((key >> 42) & ((1u << 22) - 1));
    if (geneFull)
        *geneFull = static_cast<uint32_t>((key >> 24) & ((1u << 18) - 1));
    if (umi24)
        *umi24 = static_cast<uint32_t>(key & 0xFFFFFFu);
}

static void unpackPackedSlotToShardRow(uint32_t wl, uint64_t packedSlot, ShardRow *out)
{
    uint32_t umi = 0, gene = 0, cnt = 0;
    bool ovf = false;
    unpackBridgePackedSlot(packedSlot, &umi, &gene, &cnt, &ovf);
    (void)ovf;
    out->wl = wl;
    out->gene = gene;
    out->umi24 = umi;
    out->count = cnt;
}

} // namespace

void SoloFeature::collapseUMIall_fromBridgeHash()
{
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsedSec = [&t0]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };

    const char *why = nullptr;
    if (pSolo.trackReadIdsForTags)
        why = "trackReadIdsForTags is enabled (needs per-read recordReadInfo)";
    else if (pSolo.CBmatchWL.oneExact)
        why = "CB oneExact gating is not supported on the non-Flex direct bridge path";
    else if (P.outSAMtype.empty() || P.outSAMtype[0] != "None")
        why = "--outSAMtype must be None for direct bridge collapse (or disable STAR_SOLO_NONFLEX_HASH_BRIDGE)";
    else if (pSolo.multiMap.yes.multi)
        why = "multimapper Solo output is not supported on this path (use --soloMultiMappers Unique). "
              "CR multimap rescue is applied at alignment; --soloCrMultimapRescue yes is compatible with Unique.";
    else if (!pSolo.umiFiltering.MultiGeneUMI_CR)
        why = "--soloUMIfiltering MultiGeneUMI_CR is required";
    else if (pSolo.umiFiltering.MultiGeneUMI || pSolo.umiFiltering.MultiGeneUMI_All)
        why = "MultiGeneUMI / MultiGeneUMI_All cannot be combined with this path";
    else if (!pSolo.umiDedup.yes.CR || pSolo.umiDedup.yes.N != 1u)
        why = "exactly one --soloUMIdedup mode and it must be 1MM_CR for direct bridge collapse";

    if (why != nullptr) {
        ostringstream errOut;
        errOut << "EXITING because of fatal PARAMETER error: non-Flex inline-hash bridge direct collapse is not "
                  "supported for this configuration.\n"
               << "Reason: " << why << "\n"
               << "SOLUTION: use --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloMultiMappers Unique "
                  "(optional --soloCrMultimapRescue yes; rescue is alignment-time), --outSAMtype None, or unset "
                  "STAR_SOLO_NONFLEX_HASH_BRIDGE.\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    }

    if (!readFeatSum) {
        P.inOut->logMain << "WARNING: collapseUMIall_fromBridgeHash: readFeatSum is null" << endl;
        nCB = 0;
        return;
    }

    size_t totalHashSize = 0;
    size_t threadHashCount = 0;
    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_) {
            const size_t sz = kh_size(readFeatAll[ii]->inlineHash_);
            if (sz > 0) {
                totalHashSize += sz;
                ++threadHashCount;
            }
        }
    }
    const size_t mergedAmbigHashSize = (readFeatSum->inlineHash_ != nullptr) ? kh_size(readFeatSum->inlineHash_) : 0;
    totalHashSize += mergedAmbigHashSize;

    if (totalHashSize == 0) {
        P.inOut->logMain << "WARNING: collapseUMIall_fromBridgeHash: no thread-local or merged bridge hash entries"
                         << endl;
        if (readFeatSum->inlineHash_) {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
        }
        nCB = 0;
        return;
    }

    time_t rawTime;
    time(&rawTime);

    uint64_t threadOverflowEvents = 0;
    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] != nullptr) {
            threadOverflowEvents += readFeatAll[ii]->bridgeSlotOverflowEvents_;
        }
    }
    threadOverflowEvents += readFeatSum->bridgeSlotOverflowEvents_;

    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Direct bridge-hash UMI collapse (v10: v9 + flat MultiGeneUMI_CR aggregates), "
                        "hash_entries="
                     << totalHashSize << " thread_hashes=" << threadHashCount
                     << " merged_ambiguous_hash_entries=" << mergedAmbigHashSize
                     << " thread_packed_slot_overflow_events=" << threadOverflowEvents << endl;

    khash_t(cg_agg) *globalTuple = kh_init(cg_agg);
    if (globalTuple == nullptr) {
        P.inOut->logMain << "WARNING: collapseUMIall_fromBridgeHash: kh_init(globalTuple) failed" << endl;
        nCB = 0;
        return;
    }
    const unsigned bucketHintGlobal = static_cast<unsigned>(std::max<size_t>(8, totalHashSize + 64));
    kh_resize(cg_agg, globalTuple, bucketHintGlobal);

    std::vector<uint64_t> globalPackedSlots;
    std::unordered_map<uint32_t, std::vector<uint32_t>> slotIdsByObservedWl;
    slotIdsByObservedWl.reserve(std::min(static_cast<size_t>(65536), totalHashSize));
    std::unordered_set<uint32_t> cbSeen;
    cbSeen.reserve(std::min(static_cast<size_t>(65536), totalHashSize));
    uint64_t mergeOverflowEvents = 0;

    auto drainSourceIntoGlobal = [&](SoloReadFeature *srcFeat, khash_t(cg_agg) *&hash, bool clearCompactMaps) {
        if (srcFeat == nullptr || hash == nullptr) {
            return;
        }
        for (khiter_t it = kh_begin(hash); it != kh_end(hash); ++it) {
            if (!kh_exist(hash, it)) {
                continue;
            }
            const uint64_t tkey = kh_key(hash, it);
            const uint32_t sid = kh_val(hash, it);
            if (sid >= srcFeat->bridgePackedSlots_.size()) {
                continue;
            }
            const uint64_t srcPack = srcFeat->bridgePackedSlots_[sid];

            uint32_t compactCb = 0, umi24 = 0;
            uint16_t compactGene = 0;
            unpackBridgeCgAggKey(tkey, &compactCb, &umi24, &compactGene);

            const uint32_t wlCb = srcFeat->bridgeCompactToWl(compactCb);
            if (wlCb == static_cast<uint32_t>(-1)) {
                continue;
            }
            const uint32_t geneFull = srcFeat->bridgeCompactToGene(compactGene);
            if (geneFull == static_cast<uint32_t>(-1)) {
                continue;
            }

            if (wlCb >= (1u << 22) || geneFull >= (1u << 18)) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: non-Flex bridge global collapse key overflow.\n"
                       << "wlCb=" << wlCb << " (max " << ((1u << 22) - 1) << ") or geneFull=" << geneFull << " (max "
                       << ((1u << 18) - 1) << ").\n"
                       << "SOLUTION: disable STAR_SOLO_NONFLEX_HASH_BRIDGE for this reference/whitelist or extend "
                          "packing.\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
            }

            const uint64_t gkey = packShardCollapseAggKey(wlCb, geneFull, umi24);
            int absent = 0;
            khiter_t git = kh_put(cg_agg, globalTuple, gkey, &absent);
            if (absent) {
                const uint32_t newSid = static_cast<uint32_t>(globalPackedSlots.size());
                globalPackedSlots.push_back(srcPack);
                kh_val(globalTuple, git) = newSid;
                cbSeen.insert(wlCb);
                slotIdsByObservedWl[wlCb].push_back(newSid);
            } else {
                const uint32_t gid = kh_val(globalTuple, git);
                globalPackedSlots[gid] =
                    bridgePackedSlotMerge(globalPackedSlots[gid], srcPack, &mergeOverflowEvents);
            }
        }

        kh_destroy(cg_agg, hash);
        hash = nullptr;
        std::vector<uint64_t>().swap(srcFeat->bridgePackedSlots_);
        srcFeat->bridgeSlotOverflowEvents_ = 0;
        if (clearCompactMaps) {
            decltype(srcFeat->bridgeCbCompactByWl_)().swap(srcFeat->bridgeCbCompactByWl_);
            std::vector<uint32_t>().swap(srcFeat->bridgeCbWlByCompact_);
            decltype(srcFeat->bridgeGeneCompactByFull_)().swap(srcFeat->bridgeGeneCompactByFull_);
            std::vector<uint32_t>().swap(srcFeat->bridgeGeneFullByCompact_);
        }
    };

    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_ && kh_size(readFeatAll[ii]->inlineHash_) > 0) {
            drainSourceIntoGlobal(readFeatAll[ii], readFeatAll[ii]->inlineHash_, true);
        }
    }

    if (readFeatSum->inlineHash_) {
        if (kh_size(readFeatSum->inlineHash_) > 0) {
            drainSourceIntoGlobal(readFeatSum, readFeatSum->inlineHash_, true);
        } else {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
            decltype(readFeatSum->bridgeCbCompactByWl_)().swap(readFeatSum->bridgeCbCompactByWl_);
            std::vector<uint32_t>().swap(readFeatSum->bridgeCbWlByCompact_);
            decltype(readFeatSum->bridgeGeneCompactByFull_)().swap(readFeatSum->bridgeGeneCompactByFull_);
            std::vector<uint32_t>().swap(readFeatSum->bridgeGeneFullByCompact_);
        }
    }

    const size_t totalAggEntries = kh_size(globalTuple);
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Merged into global tuple hash (unique_keys~=" << totalAggEntries
                     << ") merge_packed_slot_overflow_events=" << mergeOverflowEvents << endl;

    if (totalAggEntries == 0) {
        P.inOut->logMain << "WARNING: collapseUMIall_fromBridgeHash: no records after CB/gene mapping" << endl;
        nCB = 0;
        kh_destroy(cg_agg, globalTuple);
        return;
    }

    std::vector<uint32_t> sortedCBs(cbSeen.begin(), cbSeen.end());
    std::sort(sortedCBs.begin(), sortedCBs.end());

    const auto geneUmiLess = [](const ShardRow &a, const ShardRow &b) {
        if (a.gene != b.gene) {
            return a.gene < b.gene;
        }
        return a.umi24 < b.umi24;
    };

    size_t nCbGeneSeg = 0;
    std::vector<ShardRow> segScratch;
    segScratch.reserve(4096);
    for (uint32_t wl : sortedCBs) {
        const auto itIds = slotIdsByObservedWl.find(wl);
        if (itIds == slotIdsByObservedWl.end()) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: internal bridge missing observed-wl slot list for wlCb=" << wl
                   << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
        }
        const std::vector<uint32_t> &ids = itIds->second;
        const size_t n = ids.size();
        if (n == 0) {
            continue;
        }
        if (segScratch.size() < n) {
            segScratch.resize(n);
        }
        for (size_t i = 0; i < n; ++i) {
            const uint32_t gsid = ids[i];
            if (gsid >= globalPackedSlots.size()) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: internal bridge bad global slot id (segment count) wlCb="
                        << wl << "\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
            }
            unpackPackedSlotToShardRow(wl, globalPackedSlots[gsid], &segScratch[i]);
        }
        std::sort(segScratch.begin(), segScratch.begin() + n, geneUmiLess);
        uint32_t prevG = static_cast<uint32_t>(-1);
        for (size_t i = 0; i < n; ++i) {
            if (i == 0 || segScratch[i].gene != prevG) {
                ++nCbGeneSeg;
                prevG = segScratch[i].gene;
            }
        }
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Unique (CB,gene) segments (matrix row upper bound) ~"
                     << nCbGeneSeg << endl;

    nCB = static_cast<uint32_t>(sortedCBs.size());
    indCB.resize(nCB);
    indCBwl.assign(pSolo.cbWLsize, static_cast<uint32_t>(-1));
    for (uint32_t i = 0; i < nCB; ++i) {
        indCB[i] = sortedCBs[i];
        indCBwl[sortedCBs[i]] = i;
    }

    std::vector<ShardRow> cbScratch;
    cbScratch.reserve(4096);

    nReadPerCB.assign(nCB, 0);
    nReadPerCBmax = 0;

    countMatStride = pSolo.umiDedup.yes.N + 1;
    nUMIperCB.assign(nCB, 0);
    nGenePerCB.assign(nCB, 0);

    countCellGeneUMI.clear();
    {
        const size_t matSlots = nCbGeneSeg * countMatStride;
        const size_t minSeed = static_cast<size_t>(nCB + 1u) * countMatStride;
        countCellGeneUMI.resize(std::max(matSlots, minSeed), 0);
    }
    countCellGeneUMIindex.assign(nCB + 1, 0);

    // Minimal Unique path: no legacy multimapper redistribution matrix.
    countMatMult.s = 0;
    countMatMult.m.clear();
    countMatMult.m.shrink_to_fit();
    countMatMult.i.clear();
    countMatMult.i.shrink_to_fit();

    std::vector<uint32_t> umiArray;
    umiArray.reserve(umiArrayStride * 64);

    // Per-CB MultiGeneUMI_CR scratch (direct bridge only): flat maps + key list, no nested umi→gene maps.
    // Pack layout MSB→LSB: [umi:32][geneIdx:32] — sorts by (UMI, gene slot index) for grouping corrected UMIs.
    auto packUmiGeneIdx32 = [](uint32_t umi24, uint32_t geneIdx) -> uint64_t {
        return (uint64_t{umi24} << 32) | uint64_t{geneIdx};
    };

    std::unordered_map<uint64_t, uint32_t> corrUmiGeneCount;
    std::vector<uint64_t> corrPackedKeys;
    // Original (uncorrected) counts keyed by literal original UMI: same genes as legacy umiGeneMapCount0[origUmi].
    std::unordered_map<uintUMI, std::vector<std::pair<uint32_t, uint32_t>>> origByOrigUmi;
    unordered_map<uintUMI, uintUMI> umiCorrUnusedBridge;
    corrUmiGeneCount.reserve(256);
    origByOrigUmi.reserve(128);
    corrPackedKeys.reserve(128);

    auto collapseOneBarcodeRows = [&](const ShardRow *base, size_t nRecs, uint32_t iCB) {
        if (nRecs == 0) {
            return;
        }

        nGenePerCB[iCB] = 0;
        nUMIperCB[iCB] = 0;
        countCellGeneUMIindex[iCB + 1] = countCellGeneUMIindex[iCB];

        std::vector<uint32_t> gID;
        std::vector<size_t> gBeg, gEnd;
        for (size_t j = 0; j < nRecs;) {
            const uint32_t gid = base[j].gene;
            if (gid & geneMultMark) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: multimapper gene encountered in non-Flex bridge direct "
                          "collapse.\n"
                       << "This path requires unique-gene alignments only (--soloMultiMappers Unique).\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
            }
            gID.push_back(gid);
            gBeg.push_back(j);
            size_t k = j;
            while (k < nRecs && base[k].gene == gid)
                ++k;
            gEnd.push_back(k);
            j = k;
        }

        const uint32_t nGenes = static_cast<uint32_t>(gID.size());

        corrUmiGeneCount.clear();
        corrPackedKeys.clear();
        origByOrigUmi.clear();

        for (uint32_t iG = 0; iG < nGenes; ++iG) {
            const size_t a = gBeg[iG];
            const size_t b = gEnd[iG];
            const uint32_t nU0 = static_cast<uint32_t>(b - a);
            if (nU0 == 0)
                continue;

            if (umiArray.size() < static_cast<size_t>(nU0) * umiArrayStride)
                umiArray.resize(static_cast<size_t>(nU0) * umiArrayStride);

            for (uint32_t t = 0; t < nU0; ++t) {
                const ShardRow &rr = base[a + t];
                umiArray[static_cast<size_t>(t) * umiArrayStride + 0] = rr.umi24;
                umiArray[static_cast<size_t>(t) * umiArrayStride + 1] = rr.count;
                umiArray[static_cast<size_t>(t) * umiArrayStride + 2] = static_cast<uint32_t>(-1);
            }

            // Original (uncorrected) counts for tie-break: legacy umiGeneMapCount0[origUmi][iG].
            for (uint64_t iu = 0; iu < static_cast<uint64_t>(nU0) * umiArrayStride; iu += umiArrayStride) {
                const uintUMI ou = umiArray[iu + 0];
                const uint32_t add = umiArray[iu + 1];
                auto &vec = origByOrigUmi[ou];
                bool found = false;
                for (auto &pr : vec) {
                    if (pr.first == iG) {
                        pr.second += add;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    vec.push_back(std::pair<uint32_t, uint32_t>(iG, add));
                }
            }

            umiArrayCorrect_CR(nU0, umiArray.data(), false, false, umiCorrUnusedBridge);

            for (uint64_t iu = 0; iu < static_cast<uint64_t>(nU0) * umiArrayStride; iu += umiArrayStride) {
                const uint64_t cp = packUmiGeneIdx32(umiArray[iu + 2], iG);
                const auto ins = corrUmiGeneCount.emplace(cp, umiArray[iu + 1]);
                if (ins.second) {
                    corrPackedKeys.push_back(cp);
                } else {
                    ins.first->second += umiArray[iu + 1];
                }
            }
        }

        std::sort(corrPackedKeys.begin(), corrPackedKeys.end());

        vector<uint32> geneCounts(static_cast<size_t>(nGenes), 0);

        for (size_t p = 0; p < corrPackedKeys.size();) {
            const uint32_t corrUmi = static_cast<uint32_t>(corrPackedKeys[p] >> 32);
            size_t q = p + 1;
            while (q < corrPackedKeys.size() && static_cast<uint32_t>(corrPackedKeys[q] >> 32) == corrUmi) {
                ++q;
            }

            uint32_t maxu = 0;
            uint32_t maxg = static_cast<uint32_t>(-1);
            for (size_t k = p; k < q; ++k) {
                const uint64_t key = corrPackedKeys[k];
                const uint32_t gidx = static_cast<uint32_t>(key);
                const uint32_t ctot = corrUmiGeneCount[key];
                if (ctot > maxu) {
                    maxu = ctot;
                    maxg = gidx;
                } else if (ctot == maxu) {
                    maxg = static_cast<uint32_t>(-1);
                }
            }

            if (maxg + 1u == 0u) {
                p = q;
                continue;
            }

            uint32_t winOrig = 0;
            bool reject = false;
            const auto itOrig = origByOrigUmi.find(static_cast<uintUMI>(corrUmi));
            if (itOrig != origByOrigUmi.end()) {
                const auto &ov = itOrig->second;
                for (const auto &pr : ov) {
                    if (pr.first == maxg) {
                        winOrig = pr.second;
                        break;
                    }
                }
                for (const auto &pr : ov) {
                    if (pr.second > winOrig) {
                        reject = true;
                        break;
                    }
                }
            }

            if (!reject) {
                geneCounts[maxg]++;
            }

            p = q;
        }

        for (uint32_t ig = 0; ig < nGenes; ++ig) {
            if (geneCounts[ig] == 0)
                continue;
            nGenePerCB[iCB]++;
            nUMIperCB[iCB] += geneCounts[ig];
            countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + 0] = gID[ig];
            countCellGeneUMI[countCellGeneUMIindex[iCB + 1] + pSolo.umiDedup.countInd.CR] = geneCounts[ig];
            countCellGeneUMIindex[iCB + 1] += countMatStride;
        }
    };

    // Global whitelist order (sorted observed CBs): slot ids from merge, local sort (gene, umi).
    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        const uint32_t wl = indCB[iCB];
        const auto itIds = slotIdsByObservedWl.find(wl);
        if (itIds == slotIdsByObservedWl.end()) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: internal bridge missing observed-wl slot list for wlCb=" << wl
                   << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
        }
        const std::vector<uint32_t> &ids = itIds->second;
        const size_t nRec = ids.size();
        if (nRec == 0) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: internal bridge empty slot list for wlCb=" << wl << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
        }
        if (cbScratch.size() < nRec) {
            cbScratch.resize(nRec);
        }
        uint64_t sumReads = 0;
        for (size_t k = 0; k < nRec; ++k) {
            const uint32_t gsid = ids[k];
            if (gsid >= globalPackedSlots.size()) {
                ostringstream errOut;
                errOut << "EXITING because of fatal ERROR: internal bridge bad global slot id for wlCb=" << wl
                       << "\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
            }
            unpackPackedSlotToShardRow(wl, globalPackedSlots[gsid], &cbScratch[k]);
            sumReads += cbScratch[k].count;
        }
        std::sort(cbScratch.begin(), cbScratch.begin() + nRec, geneUmiLess);
        nReadPerCB[iCB] = static_cast<uint32_t>(sumReads > UINT32_MAX ? UINT32_MAX : sumReads);
        if (nReadPerCB[iCB] > nReadPerCBmax)
            nReadPerCBmax = nReadPerCB[iCB];
        collapseOneBarcodeRows(cbScratch.data(), nRec, iCB);
    }

    kh_destroy(cg_agg, globalTuple);
    std::vector<uint64_t>().swap(globalPackedSlots);
    std::unordered_map<uint32_t, std::vector<uint32_t>>().swap(slotIdsByObservedWl);

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished direct bridge-hash UMI collapse, nCB=" << nCB
                     << " wall=" << elapsedSec() << " s" << endl;
}
