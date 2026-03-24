// Non-Flex Solo inline-hash bridge: collapse UMIs directly from aggregated bridge hash
// (tagless key [bridgeCB24][UMI24][GENE16]) without materializeRGUFromHash() / rGeneUMI.

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
#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace {

struct BridgeHashRec {
    uint32_t wlCb;
    uint32_t gene;
    uint32_t umi24;
    uint32_t count;
};

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
        why = "multimapper Solo output is not supported on this path (use --soloMultiMappers Unique)";
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
                  "--outSAMtype None, do not request sorted-BAM read-id tracking, or unset STAR_SOLO_NONFLEX_HASH_BRIDGE.\n";
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
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Direct bridge-hash UMI collapse (no rGeneUMI materialization), hash_entries="
                     << totalHashSize
                     << " thread_hashes=" << threadHashCount
                     << " merged_ambiguous_hash_entries=" << mergedAmbigHashSize
                     << endl;

    std::vector<BridgeHashRec> recs;
    recs.reserve(totalHashSize);

    auto appendHashRecords = [&](SoloReadFeature *srcFeat, khash_t(cg_agg) *&hash, bool clearCompactMaps) {
        if (srcFeat == nullptr || hash == nullptr) {
            return;
        }
        for (khiter_t it = kh_begin(hash); it != kh_end(hash); ++it) {
            if (!kh_exist(hash, it))
                continue;
            const uint64_t key = kh_key(hash, it);
            const uint32_t val = kh_val(hash, it);

            uint32_t compactCb = 0, umi24 = 0;
            uint16_t compactGene = 0;
            unpackBridgeCgAggKey(key, &compactCb, &umi24, &compactGene);

            const uint32_t wlCb = srcFeat->bridgeCompactToWl(compactCb);
            if (wlCb == static_cast<uint32_t>(-1))
                continue;
            const uint32_t geneFull = srcFeat->bridgeCompactToGene(compactGene);
            if (geneFull == static_cast<uint32_t>(-1))
                continue;

            recs.push_back({wlCb, geneFull, umi24, val});
        }

        kh_destroy(cg_agg, hash);
        hash = nullptr;
        if (clearCompactMaps) {
            decltype(srcFeat->bridgeCbCompactByWl_)().swap(srcFeat->bridgeCbCompactByWl_);
            std::vector<uint32_t>().swap(srcFeat->bridgeCbWlByCompact_);
            decltype(srcFeat->bridgeGeneCompactByFull_)().swap(srcFeat->bridgeGeneCompactByFull_);
            std::vector<uint32_t>().swap(srcFeat->bridgeGeneFullByCompact_);
        }
    };

    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_ && kh_size(readFeatAll[ii]->inlineHash_) > 0) {
            appendHashRecords(readFeatAll[ii], readFeatAll[ii]->inlineHash_, true);
        }
    }

    if (readFeatSum->inlineHash_) {
        if (kh_size(readFeatSum->inlineHash_) > 0) {
            appendHashRecords(readFeatSum, readFeatSum->inlineHash_, true);
        } else {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
            decltype(readFeatSum->bridgeCbCompactByWl_)().swap(readFeatSum->bridgeCbCompactByWl_);
            std::vector<uint32_t>().swap(readFeatSum->bridgeCbWlByCompact_);
            decltype(readFeatSum->bridgeGeneCompactByFull_)().swap(readFeatSum->bridgeGeneCompactByFull_);
            std::vector<uint32_t>().swap(readFeatSum->bridgeGeneFullByCompact_);
        }
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Drained thread-local + merged ambiguous bridge hashes after extraction (recs="
                     << recs.size() << ")"
                     << endl;

    if (recs.empty()) {
        P.inOut->logMain << "WARNING: collapseUMIall_fromBridgeHash: no records after CB/gene mapping" << endl;
        nCB = 0;
        return;
    }

    std::sort(recs.begin(), recs.end(), [](const BridgeHashRec &a, const BridgeHashRec &b) {
        if (a.wlCb != b.wlCb)
            return a.wlCb < b.wlCb;
        if (a.gene != b.gene)
            return a.gene < b.gene;
        return a.umi24 < b.umi24;
    });

    // Upper bound on final (CB,gene) matrix rows = unique (wlCb,gene) runs (UMI multiplicity collapses to one row).
    size_t nCbGeneSeg = 1;
    for (size_t i = 1; i < recs.size(); ++i) {
        if (recs[i].wlCb != recs[i - 1].wlCb || recs[i].gene != recs[i - 1].gene)
            ++nCbGeneSeg;
    }
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Unique (CB,gene) segments (matrix row upper bound)="
                     << nCbGeneSeg << endl;

    // Unique CBs in sorted order -> indCB
    std::vector<uint32_t> sortedCBs;
    sortedCBs.reserve(4096);
    for (const auto &r : recs) {
        if (sortedCBs.empty() || sortedCBs.back() != r.wlCb)
            sortedCBs.push_back(r.wlCb);
    }

    nCB = static_cast<uint32_t>(sortedCBs.size());
    indCB.resize(nCB);
    indCBwl.assign(pSolo.cbWLsize, static_cast<uint32_t>(-1));
    for (uint32_t i = 0; i < nCB; ++i) {
        indCB[i] = sortedCBs[i];
        indCBwl[sortedCBs[i]] = i;
    }

    // Per-CB slice [cbPtr[iCB], cbPtr[iCB+1])
    std::vector<size_t> cbPtr(nCB + 1, 0);
    {
        size_t j = 0;
        for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
            cbPtr[iCB] = j;
            const uint32_t want = indCB[iCB];
            while (j < recs.size() && recs[j].wlCb == want)
                ++j;
        }
        cbPtr[nCB] = recs.size();
    }

    nReadPerCB.assign(nCB, 0);
    nReadPerCBmax = 0;
    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        uint64_t sum = 0;
        for (size_t j = cbPtr[iCB]; j < cbPtr[iCB + 1]; ++j)
            sum += recs[j].count;
        nReadPerCB[iCB] = static_cast<uint32_t>(sum > UINT32_MAX ? UINT32_MAX : sum);
        if (nReadPerCB[iCB] > nReadPerCBmax)
            nReadPerCBmax = nReadPerCB[iCB];
    }

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

    countMatMult.s = 1 + pSolo.multiMap.yes.N * pSolo.umiDedup.yes.N;
    countMatMult.m.clear();
    countMatMult.i.assign(nCB + 1, 0);

    // Per-gene UMI correction only needs nU0 * umiArrayStride; grow on demand (not nReadPerCBmax).
    std::vector<uint32_t> umiArray;
    umiArray.reserve(umiArrayStride * 64);

    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        const size_t cbBeg = cbPtr[iCB];
        const size_t cbEnd = cbPtr[iCB + 1];

        unordered_map<uintUMI, unordered_map<uint32_t, uint32_t>> umiGeneMapCount, umiGeneMapCount0;

        nGenePerCB[iCB] = 0;
        nUMIperCB[iCB] = 0;
        countCellGeneUMIindex[iCB + 1] = countCellGeneUMIindex[iCB];

        // Gene segments for this CB
        std::vector<uint32_t> gID;
        std::vector<size_t> gBeg, gEnd;
        for (size_t j = cbBeg; j < cbEnd;) {
            const uint32_t gid = recs[j].gene;
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
            while (k < cbEnd && recs[k].gene == gid)
                ++k;
            gEnd.push_back(k);
            j = k;
        }

        const uint32_t nGenes = static_cast<uint32_t>(gID.size());
        vector<unordered_map<uintUMI, uintUMI>> umiCorrected(static_cast<size_t>(nGenes));

        for (uint32_t iG = 0; iG < nGenes; ++iG) {
            const size_t a = gBeg[iG];
            const size_t b = gEnd[iG];
            const uint32_t nU0 = static_cast<uint32_t>(b - a);
            if (nU0 == 0)
                continue;

            if (umiArray.size() < static_cast<size_t>(nU0) * umiArrayStride)
                umiArray.resize(static_cast<size_t>(nU0) * umiArrayStride);

            for (uint32_t t = 0; t < nU0; ++t) {
                const BridgeHashRec &rr = recs[a + t];
                umiArray[static_cast<size_t>(t) * umiArrayStride + 0] = rr.umi24;
                umiArray[static_cast<size_t>(t) * umiArrayStride + 1] = rr.count;
                umiArray[static_cast<size_t>(t) * umiArrayStride + 2] = static_cast<uint32_t>(-1);
            }

            for (uint64_t iu = 0; iu < static_cast<uint64_t>(nU0) * umiArrayStride; iu += umiArrayStride)
                umiGeneMapCount0[umiArray[iu + 0]][iG] += umiArray[iu + 1];

            umiArrayCorrect_CR(nU0, umiArray.data(), true, false, umiCorrected[iG]);

            for (uint64_t iu = 0; iu < static_cast<uint64_t>(nU0) * umiArrayStride; iu += umiArrayStride)
                umiGeneMapCount[umiArray[iu + 2]][iG] += umiArray[iu + 1];
        }

        // MultiGeneUMI_CR global resolution (mirrors SoloFeature_collapseUMIall.cpp)
        vector<uint32> geneCounts(static_cast<size_t>(nGenes), 0);

        for (const auto &iu : umiGeneMapCount) {
            uint32 maxu = 0, maxg = static_cast<uint32>(-1);
            for (const auto &ig : iu.second) {
                if (ig.second > maxu) {
                    maxu = ig.second;
                    maxg = ig.first;
                } else if (ig.second == maxu) {
                    maxg = static_cast<uint32>(-1);
                }
            }

            if (maxg + 1 == 0)
                continue;

            for (const auto &ig : umiGeneMapCount0[iu.first]) {
                if (ig.second > umiGeneMapCount0[iu.first][maxg]) {
                    maxg = static_cast<uint32>(-1);
                    break;
                }
            }

            if (maxg + 1 != 0)
                geneCounts[maxg]++;
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

    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Finished direct bridge-hash UMI collapse, nCB=" << nCB
                     << " wall=" << elapsedSec() << " s" << endl;
}
