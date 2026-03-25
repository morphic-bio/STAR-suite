// Non-Flex Solo inline-hash bridge: collapse UMIs directly from aggregated bridge hash
// (tagless key [wlCb24][UMI24][geneFull16]) without materializeRGUFromHash() / rGeneUMI.
//
// Data flow (no global sort over the full tuple set):
// 1) Two passes over thread-local + merged hashes: count entries per whitelist CB, then CSR-fill
//    a flat (gene, umi24, count) array with per-CB slices [cbOffsets[i], cbOffsets[i+1]).
// 2) Per CB: in-place sort by (gene, umi24); per gene bucket run 1MM_CR via umiArrayCorrect_CR.
// 3) MultiGeneUMI_CR: collect flat (corrUmi, origUmi, geneIdx, count) rows for the CB, sort by
//    corrected UMI, resolve each corrected-UMI group with vector scans (no nested umi->gene maps).

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
#include <utility>
#include <vector>

namespace {

struct BridgeFlatSlot {
    uint32_t gene;
    uint32_t umi24;
    uint32_t count;
};

struct MgRow {
    uint32_t corr;
    uint32_t orig;
    uint32_t geneIdx;
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
                     << " ... Direct bridge-hash UMI collapse (CSR + flat MultiGene; no global tuple sort), hash_entries="
                     << totalHashSize
                     << " thread_hashes=" << threadHashCount
                     << " merged_ambiguous_hash_entries=" << mergedAmbigHashSize
                     << endl;

    std::unordered_map<uint32_t, size_t> cbEntryCount;
    cbEntryCount.reserve(std::min<size_t>(totalHashSize / 8 + 64, size_t{1} << 20));

    auto countHashEntriesPerWlCb = [&](SoloReadFeature *srcFeat, khash_t(cg_agg) *hash) {
        (void)srcFeat;
        if (hash == nullptr)
            return;
        for (khiter_t it = kh_begin(hash); it != kh_end(hash); ++it) {
            if (!kh_exist(hash, it))
                continue;
            const uint64_t key = kh_key(hash, it);
            uint32_t wlCb = 0, umi24 = 0;
            uint16_t gene16 = 0;
            unpackBridgeWlUmiGeneKey(key, &wlCb, &umi24, &gene16);
            (void)umi24;
            (void)gene16;
            ++cbEntryCount[wlCb];
        }
    };

    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_ && kh_size(readFeatAll[ii]->inlineHash_) > 0) {
            countHashEntriesPerWlCb(readFeatAll[ii], readFeatAll[ii]->inlineHash_);
        }
    }
    if (readFeatSum->inlineHash_ && kh_size(readFeatSum->inlineHash_) > 0) {
        countHashEntriesPerWlCb(readFeatSum, readFeatSum->inlineHash_);
    }

    std::vector<uint32_t> sortedCBs;
    sortedCBs.reserve(cbEntryCount.size());
    for (const auto &kv : cbEntryCount)
        sortedCBs.push_back(kv.first);
    std::sort(sortedCBs.begin(), sortedCBs.end());

    nCB = static_cast<uint32_t>(sortedCBs.size());
    indCB.resize(nCB);
    indCBwl.assign(pSolo.cbWLsize, static_cast<uint32_t>(-1));
    for (uint32_t i = 0; i < nCB; ++i) {
        indCB[i] = sortedCBs[i];
        indCBwl[sortedCBs[i]] = i;
    }

    std::vector<size_t> cbOffsets(nCB + 1, 0);
    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        const uint32_t wl = indCB[iCB];
        cbOffsets[iCB + 1] = cbOffsets[iCB] + cbEntryCount[wl];
    }
    if (cbOffsets[nCB] != totalHashSize) {
        ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: bridge hash CSR sizing mismatch (expected " << totalHashSize
               << " entries, got " << cbOffsets[nCB] << ").\n";
        exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
    }

    std::vector<BridgeFlatSlot> flat;
    flat.resize(totalHashSize);
    std::vector<size_t> cbWrite(nCB);
    for (uint32_t iCB = 0; iCB < nCB; ++iCB)
        cbWrite[iCB] = cbOffsets[iCB];

    auto drainHashIntoFlat = [&](SoloReadFeature *srcFeat, khash_t(cg_agg) *&hash) {
        (void)srcFeat;
        if (hash == nullptr)
            return;
        for (khiter_t it = kh_begin(hash); it != kh_end(hash); ++it) {
            if (!kh_exist(hash, it))
                continue;
            const uint64_t key = kh_key(hash, it);
            const uint32_t val = kh_val(hash, it);

            uint32_t wlCb = 0, umi24 = 0;
            uint16_t gene16 = 0;
            unpackBridgeWlUmiGeneKey(key, &wlCb, &umi24, &gene16);

            const uint32_t iCB = indCBwl[wlCb];
            const size_t pos = cbWrite[iCB]++;
            flat[pos].gene = static_cast<uint32_t>(gene16);
            flat[pos].umi24 = umi24;
            flat[pos].count = val;
        }

        kh_destroy(cg_agg, hash);
        hash = nullptr;
    };

    for (int ii = 0; ii < P.runThreadN; ++ii) {
        if (readFeatAll[ii] && readFeatAll[ii]->inlineHash_ && kh_size(readFeatAll[ii]->inlineHash_) > 0) {
            drainHashIntoFlat(readFeatAll[ii], readFeatAll[ii]->inlineHash_);
        }
    }

    if (readFeatSum->inlineHash_) {
        if (kh_size(readFeatSum->inlineHash_) > 0) {
            drainHashIntoFlat(readFeatSum, readFeatSum->inlineHash_);
        } else {
            kh_destroy(cg_agg, readFeatSum->inlineHash_);
            readFeatSum->inlineHash_ = nullptr;
        }
    }

    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        if (cbWrite[iCB] != cbOffsets[iCB + 1]) {
            ostringstream errOut;
            errOut << "EXITING because of fatal ERROR: bridge hash CSR fill mismatch for CB slot " << iCB << ".\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
        }
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... CSR flat layout built (entries=" << flat.size() << ")"
                     << endl;

    countMatStride = pSolo.umiDedup.yes.N + 1;
    nUMIperCB.assign(nCB, 0);
    nGenePerCB.assign(nCB, 0);
    nReadPerCB.assign(nCB, 0);
    nReadPerCBmax = 0;

    // Upper bound: at most one count-matrix slot per hash entry (distinct (CB,gene) <= totalHashSize).
    const size_t matSlotsUpper = totalHashSize * static_cast<size_t>(countMatStride);
    const size_t minSeed = static_cast<size_t>(nCB + 1u) * static_cast<size_t>(countMatStride);
    countCellGeneUMI.clear();
    countCellGeneUMI.resize(std::max(matSlotsUpper, minSeed), 0);
    countCellGeneUMIindex.assign(nCB + 1, 0);

    countMatMult.s = 1 + pSolo.multiMap.yes.N * pSolo.umiDedup.yes.N;
    countMatMult.m.clear();
    countMatMult.i.assign(nCB + 1, 0);

    std::vector<uint32_t> umiArray;
    umiArray.reserve(umiArrayStride * 64);

    std::vector<MgRow> mgBuf;
    std::vector<std::pair<uint32_t, uint32_t>> aggGene;
    std::vector<std::pair<uint32_t, uint32_t>> origAtCu;
    vector<uint32> geneCounts;
    static std::unordered_map<uintUMI, uintUMI> s_emptyUmiCorr;

    // Distinct (CB,gene) segments for logging only (fused into CB loop).
    std::vector<uint32_t> geneStamp(65536, 0);
    uint32_t stampVer = 1;
    size_t nCbGeneSeg = 0;

    for (uint32_t iCB = 0; iCB < nCB; ++iCB) {
        const size_t cbBeg = cbOffsets[iCB];
        const size_t cbEnd = cbOffsets[iCB + 1];
        const size_t sliceLen = cbEnd - cbBeg;
        BridgeFlatSlot *slice = flat.data() + cbBeg;

        nGenePerCB[iCB] = 0;
        nUMIperCB[iCB] = 0;
        countCellGeneUMIindex[iCB + 1] = countCellGeneUMIindex[iCB];

        if (sliceLen == 0) {
            if (++stampVer == 0u) {
                std::fill(geneStamp.begin(), geneStamp.end(), 0u);
                stampVer = 1u;
            }
            continue;
        }

        uint64_t readSum = 0;
        for (size_t j = cbBeg; j < cbEnd; ++j) {
            readSum += flat[j].count;
            const uint32_t g = flat[j].gene;
            if (geneStamp[g] != stampVer) {
                geneStamp[g] = stampVer;
                ++nCbGeneSeg;
            }
        }
        nReadPerCB[iCB] = static_cast<uint32_t>(readSum > UINT32_MAX ? UINT32_MAX : readSum);
        if (nReadPerCB[iCB] > nReadPerCBmax)
            nReadPerCBmax = nReadPerCB[iCB];
        if (++stampVer == 0u) {
            std::fill(geneStamp.begin(), geneStamp.end(), 0u);
            stampVer = 1u;
        }

        std::sort(slice, slice + sliceLen, [](const BridgeFlatSlot &a, const BridgeFlatSlot &b) {
            if (a.gene != b.gene)
                return a.gene < b.gene;
            return a.umi24 < b.umi24;
        });

        std::vector<uint32_t> gID;
        std::vector<size_t> gBeg, gEnd;
        for (size_t j = 0; j < sliceLen;) {
            const uint32_t gid = slice[j].gene;
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
            while (k < sliceLen && slice[k].gene == gid)
                ++k;
            gEnd.push_back(k);
            j = k;
        }

        const uint32_t nGenes = static_cast<uint32_t>(gID.size());
        mgBuf.clear();
        mgBuf.reserve(sliceLen);

        for (uint32_t iG = 0; iG < nGenes; ++iG) {
            const size_t a = gBeg[iG];
            const size_t b = gEnd[iG];
            const uint32_t nU0 = static_cast<uint32_t>(b - a);
            if (nU0 == 0)
                continue;

            if (umiArray.size() < static_cast<size_t>(nU0) * umiArrayStride)
                umiArray.resize(static_cast<size_t>(nU0) * umiArrayStride);

            for (uint32_t t = 0; t < nU0; ++t) {
                const BridgeFlatSlot &rr = slice[a + t];
                umiArray[static_cast<size_t>(t) * umiArrayStride + 0] = rr.umi24;
                umiArray[static_cast<size_t>(t) * umiArrayStride + 1] = rr.count;
                umiArray[static_cast<size_t>(t) * umiArrayStride + 2] = static_cast<uint32_t>(-1);
            }

            // readInfoRec=false: bridge path has no per-read recordReadInfo replay; skip umiCorr map inserts.
            umiArrayCorrect_CR(nU0, umiArray.data(), false, false, s_emptyUmiCorr);

            for (uint64_t iu = 0; iu < static_cast<uint64_t>(nU0) * umiArrayStride; iu += umiArrayStride) {
                MgRow row;
                row.orig = umiArray[iu + 0];
                row.corr = umiArray[iu + 2];
                row.geneIdx = iG;
                row.count = umiArray[iu + 1];
                mgBuf.push_back(row);
            }
        }

        std::sort(mgBuf.begin(), mgBuf.end(), [](const MgRow &x, const MgRow &y) { return x.corr < y.corr; });

        geneCounts.assign(static_cast<size_t>(nGenes), 0);

        for (size_t p = 0; p < mgBuf.size();) {
            size_t q = p + 1;
            const uint32_t cu = mgBuf[p].corr;
            while (q < mgBuf.size() && mgBuf[q].corr == cu)
                ++q;

            std::sort(mgBuf.begin() + p, mgBuf.begin() + q,
                      [](const MgRow &x, const MgRow &y) { return x.geneIdx < y.geneIdx; });

            aggGene.clear();
            for (size_t i = p; i < q;) {
                const uint32_t g0 = mgBuf[i].geneIdx;
                uint32_t s = 0;
                while (i < q && mgBuf[i].geneIdx == g0) {
                    s += mgBuf[i].count;
                    ++i;
                }
                aggGene.push_back({g0, s});
            }

            uint32_t maxu = 0;
            uint32_t maxg = static_cast<uint32_t>(-1);
            for (const auto &pr : aggGene) {
                if (pr.second > maxu) {
                    maxu = pr.second;
                    maxg = pr.first;
                } else if (pr.second == maxu) {
                    maxg = static_cast<uint32_t>(-1);
                }
            }

            if (maxg + 1u == 0u) {
                p = q;
                continue;
            }

            origAtCu.clear();
            for (size_t i = p; i < q; ++i) {
                if (mgBuf[i].orig == cu)
                    origAtCu.push_back({mgBuf[i].geneIdx, mgBuf[i].count});
            }
            std::sort(origAtCu.begin(), origAtCu.end(),
                      [](const std::pair<uint32_t, uint32_t> &x, const std::pair<uint32_t, uint32_t> &y) {
                          return x.first < y.first;
                      });
            {
                size_t o = 0;
                for (size_t i = 0; i < origAtCu.size();) {
                    const uint32_t g0 = origAtCu[i].first;
                    uint32_t s = 0;
                    while (i < origAtCu.size() && origAtCu[i].first == g0) {
                        s += origAtCu[i].second;
                        ++i;
                    }
                    origAtCu[o++] = {g0, s};
                }
                origAtCu.resize(o);
            }

            uint32_t baseOrig = 0;
            for (const auto &pr : origAtCu) {
                if (pr.first == maxg) {
                    baseOrig = pr.second;
                    break;
                }
            }

            for (const auto &pr : origAtCu) {
                if (pr.second > baseOrig) {
                    maxg = static_cast<uint32_t>(-1);
                    break;
                }
            }

            if (maxg + 1u != 0u)
                geneCounts[maxg]++;

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
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Unique (CB,gene) segments (matrix row upper bound)="
                     << nCbGeneSeg << " | Finished direct bridge-hash UMI collapse, nCB=" << nCB << " wall="
                     << elapsedSec() << " s" << endl;
}
