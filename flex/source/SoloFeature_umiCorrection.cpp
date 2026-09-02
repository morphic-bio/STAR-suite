#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "UMICorrector.h"
#include "ErrorWarning.h"
#include "hash_shims_cpp_compat.h"
#include "FlexGdna.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <iomanip>
#include <omp.h>
#include <unordered_set>

static const bool g_debugURHist = (std::getenv("STAR_DEBUG_TAG") != nullptr);

// Dump format for clique replay harness (env var: STAR_DUMP_CLIQUE_GROUPS=<path>)
static const uint32_t CLIQUE_DUMP_MAGIC   = 0x434C5155; // "CLQU"
static const uint32_t CLIQUE_DUMP_VERSION = 1;

void SoloFeature::applyCliqueCorrectionsToHash() {
    if (!readFeatSum || !readFeatSum->inlineHash_ || !umiCorrectionHash || kh_size(umiCorrectionHash) == 0) {
        return;
    }

    struct HashUpdate {
        uint64_t oldKey;
        uint64_t newKey;
        uint32_t value;
    };
    std::vector<HashUpdate> updates;

    for (khiter_t iter = kh_begin(readFeatSum->inlineHash_); iter != kh_end(readFeatSum->inlineHash_); ++iter) {
        if (!kh_exist(readFeatSum->inlineHash_, iter)) continue;

        uint64_t key = kh_key(readFeatSum->inlineHash_, iter);
        uint32_t value = kh_val(readFeatSum->inlineHash_, iter);

        khiter_t corrIt = kh_get(cg_agg, umiCorrectionHash, key);
        if (corrIt == kh_end(umiCorrectionHash)) continue;

        uint32_t correctedUmi = kh_val(umiCorrectionHash, corrIt);

        uint32_t cbIdx;
        uint16_t geneIdx;
        uint8_t tagIdx;
        unpackCgAggKey(key, &cbIdx, nullptr, &geneIdx, &tagIdx);
        uint64_t newKey = packCgAggKey(cbIdx, correctedUmi, geneIdx, tagIdx);

        updates.push_back({key, newKey, value});
    }

    for (const auto &update : updates) {
        khiter_t old_iter = kh_get(cg_agg, readFeatSum->inlineHash_, update.oldKey);
        if (old_iter != kh_end(readFeatSum->inlineHash_)) {
            kh_del(cg_agg, readFeatSum->inlineHash_, old_iter);
        }

        int absent;
        khiter_t new_iter = kh_put(cg_agg, readFeatSum->inlineHash_, update.newKey, &absent);
        if (absent) {
            kh_val(readFeatSum->inlineHash_, new_iter) = update.value;
        } else {
            kh_val(readFeatSum->inlineHash_, new_iter) = pSolo.flexMode
                ? flexGdnaMergeValue(kh_val(readFeatSum->inlineHash_, new_iter),
                                     update.value)
                : kh_val(readFeatSum->inlineHash_, new_iter) + update.value;
        }
    }
}

void SoloFeature::runCliqueCorrection() {
    if (pSolo.umiCorrectionMode == 0) return;

    if (!pSolo.inlineHashMode || !readFeatSum || !readFeatSum->inlineHash_) return;

    P.inOut->logMain << "Running clique-based UMI correction..." << endl;
    auto totalStart = std::chrono::steady_clock::now();

    // --- Phase 1: extract hash entries into a flat sortable vector ---
    // Each entry: (groupKey, umi24, count).  We sort by groupKey to process
    // groups sequentially, avoiding 234M inner-map allocations.
    struct FlatEntry {
        uint64_t groupKey;
        uint64_t cgAggKey;   // original key — needed for flat correction hash
        uint32_t umi24;
        uint32_t count;
    };

    // Pre-build allowed cbIdx set (same logic as old buildHistogramsFromHash)
    std::unordered_set<uint32_t> cbIdxAllowed;
    const bool hasAllowList = !cellsAllowSet.empty();
    if (hasAllowList) {
        for (uint32_t idx = 0; idx < static_cast<uint32_t>(pSolo.cbWLstr.size()); ++idx) {
            const std::string &wl = pSolo.cbWLstr[idx];
            if (cellsAllowSet.count(wl.substr(0, 16))) {
                cbIdxAllowed.insert(idx);
            } else if (wl.length() >= 24 && cellsAllowSet.count(wl.substr(0, 24))) {
                cbIdxAllowed.insert(idx);
            }
        }
    }

    std::vector<FlatEntry> entries;
    entries.reserve(kh_size(readFeatSum->inlineHash_));

    for (khiter_t iter = kh_begin(readFeatSum->inlineHash_); iter != kh_end(readFeatSum->inlineHash_); ++iter) {
        if (!kh_exist(readFeatSum->inlineHash_, iter)) continue;

        uint64_t key = kh_key(readFeatSum->inlineHash_, iter);
        uint32_t count = pSolo.flexMode
            ? flexGdnaValueCount(kh_val(readFeatSum->inlineHash_, iter))
            : kh_val(readFeatSum->inlineHash_, iter);

        uint32_t cbIdx, umi24;
        uint16_t geneIdx;
        uint8_t tagIdx;
        unpackCgAggKey(key, &cbIdx, &umi24, &geneIdx, &tagIdx);

        if (cbIdx >= pSolo.cbWLstr.size()) continue;
        if (hasAllowList && cbIdxAllowed.find(cbIdx) == cbIdxAllowed.end()) continue;

        // Correct within the final composite-barcode surface: CB + tag/sample + gene.
        // Do not merge UMIs across different sample tags.
        uint64_t groupKey = (static_cast<uint64_t>(cbIdx) << 24) |
                            (static_cast<uint64_t>(tagIdx) << 16) |
                            static_cast<uint64_t>(geneIdx);

        entries.push_back({groupKey, key, umi24, count});
        readsURGrouped += count;
    }

    std::sort(entries.begin(), entries.end(),
              [](const FlatEntry &a, const FlatEntry &b) { return a.groupKey < b.groupKey; });

    auto phase1End = std::chrono::steady_clock::now();
    double phase1Sec = std::chrono::duration<double>(phase1End - totalStart).count();
    P.inOut->logMain << "Clique correction phase 1 (extract+sort): "
                     << entries.size() << " entries, "
                     << std::fixed << std::setprecision(1) << phase1Sec << " sec" << endl;

    // --- Optional: dump grouped entries for external replay harness ---
    const char *dumpPath = std::getenv("STAR_DUMP_CLIQUE_GROUPS");
    if (dumpPath && dumpPath[0] != '\0') {
        FILE *df = std::fopen(dumpPath, "wb");
        if (df) {
            std::fwrite(&CLIQUE_DUMP_MAGIC, 4, 1, df);
            std::fwrite(&CLIQUE_DUMP_VERSION, 4, 1, df);
            int32_t mc = pSolo.umiMinCount;
            double rt = pSolo.umiRatioThresh;
            int32_t ms = pSolo.maxComponentSize;
            std::fwrite(&mc, sizeof(mc), 1, df);
            std::fwrite(&rt, sizeof(rt), 1, df);
            std::fwrite(&ms, sizeof(ms), 1, df);
            // Count groups
            uint64_t nGroups = 0;
            if (!entries.empty()) {
                nGroups = 1;
                for (size_t i = 1; i < entries.size(); ++i)
                    if (entries[i].groupKey != entries[i-1].groupKey) ++nGroups;
            }
            std::fwrite(&nGroups, sizeof(nGroups), 1, df);
            if (!entries.empty()) {
                size_t gStart = 0;
                for (size_t i = 1; i <= entries.size(); ++i) {
                    if (i == entries.size() || entries[i].groupKey != entries[i-1].groupKey) {
                        uint64_t gk = entries[gStart].groupKey;
                        uint32_t n = static_cast<uint32_t>(i - gStart);
                        std::fwrite(&gk, sizeof(gk), 1, df);
                        std::fwrite(&n, sizeof(n), 1, df);
                        for (size_t j = gStart; j < i; ++j) {
                            std::fwrite(&entries[j].umi24, 4, 1, df);
                            std::fwrite(&entries[j].count, 4, 1, df);
                        }
                        gStart = i;
                    }
                }
            }
            std::fclose(df);
            P.inOut->logMain << "Clique groups dumped to " << dumpPath
                             << " (" << nGroups << " groups, " << entries.size() << " entries)" << endl;
        }
    }

    // --- Phase 2: process each group's clique correction (parallel) ---
    auto phase2Start = std::chrono::steady_clock::now();

    if (umiCorrectionHash) {
        kh_destroy(cg_agg, umiCorrectionHash);
    }
    umiCorrectionHash = kh_init(cg_agg);

    // 2a. Build group boundary index: groupStarts[g] is the first entry of group g.
    //     groupStarts[nGroups] == entries.size() (sentinel).
    std::vector<size_t> groupStarts;
    if (!entries.empty()) {
        groupStarts.push_back(0);
        for (size_t i = 1; i < entries.size(); ++i) {
            if (entries[i].groupKey != entries[i-1].groupKey)
                groupStarts.push_back(i);
        }
    }
    const size_t nGroups = groupStarts.size();
    groupStarts.push_back(entries.size()); // sentinel

    UMIParams params(pSolo.umiMinCount, pSolo.umiRatioThresh, pSolo.maxComponentSize);
    const int nThreads = std::max(1, P.runThreadN);

    // Thread-local accumulators
    struct ThreadAccum {
        khash_t(cg_agg) *corrHash;
        uint64_t umisBeforeT, umisAfterT;
        uint64_t readsBeforeT, readsAfterT;
        uint32_t mergesT, componentsT, componentsCappedT, componentsBelowThresholdT;
        uint32_t maxComponentSeenT;
        uint32_t componentSizeHistT[5];
    };
    std::vector<ThreadAccum> accums(nThreads);
    for (auto &a : accums) {
        a.corrHash = kh_init(cg_agg);
        a.umisBeforeT = a.umisAfterT = 0;
        a.readsBeforeT = a.readsAfterT = 0;
        a.mergesT = a.componentsT = a.componentsCappedT = a.componentsBelowThresholdT = 0;
        a.maxComponentSeenT = 0;
        for (int k = 0; k < 5; ++k) a.componentSizeHistT[k] = 0;
    }

    // 2b. Parallel group processing
    #pragma omp parallel num_threads(nThreads)
    {
        const int tid = omp_get_thread_num();
        ThreadAccum &acc = accums[tid];
        std::vector<UMICount> localCounts;
        std::vector<uint64_t> localKeys;

        #pragma omp for schedule(dynamic, 256)
        for (size_t g = 0; g < nGroups; ++g) {
            const size_t gStart = groupStarts[g];
            const size_t gEnd   = groupStarts[g + 1];

            localCounts.clear();
            localKeys.clear();
            uint64_t groupReads = 0;
            for (size_t i = gStart; i < gEnd; ++i) {
                localCounts.emplace_back(entries[i].umi24, entries[i].count);
                localKeys.push_back(entries[i].cgAggKey);
                groupReads += entries[i].count;
            }
            acc.readsBeforeT += groupReads;
            acc.readsAfterT += groupReads;

            UMICorrectionResult result = UMICorrector::correctClique(localCounts, params);

            acc.umisBeforeT += result.uniqueUmisInput;
            acc.umisAfterT += result.uniqueUmisPostFilter - result.merges;
            acc.mergesT += result.merges;
            acc.componentsT += result.components;
            acc.componentsCappedT += result.componentsCapped;
            acc.componentsBelowThresholdT += result.componentsBelowThreshold;

            for (uint32_t compSize : result.componentSizes) {
                if (compSize == 1) acc.componentSizeHistT[0]++;
                else if (compSize == 2) acc.componentSizeHistT[1]++;
                else if (compSize == 3) acc.componentSizeHistT[2]++;
                else if (compSize == 4) acc.componentSizeHistT[3]++;
                else if (compSize > 4) acc.componentSizeHistT[4]++;
                if (compSize > acc.maxComponentSeenT) acc.maxComponentSeenT = compSize;
            }

            for (const auto &corr : result.urToUb) {
                if (corr.first == corr.second) continue;
                for (size_t i = 0; i < localCounts.size(); ++i) {
                    if (localCounts[i].ur == corr.first) {
                        int absent;
                        khiter_t ki = kh_put(cg_agg, acc.corrHash, localKeys[i], &absent);
                        kh_val(acc.corrHash, ki) = corr.second;
                    }
                }
            }
        }
    } // end omp parallel

    // 2c. Merge thread-local results
    for (const auto &acc : accums) {
        umisBeforeTotal += acc.umisBeforeT;
        umisAfterTotal += acc.umisAfterT;
        readsBeforeTotal += acc.readsBeforeT;
        readsAfterTotal += acc.readsAfterT;
        mergesTotal += acc.mergesT;
        componentsTotal += acc.componentsT;
        componentsCappedTotal += acc.componentsCappedT;
        componentsBelowThresholdTotal += acc.componentsBelowThresholdT;
        if (acc.maxComponentSeenT > maxComponentSeen) maxComponentSeen = acc.maxComponentSeenT;
        for (int k = 0; k < 5; ++k) componentSizeHist[k] += acc.componentSizeHistT[k];

        for (khiter_t it = kh_begin(acc.corrHash); it != kh_end(acc.corrHash); ++it) {
            if (!kh_exist(acc.corrHash, it)) continue;
            int absent;
            khiter_t ki = kh_put(cg_agg, umiCorrectionHash, kh_key(acc.corrHash, it), &absent);
            kh_val(umiCorrectionHash, ki) = kh_val(acc.corrHash, it);
        }
        kh_destroy(cg_agg, acc.corrHash);
    }

    auto phase2End = std::chrono::steady_clock::now();
    double phase2Sec = std::chrono::duration<double>(phase2End - phase2Start).count();
    P.inOut->logMain << "Clique correction phase 2 (group BFS): " << nGroups
                     << " groups, " << nThreads << " threads, "
                     << std::fixed << std::setprecision(1) << phase2Sec << " sec" << endl;

    // Free the flat vector now — no longer needed
    { std::vector<FlatEntry>().swap(entries); }

    P.inOut->logMain << "Clique correction complete: merges=" << mergesTotal
                     << " components=" << componentsTotal
                     << " umis_before=" << umisBeforeTotal
                     << " umis_after=" << umisAfterTotal
                     << " reads_before=" << readsBeforeTotal
                     << " reads_after=" << readsAfterTotal
                     << endl;

    // --- Phase 3: apply corrections to inline hash ---
    auto phase3Start = std::chrono::steady_clock::now();
    if (kh_size(umiCorrectionHash) > 0) {
        applyCliqueCorrectionsToHash();
    }

    // Free correction hash — no longer needed
    kh_destroy(cg_agg, umiCorrectionHash);
    umiCorrectionHash = nullptr;

    auto totalEnd = std::chrono::steady_clock::now();
    double phase3Sec = std::chrono::duration<double>(totalEnd - phase3Start).count();
    double totalSec = std::chrono::duration<double>(totalEnd - totalStart).count();
    P.inOut->logMain << "Clique correction phase 3 (apply): "
                     << std::fixed << std::setprecision(1) << phase3Sec << " sec" << endl;
    P.inOut->logMain << "Clique correction total: "
                     << std::fixed << std::setprecision(1) << totalSec << " sec" << endl;
}
