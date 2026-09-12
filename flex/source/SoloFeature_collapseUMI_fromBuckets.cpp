#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "CbBucketStore.h"
#include "ErrorWarning.h"
#include "FlexGdna.h"
#include "MexWriter.h"
#include "SampleDetector.h"
#include "TimeFunctions.h"
#include "UMICorrector.h"
#include "serviceFuns.cpp"
#include "streamFuns.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace {

struct BucketMetrics {
    uint64_t umisBefore = 0;
    uint64_t umisAfter = 0;
    uint64_t readsBefore = 0;
    uint64_t readsAfter = 0;
    uint64_t readsGrouped = 0;
    uint32_t merges = 0;
    uint32_t components = 0;
    uint32_t capped = 0;
    uint32_t below = 0;
    uint32_t maxComponent = 0;
    uint32_t componentHist[5] = {0, 0, 0, 0, 0};
};

struct BucketResult {
    std::vector<uint32_t> cbIndices;
    std::vector<std::string> barcodes;
    std::vector<uint64_t> cbTagKeys;
    std::vector<MexWriter::Triplet> triplets;
    std::vector<uint32_t> cellUmis;
    std::vector<uint32_t> cellGenes;
    std::vector<FlexGdnaCellSummary> gdnaCells;
    std::vector<FlexGdnaGeneCount> gdnaGeneCounts;
    uint64_t finalMolecules = 0;
    BucketMetrics metrics;
    uint64_t inputRecords = 0;
    uint64_t inputCounts = 0;
    uint32_t maxGene = 0;
    size_t tripletGroups = 0;
    std::string error;
};

bool sameGroup(const star::solo::PackedCbRecord &left,
               const star::solo::PackedCbRecord &right)
{
    return left.cb_index() == right.cb_index()
        && left.tag5() == right.tag5()
        && left.gene15() == right.gene15();
}

void addComponentMetrics(const UMICorrectionResult &correction,
                         BucketMetrics *metrics)
{
    metrics->umisBefore += correction.uniqueUmisInput;
    metrics->umisAfter += correction.uniqueUmisPostFilter - correction.merges;
    metrics->merges += correction.merges;
    metrics->components += correction.components;
    metrics->capped += correction.componentsCapped;
    metrics->below += correction.componentsBelowThreshold;
    for (uint32_t size : correction.componentSizes) {
        if (size == 1) ++metrics->componentHist[0];
        else if (size == 2) ++metrics->componentHist[1];
        else if (size == 3) ++metrics->componentHist[2];
        else if (size == 4) ++metrics->componentHist[3];
        else ++metrics->componentHist[4];
        metrics->maxComponent = std::max(metrics->maxComponent, size);
    }
}

} // namespace

void SoloFeature::collapseUMIall_fromBuckets()
{
    if (!readFeatSum || !readFeatSum->bucketStorageEnabled()
        || !pSolo.cbBucketStore) {
        P.inOut->logMain
            << "ERROR: collapseUMIall_fromBuckets called without a bucket store"
            << endl;
        return;
    }

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Starting bucket-parallel Flex collapse" << endl;

    // These four steps run before the parallel loop and were previously
    // untimed, which is why a ~17 s block hid between two log timestamps.
    auto preMark = std::chrono::steady_clock::now();
    auto preTick = [&preMark]() {
        const auto now = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(now - preMark).count();
        preMark = now;
        return elapsed;
    };

    // Ambiguous CBs retain their established global evidence pass. Only
    // observations chosen by that resolver enter the disjoint CB buckets.
    resolvePendingAmbiguousToHash(false);
    const double tResolveAmbig = preTick();
    readFeatSum->flushBucketSegments();
    const double tFlushSegments = preTick();

    std::string storeError;
    if (!pSolo.cbBucketStore->finalize(&storeError)) {
        exitWithError("EXITING because CB bucket finalization failed: "
                          + storeError + "\n",
                      std::cerr, P.inOut->logMain,
                      EXIT_CODE_INCONSISTENT_DATA, P);
    }

    const double tFinalizeStore = preTick();

    std::vector<std::string> geneIds;
    if (!P.pSolo.probeListPath.empty() && P.pSolo.probeListPath != "-") {
        std::ifstream probeFile(P.pSolo.probeListPath);
        std::string line;
        while (std::getline(probeFile, line)) {
            if (!line.empty() && line[0] == '#')
                continue;
            const size_t begin = line.find_first_not_of(" \t\r\n");
            const size_t end = line.find_last_not_of(" \t\r\n");
            if (begin != std::string::npos)
                geneIds.push_back(line.substr(begin, end - begin + 1));
        }
    }
    if (geneIds.empty()) {
        exitWithError("EXITING because the Flex probe list is unavailable: "
                          + P.pSolo.probeListPath + "\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    }

    const bool hasAllowList = !cellsAllowSet.empty();
    std::vector<uint8_t> cbAllowed(pSolo.cbWLsize, hasAllowList ? 0 : 1);
    if (hasAllowList) {
        for (uint32_t cb = 0; cb < pSolo.cbWLstr.size(); ++cb) {
            const std::string &wl = pSolo.cbWLstr[cb];
            if (cellsAllowSet.count(wl.substr(0, 16))
                || (wl.size() >= 24
                    && cellsAllowSet.count(wl.substr(0, 24))))
                cbAllowed[cb] = 1;
        }
    }

    const double tGeneIdsAndAllowed = preTick();
    P.inOut->logMain << "Solo timing: pre-loop: resolveAmbiguous " << tResolveAmbig
                     << " s, flushSegments " << tFlushSegments
                     << " s, finalizeStore " << tFinalizeStore
                     << " s, geneIds+allowed " << tGeneIdsAndAllowed << " s"
                     << endl << std::flush;

    const bool collectGdna = pSolo.runFlexFilter && pSolo.flexMode
        && pSolo.flexGdnaMode != ParametersSolo::FlexGdnaOff;
    const size_t gdnaGeneSlots = FlexGdnaProbeMetadata::instance().geneProbeCounts().size();
    const uint32_t bucketCount = pSolo.cbBucketStore->bucket_count();
    std::vector<BucketResult> results(bucketCount);
    pSolo.cbBucketStore->reset_bucket_claims();
    const int tailThreads = std::max(1, P.runThreadN);
    const UMIParams correctionParams(
        pSolo.umiMinCount, pSolo.umiRatioThresh, pSolo.maxComponentSize);

    // Phase timers, thread-seconds summed. The enclosing "collapse" timestamp
    // also covers the MEX write and flexfilter, so it is not a measure of this
    // loop; three optimisation attempts were aimed at the sort on the strength
    // of profile *sample* shares before anyone measured the phases directly.
    double tLoad = 0, tStream = 0, tUmi = 0;
    const auto parallelStart = std::chrono::steady_clock::now();
    auto tick = [](std::chrono::steady_clock::time_point &mark) {
        const auto now = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(now - mark).count();
        mark = now;
        return elapsed;
    };

#pragma omp parallel num_threads(tailThreads) \
    reduction(+:tLoad, tStream, tUmi)
    {
        uint32_t bucket = 0;
        while (pSolo.cbBucketStore->claim_bucket(&bucket)) {
            BucketResult &out = results[bucket];
            std::vector<std::vector<uint8_t>> segments;
            auto mark = std::chrono::steady_clock::now();
            if (!pSolo.cbBucketStore->consume_encoded_segments(
                    bucket, &segments, &out.error))
                continue;
            tLoad += tick(mark);
            size_t totalRecords = 0;
            for (const auto &segment : segments)
                totalRecords += segment.size() / star::solo::PackedCbRecord::kSerializedBytes;
            out.inputRecords = totalRecords;

            // append_segment sorted each producer-local run before publishing
            // it. Merge the runs here while other bucket threads do the same;
            // this preserves the exact global [CB, tag, gene, UMI] order while
            // avoiding another comparison sort over the complete bucket.
            const size_t sentinel = segments.size();
            std::vector<size_t> next(sentinel, 0);
            std::vector<uint64_t> currentKey(sentinel, 0);
            std::vector<star::solo::PackedCbRecord> currentRecord(sentinel);
            size_t leafCount = 1;
            while (leafCount < std::max<size_t>(1, sentinel))
                leafCount <<= 1;
            std::vector<size_t> tournament(leafCount * 2, sentinel);
            for (size_t segment = 0; segment < sentinel; ++segment) {
                if (!segments[segment].empty()) {
                    currentRecord[segment] = star::solo::PackedCbRecord::from_encoded(segments[segment].data());
                    currentKey[segment] = currentRecord[segment].group_sort_key();
                    tournament[leafCount + segment] = segment;
                }
            }
            const auto winner = [&currentKey, sentinel](size_t left,
                                                        size_t right) {
                if (left == sentinel) return right;
                if (right == sentinel) return left;
                if (currentKey[left] != currentKey[right])
                    return currentKey[left] < currentKey[right] ? left : right;
                // Exact-key values merge commutatively; a segment tie-break
                // nevertheless makes the output independent of scheduling.
                return std::min(left, right);
            };
            for (size_t node = leafCount; node-- > 1;)
                tournament[node] = winner(tournament[node * 2],
                                          tournament[node * 2 + 1]);
            // Each bucket worker still runs concurrently. Keep the complete
            // current [CB, tag, gene] group for UMI correction, then emit its
            // final counts immediately; no full merged/molecule arrays.
            std::vector<star::solo::PackedCbRecord> group;
            std::vector<UMICount> counts;
            uint32_t previousCb = UINT32_MAX;
            uint64_t previousCbTag = UINT64_MAX;
            uint32_t cell = 0;
            auto finishGroup = [&]() {
                if (group.empty()) return;
                const uint32_t cb = group.front().cb_index();
                const uint8_t tag = group.front().tag5();
                const uint16_t gene = group.front().gene15();
                const bool correctGroup = pSolo.umiCorrectionMode > 0
                    && cb < cbAllowed.size() && cbAllowed[cb] != 0;
                uint64_t groupReads = 0;
                for (const auto &record : group)
                    groupReads += flexGdnaValueCount(record.value);
                out.inputCounts += groupReads;
                if (correctGroup) {
                    counts.clear();
                    counts.reserve(group.size());
                    for (const auto &record : group)
                        counts.emplace_back(record.umi24(),
                                            flexGdnaValueCount(record.value));
                    auto umiMark = std::chrono::steady_clock::now();
                    const auto correction = UMICorrector::correctClique(
                        counts, correctionParams);
                    tUmi += tick(umiMark);
                    out.metrics.readsGrouped += groupReads;
                    out.metrics.readsBefore += groupReads;
                    out.metrics.readsAfter += groupReads;
                    addComponentMetrics(correction, &out.metrics);
                    for (auto &record : group) {
                        const auto found = correction.urToUb.find(record.umi24());
                        if (found != correction.urToUb.end())
                            record.key = packCgAggKey(cb, found->second, gene, tag);
                    }
                }
                std::sort(group.begin(), group.end(),
                          [](const star::solo::PackedCbRecord &left,
                             const star::solo::PackedCbRecord &right) {
                              return left.key < right.key;
                          });
                size_t write = 0;
                for (size_t i = 0; i < group.size();) {
                    auto merged = group[i++];
                    while (i < group.size() && group[i].key == merged.key)
                        merged.value = flexGdnaMergeValue(merged.value,
                                                          group[i++].value);
                    group[write++] = merged;
                }
                group.resize(write);
                if (cb != previousCb) {
                    out.cbIndices.push_back(cb);
                    previousCb = cb;
                }
                out.finalMolecules += group.size();
                if (tag != 0 && cb < pSolo.cbWLstr.size()
                    && tag < gCanonicalTags.size()
                    && !gCanonicalTags[tag].empty()) {
                    const uint64_t cbTag =
                        (static_cast<uint64_t>(cb) << 8) | tag;
                    if (cbTag != previousCbTag) {
                        out.barcodes.push_back(pSolo.cbWLstr[cb] + gCanonicalTags[tag]);
                        out.cbTagKeys.push_back(cbTag);
                        out.cellUmis.push_back(0);
                        out.cellGenes.push_back(0);
                        if (collectGdna) out.gdnaCells.emplace_back();
                        cell = static_cast<uint32_t>(out.barcodes.size() - 1);
                        previousCbTag = cbTag;
                    }
                    if (collectGdna) {
                        uint32_t regionCounts[4] = {0, 0, 0, 0};
                        for (const auto& record : group)
                            ++regionCounts[flexGdnaValueRegion(record.value)];
                        flexGdnaAccumulateGroup(out.gdnaCells.back(), out.gdnaGeneCounts,
                                                gene, gdnaGeneSlots, regionCounts);
                    }
                    if (gene > 0 && gene <= geneIds.size()) {
                        const uint32_t moleculeCount = static_cast<uint32_t>(group.size());
                        out.triplets.push_back(
                            {cell, static_cast<uint32_t>(gene - 1), moleculeCount});
                        out.cellUmis.back() += moleculeCount;
                        ++out.cellGenes.back();
                        out.maxGene = std::max<uint32_t>(out.maxGene, gene);
                    }
                    ++out.tripletGroups;
                }
                group.clear();
            };
            while (tournament[1] != sentinel) {
                const size_t segment = tournament[1];
                const auto record = currentRecord[segment];
                if (!group.empty() && !sameGroup(group.front(), record))
                    finishGroup();
                // Equal packed keys are adjacent in the global merge. Fold
                // them with exactly the same saturated-count/region operator.
                if (!group.empty() && group.back().key == record.key)
                    group.back().value = flexGdnaMergeValue(group.back().value,
                                                           record.value);
                else
                    group.push_back(record);
                next[segment] += star::solo::PackedCbRecord::kSerializedBytes;
                const size_t leaf = leafCount + segment;
                if (next[segment] == segments[segment].size()) {
                    tournament[leaf] = sentinel;
                    std::vector<uint8_t>().swap(segments[segment]);
                } else {
                    currentRecord[segment] = star::solo::PackedCbRecord::from_encoded(
                        segments[segment].data() + next[segment]);
                    currentKey[segment] = currentRecord[segment].group_sort_key();
                }
                for (size_t node = leaf / 2; node > 0; node /= 2)
                    tournament[node] = winner(tournament[node * 2],
                                              tournament[node * 2 + 1]);
            }
            finishGroup();
            tStream += tick(mark);
        }
    }
    P.inOut->logMain << "Solo timing: bucket loop wall "
                     << std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - parallelStart).count()
                     << " s; thread-seconds load " << tLoad
                     << ", stream_merge_correct_emit " << tStream
                     << ", umi_correct_subset " << tUmi << endl << std::flush;

    InlineMatrixBundle inlineMatrix;
    indCB.clear();
    uint32_t maxGeneIdx = 0;
    uint64_t totalInputRecords = 0;
    uint64_t totalInputCounts = 0;
    uint64_t totalFinalMolecules = 0;
    size_t nTripletGroups = 0;
    // Prefix offsets fix output order while workers fill disjoint final ranges.
    const auto fanInStart = std::chrono::steady_clock::now();
    std::vector<size_t> offCells(bucketCount + 1, 0);
    std::vector<size_t> offTriplets(bucketCount + 1, 0);
    std::vector<size_t> offGdnaCounts(bucketCount + 1, 0);
    std::vector<size_t> offCbIndices(bucketCount + 1, 0);
    for (uint32_t bucket = 0; bucket < bucketCount; ++bucket) {
        BucketResult &part = results[bucket];
        if (!part.error.empty()) {
            exitWithError("EXITING because CB bucket " + std::to_string(bucket)
                              + " could not be loaded: " + part.error + "\n",
                          std::cerr, P.inOut->logMain,
                          EXIT_CODE_INCONSISTENT_DATA, P);
        }
        totalInputRecords += part.inputRecords;
        totalInputCounts += part.inputCounts;
        totalFinalMolecules += part.finalMolecules;
        nTripletGroups += part.tripletGroups;
        maxGeneIdx = std::max(maxGeneIdx, part.maxGene);

        offCells[bucket + 1] = offCells[bucket] + part.barcodes.size();
        offTriplets[bucket + 1] = offTriplets[bucket] + part.triplets.size();
        offGdnaCounts[bucket + 1] = offGdnaCounts[bucket] + part.gdnaGeneCounts.size();
        offCbIndices[bucket + 1] = offCbIndices[bucket] + part.cbIndices.size();

        umisBeforeTotal += part.metrics.umisBefore;
        umisAfterTotal += part.metrics.umisAfter;
        readsBeforeTotal += part.metrics.readsBefore;
        readsAfterTotal += part.metrics.readsAfter;
        readsURGrouped += part.metrics.readsGrouped;
        mergesTotal += part.metrics.merges;
        componentsTotal += part.metrics.components;
        componentsCappedTotal += part.metrics.capped;
        componentsBelowThresholdTotal += part.metrics.below;
        maxComponentSeen = std::max(
            maxComponentSeen, part.metrics.maxComponent);
        for (int i = 0; i < 5; ++i)
            componentSizeHist[i] += part.metrics.componentHist[i];
    }

    SampleMatrixData &matrix = inlineMatrix.matrixData;
    matrix.countMatStride = 2;
    const size_t matrixWords = offTriplets[bucketCount] * matrix.countMatStride;
    if (matrixWords > UINT32_MAX) {
        exitWithError("EXITING because the Flex sparse matrix exceeds its offset range\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_INCONSISTENT_DATA, P);
    }
    matrix.countCellGeneUMI.resize(matrixWords);
    matrix.countCellGeneUMIindex.resize(offCells[bucketCount] + 1);
    inlineMatrix.rawMexFromCsr = true;
    inlineMatrix.matrixData.barcodes.resize(offCells[bucketCount]);
    inlineMatrix.cbTagKeys.resize(offCells[bucketCount]);
    inlineMatrix.matrixData.nUMIperCB.resize(offCells[bucketCount]);
    inlineMatrix.matrixData.nGenePerCB.resize(offCells[bucketCount]);
    inlineMatrix.gdnaCountsReady = collectGdna;
    if (collectGdna) inlineMatrix.gdnaCells.resize(offCells[bucketCount]);
    inlineMatrix.gdnaGeneCounts.resize(offGdnaCounts[bucketCount]);
    indCB.resize(offCbIndices[bucketCount]);

#pragma omp parallel for schedule(dynamic, 1) num_threads(tailThreads)
    for (int64_t bucket = 0; bucket < static_cast<int64_t>(bucketCount); ++bucket) {
        BucketResult &part = results[bucket];
        const size_t cellBase = offCells[bucket];
        std::move(part.barcodes.begin(), part.barcodes.end(),
                  inlineMatrix.matrixData.barcodes.begin() + cellBase);
        std::copy(part.cbTagKeys.begin(), part.cbTagKeys.end(),
                  inlineMatrix.cbTagKeys.begin() + cellBase);
        std::copy(part.cellUmis.begin(), part.cellUmis.end(),
                  inlineMatrix.matrixData.nUMIperCB.begin() + cellBase);
        std::copy(part.cellGenes.begin(), part.cellGenes.end(),
                  inlineMatrix.matrixData.nGenePerCB.begin() + cellBase);
        // Build the caller's storage directly in the final cell order. The
        // raw MEX writer reads this same storage, so no final COO copy exists.
        size_t at = offTriplets[bucket] * matrix.countMatStride;
        size_t entry = 0;
        for (size_t cell = 0; cell < part.barcodes.size(); ++cell) {
            matrix.countCellGeneUMIindex[cellBase + cell] = static_cast<uint32_t>(at);
            while (entry < part.triplets.size() && part.triplets[entry].cell_idx == cell) {
                const auto& triplet = part.triplets[entry++];
                matrix.countCellGeneUMI[at++] = triplet.gene_idx;
                matrix.countCellGeneUMI[at++] = triplet.count;
            }
        }
        if (collectGdna) {
            for (size_t cell = 0; cell < part.gdnaCells.size(); ++cell) {
                auto summary = part.gdnaCells[cell];
                summary.begin += offGdnaCounts[bucket];
                inlineMatrix.gdnaCells[cellBase + cell] = summary;
            }
            std::copy(part.gdnaGeneCounts.begin(), part.gdnaGeneCounts.end(),
                      inlineMatrix.gdnaGeneCounts.begin() + offGdnaCounts[bucket]);
        }
        std::copy(part.cbIndices.begin(), part.cbIndices.end(),
                  indCB.begin() + offCbIndices[bucket]);
        // Metrics were collected above; release each source as soon as its
        // disjoint output range is complete, before writer/caller allocation.
        part = BucketResult{};
    }
    matrix.countCellGeneUMIindex.back() = static_cast<uint32_t>(matrixWords);
    P.inOut->logMain << "Solo timing: bucket fan-in "
                     << std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - fanInStart).count()
                     << " s" << endl << std::flush;

    auto stepMark = std::chrono::steady_clock::now();
    auto stepTick = [&stepMark]() {
        const auto now = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(now - stepMark).count();
        stepMark = now;
        return elapsed;
    };

    nCB = static_cast<uint32_t>(indCB.size());
    indCBwl.assign(pSolo.cbWLsize, static_cast<uint32_t>(-1));
    for (uint32_t i = 0; i < nCB; ++i)
        indCBwl[indCB[i]] = i;

    if (maxGeneIdx > geneIds.size()) {
        for (uint32_t gene = static_cast<uint32_t>(geneIds.size());
             gene < maxGeneIdx; ++gene)
            geneIds.push_back("UNKNOWN_PROBE_" + std::to_string(gene + 1));
    }

    const double tIndexWl = stepTick();

    matrix.nCells = static_cast<uint32_t>(matrix.barcodes.size());
    matrix.nGenes = static_cast<uint32_t>(geneIds.size());
    matrix.features = std::move(geneIds);
    const double tMatrixSetup = stepTick();

    P.inOut->logMain << "[CB-BUCKET] backend="
                     << (pSolo.cbBucketStore->using_spill() ? "spill" : "ram")
                     << " transitioned="
                     << (pSolo.cbBucketStore->transitioned_to_spill() ? "yes" : "no")
                     << " streamed_records=" << totalInputRecords
                     << " aggregated_counts=" << totalInputCounts
                     << " final_molecules="
                     << totalFinalMolecules
                     << " buckets=" << bucketCount
                     << " tail_threads=" << tailThreads
                     << " async_merges="
                     << pSolo.cbBucketStore->async_merge_count()
                     << " async_merged_records="
                     << pSolo.cbBucketStore->async_merged_records() << endl;
    P.inOut->logMain << "Found " << nCB << " unique CBs and "
                     << nTripletGroups << " (CB, gene, tag) groups" << endl;
    P.inOut->logMain << "Found " << matrix.nCells
                     << " unique (CB, TAG) combinations" << endl;
    P.inOut->logMain << "  Genes: " << matrix.nGenes
                     << ", Entries: " << offTriplets[bucketCount] << endl;

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

    const double tPerCbVectors = stepTick();
    P.inOut->logMain << "Solo timing: post-fan-in setup: indexWl " << tIndexWl
                     << " s, matrixSetup " << tMatrixSetup
                     << " s, perCbVectors " << tPerCbVectors << " s"
                     << endl << std::flush;

    const std::string mexDir = P.outFileNamePrefix + pSolo.outFileNames[0]
        + SoloFeatureTypes::Names[featureType] + "/raw/";
    createDirectory(mexDir, P.runDirPerm, "Solo raw MEX directory", P);
    P.inOut->logMain << "Solo timing: mex dir prepared " << stepTick()
                     << " s" << endl << std::flush;
    writeMexFromInlineHashDedup(mexDir, inlineMatrix);

    if (pSolo.runFlexFilter) {
        std::string flexOutputPrefix = pSolo.flexFilterOutputPrefix;
        if (flexOutputPrefix.back() != '/')
            flexOutputPrefix += '/';
        createDirectory(flexOutputPrefix, P.runDirPerm,
                        "FlexFilter output directory", P);
        runFlexFilterInline(inlineMatrix, flexOutputPrefix);
    }

    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime)
                     << " ... Finished bucket-parallel Flex collapse" << endl;
}
