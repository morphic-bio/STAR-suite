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
    std::vector<uint64_t> moleculeKeys;
    std::vector<uint8_t> moleculeRegions;
    BucketMetrics metrics;
    uint64_t inputRecords = 0;
    uint64_t inputCounts = 0;
    uint32_t maxGene = 0;
    size_t tripletGroups = 0;
    std::string error;
};

uint64_t groupSortKey(const star::solo::PackedCbRecord &record)
{
    return (static_cast<uint64_t>(record.cb_index()) << 44)
         | (static_cast<uint64_t>(record.tag5()) << 39)
         | (static_cast<uint64_t>(record.gene15()) << 24)
         | record.umi24();
}

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

    // Ambiguous CBs retain their established global evidence pass. Only
    // observations chosen by that resolver enter the disjoint CB buckets.
    resolvePendingAmbiguousToHash(false);
    readFeatSum->flushBucketSegments();

    std::string storeError;
    if (!pSolo.cbBucketStore->finalize(&storeError)) {
        exitWithError("EXITING because CB bucket finalization failed: "
                          + storeError + "\n",
                      std::cerr, P.inOut->logMain,
                      EXIT_CODE_INCONSISTENT_DATA, P);
    }

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

    const uint32_t bucketCount = pSolo.cbBucketStore->bucket_count();
    std::vector<BucketResult> results(bucketCount);
    pSolo.cbBucketStore->reset_bucket_claims();
    const int tailThreads = std::max(1, P.runThreadN);
    const UMIParams correctionParams(
        pSolo.umiMinCount, pSolo.umiRatioThresh, pSolo.maxComponentSize);

#pragma omp parallel num_threads(tailThreads)
    {
        uint32_t bucket = 0;
        while (pSolo.cbBucketStore->claim_bucket(&bucket)) {
            BucketResult &out = results[bucket];
            std::vector<star::solo::PackedCbRecord> records;
            if (!pSolo.cbBucketStore->load_bucket(bucket, &records, &out.error))
                continue;
            out.inputRecords = records.size();

            // This ordering is injective over the packed key and also places
            // each (CB, tag, gene) group contiguously. Exact duplicate keys can
            // therefore be compacted in the same pass without a second sort.
            std::sort(records.begin(), records.end(),
                      [](const star::solo::PackedCbRecord &left,
                         const star::solo::PackedCbRecord &right) {
                          return groupSortKey(left) < groupSortKey(right);
                      });
            // Reproduce the old fused hash aggregation exactly, including
            // saturated read counts and probe-region conflict propagation.
            size_t compact = 0;
            for (size_t i = 0; i < records.size();) {
                star::solo::PackedCbRecord merged = records[i++];
                while (i < records.size() && records[i].key == merged.key)
                    merged.value = flexGdnaMergeValue(
                        merged.value, records[i++].value);
                out.inputCounts += flexGdnaValueCount(merged.value);
                records[compact++] = merged;
            }
            records.resize(compact);

            std::vector<star::solo::PackedCbRecord> molecules;
            molecules.reserve(records.size());
            std::vector<UMICount> counts;
            for (size_t begin = 0; begin < records.size();) {
                size_t end = begin + 1;
                while (end < records.size()
                       && sameGroup(records[begin], records[end]))
                    ++end;

                const uint32_t cb = records[begin].cb_index();
                const uint8_t tag = records[begin].tag5();
                const uint16_t gene = records[begin].gene15();
                const bool correctGroup = pSolo.umiCorrectionMode > 0
                    && cb < cbAllowed.size() && cbAllowed[cb] != 0;
                UMICorrectionResult correction;
                if (correctGroup) {
                    counts.clear();
                    counts.reserve(end - begin);
                    uint64_t groupReads = 0;
                    for (size_t i = begin; i < end; ++i) {
                        const uint32_t readCount =
                            flexGdnaValueCount(records[i].value);
                        counts.emplace_back(records[i].umi24(), readCount);
                        groupReads += readCount;
                    }
                    correction = UMICorrector::correctClique(
                        counts, correctionParams);
                    out.metrics.readsGrouped += groupReads;
                    out.metrics.readsBefore += groupReads;
                    out.metrics.readsAfter += groupReads;
                    addComponentMetrics(correction, &out.metrics);
                }

                const size_t moleculeBegin = molecules.size();
                for (size_t i = begin; i < end; ++i) {
                    star::solo::PackedCbRecord corrected = records[i];
                    if (correctGroup) {
                        const auto found =
                            correction.urToUb.find(corrected.umi24());
                        if (found != correction.urToUb.end())
                            corrected.key = packCgAggKey(
                                cb, found->second, gene, tag);
                    }
                    molecules.push_back(corrected);
                }
                std::sort(molecules.begin() + moleculeBegin, molecules.end(),
                          [](const star::solo::PackedCbRecord &left,
                             const star::solo::PackedCbRecord &right) {
                              return left.key < right.key;
                          });
                size_t write = moleculeBegin;
                for (size_t i = moleculeBegin; i < molecules.size();) {
                    star::solo::PackedCbRecord merged = molecules[i++];
                    while (i < molecules.size()
                           && molecules[i].key == merged.key)
                        merged.value = flexGdnaMergeValue(
                            merged.value, molecules[i++].value);
                    molecules[write++] = merged;
                }
                molecules.resize(write);
                begin = end;
            }

            uint32_t previousCb = UINT32_MAX;
            uint64_t previousCbTag = UINT64_MAX;
            uint32_t cell = 0;
            for (size_t begin = 0; begin < molecules.size();) {
                size_t end = begin + 1;
                while (end < molecules.size()
                       && sameGroup(molecules[begin], molecules[end]))
                    ++end;
                const uint32_t cb = molecules[begin].cb_index();
                const uint8_t tag = molecules[begin].tag5();
                const uint16_t gene = molecules[begin].gene15();

                if (cb != previousCb) {
                    out.cbIndices.push_back(cb);
                    previousCb = cb;
                }
                for (size_t i = begin; i < end; ++i) {
                    out.moleculeKeys.push_back(molecules[i].key);
                    out.moleculeRegions.push_back(static_cast<uint8_t>(
                        flexGdnaValueRegion(molecules[i].value)));
                }

                if (tag != 0 && cb < pSolo.cbWLstr.size()
                    && tag < gCanonicalTags.size()
                    && !gCanonicalTags[tag].empty()) {
                    const uint64_t cbTag =
                        (static_cast<uint64_t>(cb) << 8) | tag;
                    if (cbTag != previousCbTag) {
                        out.barcodes.push_back(
                            pSolo.cbWLstr[cb] + gCanonicalTags[tag]);
                        out.cbTagKeys.push_back(cbTag);
                        out.cellUmis.push_back(0);
                        out.cellGenes.push_back(0);
                        cell = static_cast<uint32_t>(out.barcodes.size() - 1);
                        previousCbTag = cbTag;
                    }
                    if (gene > 0 && gene <= geneIds.size()) {
                        const uint32_t moleculeCount =
                            static_cast<uint32_t>(end - begin);
                        out.triplets.push_back(
                            {cell, static_cast<uint32_t>(gene - 1), moleculeCount});
                        out.cellUmis.back() += moleculeCount;
                        ++out.cellGenes.back();
                        out.maxGene = std::max<uint32_t>(out.maxGene, gene);
                    }
                    ++out.tripletGroups;
                }
                begin = end;
            }
        }
    }

    InlineMatrixBundle inlineMatrix;
    indCB.clear();
    uint32_t maxGeneIdx = 0;
    uint64_t totalInputRecords = 0;
    uint64_t totalInputCounts = 0;
    size_t nTripletGroups = 0;
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
        nTripletGroups += part.tripletGroups;
        maxGeneIdx = std::max(maxGeneIdx, part.maxGene);

        const uint32_t cellBase = static_cast<uint32_t>(
            inlineMatrix.matrixData.barcodes.size());
        inlineMatrix.matrixData.barcodes.insert(
            inlineMatrix.matrixData.barcodes.end(),
            part.barcodes.begin(), part.barcodes.end());
        inlineMatrix.cbTagKeys.insert(
            inlineMatrix.cbTagKeys.end(),
            part.cbTagKeys.begin(), part.cbTagKeys.end());
        inlineMatrix.matrixData.nUMIperCB.insert(
            inlineMatrix.matrixData.nUMIperCB.end(),
            part.cellUmis.begin(), part.cellUmis.end());
        inlineMatrix.matrixData.nGenePerCB.insert(
            inlineMatrix.matrixData.nGenePerCB.end(),
            part.cellGenes.begin(), part.cellGenes.end());
        for (MexWriter::Triplet triplet : part.triplets) {
            triplet.cell_idx += cellBase;
            inlineMatrix.triplets.push_back(triplet);
        }
        inlineMatrix.gdnaMoleculeKeys.insert(
            inlineMatrix.gdnaMoleculeKeys.end(),
            part.moleculeKeys.begin(), part.moleculeKeys.end());
        inlineMatrix.gdnaMoleculeRegions.insert(
            inlineMatrix.gdnaMoleculeRegions.end(),
            part.moleculeRegions.begin(), part.moleculeRegions.end());
        indCB.insert(indCB.end(), part.cbIndices.begin(), part.cbIndices.end());

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

    nCB = static_cast<uint32_t>(indCB.size());
    indCBwl.assign(pSolo.cbWLsize, static_cast<uint32_t>(-1));
    for (uint32_t i = 0; i < nCB; ++i)
        indCBwl[indCB[i]] = i;

    if (maxGeneIdx > geneIds.size()) {
        for (uint32_t gene = static_cast<uint32_t>(geneIds.size());
             gene < maxGeneIdx; ++gene)
            geneIds.push_back("UNKNOWN_PROBE_" + std::to_string(gene + 1));
    }

    SampleMatrixData &matrix = inlineMatrix.matrixData;
    matrix.nCells = static_cast<uint32_t>(matrix.barcodes.size());
    matrix.nGenes = static_cast<uint32_t>(geneIds.size());
    matrix.countMatStride = 3;
    matrix.features = std::move(geneIds);
    matrix.countCellGeneUMIindex.assign(matrix.nCells + 1, 0);
    uint32_t matrixOffset = 0;
    size_t tripletIndex = 0;
    for (uint32_t cell = 0; cell < matrix.nCells; ++cell) {
        matrix.countCellGeneUMIindex[cell] = matrixOffset;
        while (tripletIndex < inlineMatrix.triplets.size()
               && inlineMatrix.triplets[tripletIndex].cell_idx == cell) {
            matrix.countCellGeneUMI.push_back(
                inlineMatrix.triplets[tripletIndex].gene_idx);
            matrix.countCellGeneUMI.push_back(
                inlineMatrix.triplets[tripletIndex].count);
            matrix.countCellGeneUMI.push_back(0);
            matrixOffset += matrix.countMatStride;
            ++tripletIndex;
        }
    }
    matrix.countCellGeneUMIindex[matrix.nCells] = matrixOffset;

    P.inOut->logMain << "[CB-BUCKET] backend="
                     << (pSolo.cbBucketStore->using_spill() ? "spill" : "ram")
                     << " transitioned="
                     << (pSolo.cbBucketStore->transitioned_to_spill() ? "yes" : "no")
                     << " streamed_records=" << totalInputRecords
                     << " aggregated_counts=" << totalInputCounts
                     << " final_molecules="
                     << inlineMatrix.gdnaMoleculeKeys.size()
                     << " buckets=" << bucketCount
                     << " tail_threads=" << tailThreads << endl;
    P.inOut->logMain << "Found " << nCB << " unique CBs and "
                     << nTripletGroups << " (CB, gene, tag) groups" << endl;
    P.inOut->logMain << "Found " << matrix.nCells
                     << " unique (CB, TAG) combinations" << endl;
    P.inOut->logMain << "  Genes: " << matrix.nGenes
                     << ", Entries: " << inlineMatrix.triplets.size() << endl;

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

    const std::string mexDir = P.outFileNamePrefix + pSolo.outFileNames[0]
        + SoloFeatureTypes::Names[featureType] + "/raw/";
    createDirectory(mexDir, P.runDirPerm, "Solo raw MEX directory", P);
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
