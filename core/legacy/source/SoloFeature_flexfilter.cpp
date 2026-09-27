#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "libflex/FlexFilter.h"
#include "FlexHashScreen.h"
#include "hash_shims_cpp_compat.h"
#include "Parameters.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "MexWriter.h"
#include "BorrowedBarcodeIndex.h"
#include <atomic>
#include <chrono>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <thread>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <iostream>
#include <cstdio>
#include <iomanip>

namespace {

void trim(std::string& s) {
    const char* ws = " \t\r\n";
    size_t start = s.find_first_not_of(ws);
    if (start == std::string::npos) {
        s.clear();
        return;
    }
    size_t end = s.find_last_not_of(ws);
    s = s.substr(start, end - start + 1);
}

std::vector<MexWriter::Feature> makeMexFeatures(const std::vector<std::string>& geneIds) {
    std::vector<MexWriter::Feature> features;
    features.reserve(geneIds.size());
    for (const auto& geneId : geneIds) {
        features.emplace_back(geneId, geneId, "Gene Expression");
    }
    return features;
}

} // namespace

void SoloFeature::runFlexFilterInline(
    const InlineMatrixBundle& inlineMatrix,
    const std::string& outputPrefix)
{
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Running flexfilter inline (in-memory matrix)..." << endl;

    // Load allowed tags / whitelist
    std::vector<std::string> allowedTags;
    std::vector<std::string> whitelistLabels;
    std::vector<std::string> whitelistTags;
    if (!pSolo.flexFilterAllowedTagsPath.empty()) {
        std::ifstream tagsFile(pSolo.flexFilterAllowedTagsPath);
        if (tagsFile.is_open()) {
            std::string tagLine;
            while (std::getline(tagsFile, tagLine)) {
                trim(tagLine);
                if (tagLine.empty() || tagLine[0] == '#') continue;

                std::istringstream lineStream(tagLine);
                std::string firstToken, secondToken;
                if (lineStream >> firstToken) {
                    std::string sampleLabel;
                    std::string tagSeq;
                    if (lineStream >> secondToken) {
                        sampleLabel = firstToken;
                        tagSeq = secondToken;
                    } else {
                        tagSeq = firstToken;
                    }
                    trim(sampleLabel);
                    trim(tagSeq);
                    if (!tagSeq.empty() && tagSeq.size() >= 8) {
                        tagSeq = tagSeq.substr(0, 8);
                        allowedTags.push_back(tagSeq);
                        if (!sampleLabel.empty()) {
                            whitelistLabels.push_back(sampleLabel);
                            whitelistTags.push_back(tagSeq);
                        }
                    }
                }
            }
            tagsFile.close();
        } else {
            P.inOut->logMain << "WARNING: Could not open allowed tags file: " 
                             << pSolo.flexFilterAllowedTagsPath << endl;
        }
    }

    std::vector<std::string> sampleLabels;
    std::vector<std::string> sampleTags;

    if (!whitelistLabels.empty() && whitelistLabels.size() == whitelistTags.size()) {
        sampleLabels = whitelistLabels;
        sampleTags = whitelistTags;
        P.inOut->logMain << "  Loaded " << sampleLabels.size() << " allowed tags with sample labels" << endl;
    } else if (!whitelistLabels.empty()) {
        P.inOut->logMain << "WARNING: Mismatch between sample labels and tags in "
                         << pSolo.flexFilterAllowedTagsPath << " (labels=" << whitelistLabels.size()
                         << ", tags=" << whitelistTags.size() << "); falling back to auto-derived sample names" << endl;
    }

    if (sampleLabels.empty()) {
        FlexFilter::deriveSampleWhitelist(
            inlineMatrix.matrixData.barcodes,
            sampleLabels,
            sampleTags);
    }

    if (!allowedTags.empty() && whitelistLabels.empty()) {
        std::unordered_set<std::string> allowedSet(allowedTags.begin(), allowedTags.end());
        std::vector<std::string> filteredLabels;
        std::vector<std::string> filteredTags;
        for (size_t i = 0; i < sampleTags.size(); ++i) {
            if (allowedSet.count(sampleTags[i])) {
                filteredLabels.push_back(sampleLabels[i]);
                filteredTags.push_back(sampleTags[i]);
            }
        }
        sampleLabels = std::move(filteredLabels);
        sampleTags = std::move(filteredTags);
        P.inOut->logMain << "Filtered to " << sampleLabels.size() << " allowed tags" << endl;
    }

    if (sampleLabels.empty()) {
        std::ostringstream errMsg;
        errMsg << "ERROR: No samples to process (check allowed tags filter)";
        P.inOut->logMain << errMsg.str() << endl;
        return;
    }

    // Configure FlexFilter inputs
    FlexFilter::MemoryInputs mem;
    mem.matrixData = inlineMatrix.matrixData;
    mem.observedBarcodes = inlineMatrix.matrixData.barcodes;
    mem.sampleLabels = sampleLabels;
    mem.sampleTags = sampleTags;

    FlexFilter::Config config;
    config.tagAwareCaller = pSolo.flexFilterCallerMode == "tag-aware";
    config.mitochondrialGenesPath = pSolo.cellFilterMitochondrialGenes;
    config.totalThreads = static_cast<uint32_t>(std::max(1, P.runThreadN));
    // Handle per-tag mode: multiply by number of tags
    if (pSolo.flexFilterExpectedPerTagMode) {
        config.totalExpectedCells = pSolo.flexFilterTotalExpected * static_cast<uint32_t>(sampleTags.size());
        P.inOut->logMain << "FlexFilter: Using per-tag expected cells: " << pSolo.flexFilterTotalExpected 
                         << " x " << sampleTags.size() << " tags = " << config.totalExpectedCells << " total\n";
    } else {
        config.totalExpectedCells = pSolo.flexFilterTotalExpected;
    }
    FlexFilter::populateConfigWithDefaults(config);
    if (pSolo.flexFilterEdNiters > 0) {
        config.emptydropsParams.simN = pSolo.flexFilterEdNiters;
    } else {
        config.emptydropsParams.simN = config.tagAwareCaller ? 100000 : 10000;
    }
    if (pSolo.flexFilterEdFdrThreshold > 0.0) {
        config.emptydropsParams.FDR = pSolo.flexFilterEdFdrThreshold;
    } else {
        config.emptydropsParams.FDR = config.tagAwareCaller ? 0.01 : 0.001;
    }
    // Simple EmptyDrops parameters (formerly OrdMag)
    if (pSolo.flexFilterOrdmagNsamples > 0) {
        config.simpleEmptyDropsParams.nExpectedCells = pSolo.flexFilterOrdmagNsamples;
    }
    if (pSolo.flexFilterOrdmagUmiMin > 0) {
        config.simpleEmptyDropsParams.umiMin = static_cast<uint32_t>(pSolo.flexFilterOrdmagUmiMin);
    }
    if (pSolo.flexFilterOrdmagTargetPct > 0.0) {
        config.simpleEmptyDropsParams.maxPercentile = pSolo.flexFilterOrdmagTargetPct;
    }
    
    // EmptyDrops parameters
    if (pSolo.flexFilterEdLower > 0) {
        config.emptydropsParams.indMin = pSolo.flexFilterEdLower;
    }
    if (pSolo.flexFilterEdMaxTotalBuckets > 0) {
        config.emptydropsParams.maxTotalBuckets = pSolo.flexFilterEdMaxTotalBuckets;
    }
    
    // Occupancy parameters
    if (pSolo.flexFilterTotalPartitions > 0) {
        config.totalPartitions = pSolo.flexFilterTotalPartitions;
    }
    if (pSolo.flexFilterRecoveryFactor > 0.0) {
        config.recoveryFactor = pSolo.flexFilterRecoveryFactor;
    }
    if (pSolo.flexFilterOccupancyPercentile > 0.0) {
        config.occupancyPercentile = pSolo.flexFilterOccupancyPercentile;
    }
    if (pSolo.flexFilterLowUmiThreshold > 0) {
        config.lowUMIThreshold = pSolo.flexFilterLowUmiThreshold;
    }
    
    // Simple EmptyDrops fallback configuration
    config.useSimpleEmptyDrops = pSolo.flexFilterUseSimpleED;
    config.simpleEDMinRescues = pSolo.flexFilterSimpleEDMinRescues;
    config.simpleEDMinAmbient = pSolo.flexFilterSimpleEDMinAmbient;
    config.simpleEDMinCandidates = pSolo.flexFilterSimpleEDMinCandidates;
    // If force-enabled, set disabled=false
    if (config.useSimpleEmptyDrops) {
        config.simpleEmptyDropsParams.disabled = false;
    }
    
    // Debug and testing flags
    config.debugTagLog = pSolo.flexFilterDebugTagLog;
    config.debugOutputDir = pSolo.flexFilterDebugOutputDir;
    config.disableOccupancyFilter = pSolo.flexFilterDisableOccupancy;
    config.enableInvariantChecks = pSolo.flexFilterInvariantChecks;
    
    // Output options
    config.keepCBTag = config.tagAwareCaller || pSolo.flexFilterKeepCBTag;
    if (config.tagAwareCaller) {
        config.simpleEmptyDropsParams.maxThreads = pSolo.cellFilterBootstrapThreads > 0
            ? pSolo.cellFilterBootstrapThreads : config.totalThreads;
        config.emptydropsParams.mcThreads = config.totalThreads;
        P.inOut->logMain << "Flex cell caller: tag-aware, grouped sample labels, full CB16+TAG8; "
            << "bootstrapStreams=" << config.simpleEmptyDropsParams.maxThreads
            << " totalCallerThreads=" << config.totalThreads
            << " simulations=" << config.emptydropsParams.simN
            << " FDR=" << config.emptydropsParams.FDR << "\n";
    }

    createDirectory(outputPrefix, P.runDirPerm, "FlexFilter output directory", P);

    FlexFilter filter;
    FlexFilter::Outputs outputs;
    int result = filter.runFromMemory(mem, &outputs, config);

    time(&rawTime);
    if (result != 0) {
        std::ostringstream errMsg;
        errMsg << "ERROR: FlexFilter pipeline failed with code " << result;
        P.inOut->logMain << timeMonthDayTime(rawTime) << " " << errMsg.str() << endl;
        if (pSolo.flexFilterFatalOnError) {
            exitWithError(errMsg.str(), std::cerr, P.inOut->logMain, EXIT_CODE_RUNTIME, P);
        } else {
            P.inOut->logMain << "  Continuing despite flexfilter failure (use --soloFlexFatalOnError yes to fail-fast)" << endl;
        }
        return;
    }

    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Flexfilter pipeline complete" << endl;
    P.inOut->logMain << "  Processed " << outputs.tagResults.size() << " sample groups" << endl;

    std::cout << "FlexFilter completed successfully\n";
    std::cout << "  Processed " << outputs.tagResults.size() << " sample groups\n";
    std::cout << "Writing per-sample MEX outputs...\n";

    BorrowedBarcodeIndex barcodeToIdx(inlineMatrix.matrixData.nCells);
    for (uint32_t idx = 0; idx < inlineMatrix.matrixData.nCells; ++idx) {
        barcodeToIdx.insert(inlineMatrix.matrixData.barcodes[idx], idx);
    }

    auto printTagLog = [&](const FlexFilter::Outputs::TagResults& tagResult, const std::string& label){
        P.inOut->logMain << "  [" << label << "] Retain=" << tagResult.nRetainWindow
                         << " SimpleED=" << tagResult.nSimpleCells
                         << " TailTested=" << tagResult.nTailTested
                         << " ED_Pass=" << (tagResult.nSimplePassers + tagResult.nTailPassers)
                         << " OccRemoved=" << tagResult.occupancyRemoved
                         << " Final=" << tagResult.passingBarcodes.size()
                         << " Expected=" << tagResult.expectedCells << endl;
    };

    std::vector<MexWriter::Feature> mexFeatures = makeMexFeatures(inlineMatrix.matrixData.features);
    std::vector<int32_t> exportGeneIndex(mexFeatures.size());
    std::iota(exportGeneIndex.begin(), exportGeneIndex.end(), 0);
    if (!pSolo.flexFilteredGeneList.empty() && pSolo.flexFilteredGeneList != "-") {
        std::ifstream allowFile(pSolo.flexFilteredGeneList);
        if (!allowFile) exitWithError("Cannot open --soloFlexFilteredGeneList\n", std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        std::unordered_set<std::string> allowed;
        std::string gene;
        while (std::getline(allowFile, gene)) {
            trim(gene);
            if (!gene.empty() && gene[0] != '#') allowed.insert(gene);
        }
        std::vector<MexWriter::Feature> selected;
        std::fill(exportGeneIndex.begin(), exportGeneIndex.end(), -1);
        for (size_t i = 0; i < mexFeatures.size(); ++i) {
            if (allowed.erase(inlineMatrix.matrixData.features[i])) {
                exportGeneIndex[i] = selected.size();
                selected.push_back(mexFeatures[i]);
            }
        }
        if (selected.empty() || !allowed.empty())
            exitWithError("Empty or unmatched gene IDs in --soloFlexFilteredGeneList\n", std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
        mexFeatures.swap(selected);
        P.inOut->logMain << "Flex feature universes: calling=" << exportGeneIndex.size()
            << " filteredExport=" << mexFeatures.size() << "\n";
    }

    std::string summaryPath = outputPrefix;
    if (!summaryPath.empty() && summaryPath.back() != '/') {
        summaryPath += '/';
    }
    summaryPath += "flexfilter_summary.tsv";
    std::ofstream summaryFile(summaryPath);
    if (summaryFile.is_open()) {
        summaryFile << "Sample\tExpected\tRetain\tSimple_ED\tTail_Tested\tED_Pass\tOcc_Rem\tFinal\tTotal_UMIs\n";
    }

    std::cout << "\nSummary (saved to " << summaryPath << "):\n";
    std::printf("%-15s %10s %8s %10s %12s %10s %10s %8s %14s\n",
           "Sample", "Expected", "Retain", "Simple_ED",
           "Tail_Tested", "ED_Pass", "Occ_Rem", "Final", "Total_UMIs");

    uint32_t totalExpected = 0;
    uint32_t totalRetain = 0;
    uint32_t totalSimpleED = 0;
    uint32_t totalTailTested = 0;
    uint32_t totalEDPass = 0;
    uint32_t totalOccRemoved = 0;
    uint32_t totalFinal = 0;
    uint64_t totalUMIs = 0;

    uint32_t stride = inlineMatrix.matrixData.countMatStride;

    struct SampleMexResult {
        bool skipped = true;
        int writeResult = -1;
        uint32_t finalCells = 0;
        size_t entries = 0;
        uint64_t sampleUMI = 0;
    };
    const size_t nSamples = outputs.tagResults.size();
    std::vector<SampleMexResult> sampleResults(nSamples);
    std::vector<std::string> samplePrefixes(nSamples);
    std::vector<uint8_t> hasMappedBarcode(nSamples, 0);
    const auto sampleMexStart = std::chrono::steady_clock::now();

    // Directory creation writes to the shared STAR log, so keep it ordered and
    // outside the worker pool. The expensive matrix construction and output
    // below touch only per-sample state and disjoint files.
    for (size_t sample = 0; sample < nSamples; ++sample) {
        const auto& tagResult = outputs.tagResults[sample];
        for (const auto& bc : tagResult.passingBarcodes) {
            if (barcodeToIdx.find(bc) != UINT32_MAX) {
                hasMappedBarcode[sample] = 1;
                break;
            }
        }
        if (!hasMappedBarcode[sample])
            continue;
        std::string& samplePrefix = samplePrefixes[sample];
        samplePrefix = outputPrefix;
        if (!samplePrefix.empty() && samplePrefix.back() != '/')
            samplePrefix += '/';
        samplePrefix += tagResult.sampleLabel + "/Gene/filtered/";
        createDirectory(samplePrefix, P.runDirPerm,
                        "FlexFilter filtered MEX directory", P);
    }

    const unsigned int outputThreads =
        static_cast<unsigned int>(std::max(1, P.runThreadN));
    // A few concurrent matrices expose sample-level parallelism, while a small
    // worker group inside each writer keeps large samples from becoming a new
    // serial tail. The product never exceeds the STAR thread budget.
    const unsigned int desiredThreadsPerMatrix = 4u;
    const unsigned int sampleWorkers = nSamples == 0u ? 0u : std::min(
        static_cast<unsigned int>(nSamples),
        std::max(1u, outputThreads / desiredThreadsPerMatrix));
    const unsigned int matrixThreads = sampleWorkers == 0u
        ? 1u : std::max(1u, outputThreads / sampleWorkers);
    std::atomic<size_t> nextSample(0);
    auto writeSample = [&]() {
        for (;;) {
            const size_t sample = nextSample.fetch_add(1, std::memory_order_relaxed);
            if (sample >= nSamples)
                break;
            if (!hasMappedBarcode[sample])
                continue;
            const auto& tagResult = outputs.tagResults[sample];
            SampleMexResult& sampleResult = sampleResults[sample];

            std::unordered_map<uint32_t, uint32_t> oldToNew;
            std::vector<std::string> filteredBarcodes;
            filteredBarcodes.reserve(tagResult.passingBarcodes.size());
            for (const auto& bc : tagResult.passingBarcodes) {
                const uint32_t oldIdx = barcodeToIdx.find(bc);
                if (oldIdx == UINT32_MAX)
                    continue;
                if (oldToNew.find(oldIdx) != oldToNew.end())
                    continue;
                uint32_t newIdx = static_cast<uint32_t>(filteredBarcodes.size());
                oldToNew[oldIdx] = newIdx;
                filteredBarcodes.push_back(bc);
            }

            if (filteredBarcodes.empty())
                continue;

            std::vector<MexWriter::Triplet> filteredTriplets;
            filteredTriplets.reserve(filteredBarcodes.size() * 8);
            for (const auto& kv : oldToNew) {
                uint32_t oldIdx = kv.first;
                uint32_t newIdx = kv.second;
                uint32_t start = inlineMatrix.matrixData.countCellGeneUMIindex[oldIdx];
                uint32_t end = inlineMatrix.matrixData.countCellGeneUMIindex[oldIdx + 1];
                for (uint32_t ptr = start; ptr < end; ptr += stride) {
                    uint32_t geneIdx = inlineMatrix.matrixData.countCellGeneUMI[ptr];
                    uint32_t count = inlineMatrix.matrixData.countCellGeneUMI[ptr + 1];
                    if (count == 0)
                        continue;
                    if (exportGeneIndex[geneIdx] < 0) continue;
                    filteredTriplets.push_back({newIdx, static_cast<uint32_t>(exportGeneIndex[geneIdx]), count});
                }
            }

            const int cbLen = config.keepCBTag ? -1 : 16;
            sampleResult.writeResult = MexWriter::writeMex(
                samplePrefixes[sample], filteredBarcodes, mexFeatures,
                filteredTriplets, cbLen, matrixThreads);
            sampleResult.skipped = false;
            sampleResult.finalCells =
                static_cast<uint32_t>(filteredBarcodes.size());
            sampleResult.entries = filteredTriplets.size();

            // Preserve the established summary accounting, including any
            // duplicate barcode entries in the filter result.
            for (const auto& bc : tagResult.passingBarcodes) {
                const uint32_t cell = barcodeToIdx.find(bc);
                if (cell != UINT32_MAX)
                    sampleResult.sampleUMI += inlineMatrix.matrixData.nUMIperCB[cell];
            }
        }
    };

    std::vector<std::thread> outputWorkers;
    if (sampleWorkers > 0u) {
        outputWorkers.reserve(sampleWorkers - 1u);
        for (unsigned int worker = 1; worker < sampleWorkers; ++worker)
            outputWorkers.emplace_back(writeSample);
        writeSample();
        for (std::thread& worker : outputWorkers)
            worker.join();
    }
    P.inOut->logMain << "Solo timing: per-sample MEX "
                     << std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - sampleMexStart).count()
                     << " s, sample_workers=" << sampleWorkers
                     << ", matrix_threads=" << matrixThreads
                     << endl << std::flush;

    // Logs and summaries remain in whitelist order, independent of scheduling.
    for (size_t sample = 0; sample < nSamples; ++sample) {
        const auto& tagResult = outputs.tagResults[sample];
        const SampleMexResult& sampleResult = sampleResults[sample];
        printTagLog(tagResult, tagResult.sampleLabel);
        if (sampleResult.skipped) {
            P.inOut->logMain << "  Skipping " << tagResult.sampleLabel
                             << " (no passing barcodes mapped)\n";
            continue;
        }
        if (sampleResult.writeResult != 0) {
            std::cerr << "  ERROR: MexWriter failed for " << tagResult.sampleLabel
                      << " (barcodes=" << sampleResult.finalCells
                      << ", entries=" << sampleResult.entries << ")\n";
        } else {
            P.inOut->logMain << "  " << tagResult.sampleLabel << " ("
                             << tagResult.tag << "): "
                             << sampleResult.finalCells << " cells, "
                             << sampleResult.entries << " entries" << endl;
        }

        // Summary statistics
        uint32_t retainWindow = tagResult.nRetainWindow;
        uint32_t simpleED = tagResult.nSimpleCells;
        uint32_t tailTested = tagResult.nTailTested;
        uint32_t edPass = tagResult.nSimplePassers + tagResult.nTailPassers;
        uint32_t occRemoved = tagResult.occupancyRemoved;
        uint32_t finalCells = sampleResult.finalCells;
        uint64_t sampleUMI = sampleResult.sampleUMI;

        std::printf("%-15s %10u %8u %10u %12u %10u %10u %8u %14lu\n",
               tagResult.sampleLabel.c_str(),
               tagResult.expectedCells,
               retainWindow,
               simpleED,
               tailTested,
               edPass,
               occRemoved,
               finalCells,
               sampleUMI);

        if (summaryFile.is_open()) {
            summaryFile << tagResult.sampleLabel << '\t'
                        << tagResult.expectedCells << '\t'
                        << retainWindow << '\t'
                        << simpleED << '\t'
                        << tailTested << '\t'
                        << edPass << '\t'
                        << occRemoved << '\t'
                        << finalCells << '\t'
                        << sampleUMI << '\n';
        }

        totalExpected += tagResult.expectedCells;
        totalRetain += retainWindow;
        totalSimpleED += simpleED;
        totalTailTested += tailTested;
        totalEDPass += edPass;
        totalOccRemoved += occRemoved;
        totalFinal += finalCells;
        totalUMIs += sampleUMI;
    }

    std::printf("%-15s %10u %8u %10u %12u %10u %10u %8u %14lu\n",
           "TOTAL",
           totalExpected,
           totalRetain,
           totalSimpleED,
           totalTailTested,
           totalEDPass,
           totalOccRemoved,
           totalFinal,
           totalUMIs);

    if (summaryFile.is_open()) {
        summaryFile << "TOTAL\t"
                    << totalExpected << '\t'
                    << totalRetain << '\t'
                    << totalSimpleED << '\t'
                    << totalTailTested << '\t'
                    << totalEDPass << '\t'
                    << totalOccRemoved << '\t'
                    << totalFinal << '\t'
                    << totalUMIs << '\n';
        summaryFile.close();
    }

}
