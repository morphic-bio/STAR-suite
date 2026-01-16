#include "SoloFeature.h"
#include "libflex/FlexFilter.h"
#include "Parameters.h"
#include "TimeFunctions.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include "MexWriter.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <iostream>
#include <cstdio>

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
        config.emptydropsParams.simN = 10000;
    }
    if (pSolo.flexFilterEdFdrThreshold > 0.0) {
        config.emptydropsParams.FDR = pSolo.flexFilterEdFdrThreshold;
    } else {
        config.emptydropsParams.FDR = 0.001;
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
    config.keepCBTag = pSolo.flexFilterKeepCBTag;

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
    P.inOut->logMain << "  Processed " << outputs.tagResults.size() << " tags" << endl;

    std::cout << "FlexFilter completed successfully\n";
    std::cout << "  Processed " << outputs.tagResults.size() << " tags\n";
    std::cout << "Writing per-sample MEX outputs...\n";

    std::unordered_map<std::string, uint32_t> barcodeToIdx;
    barcodeToIdx.reserve(inlineMatrix.matrixData.nCells * 2);
    for (uint32_t idx = 0; idx < inlineMatrix.matrixData.nCells; ++idx) {
        barcodeToIdx[inlineMatrix.matrixData.barcodes[idx]] = idx;
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

    for (const auto& tagResult : outputs.tagResults) {
        printTagLog(tagResult, tagResult.sampleLabel);
        std::unordered_map<uint32_t, uint32_t> oldToNew;
        std::vector<std::string> filteredBarcodes;
        filteredBarcodes.reserve(tagResult.passingBarcodes.size());
        for (const auto& bc : tagResult.passingBarcodes) {
            auto it = barcodeToIdx.find(bc);
            if (it == barcodeToIdx.end()) {
                continue;
            }
            uint32_t oldIdx = it->second;
            if (oldToNew.find(oldIdx) != oldToNew.end()) {
                continue;
            }
            uint32_t newIdx = filteredBarcodes.size();
            oldToNew[oldIdx] = newIdx;
            filteredBarcodes.push_back(bc);
        }

        if (filteredBarcodes.empty()) {
            P.inOut->logMain << "  Skipping " << tagResult.sampleLabel << " (no passing barcodes mapped)\n";
            continue;
        }

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
                if (count == 0) continue;
                filteredTriplets.push_back({newIdx, geneIdx, count});
            }
        }

        std::string samplePrefix = outputPrefix;
        if (!samplePrefix.empty() && samplePrefix.back() != '/') {
            samplePrefix += '/';
        }
        samplePrefix += tagResult.sampleLabel + "/Gene/filtered/";
        createDirectory(samplePrefix, P.runDirPerm, "FlexFilter filtered MEX directory", P);

        // Per-sample MEX: strip sample tag from barcodes (16bp output) unless keepCBTag is set
        int cb_len = config.keepCBTag ? -1 : 16;
        int writeResult = MexWriter::writeMex(samplePrefix, filteredBarcodes, mexFeatures, filteredTriplets, cb_len);
        if (writeResult != 0) {
            std::cerr << "  ERROR: MexWriter failed for " << tagResult.sampleLabel
                      << " (barcodes=" << filteredBarcodes.size()
                      << ", entries=" << filteredTriplets.size() << ")\n";
        } else {
            P.inOut->logMain << "  " << tagResult.sampleLabel << " (" << tagResult.tag << "): "
                             << filteredBarcodes.size() << " cells, "
                             << filteredTriplets.size() << " entries" << endl;
        }

        // Summary statistics
        uint32_t retainWindow = tagResult.nRetainWindow;
        uint32_t simpleED = tagResult.nSimpleCells;  // Simple EmptyDrops cells (fallback)
        uint32_t tailTested = tagResult.nTailTested;
        uint32_t edPass = tagResult.nSimplePassers + tagResult.nTailPassers;
        uint32_t occRemoved = tagResult.occupancyRemoved;
        uint32_t finalCells = static_cast<uint32_t>(filteredBarcodes.size());

        uint64_t sampleUMI = 0;
        for (const auto& bc : tagResult.passingBarcodes) {
            auto it = barcodeToIdx.find(bc);
            if (it != barcodeToIdx.end()) {
                sampleUMI += inlineMatrix.matrixData.nUMIperCB[it->second];
            }
        }

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
