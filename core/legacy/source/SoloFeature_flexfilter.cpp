#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "libflex/FlexFilter.h"
#include "FlexGdna.h"
#include "FlexHashScreen.h"
#include "hash_shims_cpp_compat.h"
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

std::string jsonEscape(const std::string& value) {
    std::ostringstream out;
    for (const unsigned char c : value) {
        switch (c) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                        << static_cast<unsigned int>(c) << std::dec;
                } else {
                    out << static_cast<char>(c);
                }
        }
    }
    return out.str();
}

void writeGdnaJson(const std::string& path,
                   const std::string& scope,
                   const FlexGdnaEstimate& estimate) {
    std::ofstream out(path.c_str());
    if (!out.is_open())
        return;

    out << std::setprecision(17)
        << "{\n"
        << "  \"scope\": \"" << jsonEscape(scope) << "\",\n"
        << "  \"valid\": " << (estimate.valid ? "true" : "false") << ",\n"
        << "  \"status\": \"" << jsonEscape(estimate.status) << "\",\n";
    if (estimate.valid) {
        out << "  \"estimated_gdna_content\": " << estimate.estimatedGdnaFraction << ",\n"
            << "  \"estimated_gdna_unspliced_threshold\": " << estimate.threshold << ",\n"
            << "  \"estimated_gdna_per_probe\": " << estimate.estimatedGdnaPerProbe << ",\n"
            << "  \"model_constant\": " << estimate.modelConstant << ",\n"
            << "  \"model_slope\": " << estimate.modelSlope << ",\n"
            << "  \"model_critical_point\": " << estimate.modelCriticalPoint << ",\n"
            << "  \"model_rss\": " << estimate.modelRss << ",\n";
    } else {
        out << "  \"estimated_gdna_content\": null,\n"
            << "  \"estimated_gdna_unspliced_threshold\": null,\n"
            << "  \"estimated_gdna_per_probe\": null,\n"
            << "  \"model_constant\": null,\n"
            << "  \"model_slope\": null,\n"
            << "  \"model_critical_point\": null,\n"
            << "  \"model_rss\": null,\n";
    }
    out << "  \"control_genes\": " << estimate.controlGenes << ",\n"
        << "  \"gene_assigned_filtered_molecules\": "
        << estimate.totalFilteredMolecules << ",\n"
        << "  \"classified_filtered_molecules\": " << estimate.classifiedMolecules << ",\n"
        << "  \"unknown_region_filtered_molecules\": " << estimate.unknownMolecules << ",\n"
        << "  \"conflicting_region_filtered_molecules\": " << estimate.conflictingMolecules
        << ",\n"
        << "  \"gene_unassigned_filtered_molecules\": " << estimate.unassignedMolecules
        << "\n}\n";
}

struct GdnaMoleculeBucket {
    std::vector<FlexGdnaGeneMoleculeCounts> genes;
    uint64_t classified = 0;
    uint64_t unknown = 0;
    uint64_t conflicting = 0;
    uint64_t unassigned = 0;
};

void addGdnaMolecule(GdnaMoleculeBucket& bucket,
                     uint16_t gene,
    FlexGdnaRegion region) {
    if (gene == 0 || gene >= bucket.genes.size()) {
        ++bucket.unassigned;
        return;
    }
    if (region == FlexGdnaSpliced) {
        ++bucket.genes[gene].spliced;
        ++bucket.classified;
    } else if (region == FlexGdnaUnspliced) {
        ++bucket.genes[gene].unspliced;
        ++bucket.classified;
    } else if (region == FlexGdnaConflicting) {
        ++bucket.conflicting;
    } else {
        ++bucket.unknown;
    }
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

    // Cell Ranger's Flex gDNA diagnostic is defined on final molecule
    // families in filtered cells. The aggregate hash is still live here, so
    // scan it once and avoid requiring a BAM or materialized molecule table.
    if (pSolo.flexMode
        && pSolo.flexGdnaMode != ParametersSolo::FlexGdnaOff) {
        const FlexGdnaProbeMetadata& metadata = FlexGdnaProbeMetadata::instance();
        const size_t nGeneSlots = metadata.geneProbeCounts().size();
        std::vector<GdnaMoleculeBucket> buckets(outputs.tagResults.size());
        GdnaMoleculeBucket libraryBucket;
        for (auto& bucket : buckets)
            bucket.genes.resize(nGeneSlots);
        libraryBucket.genes.resize(nGeneSlots);

        bool identityComplete =
            inlineMatrix.cbTagKeys.size() == inlineMatrix.matrixData.nCells;
        std::vector<int32_t> sampleByCell(inlineMatrix.matrixData.nCells, -1);
        for (size_t sample = 0; sample < outputs.tagResults.size(); ++sample) {
            for (const std::string& barcode : outputs.tagResults[sample].passingBarcodes) {
                const auto it = barcodeToIdx.find(barcode);
                if (it == barcodeToIdx.end())
                    continue;
                int32_t& assignment = sampleByCell[it->second];
                if (assignment >= 0 && assignment != static_cast<int32_t>(sample)) {
                    assignment = -2;
                    identityComplete = false;
                } else if (assignment != -2) {
                    assignment = static_cast<int32_t>(sample);
                }
            }
        }

        std::unordered_map<uint64_t, uint32_t> cellByCbTag;
        if (identityComplete) {
            cellByCbTag.reserve(inlineMatrix.cbTagKeys.size() * 2u);
            for (uint32_t cell = 0; cell < inlineMatrix.cbTagKeys.size(); ++cell) {
                if (!cellByCbTag.emplace(inlineMatrix.cbTagKeys[cell], cell).second) {
                    identityComplete = false;
                    break;
                }
            }
        }

        if (identityComplete && readFeatSum != nullptr && readFeatSum->inlineHash_ != nullptr) {
            khash_t(cg_agg)* hash = readFeatSum->inlineHash_;
            for (khiter_t iter = kh_begin(hash); iter != kh_end(hash); ++iter) {
                if (!kh_exist(hash, iter))
                    continue;
                const uint64_t key = kh_key(hash, iter);
                const uint32_t cb = static_cast<uint32_t>((key >> 44) & 0xFFFFFu);
                const uint8_t tag = static_cast<uint8_t>(key & 0x1Fu);
                const auto cellIt =
                    cellByCbTag.find((static_cast<uint64_t>(cb) << 8) | tag);
                if (cellIt == cellByCbTag.end())
                    continue;
                const int32_t sample = sampleByCell[cellIt->second];
                if (sample < 0)
                    continue;

                const uint16_t gene =
                    static_cast<uint16_t>((key >> 5) & 0x7FFFu);
                const FlexGdnaRegion region =
                    flexGdnaValueRegion(kh_val(hash, iter));
                addGdnaMolecule(buckets[static_cast<size_t>(sample)], gene, region);
                addGdnaMolecule(libraryBucket, gene, region);
            }
        } else if (identityComplete
                   && inlineMatrix.gdnaMoleculeKeys.size()
                        == inlineMatrix.gdnaMoleculeRegions.size()
                   && !inlineMatrix.gdnaMoleculeKeys.empty()) {
            for (size_t molecule = 0;
                 molecule < inlineMatrix.gdnaMoleculeKeys.size(); ++molecule) {
                const uint64_t key = inlineMatrix.gdnaMoleculeKeys[molecule];
                const uint32_t cb = static_cast<uint32_t>((key >> 44) & 0xFFFFFu);
                const uint8_t tag = static_cast<uint8_t>(key & 0x1Fu);
                const auto cellIt = cellByCbTag.find(
                    (static_cast<uint64_t>(cb) << 8) | tag);
                if (cellIt == cellByCbTag.end())
                    continue;
                const int32_t sample = sampleByCell[cellIt->second];
                if (sample < 0)
                    continue;
                const uint16_t gene =
                    static_cast<uint16_t>((key >> 5) & 0x7FFFu);
                const FlexGdnaRegion region = static_cast<FlexGdnaRegion>(
                    inlineMatrix.gdnaMoleculeRegions[molecule]);
                addGdnaMolecule(buckets[static_cast<size_t>(sample)], gene, region);
                addGdnaMolecule(libraryBucket, gene, region);
            }
        } else {
            identityComplete = false;
        }

        const bool cacheComplete =
            !pSolo.hashScreenEnabled
            || FlexHashScreenCache::instance().hasRegionMetadata();
        const bool diagnosticComplete =
            pSolo.flexGdnaReady && cacheComplete && identityComplete;
        const char* unavailableStatus = !pSolo.flexGdnaReady
            ? "probe_region_metadata_unavailable"
            : (!cacheComplete ? "cache_region_metadata_unavailable"
                              : "filtered_barcode_identity_unavailable");

        std::vector<FlexGdnaEstimate> estimates;
        estimates.reserve(buckets.size() + 1u);
        for (const GdnaMoleculeBucket& bucket : buckets) {
            FlexGdnaEstimate estimate = flexGdnaEstimate(
                metadata, bucket.genes, bucket.classified,
                bucket.unknown, bucket.conflicting, bucket.unassigned);
            if (!diagnosticComplete) {
                estimate.valid = false;
                estimate.status = unavailableStatus;
            }
            estimates.push_back(estimate);
        }
        FlexGdnaEstimate libraryEstimate = flexGdnaEstimate(
            metadata, libraryBucket.genes, libraryBucket.classified,
            libraryBucket.unknown, libraryBucket.conflicting,
            libraryBucket.unassigned);
        if (!diagnosticComplete) {
            libraryEstimate.valid = false;
            libraryEstimate.status = unavailableStatus;
        }

        std::string gdnaSummaryPath = outputPrefix;
        if (!gdnaSummaryPath.empty() && gdnaSummaryPath.back() != '/')
            gdnaSummaryPath += '/';
        gdnaSummaryPath += "flex_gdna_summary.tsv";
        std::ofstream gdnaSummary(gdnaSummaryPath.c_str());
        if (gdnaSummary.is_open()) {
            gdnaSummary
                << "scope\tstatus\testimated_gdna_content"
                   "\testimated_gdna_percent\testimated_gdna_unspliced_threshold"
                   "\testimated_gdna_per_probe\tgene_assigned_filtered_molecules"
                   "\tclassified_filtered_molecules"
                   "\tunknown_region_filtered_molecules"
                   "\tconflicting_region_filtered_molecules"
                   "\tgene_unassigned_filtered_molecules\tcontrol_genes"
                   "\tmodel_constant\tmodel_slope\tmodel_critical_point\tmodel_rss\n";
        }

        auto writeSummaryRow = [&](const std::string& scope,
                                   const FlexGdnaEstimate& estimate) {
            if (!gdnaSummary.is_open())
                return;
            gdnaSummary << std::setprecision(17)
                        << scope << '\t' << estimate.status << '\t';
            if (estimate.valid) {
                gdnaSummary
                    << estimate.estimatedGdnaFraction << '\t'
                    << estimate.estimatedGdnaFraction * 100.0 << '\t'
                    << estimate.threshold << '\t'
                    << estimate.estimatedGdnaPerProbe << '\t';
            } else {
                gdnaSummary << "NA\tNA\tNA\tNA\t";
            }
            gdnaSummary
                << estimate.totalFilteredMolecules << '\t'
                << estimate.classifiedMolecules << '\t'
                << estimate.unknownMolecules << '\t'
                << estimate.conflictingMolecules << '\t'
                << estimate.unassignedMolecules << '\t'
                << estimate.controlGenes << '\t';
            if (estimate.valid) {
                gdnaSummary
                    << estimate.modelConstant << '\t'
                    << estimate.modelSlope << '\t'
                    << estimate.modelCriticalPoint << '\t'
                    << estimate.modelRss;
            } else {
                gdnaSummary << "NA\tNA\tNA\tNA";
            }
            gdnaSummary << '\n';
        };

        for (size_t sample = 0; sample < outputs.tagResults.size(); ++sample) {
            const std::string& label = outputs.tagResults[sample].sampleLabel;
            writeSummaryRow(label, estimates[sample]);
            std::string samplePrefix = outputPrefix;
            if (!samplePrefix.empty() && samplePrefix.back() != '/')
                samplePrefix += '/';
            samplePrefix += label + "/Gene/filtered/";
            createDirectory(samplePrefix, P.runDirPerm,
                            "Flex gDNA per-sample metrics directory", P);
            writeGdnaJson(samplePrefix + "gdna_metrics.json", label, estimates[sample]);
        }
        writeSummaryRow("library", libraryEstimate);
        writeGdnaJson(outputPrefix + (outputPrefix.empty() || outputPrefix.back() == '/'
                         ? "" : "/") + "flex_gdna_library.json",
                      "library", libraryEstimate);

        P.inOut->logMain
            << "Flex gDNA diagnostic: " << libraryEstimate.status;
        if (libraryEstimate.valid) {
            P.inOut->logMain
                << " (estimated_content=" << libraryEstimate.estimatedGdnaFraction
                << ", threshold=" << libraryEstimate.threshold
                << ", gene_assigned_molecules="
                << libraryEstimate.totalFilteredMolecules
                << ")";
        }
        P.inOut->logMain << "\n";

        if (pSolo.flexGdnaMode == ParametersSolo::FlexGdnaRequired
            && !libraryEstimate.valid) {
            std::ostringstream errMsg;
            errMsg << "ERROR: required Flex gDNA diagnostic failed: "
                   << libraryEstimate.status;
            exitWithError(errMsg.str(), std::cerr, P.inOut->logMain,
                          EXIT_CODE_RUNTIME, P);
        }
    }
}
