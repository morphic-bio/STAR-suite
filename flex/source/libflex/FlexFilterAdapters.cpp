#include "FlexFilterAdapters.h"
#include "OrdMagStage.h"
#include "ErrorWarning.h"
#include "streamFuns.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>

using namespace std;

// STAR adapter implementation

bool STARFlexFilterAdapter::ensureDirectory(const string& path, mode_t mode) {
    // Use mkdir -p via system call (simple approach)
    string cmd = "mkdir -p " + path;
    int result = system(cmd.c_str());
    if (result == 0 && mode != 0755) {
        // Try to set permissions if different from default
        chmod(path.c_str(), mode);
    }
    return result == 0;
}

void STARFlexFilterAdapter::writeOrdMagOutputs(
    const FlexFilter::Outputs::TagResults& tagResult,
    const string& outputPrefix) {
    
    // Write summary stats (OrdMag is no longer used for ED, but retain directory for compatibility)
    string ordmagDir = outputPrefix + "OrdMag/";
    if (!ensureDirectory(ordmagDir, 0755)) {
        return;
    }
    
    // Write summary file
    string summaryFile = ordmagDir + "filter_summary.json";
    ofstream summaryOut(summaryFile);
    if (summaryOut.is_open()) {
        summaryOut << "{\n";
        summaryOut << "  \"expected_cells\": " << tagResult.expectedCells << ",\n";
        summaryOut << "  \"retain_window\": " << tagResult.nRetainWindow << ",\n";
        summaryOut << "  \"simple_cells\": " << tagResult.nSimpleCells << ",\n";
        summaryOut << "  \"tail_tested\": " << tagResult.nTailTested << ",\n";
        summaryOut << "  \"simple_passers\": " << tagResult.nSimplePassers << ",\n";
        summaryOut << "  \"tail_passers\": " << tagResult.nTailPassers << ",\n";
        summaryOut << "  \"occupancy_removed\": " << tagResult.occupancyRemoved << ",\n";
        summaryOut << "  \"final_cells\": " << tagResult.passingBarcodes.size() << "\n";
        summaryOut << "}\n";
        summaryOut.close();
    }
}

void STARFlexFilterAdapter::writeEmptyDropsOutputs(
    const FlexFilter::Outputs::TagResults& tagResult,
    const string& outputPrefix) {
    
    string edDir = outputPrefix + "InlineDrops/";
    if (!ensureDirectory(edDir, 0755)) {
        return;
    }
    
    // Write passing_barcodes.txt
    string passingFile = edDir + "passing_barcodes.txt";
    ofstream passingOut(passingFile);
    if (passingOut.is_open()) {
        for (const string& barcode : tagResult.passingBarcodes) {
            passingOut << barcode << "\n";
        }
        passingOut.close();
    }
    
    // Write pvalues.csv
    string pvaluesFile = edDir + "pvalues.csv";
    ofstream pvaluesOut(pvaluesFile);
    if (pvaluesOut.is_open()) {
        pvaluesOut << "barcode,cell_index,umi_count,obs_log_prob,p_value,p_adjusted,rank\n";
        
        // Map cellIndex to barcode using tagBarcodes
        vector<string> barcodes = tagResult.tagBarcodes;
        if (barcodes.empty()) {
            barcodes = tagResult.filteredBarcodes; // Fallback
        }
        
        for (const auto& result : tagResult.emptydropsResults) {
            string barcode = (result.cellIndex < barcodes.size()) ? barcodes[result.cellIndex] : "UNKNOWN";
            // Remove -1 suffix if present
            if (barcode.size() > 2 && barcode.substr(barcode.size() - 2) == "-1") {
                barcode = barcode.substr(0, barcode.size() - 2);
            }
            
            // Get UMI count if available (would need to pass through from FlexFilter)
            uint32_t umiCount = 0; // TODO: Store UMI counts in TagResults
            
            pvaluesOut << barcode << "," << result.cellIndex << "," << umiCount << ","
                      << result.obsLogProb << "," << result.pValue << ","
                      << result.pAdjusted << "," << result.monteCarloRank << "\n";
        }
        pvaluesOut.close();
    }
}

void STARFlexFilterAdapter::writeMEXFiles(
    const FlexFilter::Outputs::TagResults& /* tagResult */,
    const SampleMatrixData& /* matrixData */,
    const string& /* outputPrefix */) {
    
    // This would need to reconstruct MEX from filtered results
    // For now, this is a placeholder - actual implementation would filter matrixData
    // to only include passing barcodes and write MEX format
    
    // TODO: Implement MEX writing from filtered matrix
    // Similar to PerSamplePipeline::writeMEXFiles but using tagResult.passingBarcodes
}

void STARFlexFilterAdapter::writeAggregateFilterSummary(
    const FlexFilter::Outputs& outputs,
    const string& outputPrefix) {
    
    string summaryFile = outputPrefix + "filter_summary.json";
    ofstream summaryOut(summaryFile);
    if (!summaryOut.is_open()) {
        return;
    }
    
    // Compute totals
    uint32_t totalFinal = 0;
    uint32_t totalOccRemoved = 0;
    for (const auto& tagResult : outputs.tagResults) {
        totalFinal += static_cast<uint32_t>(tagResult.passingBarcodes.size());
        totalOccRemoved += tagResult.occupancyRemoved;
    }
    
    summaryOut << "{\n";
    summaryOut << "  \"total_samples\": " << outputs.tagResults.size() << ",\n";
    summaryOut << "  \"total_final_cells\": " << totalFinal << ",\n";
    summaryOut << "  \"total_occupancy_removed\": " << totalOccRemoved << ",\n";
    summaryOut << "  \"tag_breakdown\": [\n";
    
    for (size_t i = 0; i < outputs.tagResults.size(); i++) {
        const auto& tagResult = outputs.tagResults[i];
        summaryOut << "    {\n";
        summaryOut << "      \"tag\": \"" << tagResult.tag << "\",\n";
        summaryOut << "      \"expected_cells\": " << tagResult.expectedCells << ",\n";
        summaryOut << "      \"simple_cells\": " << tagResult.nSimpleCells << ",\n";
        summaryOut << "      \"tail_tested\": " << tagResult.nTailTested << ",\n";
        summaryOut << "      \"simple_passers\": " << tagResult.nSimplePassers << ",\n";
        summaryOut << "      \"tail_passers\": " << tagResult.nTailPassers << ",\n";
        summaryOut << "      \"occupancy_removed\": " << tagResult.occupancyRemoved << ",\n";
        summaryOut << "      \"final_cells\": " << tagResult.passingBarcodes.size() << "\n";
        summaryOut << "    }";
        if (i < outputs.tagResults.size() - 1) {
            summaryOut << ",";
        }
        summaryOut << "\n";
    }
    
    summaryOut << "  ]\n";
    summaryOut << "}\n";
    summaryOut.close();
}

void STARFlexFilterAdapter::writeOutputs(
    const FlexFilter::Outputs& outputs,
    const string& outputPrefix,
    mode_t /* runDirPerm */) {
    
    // Write outputs for each tag/sample
    // outputPrefix should be like "Solo.out/sample_outs/" - we append sampleLabel/Gene/filtered/
    for (const auto& tagResult : outputs.tagResults) {
        string samplePrefix = outputPrefix + tagResult.sampleLabel + "/Gene/filtered/";
        
        writeOrdMagOutputs(tagResult, samplePrefix);
        writeEmptyDropsOutputs(tagResult, samplePrefix);
        // writeMEXFiles would be called here if matrix data available
    }
    
    // Write aggregate filter_summary.json at outputPrefix level
    // Note: Per-sample filter_summary.json is written by writeOrdMagOutputs
    // This aggregate one is optional but can be useful for overall stats
}

// Standalone adapter implementation

bool StandaloneFlexFilterAdapter::ensureDirectory(const string& path) {
    string cmd = "mkdir -p " + path;
    return system(cmd.c_str()) == 0;
}

void StandaloneFlexFilterAdapter::writePerTagFilteredBarcodes(
    const FlexFilter::Outputs::TagResults& tagResult,
    const string& outputPrefix) {
    
    string byTagDir = outputPrefix + "by_tag/";
    if (!ensureDirectory(byTagDir)) {
        return;
    }
    
    string filePath = byTagDir + tagResult.tag + "_filtered_barcodes.tsv";
    ofstream out(filePath);
    if (!out.is_open()) {
        return;
    }
    
    out << "barcode\n";
    // Barcodes are already sorted by index (matching golden file order)
    for (const string& barcode : tagResult.filteredBarcodes) {
        out << barcode << "\n";
    }
    out.close();
}

void StandaloneFlexFilterAdapter::writeEmptyDropsPassing(
    const FlexFilter::Outputs::TagResults& tagResult,
    const string& outputPrefix) {
    
    string edDir = outputPrefix + "emptydrops_" + tagResult.tag + "/";
    if (!ensureDirectory(edDir)) {
        return;
    }
    
    string passingFile = edDir + "passing_barcodes.txt";
    ofstream out(passingFile);
    if (!out.is_open()) {
        return;
    }
    
    // Barcodes are already sorted by index (matching golden file order)
    for (const string& barcode : tagResult.passingBarcodes) {
        out << barcode << "\n";
    }
    out.close();
}

void StandaloneFlexFilterAdapter::writePerTagExpected(
    const FlexFilter::Outputs& outputs,
    const string& outputPrefix) {
    
    string filePath = outputPrefix + "per_tag_expected.tsv";
    ofstream out(filePath);
    if (!out.is_open()) {
        return;
    }
    
    out << "tag\texpected_cells\n";
    for (const auto& tagResult : outputs.tagResults) {
        out << tagResult.tag << "\t" << tagResult.expectedCells << "\n";
    }
    out.close();
}

void StandaloneFlexFilterAdapter::writeAggregateFilterSummary(
    const FlexFilter::Outputs& outputs,
    const string& outputPrefix) {
    
    string summaryFile = outputPrefix + "filter_summary.json";
    ofstream summaryOut(summaryFile);
    if (!summaryOut.is_open()) {
        return;
    }
    
    // Compute aggregate statistics
    uint32_t totalBarcodesInput = 0;
    uint32_t totalRetain = 0;
    uint32_t totalSimple = 0;
    uint32_t totalTailTested = 0;
    uint32_t totalSimplePass = 0;
    uint32_t totalTailPass = 0;
    uint32_t totalOccRemoved = 0;
    uint32_t totalFinal = 0;
    uint32_t totalExpectedCells = 0;
    
    for (const auto& tagResult : outputs.tagResults) {
        totalBarcodesInput += static_cast<uint32_t>(tagResult.tagBarcodes.size());
        totalRetain += tagResult.nRetainWindow;
        totalSimple += tagResult.nSimpleCells;
        totalTailTested += tagResult.nTailTested;
        totalSimplePass += tagResult.nSimplePassers;
        totalTailPass += tagResult.nTailPassers;
        totalOccRemoved += tagResult.occupancyRemoved;
        totalFinal += static_cast<uint32_t>(tagResult.passingBarcodes.size());
        totalExpectedCells += tagResult.expectedCells;
    }
    
    summaryOut << "{\n";
    summaryOut << "  \"total_barcodes_input\": " << totalBarcodesInput << ",\n";
    summaryOut << "  \"total_retain_window\": " << totalRetain << ",\n";
    summaryOut << "  \"total_simple_cells\": " << totalSimple << ",\n";
    summaryOut << "  \"total_tail_tested\": " << totalTailTested << ",\n";
    summaryOut << "  \"total_simple_passers\": " << totalSimplePass << ",\n";
    summaryOut << "  \"total_tail_passers\": " << totalTailPass << ",\n";
    summaryOut << "  \"total_occupancy_removed\": " << totalOccRemoved << ",\n";
    summaryOut << "  \"total_final_cells\": " << totalFinal << ",\n";
    summaryOut << "  \"total_expected_cells\": " << totalExpectedCells << ",\n";
    summaryOut << "  \"tag_breakdown\": [\n";
    
    for (size_t i = 0; i < outputs.tagResults.size(); i++) {
        const auto& tagResult = outputs.tagResults[i];
        summaryOut << "    {\n";
        summaryOut << "      \"tag\": \"" << tagResult.tag << "\",\n";
        summaryOut << "      \"expected_cells\": " << tagResult.expectedCells << ",\n";
        summaryOut << "      \"retain_window\": " << tagResult.nRetainWindow << ",\n";
        summaryOut << "      \"simple_cells\": " << tagResult.nSimpleCells << ",\n";
        summaryOut << "      \"tail_tested\": " << tagResult.nTailTested << ",\n";
        summaryOut << "      \"simple_passers\": " << tagResult.nSimplePassers << ",\n";
        summaryOut << "      \"tail_passers\": " << tagResult.nTailPassers << ",\n";
        summaryOut << "      \"occupancy_removed\": " << tagResult.occupancyRemoved << ",\n";
        summaryOut << "      \"final_cells\": " << tagResult.passingBarcodes.size() << "\n";
        summaryOut << "    }";
        if (i < outputs.tagResults.size() - 1) {
            summaryOut << ",";
        }
        summaryOut << "\n";
    }
    
    summaryOut << "  ]\n";
    summaryOut << "}\n";
    summaryOut.close();
}

void StandaloneFlexFilterAdapter::writeOutputs(
    const FlexFilter::Outputs& outputs,
    const string& outputPrefix) {
    
    // Write per-tag outputs
    for (const auto& tagResult : outputs.tagResults) {
        writePerTagFilteredBarcodes(tagResult, outputPrefix);
        writeEmptyDropsPassing(tagResult, outputPrefix);
    }
    
    // Write aggregate files
    writePerTagExpected(outputs, outputPrefix);
    writeAggregateFilterSummary(outputs, outputPrefix);
}

