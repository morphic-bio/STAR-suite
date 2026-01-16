#include "FlexFilterIO.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cstdio>
#include <sys/stat.h>
#include <sys/types.h>

namespace FlexFilterIO {

//-----------------------------------------------------------------------------
// Utility Functions
//-----------------------------------------------------------------------------

void trimString(std::string& s) {
    const char* ws = " \t\r\n";
    size_t start = s.find_first_not_of(ws);
    if (start == std::string::npos) {
        s.clear();
        return;
    }
    size_t end = s.find_last_not_of(ws);
    s = s.substr(start, end - start + 1);
}

uint64_t sumUMIFromMatrix(const std::string& matrixPath) {
    std::ifstream mtx(matrixPath);
    if (!mtx.is_open()) {
        return 0;
    }
    std::string line;
    // Skip comments
    while (std::getline(mtx, line)) {
        if (!line.empty() && line[0] != '%') {
            break; // dimension line read into 'line'
        }
    }
    // Sum counts
    uint64_t total = 0;
    uint32_t r, c, v;
    while (mtx >> r >> c >> v) {
        total += v;
    }
    return total;
}

//-----------------------------------------------------------------------------
// Directory Creation
//-----------------------------------------------------------------------------

bool defaultCreateDirectory(const std::string& path, mode_t perm) {
    // Create directory hierarchy (mkdir -p equivalent)
    std::string currentPath;
    std::istringstream pathStream(path);
    std::string segment;
    
    // Handle absolute paths
    if (!path.empty() && path[0] == '/') {
        currentPath = "/";
    }
    
    while (std::getline(pathStream, segment, '/')) {
        if (segment.empty()) continue;
        
        if (!currentPath.empty() && currentPath.back() != '/') {
            currentPath += '/';
        }
        currentPath += segment;
        
        struct stat st;
        if (stat(currentPath.c_str(), &st) != 0) {
            // Directory doesn't exist, create it
            if (mkdir(currentPath.c_str(), perm) != 0) {
                std::cerr << "WARNING: Failed to create directory: " << currentPath << std::endl;
                return false;
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------
// Whitelist Loading
//-----------------------------------------------------------------------------

bool loadSampleWhitelist(
    const std::string& path,
    SampleWhitelist& whitelist,
    std::string& errorMsg)
{
    whitelist.sampleLabels.clear();
    whitelist.sampleTags.clear();
    whitelist.labeledCount = 0;
    whitelist.tagOnlyCount = 0;
    
    if (path.empty()) {
        errorMsg = "Empty whitelist path";
        return false;
    }
    
    std::ifstream file(path);
    if (!file.is_open()) {
        errorMsg = "Could not open sample whitelist file: " + path;
        return false;
    }
    
    std::vector<std::string> tempLabels;
    std::vector<std::string> tempTags;
    
    std::string line;
    uint32_t lineNum = 0;
    while (std::getline(file, line)) {
        lineNum++;
        trimString(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Try TAB-delimited first: sample_name<TAB>tag_sequence
        size_t tabPos = line.find('\t');
        if (tabPos != std::string::npos) {
            std::string label = line.substr(0, tabPos);
            std::string tag = line.substr(tabPos + 1);
            
            // Trim fields
            trimString(label);
            trimString(tag);
            
            if (!tag.empty()) {
                // Limit tag to first 8bp
                if (tag.size() > 8) {
                    tag = tag.substr(0, 8);
                }
                if (label.empty()) {
                    // Fall back to auto-name when the label column is blank
                    label = "TAG_" + tag;
                    whitelist.tagOnlyCount++;
                } else {
                    whitelist.labeledCount++;
                }
                tempLabels.push_back(label);
                tempTags.push_back(tag);
                continue;
            }
        }
        
        // Single-token: tag sequence only (limit to 8bp)
        std::string token = line;
        trimString(token);
        
        if (!token.empty() && token.size() >= 8) {
            std::string tag = token.substr(0, 8);  // Limit to first 8bp
            std::string label = "TAG_" + tag;  // Auto-generate label
            tempLabels.push_back(label);
            tempTags.push_back(tag);
            whitelist.tagOnlyCount++;
        } else if (!token.empty()) {
            std::cerr << "WARNING: Line " << lineNum << ": Tag too short (<8bp): " << token << std::endl;
        }
    }
    file.close();
    
    if (tempLabels.empty()) {
        errorMsg = "No valid samples found in whitelist file: " + path;
        return false;
    }
    
    // Preserve input order from the whitelist file to keep stable sample ordering
    whitelist.sampleLabels = std::move(tempLabels);
    whitelist.sampleTags = std::move(tempTags);
    
    // Log counts
    std::cout << "Loaded " << whitelist.sampleLabels.size() << " samples ("
              << whitelist.labeledCount << " labeled, "
              << whitelist.tagOnlyCount << " tag-only)" << std::endl;
    
    return true;
}

//-----------------------------------------------------------------------------
// MEX Loading
//-----------------------------------------------------------------------------

bool loadCompositeMex(
    const std::string& mexDir,
    CompositeMexData& data,
    std::string& errorMsg)
{
    data.barcodes.clear();
    data.features.clear();
    data.cellGeneMap.clear();
    data.nCells = 0;
    data.nGenes = 0;
    
    // Determine file paths (with InlineHashDedup_ prefix fallback)
    std::string barcodesPath = mexDir + "/barcodes.tsv";
    std::string featuresPath = mexDir + "/features.tsv";
    std::string matrixPath = mexDir + "/matrix.mtx";
    
    // Check if standard files exist, otherwise try prefix
    std::ifstream testBc(barcodesPath);
    if (!testBc.is_open()) {
        barcodesPath = mexDir + "/InlineHashDedup_barcodes.tsv";
        featuresPath = mexDir + "/InlineHashDedup_features.tsv";
        matrixPath = mexDir + "/InlineHashDedup_matrix.mtx";
    }
    testBc.close();
    
    // Load barcodes
    std::ifstream bcFile(barcodesPath);
    if (!bcFile.is_open()) {
        errorMsg = "Cannot open barcodes file: " + barcodesPath;
        return false;
    }
    std::string line;
    while (std::getline(bcFile, line)) {
        trimString(line);
        if (!line.empty()) {
            data.barcodes.push_back(line);
        }
    }
    bcFile.close();
    
    if (data.barcodes.empty()) {
        errorMsg = "No barcodes found in: " + barcodesPath;
        return false;
    }
    
    // Load features (id, name, type) - matching external column handling exactly
    std::ifstream ftFile(featuresPath);
    if (!ftFile.is_open()) {
        errorMsg = "Cannot open features file: " + featuresPath;
        return false;
    }
    uint32_t featureLineNum = 0;
    while (std::getline(ftFile, line)) {
        featureLineNum++;
        trimString(line);
        if (line.empty()) continue;
        
        std::string geneId, geneName, featureType;
        std::istringstream ls(line);
        std::getline(ls, geneId, '\t');
        std::getline(ls, geneName, '\t');
        std::getline(ls, featureType, '\t');
        
        trimString(geneId);
        trimString(geneName);
        trimString(featureType);
        
        if (geneId.empty()) {
            std::cerr << "WARNING: Malformed feature line " << featureLineNum 
                      << " (missing gene ID): " << line << std::endl;
            continue;
        }
        
        // Default name and type if missing
        if (geneName.empty()) {
            geneName = geneId;
            std::cerr << "WARNING: Feature line " << featureLineNum 
                      << " missing gene name, using ID: " << geneId << std::endl;
        }
        if (featureType.empty()) {
            featureType = "Gene Expression";
        }
        
        data.features.emplace_back(geneId, geneName, featureType);
    }
    ftFile.close();
    
    if (data.features.empty()) {
        errorMsg = "No features found in: " + featuresPath;
        return false;
    }
    
    // Load matrix (Matrix Market COO format)
    std::ifstream mtxFile(matrixPath);
    if (!mtxFile.is_open()) {
        errorMsg = "Cannot open matrix file: " + matrixPath;
        return false;
    }
    
    // Skip comments
    while (std::getline(mtxFile, line)) {
        if (line.empty() || line[0] != '%') break;
    }
    
    // Parse dimensions
    uint32_t numGenes = 0, numCells = 0, numNonZero = 0;
    {
        std::istringstream dimStream(line);
        dimStream >> numGenes >> numCells >> numNonZero;
    }
    
    if (numCells == 0 || numGenes == 0) {
        errorMsg = "Invalid matrix dimensions in: " + matrixPath;
        return false;
    }
    
    data.nCells = numCells;
    data.nGenes = numGenes;
    data.cellGeneMap.assign(numCells, {});
    
    // Read triplets: gene_row cell_col count (1-based)
    uint32_t geneRow, cellCol, count;
    while (mtxFile >> geneRow >> cellCol >> count) {
        if (geneRow == 0 || cellCol == 0) continue;
        uint32_t geneIdx = geneRow - 1;
        uint32_t cellIdx = cellCol - 1;
        if (cellIdx < numCells && geneIdx < numGenes) {
            data.cellGeneMap[cellIdx].emplace_back(geneIdx, count);
        }
    }
    mtxFile.close();
    
    std::cout << "Loaded MEX: " << data.barcodes.size() << " barcodes, "
              << data.features.size() << " features, "
              << numNonZero << " entries" << std::endl;
    
    return true;
}

//-----------------------------------------------------------------------------
// MEX Writing
//-----------------------------------------------------------------------------

void writeFilteredMexForTag(
    const std::string& outputPrefix,
    const std::string& sampleLabel,
    const std::vector<std::string>& passingBarcodes,
    const CompositeMexData& mexData,
    CreateDirectoryFunc createDir,
    mode_t runDirPerm,
    MexWriteResult& result,
    int cb_len)
{
    result.success = false;
    result.cellsWritten = 0;
    result.entriesWritten = 0;
    result.missingBarcodes = 0;
    result.errorMessage.clear();
    
    if (passingBarcodes.empty()) {
        std::cout << "    Skipping " << sampleLabel << " (no passing barcodes)" << std::endl;
        result.success = true;  // Not an error, just empty
        return;
    }
    
    // Build output directory path
    std::string sampleDir = outputPrefix;
    if (!sampleDir.empty() && sampleDir.back() != '/') {
        sampleDir += '/';
    }
    sampleDir += sampleLabel + "/Gene/filtered/";
    
    // Create directory
    CreateDirectoryFunc dirFunc = createDir ? createDir : defaultCreateDirectory;
    if (!dirFunc(sampleDir, runDirPerm)) {
        result.errorMessage = "Failed to create directory: " + sampleDir;
        std::cerr << "WARNING: " << result.errorMessage << std::endl;
        return;
    }
    
    std::cout << "    Writing MEX to: " << sampleDir << std::endl;
    
    // Build barcode index map
    std::unordered_map<std::string, uint32_t> barcodeToIdx;
    barcodeToIdx.reserve(mexData.barcodes.size() * 2);
    for (uint32_t i = 0; i < mexData.barcodes.size(); i++) {
        barcodeToIdx[mexData.barcodes[i]] = i;
    }
    
    // Map passing barcodes to old indices, build filtered barcode list
    std::unordered_map<uint32_t, uint32_t> oldToNew;
    std::vector<std::string> filteredBarcodes;
    filteredBarcodes.reserve(passingBarcodes.size());
    
    for (const auto& bc : passingBarcodes) {
        auto it = barcodeToIdx.find(bc);
        if (it == barcodeToIdx.end()) {
            std::cerr << "WARNING: Passing barcode not found in MEX: " << bc << std::endl;
            result.missingBarcodes++;
            continue;
        }
        uint32_t oldIdx = it->second;
        if (oldToNew.find(oldIdx) != oldToNew.end()) {
            continue;  // Skip duplicate
        }
        uint32_t newIdx = static_cast<uint32_t>(filteredBarcodes.size());
        oldToNew[oldIdx] = newIdx;
        filteredBarcodes.push_back(bc);
    }
    
    if (filteredBarcodes.empty()) {
        std::cout << "    Skipping " << sampleLabel << " (no passing barcodes mapped)" << std::endl;
        result.success = true;
        return;
    }
    
    // Build filtered triplets - iterate in cell index order for deterministic output
    std::vector<MexWriter::Triplet> filteredTriplets;
    filteredTriplets.reserve(mexData.cellGeneMap.size() * 2);  // rough estimate
    
    for (uint32_t oldIdx = 0; oldIdx < mexData.cellGeneMap.size(); oldIdx++) {
        auto it = oldToNew.find(oldIdx);
        if (it == oldToNew.end()) continue;  // Skip cells not in passing set
        
        uint32_t newIdx = it->second;
        for (const auto& gc : mexData.cellGeneMap[oldIdx]) {
            filteredTriplets.push_back({newIdx, gc.first, gc.second});
        }
    }
    
    // Ensure directory path ends with /
    std::string samplePrefix = sampleDir;
    if (!samplePrefix.empty() && samplePrefix.back() != '/') {
        samplePrefix.push_back('/');
    }
    
    // Write MEX (cb_len controls barcode truncation: 16 = strip tag, -1 = keep full)
    int writeResult = MexWriter::writeMex(samplePrefix, filteredBarcodes, mexData.features, filteredTriplets, cb_len);
    
    if (writeResult != 0) {
        result.errorMessage = "MexWriter failed for " + sampleLabel;
        std::cerr << "  ERROR: " << result.errorMessage 
                  << " (barcodes=" << filteredBarcodes.size()
                  << ", entries=" << filteredTriplets.size() << ")" << std::endl;
        return;
    }
    
    result.success = true;
    result.cellsWritten = static_cast<uint32_t>(filteredBarcodes.size());
    result.entriesWritten = static_cast<uint32_t>(filteredTriplets.size());
    
    std::cout << "  " << sampleLabel << ": "
              << filteredBarcodes.size() << " cells, "
              << filteredTriplets.size() << " entries" << std::endl;
}

//-----------------------------------------------------------------------------
// Summary Emission
//-----------------------------------------------------------------------------

void emitSummary(
    const std::string& outputPrefix,
    const FlexFilter::Outputs& outputs,
    const FlexFilter::Config& config)
{
    // Build summary path
    std::string summaryPath = outputPrefix;
    if (!summaryPath.empty() && summaryPath.back() != '/') {
        summaryPath.push_back('/');
    }
    summaryPath += "flexfilter_summary.tsv";
    
    // Open summary file
    std::ofstream summaryFile(summaryPath);
    if (summaryFile.is_open()) {
        summaryFile << "Sample\tExpected\tRetain\tSimple\tTail_Tested\tSimple_Pass\tTail_Pass\tOcc_Removed\tFinal\tTotal_UMIs\n";
    }
    
    // Print header to stdout
    std::cout << "\nSummary (saved to " << summaryPath << "):" << std::endl;
    std::printf("%-12s %8s %8s %8s %10s %10s %10s %10s %8s %14s\n",
           "Sample", "Expected", "Retain", "Simple", "Tail_Test",
           "Simp_Pass", "Tail_Pass", "Occ_Rem", "Final", "Total_UMIs");
    
    // Totals
    uint32_t totalExpected = 0;
    uint32_t totalRetain = 0;
    uint32_t totalSimple = 0;
    uint32_t totalTailTested = 0;
    uint32_t totalSimplePass = 0;
    uint32_t totalTailPass = 0;
    uint32_t totalOccRemoved = 0;
    uint32_t totalFinal = 0;
    uint64_t totalUMIs = 0;
    
    (void)config;  // Unused in new summary format
    
    for (const auto& tagResult : outputs.tagResults) {
        // Sum UMIs from output matrix
        std::string mtxPath = outputPrefix;
        if (!mtxPath.empty() && mtxPath.back() != '/') {
            mtxPath.push_back('/');
        }
        mtxPath += tagResult.sampleLabel + "/Gene/filtered/matrix.mtx";
        uint64_t totalUMI = sumUMIFromMatrix(mtxPath);
        uint32_t finalCells = static_cast<uint32_t>(tagResult.passingBarcodes.size());
        
        // Print row
        std::printf("%-12s %8u %8u %8u %10u %10u %10u %10u %8u %14lu\n",
               tagResult.sampleLabel.c_str(),
               tagResult.expectedCells,
               tagResult.nRetainWindow,
               tagResult.nSimpleCells,
               tagResult.nTailTested,
               tagResult.nSimplePassers,
               tagResult.nTailPassers,
               tagResult.occupancyRemoved,
               finalCells,
               totalUMI);
        
        // Write to TSV
        if (summaryFile.is_open()) {
            summaryFile << tagResult.sampleLabel << '\t'
                        << tagResult.expectedCells << '\t'
                        << tagResult.nRetainWindow << '\t'
                        << tagResult.nSimpleCells << '\t'
                        << tagResult.nTailTested << '\t'
                        << tagResult.nSimplePassers << '\t'
                        << tagResult.nTailPassers << '\t'
                        << tagResult.occupancyRemoved << '\t'
                        << finalCells << '\t'
                        << totalUMI << '\n';
        }
        
        // Accumulate totals
        totalExpected += tagResult.expectedCells;
        totalRetain += tagResult.nRetainWindow;
        totalSimple += tagResult.nSimpleCells;
        totalTailTested += tagResult.nTailTested;
        totalSimplePass += tagResult.nSimplePassers;
        totalTailPass += tagResult.nTailPassers;
        totalOccRemoved += tagResult.occupancyRemoved;
        totalFinal += finalCells;
        totalUMIs += totalUMI;
    }
    
    // Print TOTAL row
    std::printf("%-12s %8u %8u %8u %10u %10u %10u %10u %8u %14lu\n",
           "TOTAL",
           totalExpected,
           totalRetain,
           totalSimple,
           totalTailTested,
           totalSimplePass,
           totalTailPass,
           totalOccRemoved,
           totalFinal,
           totalUMIs);
    
    // Write TOTAL to TSV
    if (summaryFile.is_open()) {
        summaryFile << "TOTAL\t"
                    << totalExpected << '\t'
                    << totalRetain << '\t'
                    << totalSimple << '\t'
                    << totalTailTested << '\t'
                    << totalSimplePass << '\t'
                    << totalTailPass << '\t'
                    << totalOccRemoved << '\t'
                    << totalFinal << '\t'
                    << totalUMIs << '\n';
        summaryFile.close();
    }
}

} // namespace FlexFilterIO

