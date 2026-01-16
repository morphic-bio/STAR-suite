#include "FlexFilter.h"
// Inlined barcode utilities (formerly in TagDominanceCheck)
#include "ErrorWarning.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <chrono>
#include <thread>
#include <atomic>
#include <mutex>

using namespace std;

// ============================================================================
// Inlined barcode utilities (formerly in TagDominanceCheck)
// ============================================================================

// Canonicalize composite barcode label before extracting TAG8
// Handles: strips "-1" suffix, handles "S{sample}:CB" format
static string canonicalizeBarcode(const string& barcode) {
    string result = barcode;
    // Strip "-1" suffix if present
    if (result.size() > 2 && result.substr(result.size() - 2) == "-1") {
        result = result.substr(0, result.size() - 2);
    }
    // Handle "S{sample}:CB" format (extract CB portion after ':')
    size_t colonPos = result.find(':');
    if (colonPos != string::npos && colonPos + 1 < result.size()) {
        result = result.substr(colonPos + 1);
    }
    return result;
}

// Extract TAG8 (last 8 chars) from canonicalized composite barcode
static string extractTAG8(const string& canonicalBarcode, uint32_t tagLength = 8) {
    if (canonicalBarcode.size() < tagLength) return "";
    return canonicalBarcode.substr(canonicalBarcode.size() - tagLength);
}

// Extract CB16 and TAG8 in one call (fewer allocations)
static pair<string, string> extractCB16AndTAG8(const string& barcode, uint32_t tagLength = 8) {
    string canonical = canonicalizeBarcode(barcode);
    if (canonical.size() < 16 + tagLength) {
        return make_pair(string(), string());
    }
    string cb16 = canonical.substr(0, canonical.size() - tagLength);
    string tag8 = canonical.substr(canonical.size() - tagLength);
    return make_pair(cb16, tag8);
}

// ============================================================================

// Helper: Filter matrix by keeping a specific set of cell indices
// OPTIMIZATION: Added vector pre-reservation
static SampleMatrixData filterMatrixByKeepIndices(
    const SampleMatrixData& fullMatrix,
    const vector<uint32_t>& keepIndices)
{
    SampleMatrixData out;
    out.countMatStride = fullMatrix.countMatStride;
    out.nGenes = fullMatrix.nGenes;
    out.features = fullMatrix.features;
    
    if (keepIndices.empty()) {
        out.nCells = 0;
        return out;
    }
    
    // OPTIMIZATION: Count total gene entries first for pre-allocation
    uint64_t totalGeneEntries = 0;
    vector<uint32_t> keepMap(fullMatrix.nCells, static_cast<uint32_t>(-1));
    out.barcodes.reserve(keepIndices.size());
    for (uint32_t idx : keepIndices) {
        if (idx >= fullMatrix.barcodes.size()) continue;
        keepMap[idx] = static_cast<uint32_t>(out.barcodes.size());
        out.barcodes.push_back(fullMatrix.barcodes[idx]);
        totalGeneEntries += fullMatrix.nGenePerCB[idx];
    }
    
    out.nCells = static_cast<uint32_t>(out.barcodes.size());
    out.nUMIperCB.assign(out.nCells, 0);
    out.nGenePerCB.assign(out.nCells, 0);
    out.countCellGeneUMIindex.assign(out.nCells + 1, 0);  // +1 for sentinel
    out.countCellGeneUMI.reserve(totalGeneEntries * out.countMatStride);  // Pre-allocate
    
    uint32_t currentCell = static_cast<uint32_t>(-1);
    uint64_t matrixIdx = 0;
    for (uint32_t i = 0; i < fullMatrix.nCells; i++) {
        uint32_t newIdx = keepMap[i];
        if (newIdx == static_cast<uint32_t>(-1)) continue;
        
        if (newIdx != currentCell) {
            currentCell = newIdx;
            out.countCellGeneUMIindex[newIdx] = matrixIdx * out.countMatStride;
        }
        
        uint64_t startIdx = fullMatrix.countCellGeneUMIindex[i];
        uint32_t nGenes = fullMatrix.nGenePerCB[i];
        for (uint32_t g = 0; g < nGenes; g++) {
            uint32_t geneIdx = fullMatrix.countCellGeneUMI[startIdx + g * fullMatrix.countMatStride];
            uint32_t count = fullMatrix.countCellGeneUMI[startIdx + g * fullMatrix.countMatStride + 1];
            out.countCellGeneUMI.push_back(geneIdx);
            out.countCellGeneUMI.push_back(count);
            out.countCellGeneUMI.push_back(0);
            out.nGenePerCB[newIdx]++;
            out.nUMIperCB[newIdx] += count;
            matrixIdx++;
        }
    }
    // Set sentinel element (final index = total matrix size)
    out.countCellGeneUMIindex[out.nCells] = matrixIdx * out.countMatStride;
    return out;
}

// Helper: Get external defaults for Simple EmptyDrops parameters
// These match the defaults from /mnt/pikachu/flex/scripts/run_ordmag_nonambient.py
SimpleEmptyDropsParams FlexFilter::getExternalSimpleEmptyDropsDefaults() {
    SimpleEmptyDropsParams params;
    params.nExpectedCells = 3000;  // DEFAULT_CR_EXPECTED_CELLS
    params.maxPercentile = 0.99;   // DEFAULT_CR_MAX_PERCENTILE
    params.maxMinRatio = 10.0;      // DEFAULT_CR_MAX_MIN_RATIO
    params.umiMin = 500;            // DEFAULT_CR_UMI_MIN
    params.umiMinFracMedian = 0.01; // DEFAULT_CR_UMI_MIN_FRAC_MEDIAN
    params.candMaxN = 20000;        // DEFAULT_CR_CAND_MAX_N
    params.indMin = 45000;           // DEFAULT_CR_IND_MIN
    params.indMax = 90000;           // DEFAULT_CR_IND_MAX
    params.disabled = true;         // Simple EmptyDrops is disabled by default (fallback only)
    return params;
}

// Helper: Get external defaults for EmptyDrops parameters
// These match the defaults from /mnt/pikachu/flex/scripts/run_ordmag_nonambient.py and run_emptydrops_cr_canonical.R
EmptyDropsParams FlexFilter::getExternalEmptyDropsDefaults() {
    EmptyDropsParams params;
    params.indMin = 45000;              // DEFAULT_CR_IND_MIN (ambient range)
    params.indMax = 90000;               // DEFAULT_CR_IND_MAX (ambient range)
    params.umiMin = 500;                 // DEFAULT_CR_UMI_MIN
    params.umiMinFracMedian = 0.01;      // DEFAULT_CR_UMI_MIN_FRAC_MEDIAN
    params.candMaxN = 20000;             // DEFAULT_CR_CAND_MAX_N
    params.FDR = 0.01;                   // Default FDR threshold (CR parity)
    params.rawPvalueThreshold = 0.05;   // Default p-value threshold (calibrated to Cell Ranger output)
    params.useFDRGate = true;            // Use FDR gate by default (matches CR/R)
    params.simN = 10000;                 // Monte Carlo iterations (CR/R default)
    params.seed = 19760110LLU;          // Seed used in test harness (matches golden generation)
    params.lowerTestingBound = 500;     // R's emptyDropsCellRanger umi.min parameter: cells with UMI <= 500 excluded from candidates
    params.ambientUmiMax = 100;         // R's emptyDropsCellRanger lower parameter: cells with UMI <= 100 used for ambient
    return params;
}

// Helper: Populate config with external defaults for any zero/unset fields
// Only fills fields that are zero or unset, preserving explicit overrides
void FlexFilter::populateConfigWithDefaults(Config& config) {
    SimpleEmptyDropsParams simpleEDDefaults = getExternalSimpleEmptyDropsDefaults();
    EmptyDropsParams emptydropsDefaults = getExternalEmptyDropsDefaults();
    
    // Populate Simple EmptyDrops params (only if zero/unset)
    if (config.simpleEmptyDropsParams.nExpectedCells == 0) {
        config.simpleEmptyDropsParams.nExpectedCells = simpleEDDefaults.nExpectedCells;
    }
    if (config.simpleEmptyDropsParams.maxPercentile == 0.0) {
        config.simpleEmptyDropsParams.maxPercentile = simpleEDDefaults.maxPercentile;
    }
    if (config.simpleEmptyDropsParams.maxMinRatio == 0.0) {
        config.simpleEmptyDropsParams.maxMinRatio = simpleEDDefaults.maxMinRatio;
    }
    if (config.simpleEmptyDropsParams.umiMin == 0) {
        config.simpleEmptyDropsParams.umiMin = simpleEDDefaults.umiMin;
    }
    if (config.simpleEmptyDropsParams.umiMinFracMedian == 0.0) {
        config.simpleEmptyDropsParams.umiMinFracMedian = simpleEDDefaults.umiMinFracMedian;
    }
    if (config.simpleEmptyDropsParams.candMaxN == 0) {
        config.simpleEmptyDropsParams.candMaxN = simpleEDDefaults.candMaxN;
    }
    if (config.simpleEmptyDropsParams.indMin == 0) {
        config.simpleEmptyDropsParams.indMin = simpleEDDefaults.indMin;
    }
    if (config.simpleEmptyDropsParams.indMax == 0) {
        config.simpleEmptyDropsParams.indMax = simpleEDDefaults.indMax;
    }
    // Note: disabled defaults to true (fallback only)
    
    // Populate EmptyDrops params (only if zero/unset)
    if (config.emptydropsParams.indMin == 0) {
        config.emptydropsParams.indMin = emptydropsDefaults.indMin;
    }
    if (config.emptydropsParams.indMax == 0) {
        config.emptydropsParams.indMax = emptydropsDefaults.indMax;
    }
    if (config.emptydropsParams.umiMin == 0) {
        config.emptydropsParams.umiMin = emptydropsDefaults.umiMin;
    }
    if (config.emptydropsParams.umiMinFracMedian == 0.0) {
        config.emptydropsParams.umiMinFracMedian = emptydropsDefaults.umiMinFracMedian;
    }
    if (config.emptydropsParams.candMaxN == 0) {
        config.emptydropsParams.candMaxN = emptydropsDefaults.candMaxN;
    }
    if (config.emptydropsParams.FDR == 0.0) {
        config.emptydropsParams.FDR = emptydropsDefaults.FDR;
    }
    if (config.emptydropsParams.rawPvalueThreshold == 0.0) {
        config.emptydropsParams.rawPvalueThreshold = emptydropsDefaults.rawPvalueThreshold;
    }
    // useFDRGate defaults to false, so we only set if it's explicitly true (no change needed)
    // But we can set it to false if it's uninitialized (but bool can't be checked for zero)
    // So we'll leave useFDRGate as-is (defaults to false in struct)
    if (config.emptydropsParams.simN == 0) {
        config.emptydropsParams.simN = emptydropsDefaults.simN;
    }
    if (config.emptydropsParams.seed == 0) {
        config.emptydropsParams.seed = emptydropsDefaults.seed;
    }
    if (config.emptydropsParams.lowerTestingBound == 0) {
        config.emptydropsParams.lowerTestingBound = emptydropsDefaults.lowerTestingBound;
    }
    if (config.emptydropsParams.ambientUmiMax == 0) {
        config.emptydropsParams.ambientUmiMax = emptydropsDefaults.ambientUmiMax;
    }
}

// Helper: Extract tag from composite barcode
string FlexFilter::extractTag(const string& compositeBarcode) {
    string canonical = canonicalizeBarcode(compositeBarcode);
    return extractTAG8(canonical, 8);
}

// Helper: Compute per-tag expected cells allocation
// Matches Python's allocate_expected() logic: uses ncells_simple as weights, floor+remainder distribution
map<string, uint32_t> FlexFilter::computeTagAllocations(
    const SampleMatrixData& matrixData,
    const vector<string>& sampleTags,
    const map<string, uint32_t>& override,
    const SimpleEmptyDropsParams& simpleEDParams,
    bool useNcellsSimpleWeights,
    uint32_t totalExpectedCells) {
    
    map<string, uint32_t> allocations;
    
    // If override provided, use that
    if (!override.empty()) {
        for (const auto& tag : sampleTags) {
            auto it = override.find(tag);
            if (it != override.end()) {
                allocations[tag] = it->second;
            } else {
                allocations[tag] = 0;
            }
        }
        return allocations;
    }
    
    // OPTIMIZATION: Single-pass bucketing of cells by tag
    // Build tag index for O(1) lookup
    unordered_map<string, size_t> tagToIdx;
    tagToIdx.reserve(sampleTags.size());
    for (size_t i = 0; i < sampleTags.size(); i++) {
        tagToIdx[sampleTags[i]] = i;
    }
    
    // Single pass: bucket UMI counts by tag
    vector<vector<uint32_t>> tagUMIsBucketed(sampleTags.size());
    vector<uint64_t> tagUMISums(sampleTags.size(), 0);
    for (size_t i = 0; i < tagUMIsBucketed.size(); i++) {
        tagUMIsBucketed[i].reserve(matrixData.barcodes.size() / sampleTags.size());  // estimate
    }
    
    for (size_t i = 0; i < matrixData.barcodes.size(); i++) {
        string cellTag = extractTag(matrixData.barcodes[i]);
        auto it = tagToIdx.find(cellTag);
        if (it != tagToIdx.end()) {
            size_t tagIdx = it->second;
            tagUMIsBucketed[tagIdx].push_back(matrixData.nUMIperCB[i]);
            tagUMISums[tagIdx] += matrixData.nUMIperCB[i];
        }
    }
    
    // Compute weights for tag allocation
    vector<double> weights;
    vector<string> tagOrder; // Track tag order for allocation
    weights.reserve(sampleTags.size());
    tagOrder.reserve(sampleTags.size());
    
    // Suppress unused parameter warning (kept for API compatibility)
    (void)simpleEDParams;
    
    // useNcellsSimpleWeights=false means use UMI weighting (the new default)
    // useNcellsSimpleWeights=true means equal weighting (user disabled UMI weights)
    bool useUMIWeighting = !useNcellsSimpleWeights;
    
    if (useUMIWeighting) {
        // Winsorized UMI sums (fast O(n) per tag)
        // Cap UMIs at 95th percentile to prevent outlier blobs from skewing allocation
        for (size_t ti = 0; ti < sampleTags.size(); ti++) {
            const auto& tag = sampleTags[ti];
            tagOrder.push_back(tag);
            
            vector<uint32_t> tagUMIs = tagUMIsBucketed[ti];  // Copy for nth_element (modifies in place)
            
            if (tagUMIs.empty()) {
                weights.push_back(1.0);  // Fallback for empty tags
                continue;
            }
            
            // Find 95th percentile cap using O(n) nth_element
            size_t p95_idx = max(static_cast<size_t>(1), tagUMIs.size() * 95 / 100);
            if (p95_idx >= tagUMIs.size()) p95_idx = tagUMIs.size() - 1;
            
            nth_element(tagUMIs.begin(), tagUMIs.begin() + p95_idx, tagUMIs.end());
            uint32_t cap = tagUMIs[p95_idx];
            
            // Winsorized sum: cap values above 95th percentile
            uint64_t winsorizedSum = 0;
            for (uint32_t u : tagUMIsBucketed[ti]) {  // Use original (unsorted) data
                winsorizedSum += min(u, cap);
            }
            
            weights.push_back(static_cast<double>(max(winsorizedSum, static_cast<uint64_t>(1))));
        }
    } else {
        // Equal weighting: user disabled UMI weights, allocate equally
        for (size_t ti = 0; ti < sampleTags.size(); ti++) {
            tagOrder.push_back(sampleTags[ti]);
            weights.push_back(1.0);  // Equal weight for all tags
        }
    }
    
    // Compute total expected
    // If totalExpectedCells is provided (non-zero), use it; otherwise default to 3000 * numTags
    uint32_t totalExpected = (totalExpectedCells > 0) 
        ? totalExpectedCells 
        : 3000 * max(static_cast<uint32_t>(1), static_cast<uint32_t>(sampleTags.size()));
    
    // Allocate using Python's allocate_expected() algorithm (lines 286-320)
    if (sampleTags.empty()) {
        return allocations;
    }
    
    size_t n = weights.size();
    vector<double> weightsArray = weights;
    
    // Handle all-zero or negative weights
    bool allNonPositive = true;
    for (double w : weightsArray) {
        if (w > 0.0) {
            allNonPositive = false;
            break;
        }
    }
    
    if (allNonPositive) {
        // Equal distribution
        uint32_t base = totalExpected / static_cast<uint32_t>(n);
        uint32_t remainder = totalExpected - base * static_cast<uint32_t>(n);
        for (size_t i = 0; i < n; i++) {
            allocations[tagOrder[i]] = base + (i < remainder ? 1 : 0);
        }
        return allocations;
    }
    
    // Normalize weights to non-negative
    for (size_t i = 0; i < n; i++) {
        weightsArray[i] = max(weightsArray[i], 0.0);
    }
    
    // Compute weight sum
    double weightSum = 0.0;
    for (double w : weightsArray) {
        weightSum += w;
    }
    
    if (weightSum == 0.0) {
        // Fallback: equal weights
        for (size_t i = 0; i < n; i++) {
            weightsArray[i] = 1.0;
        }
        weightSum = static_cast<double>(n);
    }
    
    // Compute raw allocations and floor them
    vector<double> raw(n);
    vector<uint32_t> floored(n);
    for (size_t i = 0; i < n; i++) {
        raw[i] = static_cast<double>(totalExpected) * weightsArray[i] / weightSum;
        floored[i] = static_cast<uint32_t>(floor(raw[i]));
    }
    
    // Compute remainder
    uint32_t allocatedSum = 0;
    for (uint32_t f : floored) {
        allocatedSum += f;
    }
    int32_t remainder = static_cast<int32_t>(totalExpected) - static_cast<int32_t>(allocatedSum);
    
    // Distribute remainder using argsort of fractional parts
    if (remainder > 0) {
        // Sort by fractional part (descending)
        vector<pair<double, size_t>> fractional;
        for (size_t i = 0; i < n; i++) {
            fractional.push_back({raw[i] - static_cast<double>(floored[i]), i});
        }
        sort(fractional.rbegin(), fractional.rend());
        
        // Distribute remainder
        for (int32_t i = 0; i < remainder; i++) {
            size_t idx = fractional[i % n].second;
            floored[idx]++;
        }
    } else if (remainder < 0) {
        // Sort by fractional part (ascending)
        vector<pair<double, size_t>> fractional;
        for (size_t i = 0; i < n; i++) {
            fractional.push_back({raw[i] - static_cast<double>(floored[i]), i});
        }
        sort(fractional.begin(), fractional.end());
        
        // Subtract from lowest fractional parts
        for (int32_t i = 0; i < -remainder; i++) {
            size_t idx = fractional[i % n].second;
            if (floored[idx] > 0) {
                floored[idx]--;
            }
        }
    }
    
    // Clip to non-negative
    for (size_t i = 0; i < n; i++) {
        floored[i] = max(floored[i], static_cast<uint32_t>(0));
    }
    
    // Handle deficit (matching Python lines 317-319)
    uint32_t finalSum = 0;
    for (uint32_t f : floored) {
        finalSum += f;
    }
    int32_t deficit = static_cast<int32_t>(totalExpected) - static_cast<int32_t>(finalSum);
    if (deficit > 0) {
        for (int32_t i = 0; i < deficit && i < static_cast<int32_t>(n); i++) {
            floored[i]++;
        }
    }
    
    // Assign allocations
    for (size_t i = 0; i < n; i++) {
        allocations[tagOrder[i]] = floored[i];
    }
    
    return allocations;
}

// Helper: Filter matrix to tag-specific barcodes
// OPTIMIZATION: Added vector pre-reservation
SampleMatrixData FlexFilter::filterMatrixToTag(
    const SampleMatrixData& fullMatrix,
    const string& tag) {
    
    SampleMatrixData tagMatrix;
    tagMatrix.countMatStride = fullMatrix.countMatStride;
    tagMatrix.nGenes = fullMatrix.nGenes;
    tagMatrix.features = fullMatrix.features;
    
    // First pass: identify matching cells and count entries
    vector<uint32_t> matchingIndices;
    matchingIndices.reserve(fullMatrix.nCells / 4);  // estimate ~25% cells per tag
    uint64_t totalGeneEntries = 0;
    
    for (uint32_t i = 0; i < fullMatrix.nCells; i++) {
        string cellTag = extractTag(fullMatrix.barcodes[i]);
        if (cellTag == tag) {
            matchingIndices.push_back(i);
            totalGeneEntries += fullMatrix.nGenePerCB[i];
        }
    }
    
    tagMatrix.nCells = static_cast<uint32_t>(matchingIndices.size());
    
    if (tagMatrix.nCells == 0) {
        return tagMatrix;
    }
    
    // OPTIMIZATION: Pre-allocate all vectors
    tagMatrix.barcodes.reserve(tagMatrix.nCells);
    tagMatrix.nUMIperCB.resize(tagMatrix.nCells, 0);
    tagMatrix.nGenePerCB.resize(tagMatrix.nCells, 0);
    tagMatrix.countCellGeneUMIindex.resize(tagMatrix.nCells + 1, 0);  // +1 for sentinel
    tagMatrix.countCellGeneUMI.reserve(totalGeneEntries * tagMatrix.countMatStride);
    
    // Second pass: copy data (using pre-identified indices)
    uint64_t matrixOffset = 0;
    for (uint32_t newIdx = 0; newIdx < matchingIndices.size(); newIdx++) {
        uint32_t oldIdx = matchingIndices[newIdx];
        
        tagMatrix.barcodes.push_back(fullMatrix.barcodes[oldIdx]);
        tagMatrix.countCellGeneUMIindex[newIdx] = matrixOffset * tagMatrix.countMatStride;
        
        // Copy matrix entries for this cell
        uint64_t startIdx = fullMatrix.countCellGeneUMIindex[oldIdx];
        uint32_t nGenes = fullMatrix.nGenePerCB[oldIdx];
        
        for (uint32_t g = 0; g < nGenes; g++) {
            uint32_t geneIdx = fullMatrix.countCellGeneUMI[startIdx + g * fullMatrix.countMatStride];
            uint32_t count = fullMatrix.countCellGeneUMI[startIdx + g * fullMatrix.countMatStride + 1];
            
            tagMatrix.countCellGeneUMI.push_back(geneIdx);
            tagMatrix.countCellGeneUMI.push_back(count);
            tagMatrix.countCellGeneUMI.push_back(0); // Padding
            
            tagMatrix.nGenePerCB[newIdx]++;
            tagMatrix.nUMIperCB[newIdx] += count;
            matrixOffset++;
        }
    }
    
    // Set sentinel element (final index = total matrix size)
    tagMatrix.countCellGeneUMIindex[tagMatrix.nCells] = matrixOffset * tagMatrix.countMatStride;
    
    return tagMatrix;
}

// Internal implementation
int FlexFilter::runInternal(
    const SampleMatrixData& matrixData,
    const vector<string>& sampleLabels,
    const vector<string>& sampleTags,
        const map<string, uint32_t>& tagExpectedCells,
        Outputs* outputs,
        const Config& config_in) {
    
    // Copy config (unified pipeline: retain window -> ED -> occupancy post-filter)
    Config config = config_in;
    
    if (outputs == nullptr) {
        return 1;
    }

    // Ensure debug output directory exists if debug logging enabled
    if (config.debugTagLog && !config.debugOutputDir.empty()) {
        mkdir(config.debugOutputDir.c_str(), 0775);
    }

    outputs->tagResults.clear();
    outputs->tagResults.resize(sampleTags.size());
    
    // Step 1: Compute per-tag expected cells allocation
    map<string, uint32_t> allocations;
    try {
        cerr << "DEBUG: Computing tag allocations: " << matrixData.barcodes.size() << " cells, " << sampleTags.size() << " tags\n";
        // useUMIWeights=true (default): use Winsorized UMI weighting
        // useUMIWeights=false: equal allocation (user disabled via --no-tag-weights-umi)
        bool skipUMIWeights = !config.useUMIWeights;
        allocations = computeTagAllocations(matrixData, sampleTags, tagExpectedCells, config.simpleEmptyDropsParams, skipUMIWeights, config.totalExpectedCells);
        cerr << "DEBUG: Tag allocations computed " << (config.useUMIWeights ? "(UMI-weighted)" : "(equal)") << "\n";
    } catch (const std::bad_alloc& e) {
        cerr << "ERROR: bad_alloc in computeTagAllocations: " << e.what() << "\n";
        cerr << "  Matrix size: " << matrixData.barcodes.size() << " cells, " << sampleTags.size() << " tags\n";
        throw;
    } catch (const std::exception& e) {
        cerr << "ERROR: Exception in computeTagAllocations: " << e.what() << "\n";
        throw;
    }
    
    // Debug: Open log file if debug enabled
    ofstream debugLog;
    string debugLogPath;
    if (config.debugTagLog) {
        debugLogPath = config.debugOutputDir.empty() ? "" : (config.debugOutputDir + "/flexfilter_debug.log");
        if (!debugLogPath.empty()) {
            debugLog.open(debugLogPath, ios::app);
        }
    }
    auto debugOut = [&](const string& msg) {
        if (config.debugTagLog) {
            if (debugLog.is_open()) {
                debugLog << msg << flush;
            } else {
                cerr << msg << flush;
            }
        }
    };
    
    // Step 3: Process each tag in parallel using std::thread
    // Note: tagResults is pre-sized; each thread writes to a different index (thread-safe).
    // Mutex protects shared debug output (debugLog, cerr).
    
    // Determine total available cores
    unsigned int totalCores = std::thread::hardware_concurrency();
    if (totalCores == 0) totalCores = 4;  // Fallback
    // Respect OMP_NUM_THREADS for compatibility
    const char* ompNumThreadsEnv = getenv("OMP_NUM_THREADS");
    if (ompNumThreadsEnv) {
        int envThreads = atoi(ompNumThreadsEnv);
        if (envThreads > 0) totalCores = static_cast<unsigned int>(envThreads);
    }
    
    // Thread allocation strategy:
    // Process tags SEQUENTIALLY to maximize MC parallelism
    // All cores go to EmptyDrops MC for each tag
    unsigned int nTags = static_cast<unsigned int>(sampleTags.size());
    unsigned int numThreads = 1;  // Sequential tag processing
    unsigned int mcThreadsPerTag = totalCores;  // All cores for MC
    
    // Set MC threads in EmptyDrops params
    config.emptydropsParams.mcThreads = mcThreadsPerTag;
    
    cerr << "[FlexFilter] Processing " << nTags << " tags with " << numThreads << " tag threads" << endl;
    cerr << "[FlexFilter] EmptyDrops MC: " << mcThreadsPerTag << " threads per tag (total cores: " << totalCores << ")" << endl;
    
    // Mutex for thread-safe debug output
    mutex outputMutex;
    
    // Atomic counter for work distribution
    atomic<size_t> nextTagIdx(0);
    
    // Thread-safe debug output lambda
    auto threadSafeDebugOut = [&](const string& msg) {
        if (config.debugTagLog) {
            lock_guard<mutex> lock(outputMutex);
            if (debugLog.is_open()) {
                debugLog << msg << flush;
            } else {
                cerr << msg << flush;
            }
        }
    };
    
    // Emit one-time message about debug log location
    if (config.debugTagLog && !debugLogPath.empty()) {
        threadSafeDebugOut("[DEBUG] flexfilter_debug.log=" + debugLogPath + "\n");
    }
    
    // Worker function for processing a single tag
    auto processTag = [&](size_t i) {
        const string& tag = sampleTags[i];
        const string& sampleLabel = (i < sampleLabels.size()) ? sampleLabels[i] : tag;

        auto t0 = chrono::steady_clock::now();
        auto last = t0;
        {
            lock_guard<mutex> lock(outputMutex);
            cerr << "[TIMING] tag=" << tag << " start\n";
        }
        auto logTiming = [&](const char* stage) {
            auto now = chrono::steady_clock::now();
            auto delta = chrono::duration_cast<chrono::milliseconds>(now - last).count();
            auto total = chrono::duration_cast<chrono::milliseconds>(now - t0).count();
            {
                lock_guard<mutex> lock(outputMutex);
                cerr << "[TIMING] tag=" << tag << " stage=" << stage
                     << " delta_ms=" << delta << " total_ms=" << total << endl;
            }
            last = now;
        };
        
        // Use thread-safe debug output
        auto debugOut = threadSafeDebugOut;
        
        debugOut("\n=== Processing tag: " + tag + " (sample: " + sampleLabel + ") ===\n");
        
        // Filter matrix to this tag only
        SampleMatrixData rawTagMatrix = filterMatrixToTag(matrixData, tag);
        
        debugOut("  After tag filtering: " + to_string(rawTagMatrix.nCells) + " cells\n");
        logTiming("filterMatrixToTag");
        
        if (rawTagMatrix.nCells == 0) {
            debugOut("  SKIP: No cells for this tag\n");
            return;  // Exit processTag lambda (was continue in original loop)
        }
        
        Outputs::TagResults tagResult;
        tagResult.tag = tag;
        tagResult.sampleLabel = sampleLabel;
        tagResult.expectedCells = (allocations.find(tag) != allocations.end()) ? allocations[tag] : 0;
        tagResult.tagBarcodes = rawTagMatrix.barcodes;  // Store full barcode list for index mapping
        
        // EmptyDrops pipeline: retain window -> ED -> occupancy post-filter
        debugOut("  [EmptyDrops] Expected cells: " + to_string(tagResult.expectedCells) + "\n");

        // Build compact matrix from tag-filtered matrix for EmptyDrops
        // This matches R's approach: retain top N by UMI, ambient = those with UMI <= 100 within retain
        
        // Step 1: Sort raw tag-filtered cells by UMI descending and apply retain window
        vector<pair<uint32_t, uint32_t>> rawUmiIdx;  // (UMI, original index in rawTagMatrix)
        rawUmiIdx.reserve(rawTagMatrix.nCells);
        for (uint32_t i = 0; i < rawTagMatrix.nCells; i++) {
            uint32_t umi = (i < rawTagMatrix.nUMIperCB.size()) ? rawTagMatrix.nUMIperCB[i] : 0;
            rawUmiIdx.push_back({umi, i});
        }
        // Sort by UMI descending (stable sort to match R's tie-breaking)
        stable_sort(rawUmiIdx.begin(), rawUmiIdx.end(), [](const pair<uint32_t,uint32_t>& a, const pair<uint32_t,uint32_t>& b) {
            return a.first > b.first;
        });
        
        // Apply retain window (top N by UMI) - use edRetainCount or default to all
        uint32_t retainCount = config.edRetainCount > 0 ? min(config.edRetainCount, rawTagMatrix.nCells) : rawTagMatrix.nCells;
        vector<uint32_t> retainIndices;  // indices into rawTagMatrix that are in the retain window
        retainIndices.reserve(retainCount);
        for (uint32_t i = 0; i < retainCount && i < rawUmiIdx.size(); i++) {
            retainIndices.push_back(rawUmiIdx[i].second);
        }
        debugOut("  [EmptyDrops] Retain window: " + to_string(retainIndices.size()) + " cells (from " + to_string(rawTagMatrix.nCells) + " raw tag-filtered)\n");
        
        // Step 2: Compute gene sums across retain window to identify nonzero genes
        vector<uint32_t> geneSums(rawTagMatrix.nGenes, 0);
        for (uint32_t rawIdx : retainIndices) {
            uint64_t startIdx = rawTagMatrix.countCellGeneUMIindex[rawIdx];
            uint32_t nGenes = rawTagMatrix.nGenePerCB[rawIdx];
            for (uint32_t g = 0; g < nGenes; g++) {
                uint32_t geneIdx = rawTagMatrix.countCellGeneUMI[startIdx + g * rawTagMatrix.countMatStride];
                uint32_t count = rawTagMatrix.countCellGeneUMI[startIdx + g * rawTagMatrix.countMatStride + 1];
                if (geneIdx < rawTagMatrix.nGenes) {
                    geneSums[geneIdx] += count;
                }
            }
        }
        
        // Step 3: Build compact feature list and remap vectors
        vector<uint32_t> compactFeatures;  // compact index -> original gene ID
        vector<uint32_t> geneIdToCompactIdx(rawTagMatrix.nGenes, UINT32_MAX);  // original gene ID -> compact index
        for (uint32_t i = 0; i < rawTagMatrix.nGenes; i++) {
            if (geneSums[i] > 0) {
                uint32_t compactIdx = static_cast<uint32_t>(compactFeatures.size());
                compactFeatures.push_back(i);
                geneIdToCompactIdx[i] = compactIdx;
            }
        }
        debugOut("  [EmptyDrops] Compact feature set: " + to_string(compactFeatures.size()) + " nonzero genes (of " + to_string(rawTagMatrix.nGenes) + " total)\n");
        
        // Step 4: Build compact matrix for cells in retain window
        // Map from rawTagMatrix index to compact matrix index
        vector<uint32_t> rawToCompactCellIdx(rawTagMatrix.nCells, UINT32_MAX);
        vector<uint32_t> countCellGeneUMI_compact;
        vector<uint32_t> countCellGeneUMIindex_compact;
        vector<uint32_t> nGenePerCB_compact;
        vector<uint32_t> nUMIperCB_compact;
        
        countCellGeneUMIindex_compact.reserve(retainIndices.size());
        nGenePerCB_compact.reserve(retainIndices.size());
        nUMIperCB_compact.reserve(retainIndices.size());
        
        for (size_t ci = 0; ci < retainIndices.size(); ci++) {
            uint32_t rawIdx = retainIndices[ci];
            rawToCompactCellIdx[rawIdx] = static_cast<uint32_t>(ci);
            
            countCellGeneUMIindex_compact.push_back(static_cast<uint32_t>(countCellGeneUMI_compact.size()));
            uint64_t startIdx = rawTagMatrix.countCellGeneUMIindex[rawIdx];
            uint32_t nGenes = rawTagMatrix.nGenePerCB[rawIdx];
            uint32_t nCompactGenes = 0;
            for (uint32_t g = 0; g < nGenes; g++) {
                uint32_t origGeneId = rawTagMatrix.countCellGeneUMI[startIdx + g * rawTagMatrix.countMatStride];
                uint32_t count = rawTagMatrix.countCellGeneUMI[startIdx + g * rawTagMatrix.countMatStride + 1];
                if (origGeneId < geneIdToCompactIdx.size() && geneIdToCompactIdx[origGeneId] != UINT32_MAX) {
                    uint32_t compactGeneId = geneIdToCompactIdx[origGeneId];
                    countCellGeneUMI_compact.push_back(compactGeneId);
                    countCellGeneUMI_compact.push_back(count);
                    nCompactGenes++;
                }
            }
            nGenePerCB_compact.push_back(nCompactGenes);
            nUMIperCB_compact.push_back(rawTagMatrix.nUMIperCB[rawIdx]);
        }
        // Add final element to countCellGeneUMIindex_compact (cumulative index array needs N+1 elements)
        countCellGeneUMIindex_compact.push_back(static_cast<uint32_t>(countCellGeneUMI_compact.size()));
        debugOut("  [EmptyDrops] Compact matrix: " + to_string(countCellGeneUMI_compact.size() / 2) + " gene-count pairs\n");
        
        // Step 5: Build compact ambient counts from cells with UMI <= ambientUmiMax WITHIN the retain window
        vector<uint32_t> compactAmbCount(compactFeatures.size(), 0);
        size_t ambientCellsUsed = 0;
        uint32_t ambientUmiMax = config.emptydropsParams.ambientUmiMax;  // R's lower parameter (default 100)
        for (size_t ci = 0; ci < retainIndices.size(); ci++) {
            uint32_t umi = nUMIperCB_compact[ci];
            if (umi > ambientUmiMax) continue;  // only include cells with UMI <= ambientUmiMax
            // Use compact matrix for ambient
            uint32_t startIdx = countCellGeneUMIindex_compact[ci];
            uint32_t nGenes = nGenePerCB_compact[ci];
            for (uint32_t g = 0; g < nGenes; g++) {
                uint32_t compactGeneId = countCellGeneUMI_compact[startIdx + g * 2];
                uint32_t count = countCellGeneUMI_compact[startIdx + g * 2 + 1];
                if (compactGeneId < compactAmbCount.size()) {
                    compactAmbCount[compactGeneId] += count;
                }
            }
            ambientCellsUsed++;
        }
        debugOut("  [EmptyDrops] Compact ambient: built from " + to_string(ambientCellsUsed) + " cells (UMI <= " + to_string(ambientUmiMax) + " within retain)\n");
        
        // Step 5: Run Good-Turing on compact ambient counts
        AmbientProfile compactAmbProfile = EmptyDropsMultinomial::computeAmbientProfile(
            compactAmbCount,
            static_cast<uint32_t>(compactFeatures.size()),
            vector<uint32_t>(),  // featDetVec (empty = use all)
            0  // featDetN
        );
        debugOut("  [EmptyDrops] Compact ambient profile: featuresDetected=" + to_string(compactAmbProfile.featuresDetected) + ", ambProfilePnon0.size()=" + to_string(compactAmbProfile.ambProfilePnon0.size()) + "\n");
        
        // Step 6: Write ambient_compact.tsv diagnostic
        if (config.debugTagLog && !config.debugOutputDir.empty()) {
            string ambCompactPath = config.debugOutputDir + "/ambient_compact.tsv";
            ofstream ambCompactOut(ambCompactPath, ios::out | ios::trunc);
            if (ambCompactOut.is_open()) {
                ambCompactOut.setf(ios::fixed);
                ambCompactOut.precision(15);
                double checkSum = 0.0;
                for (size_t i = 0; i < compactFeatures.size(); i++) {
                    uint32_t origGeneId = compactFeatures[i];
                    double p = 0.0;
                    if (i < compactAmbProfile.ambProfileLogP.size() && compactAmbProfile.ambProfileLogP[i] > -1e20) {
                        p = exp(compactAmbProfile.ambProfileLogP[i]);
                    }
                    checkSum += p;
                    if (origGeneId < rawTagMatrix.features.size()) {
                        ambCompactOut << rawTagMatrix.features[origGeneId] << "\t" << p << "\n";
                    }
                }
                ambCompactOut.close();
                debugOut("  [EmptyDrops] Wrote compact ambient to " + ambCompactPath + " (sum=" + to_string(checkSum) + ")\n");
            }
        }
        
        // Get lower testing bound from params (default 500, R's umi.min)
        uint32_t lowerTestingBound = config.emptydropsParams.lowerTestingBound;
        // Note: ambientUmiMax already declared above (line ~948)
        
        // ============================================================================
        // NEW LOGIC: EmptyDrops FIRST, Simple EmptyDrops as FALLBACK
        // Simple EmptyDrops (formerly OrdMag) is disabled by default and only runs when:
        //   1. Force-enabled via --use-simple-empty-drops
        //   2. Too few ambient cells (< simpleEDMinAmbient)
        //   3. Too few candidates above lower bound (< simpleEDMinCandidates)
        //   4. Too few ED rescues (< simpleEDMinRescues)
        // ============================================================================
        
        // Step 1: Count ambient cells and candidates to check pre-conditions
        uint32_t nAmbientCells = 0;
        uint32_t nCandidatesAboveLower = 0;
        for (size_t ci = 0; ci < nUMIperCB_compact.size(); ci++) {
            uint32_t umi = nUMIperCB_compact[ci];
            if (umi <= ambientUmiMax) {
                nAmbientCells++;
            }
            if (umi > lowerTestingBound) {
                nCandidatesAboveLower++;
            }
        }
        
        debugOut("  [Pre-check] Ambient cells (UMI <= " + to_string(ambientUmiMax) + "): " + to_string(nAmbientCells) + "\n");
        debugOut("  [Pre-check] Candidates above lower bound (UMI > " + to_string(lowerTestingBound) + "): " + to_string(nCandidatesAboveLower) + "\n");
        
        // Check fallback conditions based on pre-conditions
        bool needsFallback = false;
        vector<string> fallbackReasons;
        
        if (nAmbientCells < config.simpleEDMinAmbient) {
            needsFallback = true;
            fallbackReasons.push_back("too few ambient cells (" + to_string(nAmbientCells) + " < " + to_string(config.simpleEDMinAmbient) + ")");
        }
        if (nCandidatesAboveLower < config.simpleEDMinCandidates) {
            needsFallback = true;
            fallbackReasons.push_back("too few candidates (" + to_string(nCandidatesAboveLower) + " < " + to_string(config.simpleEDMinCandidates) + ")");
        }
        
        // Step 2: Build candidate indices for EmptyDrops (ALL candidates above lower bound, no simple cells)
        vector<uint32_t> edCandidateIndices;
        vector<uint32_t> edCandidateCounts;
        uint32_t nExcluded = 0;
        
        for (size_t ci = 0; ci < retainIndices.size(); ci++) {
            uint32_t umi = nUMIperCB_compact[ci];
            if (umi > lowerTestingBound) {
                edCandidateIndices.push_back(static_cast<uint32_t>(ci));
                edCandidateCounts.push_back(umi);
            } else {
                nExcluded++;
            }
        }
        
        // Use seed=1 for reproducibility
        config.emptydropsParams.seed = 1;
        if (config.emptydropsParams.rawPvalueThreshold == 0.0) {
            config.emptydropsParams.rawPvalueThreshold = 0.05;
        }
        
        // Populate summary stats (nSimpleCells will be updated later if fallback runs)
        tagResult.nRetainWindow = static_cast<uint32_t>(retainIndices.size());
        tagResult.nSimpleCells = 0;  // Initially 0 - may be updated by fallback
        tagResult.nTailTested = static_cast<uint32_t>(edCandidateIndices.size());
        
        debugOut("  [EmptyDrops] Lower testing bound: " + to_string(lowerTestingBound) + " UMI\n");
        debugOut("  [EmptyDrops] Candidates: " + to_string(edCandidateIndices.size()) + " (excluded " + to_string(nExcluded) + " cells with UMI <= " + to_string(lowerTestingBound) + ")\n");
        
        // Store mapping from compact index back to raw barcode for output
        vector<string> retainBarcodes;
        retainBarcodes.reserve(retainIndices.size());
        for (uint32_t rawIdx : retainIndices) {
            retainBarcodes.push_back(rawTagMatrix.barcodes[rawIdx]);
        }
        tagResult.retainBarcodes = retainBarcodes;  // Store for ED results output
        
        // Step 3: Run EmptyDrops on ALL candidates (no simple cells - they all get tested)
        if (compactAmbProfile.ambProfilePnon0.empty()) {
            debugOut("  [EmptyDrops] WARNING: Empty compact ambient profile, skipping EmptyDrops\n");
            tagResult.emptydropsResults.clear();
            needsFallback = true;
            fallbackReasons.push_back("empty ambient profile");
        } else {
            // Invariant checks (debug only)
            if (config.debugTagLog) {
                bool invariantFailed = false;
                string invariantMsg;
                
                // Check: countCellGeneUMIindex_compact.size() == retainIndices.size() + 1
                if (countCellGeneUMIindex_compact.size() != retainIndices.size() + 1) {
                    invariantFailed = true;
                    invariantMsg = "countCellGeneUMIindex_compact.size()=" + to_string(countCellGeneUMIindex_compact.size()) + 
                                   " != retainIndices.size()+1=" + to_string(retainIndices.size() + 1);
                }
                
                // Check: countCellGeneUMI_compact.size() % countMatStride == 0
                uint32_t countMatStride = 2;
                if (countCellGeneUMI_compact.size() % countMatStride != 0) {
                    invariantFailed = true;
                    if (!invariantMsg.empty()) invariantMsg += "; ";
                    invariantMsg += "countCellGeneUMI_compact.size()=" + to_string(countCellGeneUMI_compact.size()) + 
                                    " not divisible by countMatStride=" + to_string(countMatStride);
                }
                
                // Check: ambientForSampler.size() > 0 (using compactAmbProfile.ambProfileLogP as proxy)
                if (compactAmbProfile.ambProfileLogP.empty()) {
                    invariantFailed = true;
                    if (!invariantMsg.empty()) invariantMsg += "; ";
                    invariantMsg += "compactAmbProfile.ambProfileLogP.size()=0";
                }
                
                if (invariantFailed) {
                    threadSafeDebugOut("[ED_INVARIANT_FAIL] tag=" + tag + " " + invariantMsg + "\n");
                    throw runtime_error("[ED_INVARIANT_FAIL] tag=" + tag + " " + invariantMsg);
                }
            }
            
            // Call EmptyDrops with nSimpleCells=0 (all candidates get tested)
            tagResult.emptydropsResults = EmptyDropsMultinomial::computePValues(
                compactAmbProfile,
                edCandidateIndices,
                edCandidateCounts,
                countCellGeneUMI_compact,
                countCellGeneUMIindex_compact,
                nGenePerCB_compact,
                2,  // Compact stride
                1,  // umiDedupCountIndMain
                config.emptydropsParams,
                0,  // nSimpleCells = 0 (no auto-pass, all tested)
                rawTagMatrix.features,
                static_cast<uint32_t>(retainIndices.size()),
                config.debugOutputDir,
                tag,  // Pass tag name for diagnostics
                config.enableInvariantChecks  // Enable invariant checks if flag is set
            );
        }
        
        debugOut("  [EmptyDrops] Computed p-values for " + to_string(tagResult.emptydropsResults.size()) + " candidates\n");
        
        // Debug dump (gated by --debug-tag-log)
        if (config.debugTagLog && !tagResult.emptydropsResults.empty()) {
            cout << "[ED_DEBUG] " << tag << " candidates=" << tagResult.emptydropsResults.size() << endl;
            cout << "cellIndex\tbarcode\tumi\tobsLogProb\tmonteCarloRank\tpValue" << endl;
            for (const auto& res : tagResult.emptydropsResults) {
                string barcode = (res.cellIndex < retainBarcodes.size()) ? retainBarcodes[res.cellIndex] : "UNKNOWN";
                uint32_t umi = (res.cellIndex < nUMIperCB_compact.size()) ? nUMIperCB_compact[res.cellIndex] : 0;
                cout << res.cellIndex << '\t' << barcode << '\t' << umi << '\t'
                     << res.obsLogProb << '\t' << res.monteCarloRank << '\t' << res.pValue << endl;
            }
            cout << "[ED_DEBUG_END]" << endl;
        }
        
        // Step 4: Count ED rescues and check if fallback needed
        uint32_t nEDRescues = 0;
        for (const auto& result : tagResult.emptydropsResults) {
            if (result.passesRawP) {
                nEDRescues++;
            }
        }
        
        debugOut("  [EmptyDrops] ED rescues: " + to_string(nEDRescues) + "\n");
        
        if (nEDRescues < config.simpleEDMinRescues) {
            needsFallback = true;
            fallbackReasons.push_back("too few ED rescues (" + to_string(nEDRescues) + " < " + to_string(config.simpleEDMinRescues) + ")");
        }
        
        // Step 5: Run Simple EmptyDrops as fallback if needed OR if force-enabled
        // NOTE: When fallback is needed, Simple ED ALWAYS runs (disabled flag only prevents voluntary use)
        uint32_t retainThreshold = UINT32_MAX;
        unordered_set<uint32_t> simplePasserIndices;
        
        bool runSimpleED = config.useSimpleEmptyDrops || needsFallback;
        
        if (runSimpleED) {
            if (config.useSimpleEmptyDrops) {
                cerr << "  [Simple EmptyDrops] Running (force-enabled via --use-simple-empty-drops)" << endl;
            } else {
                // Join all fallback reasons
                string reasonStr;
                for (size_t r = 0; r < fallbackReasons.size(); r++) {
                    if (r > 0) reasonStr += "; ";
                    reasonStr += fallbackReasons[r];
                }
                cerr << "  [Simple EmptyDrops] Running as fallback: " << reasonStr << endl;
            }
            
            // Run Simple EmptyDrops bootstrap with all available cores
            // (tags are processed sequentially, so we can use all cores for bootstrap)
            SimpleEmptyDropsParams simpleParams = config.simpleEmptyDropsParams;
            simpleParams.nExpectedCells = 0;  // Let bootstrap estimate
            simpleParams.useBootstrap = true;
            simpleParams.umiMin = lowerTestingBound;
            simpleParams.maxThreads = 0;  // 0 = auto: use all available cores
            
            SimpleEmptyDropsResult simpleResult = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
                nUMIperCB_compact,
                static_cast<uint32_t>(nUMIperCB_compact.size()),
                simpleParams
            );
            
            retainThreshold = simpleResult.retainThreshold;
            
            cerr << "  [Simple EmptyDrops] Estimated recovered_cells: " << simpleParams.nExpectedCells << endl;
            cerr << "  [Simple EmptyDrops] Retain threshold: " << retainThreshold << " UMI" << endl;
            
            // Collect simple cell indices (cells with UMI >= retainThreshold)
            for (size_t ci = 0; ci < nUMIperCB_compact.size(); ci++) {
                if (nUMIperCB_compact[ci] >= retainThreshold) {
                    simplePasserIndices.insert(static_cast<uint32_t>(ci));
                }
            }
            
            cerr << "  [Simple EmptyDrops] Found " << simplePasserIndices.size() << " cells with UMI >= " << retainThreshold << endl;
            
            tagResult.nSimpleCells = static_cast<uint32_t>(simplePasserIndices.size());
        } else {
            // This branch only reached if !useSimpleEmptyDrops && !needsFallback
            debugOut("  [Simple EmptyDrops] Skipped (not needed as fallback)\n");
        }
        
        // Step 6: Count final passers (ED passers + Simple ED passers, avoiding duplicates)
        tagResult.nSimplePassers = 0;
        tagResult.nTailPassers = 0;
        
        for (const auto& result : tagResult.emptydropsResults) {
            if (result.passesRawP) {
                // Check if this cell is also a simple cell
                if (simplePasserIndices.count(result.cellIndex) > 0) {
                    tagResult.nSimplePassers++;  // Counted as simple passer
                } else {
                    tagResult.nTailPassers++;    // Pure ED rescue
                }
            }
        }
        
        // Add simple cells that weren't in ED candidates (UMI > retain but excluded from ED for some reason)
        // This shouldn't happen normally but handles edge cases
        uint32_t simpleOnlyPassers = 0;
        for (uint32_t simpleIdx : simplePasserIndices) {
            bool inEDResults = false;
            for (const auto& result : tagResult.emptydropsResults) {
                if (result.cellIndex == simpleIdx && result.passesRawP) {
                    inEDResults = true;
                    break;
                }
            }
            if (!inEDResults) {
                simpleOnlyPassers++;
            }
        }
        tagResult.nSimplePassers += simpleOnlyPassers;
        
        uint32_t totalPassers = tagResult.nSimplePassers + tagResult.nTailPassers;
        debugOut("  [Final] Total passers: " + to_string(totalPassers) + " (simple=" + to_string(tagResult.nSimplePassers) + ", tail=" + to_string(tagResult.nTailPassers) + ")\n");
        
        logTiming("EmptyDrops");
        
        // Collect all passers (ED passers + Simple ED passers, avoiding duplicates)
        unordered_set<uint32_t> passerSet;
        
        // Add ED passers
        for (const auto& result : tagResult.emptydropsResults) {
            if (result.passesRawP) {
                passerSet.insert(result.cellIndex);
            }
        }
        
        // Add Simple ED passers (if fallback ran)
        for (uint32_t simpleIdx : simplePasserIndices) {
            passerSet.insert(simpleIdx);
        }
        
        // Convert to vector
        vector<uint32_t> edPasserCompactIndices(passerSet.begin(), passerSet.end());
        sort(edPasserCompactIndices.begin(), edPasserCompactIndices.end());
        
        debugOut("  [Final] Total passers: " + to_string(edPasserCompactIndices.size()) + " (ED=" + to_string(nEDRescues) + ", SimpleED=" + to_string(simplePasserIndices.size()) + ")\n");
        
        // Store ED passers for combined occupancy filter (run after all tags complete)
        // Per-tag occupancy is DISABLED - combined filter runs at the end
        vector<uint32_t> finalPasserCompactIndices = edPasserCompactIndices;
        tagResult.occupancyRemoved = 0;  // Will be updated by combined filter
        
        // Store ED passer barcodes for combined occupancy filter
        tagResult.edPasserBarcodes.clear();
        tagResult.edPasserBarcodes.reserve(edPasserCompactIndices.size());
        for (uint32_t compactIdx : edPasserCompactIndices) {
            if (compactIdx < retainBarcodes.size()) {
                tagResult.edPasserBarcodes.push_back(retainBarcodes[compactIdx]);
            }
        }
        
        debugOut("  [Occupancy] Per-tag occupancy SKIPPED (combined filter runs after all tags)\n");
        debugOut("  [Occupancy] ED passers stored: " + to_string(tagResult.edPasserBarcodes.size()) + "\n");
        
        logTiming("OccupancyPostFilter");
        
        // Build final passing barcodes
        tagResult.passingBarcodes.clear();
        tagResult.filteredBarcodes.clear();
        for (uint32_t compactIdx : finalPasserCompactIndices) {
            if (compactIdx < retainBarcodes.size()) {
                tagResult.passingBarcodes.push_back(retainBarcodes[compactIdx]);
                tagResult.filteredBarcodes.push_back(retainBarcodes[compactIdx]);
            }
        }
        
        debugOut("  [Final] Passing barcodes: " + to_string(tagResult.passingBarcodes.size()) + "\n");
        
        outputs->tagResults[i] = std::move(tagResult);
    };  // End of processTag lambda
    
    // Launch worker threads
    vector<thread> workers;
    for (unsigned int t = 0; t < numThreads; t++) {
        workers.emplace_back([&]() {
            for (;;) {
                size_t tagIdx = nextTagIdx.fetch_add(1, memory_order_relaxed);
                if (tagIdx >= sampleTags.size()) break;
                processTag(tagIdx);
            }
        });
    }
    
    // Wait for all threads to complete
    for (auto& t : workers) {
        t.join();
    }
    
    cerr << "[FlexFilter] All " << sampleTags.size() << " tags processed" << endl;
    
    // ===== COMBINED OCCUPANCY FILTER (after all tags) =====
    // This matches CR's approach: run occupancy on ALL filtered cells from ALL tags
    if (!config.disableOccupancyFilter) {
        cerr << "[FlexFilter] Running combined occupancy filter..." << endl;
        
        // Collect ALL ED passers from ALL tags
        vector<string> allEdPasserBarcodes;
        vector<uint32_t> allEdPasserUMIs;
        
        // First pass: count total for pre-allocation
        size_t totalEdPassers = 0;
        for (const auto& tagResult : outputs->tagResults) {
            totalEdPassers += tagResult.edPasserBarcodes.size();
        }
        allEdPasserBarcodes.reserve(totalEdPassers);
        allEdPasserUMIs.reserve(totalEdPassers);
        
        // Build combined list
        for (const auto& tagResult : outputs->tagResults) {
            for (const auto& barcode : tagResult.edPasserBarcodes) {
                allEdPasserBarcodes.push_back(barcode);
                allEdPasserUMIs.push_back(1);  // UMI not used for unique tag counting
            }
        }
        
        cerr << "[FlexFilter] Combined ED passers: " << allEdPasserBarcodes.size() << endl;
        
        // Compute lambda estimate for Monte Carlo (from combined passers)
        double lambdaEstimate = 0.0;
        if (config.occupancyMode == OccupancyMode::MonteCarlo) {
            // Count cells per GEM (CB16) from all ED passers
            unordered_map<string, uint32_t> cellsPerGem;
            size_t cb16Len = 16;
            for (const auto& barcode : allEdPasserBarcodes) {
                if (barcode.size() >= cb16Len) {
                    string gem = barcode.substr(0, cb16Len);
                    cellsPerGem[gem]++;
                }
            }
            
            // Build distribution
            unordered_map<uint32_t, uint32_t> countDistribution;
            for (const auto& entry : cellsPerGem) {
                countDistribution[entry.second]++;
            }
            
            // Add empty GEMs
            uint32_t gemsWithCells = static_cast<uint32_t>(cellsPerGem.size());
            uint32_t totalExpectedGems = static_cast<uint32_t>(config.totalPartitions * config.recoveryFactor);
            uint32_t emptyGems = (totalExpectedGems > gemsWithCells) ? (totalExpectedGems - gemsWithCells) : 0;
            countDistribution[0] = emptyGems;
            
            // Compute weighted average
            double sumWeighted = 0.0;
            uint64_t sumWeights = 0;
            for (const auto& entry : countDistribution) {
                sumWeighted += static_cast<double>(entry.first) * entry.second;
                sumWeights += entry.second;
            }
            lambdaEstimate = (sumWeights > 0) ? (sumWeighted / sumWeights) : 0.0;
            
            cerr << "[FlexFilter] Combined occupancy lambda: " << lambdaEstimate 
                 << " (GEMs: " << gemsWithCells << ", empty: " << emptyGems << ")" << endl;
        }
        
        // Run occupancy filter on combined set
        OccupancyGuard::Config occConfig;
        occConfig.totalPartitions = config.totalPartitions;
        occConfig.recoveryFactor = config.recoveryFactor;
        occConfig.percentile = config.occupancyPercentile;
        
        vector<uint32_t> occupancyFlagged = OccupancyGuard::filterHighOccupancy(
            allEdPasserUMIs,
            allEdPasserBarcodes,
            occConfig,
            config.occupancyMode,
            lambdaEstimate,
            8,  // tagLength
            config.occupancySimulatedGems,
            config.emptydropsParams.seed
        );
        
        // Build set of flagged barcodes
        unordered_set<string> occupancyFlaggedBarcodes;
        for (uint32_t idx : occupancyFlagged) {
            if (idx < allEdPasserBarcodes.size()) {
                occupancyFlaggedBarcodes.insert(allEdPasserBarcodes[idx]);
            }
        }
        
        cerr << "[FlexFilter] Combined occupancy flagged " << occupancyFlaggedBarcodes.size() << " barcodes" << endl;
        
        // Update each tag's results: remove flagged barcodes
        uint32_t totalRemoved = 0;
        for (auto& tagResult : outputs->tagResults) {
            vector<string> filteredPassers;
            filteredPassers.reserve(tagResult.edPasserBarcodes.size());
            
            uint32_t removedFromTag = 0;
            for (const auto& barcode : tagResult.edPasserBarcodes) {
                if (occupancyFlaggedBarcodes.find(barcode) == occupancyFlaggedBarcodes.end()) {
                    filteredPassers.push_back(barcode);
                } else {
                    removedFromTag++;
                }
            }
            
            tagResult.occupancyRemoved = removedFromTag;
            tagResult.passingBarcodes = filteredPassers;
            tagResult.filteredBarcodes = filteredPassers;
            totalRemoved += removedFromTag;
            
            if (removedFromTag > 0) {
                cerr << "[FlexFilter] Tag " << tagResult.tag << ": removed " << removedFromTag 
                     << " by occupancy, final: " << filteredPassers.size() << endl;
            }
        }
        
        cerr << "[FlexFilter] Combined occupancy removed " << totalRemoved << " cells total" << endl;
    }
    
    // Debug: For BC006/BC008, compare against golden files if available
    if (config.debugTagLog) {
        // This will be handled by test harness comparing outputs
        debugOut("\n=== FlexFilter processing complete ===\n");
    }
    
    if (debugLog.is_open()) {
        debugLog.close();
    }
    
    return 0;
}

// Helper: Load MEX files from directory
static bool loadMEXFiles(
    const string& matrixDir,
    vector<string>& barcodes,
    vector<string>& features,
    SampleMatrixData& matrixData)
{
    // Load barcodes.tsv
    string barcodesPath = matrixDir + "/barcodes.tsv";
    ifstream barcodesFile(barcodesPath);
    if (!barcodesFile.is_open()) {
        cerr << "ERROR: Cannot open barcodes file: " << barcodesPath << endl;
        return false;
    }
    
    string line;
    while (getline(barcodesFile, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty()) {
            barcodes.push_back(line);
        }
    }
    barcodesFile.close();
    
    if (barcodes.empty()) {
        cerr << "ERROR: No barcodes found in " << barcodesPath << endl;
        return false;
    }
    
    // Load features.tsv
    string featuresPath = matrixDir + "/features.tsv";
    ifstream featuresFile(featuresPath);
    if (!featuresFile.is_open()) {
        cerr << "ERROR: Cannot open features file: " << featuresPath << endl;
        return false;
    }
    
    while (getline(featuresFile, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (!line.empty()) {
            // Extract first column (gene ID) - tab-separated
            size_t tabPos = line.find('\t');
            string geneId = (tabPos != string::npos) ? line.substr(0, tabPos) : line;
            features.push_back(geneId);
        }
    }
    featuresFile.close();
    
    if (features.empty()) {
        cerr << "ERROR: No features found in " << featuresPath << endl;
        return false;
    }
    
    // Load matrix.mtx (Matrix Market format)
    string matrixPath = matrixDir + "/matrix.mtx";
    ifstream matrixFile(matrixPath);
    if (!matrixFile.is_open()) {
        cerr << "ERROR: Cannot open matrix file: " << matrixPath << endl;
        return false;
    }
    
    // Skip header lines (starting with %)
    while (getline(matrixFile, line)) {
        if (line.empty() || line[0] != '%') {
            break;
        }
    }
    
    // Parse dimensions line: num_genes num_cells num_nonzero
    uint32_t numGenes, numCells, numNonzero;
    stringstream dimStream(line);
    dimStream >> numGenes >> numCells >> numNonzero;
    
    if (numGenes != features.size()) {
        cerr << "WARNING: Matrix num_genes (" << numGenes << ") != features size (" << features.size() << ")" << endl;
    }
    if (numCells != barcodes.size()) {
        cerr << "WARNING: Matrix num_cells (" << numCells << ") != barcodes size (" << barcodes.size() << ")" << endl;
    }
    
    // Initialize matrix data structures
    matrixData.nCells = numCells;
    matrixData.nGenes = numGenes;
    matrixData.countMatStride = 3;
    matrixData.barcodes = barcodes;
    matrixData.features = features;
    
    // Initialize per-cell data
    matrixData.nUMIperCB.resize(numCells, 0);
    matrixData.nGenePerCB.resize(numCells, 0);
    matrixData.countCellGeneUMIindex.resize(numCells + 1, 0);
    
    // Temporary storage: map cell_idx -> [(gene_idx, count)]
    vector<vector<pair<uint32_t, uint32_t>>> cellGeneMap(numCells);
    
    // Read triplets: gene_row cell_col count (1-based indices in MTX)
    uint32_t geneRow, cellCol, count;
    while (matrixFile >> geneRow >> cellCol >> count) {
        if (geneRow == 0 || cellCol == 0) {
            cerr << "ERROR: Matrix contains 0-based indices (expected 1-based)" << endl;
            return false;
        }
        
        // Convert to 0-based
        uint32_t geneIdx = geneRow - 1;
        uint32_t cellIdx = cellCol - 1;
        
        if (cellIdx >= numCells || geneIdx >= numGenes) {
            cerr << "ERROR: Matrix indices out of bounds: cell=" << cellIdx << " gene=" << geneIdx << endl;
            return false;
        }
        
        // Store in temporary map
        cellGeneMap[cellIdx].push_back({geneIdx, count});
        matrixData.nUMIperCB[cellIdx] += count;
    }
    matrixFile.close();
    
    // Build sparse matrix in countCellGeneUMI format
    uint32_t offset = 0;
    for (uint32_t cellIdx = 0; cellIdx < numCells; cellIdx++) {
        matrixData.countCellGeneUMIindex[cellIdx] = offset;
        matrixData.nGenePerCB[cellIdx] = cellGeneMap[cellIdx].size();
        
        for (const auto& geneCount : cellGeneMap[cellIdx]) {
            matrixData.countCellGeneUMI.push_back(geneCount.first);   // geneIdx
            matrixData.countCellGeneUMI.push_back(geneCount.second);  // count
            matrixData.countCellGeneUMI.push_back(0);                 // reserved/multi
            offset += 3;
        }
    }
    matrixData.countCellGeneUMIindex[numCells] = offset;
    
    cout << "Loaded MEX files:\n";
    cout << "  Barcodes: " << barcodes.size() << "\n";
    cout << "  Features: " << features.size() << "\n";
    cout << "  Non-zero entries: " << numNonzero << "\n";
    
    return true;
}

void FlexFilter::deriveSampleWhitelist(
    const vector<string>& compositeBarcodes,
    vector<string>& sampleLabels,
    vector<string>& sampleTags)
{
    unordered_map<string, int> tagCounts;
    
    for (const auto& barcode : compositeBarcodes) {
        if (barcode.size() >= 24) {
            string tag = barcode.substr(16, 8);
            tagCounts[tag]++;
        } else if (barcode.size() == 16) {
            continue;
        } else {
            cerr << "WARNING: Unexpected barcode length: " << barcode << " (" << barcode.size() << " chars)" << endl;
        }
    }
    
    vector<pair<string, int>> tagVec(tagCounts.begin(), tagCounts.end());
    sort(tagVec.begin(), tagVec.end(), [](const pair<string, int>& a, const pair<string, int>& b) {
        return a.second > b.second || (a.second == b.second && a.first < b.first);
    });
    
    for (size_t i = 0; i < tagVec.size(); i++) {
        const string& tag = tagVec[i].first;
        stringstream labelStream;
        labelStream << "BC" << setw(3) << setfill('0') << (i + 1);
        sampleLabels.push_back(labelStream.str());
        sampleTags.push_back(tag);
    }
    
    cout << "Auto-derived sample whitelist:\n";
    for (size_t i = 0; i < sampleLabels.size(); i++) {
        cout << "  " << sampleLabels[i] << ": " << sampleTags[i]
             << " (" << tagVec[i].second << " barcodes)" << endl;
    }
}

// File-based entry point
int FlexFilter::runFromFiles(const FileInputs& inputs, Outputs* outputs, const Config& config) {
    // Auto-populate config with external defaults for any zero/unset fields
    Config populatedConfig = config;
    populateConfigWithDefaults(populatedConfig);
    
    // Load MEX files
    vector<string> barcodes;
    vector<string> features;
    SampleMatrixData matrixData;
    
    if (!loadMEXFiles(inputs.matrixDir, barcodes, features, matrixData)) {
        cerr << "ERROR: Failed to load MEX files from " << inputs.matrixDir << endl;
        return 1;
    }
    
    // Use explicit sample whitelist if provided, otherwise auto-derive
    vector<string> sampleLabels;
    vector<string> sampleTags;
    
    if (!inputs.sampleLabels.empty() && !inputs.sampleTags.empty()) {
        // Use explicit whitelist
        sampleLabels = inputs.sampleLabels;
        sampleTags = inputs.sampleTags;
        cout << "Using explicit sample whitelist:\n";
        for (size_t i = 0; i < sampleLabels.size(); i++) {
            cout << "  " << sampleLabels[i] << ": " << sampleTags[i] << "\n";
        }

        // For debugging: derive tags from barcodes and print for comparison
        vector<string> autoLabels;
        vector<string> autoTags;
        deriveSampleWhitelist(barcodes, autoLabels, autoTags);
        cout << "Auto-derived sample whitelist (for comparison):\n";
        for (size_t i = 0; i < autoLabels.size(); i++) {
            cout << "  " << autoLabels[i] << ": " << autoTags[i];
            if (i < sampleTags.size() && i < sampleLabels.size()) {
                cout << " (explicit: " << sampleLabels[i] << ": " << sampleTags[i] << ")";
            }
            cout << "\n";
        }
    } else {
        // Auto-derive from barcodes
        deriveSampleWhitelist(barcodes, sampleLabels, sampleTags);
    }
    
    // Filter by allowed tags if specified (only applies if auto-derived)
    if (!inputs.allowedTags.empty() && inputs.sampleLabels.empty()) {
        unordered_set<string> allowedSet(inputs.allowedTags.begin(), inputs.allowedTags.end());
        vector<string> filteredLabels;
        vector<string> filteredTags;
        
        for (size_t i = 0; i < sampleTags.size(); i++) {
            if (allowedSet.count(sampleTags[i])) {
                filteredLabels.push_back(sampleLabels[i]);
                filteredTags.push_back(sampleTags[i]);
            }
        }
        
        sampleLabels = filteredLabels;
        sampleTags = filteredTags;
        
        cout << "Filtered to " << sampleLabels.size() << " allowed tags\n";
    }
    
    if (sampleLabels.empty()) {
        cerr << "ERROR: No samples to process (check allowed tags filter)" << endl;
        return 1;
    }
    
    // Call in-memory entry point
    return runInternal(
        matrixData,
        sampleLabels,
        sampleTags,
        inputs.tagExpectedCells,
        outputs,
        populatedConfig
    );
}

// In-memory entry point
int FlexFilter::runFromMemory(const MemoryInputs& inputs, Outputs* outputs, const Config& config) {
    // Auto-populate config with external defaults for any zero/unset fields
    Config populatedConfig = config;
    populateConfigWithDefaults(populatedConfig);
    
    return runInternal(
        inputs.matrixData,
        inputs.sampleLabels,
        inputs.sampleTags,
        inputs.tagExpectedCells,
        outputs,
        populatedConfig
    );
}
