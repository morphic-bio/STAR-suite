#include "OrdMagStage.h"
#include "ParametersSolo.h"
#include "streamFuns.h"
#include "pcg_random.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <numeric>
#include <iostream>
#include <thread>
#include <vector>
using namespace std;

// Find number of cells within order of magnitude of baseline
// Matches Python: find_within_ordmag(x, baseline_idx)
uint32 SimpleEmptyDropsStage::findWithinOrdmag(
    const vector<uint32>& counts,
    uint32 baselineIdx
) {
    if (counts.empty()) return 0;
    
    // Sort ascending (like Python np.sort)
    vector<uint32> sorted = counts;
    sort(sorted.begin(), sorted.end());
    
    // Get baseline from the high end (like Python x_ascending[-(baseline_idx + 1)])
    uint32 n = sorted.size();
    if (baselineIdx >= n) baselineIdx = n - 1;
    uint32 baselinePos = n - baselineIdx - 1;
    uint32 baseline = sorted[baselinePos];
    
    // Cutoff = max(1, round(0.1 * baseline))
    uint32 cutoff = max((uint32)1, (uint32)round(0.1 * baseline));
    
    // Count cells >= cutoff (like Python: len(x) - np.searchsorted(x_ascending, cutoff))
    // searchsorted finds first position where cutoff would be inserted
    auto it = lower_bound(sorted.begin(), sorted.end(), cutoff);
    uint32 countAbove = n - (it - sorted.begin());
    
    return countAbove;
}

// Estimate recovered cells by minimizing loss
// Matches Python: estimate_recovered_cells_ordmag(nonzero_bc_counts, max_expected_cells)
// Using float precision for speed - sufficient for loss calculations
pair<uint32, double> SimpleEmptyDropsStage::estimateRecoveredCellsOrdmag(
    const vector<uint32>& nonzeroCounts,
    uint32 maxExpectedCells,
    double recoveredCellsQuantile
) {
    if (nonzeroCounts.empty()) return make_pair(0, 0.0);
    
    uint32 n = nonzeroCounts.size();
    
    // Generate log2-spaced range of recovered_cells values (1 to maxExpectedCells)
    // Python: recovered_cells = np.linspace(1, np.log2(max_expected_cells), 2000)
    // Python: recovered_cells = np.unique(np.round(np.power(2, recovered_cells)).astype(int))
    vector<uint32> recoveredCellsOptions;
    float log2Max = log2f((float)maxExpectedCells);
    float quantileF = (float)recoveredCellsQuantile;
    for (int i = 0; i < 2000; i++) {
        float log2Val = 1.0f + (log2Max - 1.0f) * i / 1999.0f;
        uint32 val = (uint32)roundf(powf(2.0f, log2Val));
        if (recoveredCellsOptions.empty() || val != recoveredCellsOptions.back()) {
            recoveredCellsOptions.push_back(val);
        }
    }
    
    // Search for best recovered_cells (float precision sufficient for loss)
    uint32 bestRecoveredCells = 1;
    float bestLoss = 1e30f;
    
    for (uint32 recoveredCells : recoveredCellsOptions) {
        // baseline_bc_idx = round(recovered_cells * (1 - quantile))
        uint32 baselineIdx = (uint32)roundf(recoveredCells * (1.0f - quantileF));
        if (baselineIdx >= n) baselineIdx = n - 1;
        
        // Get filtered cells count
        uint32 filteredCells = findWithinOrdmag(nonzeroCounts, baselineIdx);
        
        // Loss = (filtered - recovered)^2 / recovered (float precision)
        float diff = (float)filteredCells - (float)recoveredCells;
        float loss = (diff * diff) / (float)recoveredCells;
        
        if (loss < bestLoss) {
            bestLoss = loss;
            bestRecoveredCells = recoveredCells;
        }
    }
    
    return make_pair(bestRecoveredCells, (double)bestLoss);
}

// Run Cell Ranger-style filtering with bootstrap
// Matches Python: filter_cellular_barcodes_ordmag()
OrdMagResult SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
    const vector<uint32>& nUMIperCB,
    uint32 nCB,
    OrdMagParams& params
) {
    OrdMagResult result;
    
    if (nCB == 0) {
        result.retainThreshold = 0;
        result.nCellsSimple = 0;
        result.minUMI = 0;
        result.medianVal = 0;
        result.candidateLastRank = 0;
        result.ambientRange = make_pair(0, 0);
        return result;
    }
    
    // Random state (fixed seed for reproducibility, like Python np.random.RandomState(0))
    // Using PCG32 for speed (much faster than mt19937)
    pcg32 rng(0);
    
    // Get non-zero counts
    vector<uint32> nonzeroCounts;
    for (uint32 i = 0; i < nCB; i++) {
        if (nUMIperCB[i] > 0) {
            nonzeroCounts.push_back(nUMIperCB[i]);
        }
    }
    
    if (nonzeroCounts.empty()) {
        cerr << "WARNING: All barcodes have zero counts for ordmag" << endl;
        result.retainThreshold = 0;
        result.nCellsSimple = 0;
        return result;
    }
    
    uint32 nNonzero = nonzeroCounts.size();
    
    // Determine maxExpectedCells if not set
    uint32 maxExpectedCells = params.maxExpectedCells;
    if (maxExpectedCells == 0) {
        // Default: use indMin/2 as reasonable max (like CR uses empty_drops_range[0])
        maxExpectedCells = min(params.indMin / 2, (uint32)262144);
        if (maxExpectedCells < 1000) maxExpectedCells = 90000;  // Fallback
    }
    
    // Step 1: Estimate recovered_cells if not provided
    uint32 recoveredCells = params.nExpectedCells;
    
    if (recoveredCells == 0 && params.useBootstrap) {
        // Bootstrap to estimate recovered_cells (parallel with std::thread)
        uint32 nBoot = params.nBootstrapSamples;
        vector<uint32> bootRecovered(nBoot);
        vector<double> bootLoss(nBoot);
        
        // Determine thread count: use maxThreads if set, else auto-detect
        uint32 nThreads = params.maxThreads;
        if (nThreads == 0) {
            nThreads = thread::hardware_concurrency();
            if (nThreads == 0) nThreads = 4;  // Fallback
            // Respect OMP_NUM_THREADS if set
            const char* ompEnv = getenv("OMP_NUM_THREADS");
            if (ompEnv) {
                int envThreads = atoi(ompEnv);
                if (envThreads > 0) nThreads = (uint32)envThreads;
            }
        }
        nThreads = min(nThreads, nBoot);
        vector<thread> threads;
        
        uint32 baseSeed = (params.bootstrapSeed > 0) ? params.bootstrapSeed : 1;
        
        auto bootstrapWorker = [&nonzeroCounts, nNonzero, maxExpectedCells, &params, &bootRecovered, &bootLoss](uint32 startIdx, uint32 endIdx, uint32 seed) {
            // PCG32 is ~3x faster than mt19937 for uniform int sampling
            pcg32 localRng(seed, startIdx);  // seed + stream for uniqueness
            for (uint32 b = startIdx; b < endIdx; b++) {
                vector<uint32> bootstrap(nNonzero);
                for (uint32 i = 0; i < nNonzero; i++) {
                    // Fast bounded random using PCG's optimized method
                    bootstrap[i] = nonzeroCounts[localRng(nNonzero)];
                }
                pair<uint32, double> estResult = estimateRecoveredCellsOrdmag(
                    bootstrap, maxExpectedCells, params.recoveredCellsQuantile);
                bootRecovered[b] = estResult.first;
                bootLoss[b] = estResult.second;
            }
        };
        
        uint32 chunkSize = (nBoot + nThreads - 1) / nThreads;
        for (uint32 t = 0; t < nThreads; t++) {
            uint32 startIdx = t * chunkSize;
            uint32 endIdx = min(startIdx + chunkSize, nBoot);
            if (startIdx < endIdx) {
                threads.emplace_back(bootstrapWorker, startIdx, endIdx, baseSeed + t * 1000);
            }
        }
        for (auto& th : threads) th.join();
        
        // Sum results
        double sumRecovered = 0.0, sumLoss = 0.0;
        for (uint32 b = 0; b < nBoot; b++) {
            sumRecovered += bootRecovered[b];
            sumLoss += bootLoss[b];
        }
        
        recoveredCells = (uint32)round(sumRecovered / nBoot);
        recoveredCells = max(recoveredCells, params.minRecoveredCells);
        double avgLoss = sumLoss / nBoot;
        
        cout << "Found recovered_cells = " << recoveredCells << " with loss = " << avgLoss << endl;
        
        // Update params for output
        params.nExpectedCells = recoveredCells;
    } else if (recoveredCells == 0) {
        // No bootstrap, use default
        recoveredCells = 3000;
        params.nExpectedCells = recoveredCells;
        cout << "Using default recovered_cells = " << recoveredCells << endl;
    } else {
        recoveredCells = max(recoveredCells, params.minRecoveredCells);
        cout << "Using provided recovered_cells = " << recoveredCells << endl;
    }
    
    // Step 2: Compute baseline index
    uint32 baselineIdx = (uint32)round((double)recoveredCells * (1.0 - params.recoveredCellsQuantile));
    if (baselineIdx >= nNonzero) baselineIdx = nNonzero - 1;
    
    // Step 3: Bootstrap to get top_n with variance (parallel with std::thread)
    vector<uint32> topNBoot(params.nBootstrapSamples);
    
    if (params.useBootstrap) {
        uint32 nBoot = params.nBootstrapSamples;
        // Determine thread count: use maxThreads if set, else auto-detect
        uint32 nThreads = params.maxThreads;
        if (nThreads == 0) {
            nThreads = thread::hardware_concurrency();
            if (nThreads == 0) nThreads = 4;
            const char* ompEnv = getenv("OMP_NUM_THREADS");
            if (ompEnv) {
                int envThreads = atoi(ompEnv);
                if (envThreads > 0) nThreads = (uint32)envThreads;
            }
        }
        nThreads = min(nThreads, nBoot);
        vector<thread> threads;
        uint32 baseSeed2 = (params.bootstrapSeed > 0) ? (params.bootstrapSeed + 10000) : 100;
        
        auto topNWorker = [&nonzeroCounts, nNonzero, baselineIdx, &topNBoot](uint32 startIdx, uint32 endIdx, uint32 seed) {
            // PCG32 is ~3x faster than mt19937 for uniform int sampling
            pcg32 localRng(seed, startIdx + 10000);  // seed + stream for uniqueness
            for (uint32 b = startIdx; b < endIdx; b++) {
                vector<uint32> bootstrap(nNonzero);
                for (uint32 i = 0; i < nNonzero; i++) {
                    // Fast bounded random using PCG's optimized method
                    bootstrap[i] = nonzeroCounts[localRng(nNonzero)];
                }
                topNBoot[b] = findWithinOrdmag(bootstrap, baselineIdx);
            }
        };
        
        uint32 chunkSize = (nBoot + nThreads - 1) / nThreads;
        for (uint32 t = 0; t < nThreads; t++) {
            uint32 startIdx = t * chunkSize;
            uint32 endIdx = min(startIdx + chunkSize, nBoot);
            if (startIdx < endIdx) {
                threads.emplace_back(topNWorker, startIdx, endIdx, baseSeed2 + t * 1000);
            }
        }
        for (auto& th : threads) th.join();
    } else {
        // No bootstrap - just run once on actual data
        uint32 topN = findWithinOrdmag(nonzeroCounts, baselineIdx);
        for (uint32 b = 0; b < params.nBootstrapSamples; b++) {
            topNBoot[b] = topN;
        }
    }
    
    // Step 4: Summarize bootstrap results (matches summarize_bootstrapped_top_n)
    double sumTopN = 0.0;
    for (uint32 t : topNBoot) sumTopN += t;
    double meanTopN = sumTopN / params.nBootstrapSamples;
    
    double sumSqDiff = 0.0;
    for (uint32 t : topNBoot) {
        double diff = t - meanTopN;
        sumSqDiff += diff * diff;
    }
    double varTopN = sumSqDiff / params.nBootstrapSamples;
    double sdTopN = sqrt(varTopN);
    
    // Round to get number of cells
    uint32 nCellsSimple = (uint32)round(meanTopN);
    if (nCellsSimple > nNonzero) nCellsSimple = nNonzero;
    
    cout << "Bootstrap mean = " << meanTopN << ", sd = " << sdTopN << ", nCellsSimple = " << nCellsSimple << endl;
    
    // Step 5: Get actual filtered barcodes (top nCellsSimple by UMI)
    // Create index array sorted by UMI count (descending)
    typedef struct {uint32 index; uint32 count;} IndCount;
    vector<IndCount> indCount(nCB);
    for (uint32 ii = 0; ii < nCB; ii++) {
        indCount[ii].index = ii;
        indCount[ii].count = nUMIperCB[ii];
    }
    
    sort(indCount.begin(), indCount.end(), [](const IndCount& ic1, const IndCount& ic2) {
        return (ic1.count > ic2.count) || (ic1.count == ic2.count && ic1.index < ic2.index);
    });
    
    // Ensure we select all barcodes with count >= cutoff (handle ties)
    // Python: Make sure all barcodes with count x are selected if any with count >= x is selected
    uint32 cutoffCount = 0;
    if (nCellsSimple > 0 && nCellsSimple <= nCB) {
        cutoffCount = indCount[nCellsSimple - 1].count;
        // Extend to include ties
        while (nCellsSimple < nCB && indCount[nCellsSimple].count == cutoffCount) {
            // Check if we're grabbing too many (> 20% more than original)
            if (nCellsSimple > (uint32)(meanTopN * 1.2)) {
                nCellsSimple = (uint32)round(meanTopN);  // Revert to original
                break;
            }
            nCellsSimple++;
        }
    }
    
    // Compute retain threshold (UMI of last passing cell)
    uint32 retainThreshold = (nCellsSimple > 0 && nCellsSimple <= nCB) 
        ? indCount[nCellsSimple - 1].count : 0;
    
    // Extract passing indices
    for (uint32 ii = 0; ii < nCellsSimple; ii++) {
        result.passingIndices.push_back(indCount[ii].index);
    }
    
    // Compute median value
    uint32 medianVal;
    if (nCellsSimple == 0) {
        medianVal = (nCB > 0) ? indCount[0].count : 0;
    } else {
        uint32 medianRank = max((uint32)floor((double)nCellsSimple / 2.0), (uint32)1);
        medianRank = min(medianRank, nCB);
        medianVal = indCount[medianRank - 1].count;
    }
    
    // Compute minUMI for candidates (matches EmptyDrops lower bound)
    uint32 minUMIFromFrac = (uint32)round(params.umiMinFracMedian * medianVal);
    uint32 minUMI = max(params.umiMin, minUMIFromFrac);
    
    // Compute candidate limit
    uint32 candLimit = 0;
    for (uint32 ii = 0; ii < nCB; ii++) {
        if (nUMIperCB[ii] > minUMI) {
            candLimit++;
        }
    }
    candLimit = max(candLimit, nCellsSimple);
    uint32 iCandLast = min(nCellsSimple + params.candMaxN, min(candLimit, nCB));
    
    // Extract candidate indices
    for (uint32 ii = 0; ii < iCandLast; ii++) {
        result.candidateIndices.push_back(indCount[ii].index);
    }
    
    // Extract ambient indices (same logic as non-bootstrap version)
    uint32 scaledIndMin = min(params.indMin, nCB);
    uint32 scaledIndMax = min(params.indMax, nCB);
    uint32 minAmbientCells = (nCB >= 1000) ? (nCB / 10) : min((uint32)100, nCB);
    uint32 ambientWindowSize = (scaledIndMax >= scaledIndMin) ? (scaledIndMax - scaledIndMin + 1) : 0;
    
    if (scaledIndMax < scaledIndMin || ambientWindowSize < minAmbientCells || scaledIndMin == scaledIndMax) {
        uint32 fallbackSize = min(minAmbientCells, nCB);
        uint32 fallbackStart = (nCB >= fallbackSize) ? (nCB - fallbackSize) : 0;
        result.ambientIndices.clear();
        for (uint32 ii = fallbackStart; ii < nCB; ii++) {
            result.ambientIndices.push_back(indCount[ii].index);
        }
        result.ambientRange = make_pair(fallbackStart + 1, nCB);
    } else {
        uint32 ambientStart = (scaledIndMin > 0) ? scaledIndMin - 1 : 0;
        uint32 ambientEnd = scaledIndMax;
        if (ambientEnd > nCB) ambientEnd = nCB;
        result.ambientIndices.clear();
        for (uint32 ii = ambientStart; ii < ambientEnd; ii++) {
            if (ii < nCB) {
                result.ambientIndices.push_back(indCount[ii].index);
            }
        }
        result.ambientRange = make_pair(scaledIndMin, ambientEnd);
    }
    
    // Store results
    result.retainThreshold = retainThreshold;
    result.nCellsSimple = nCellsSimple;
    result.minUMI = minUMI;
    result.medianVal = medianVal;
    result.candidateLastRank = iCandLast;
    
    cout << "Median UMIs of initial cell calls: " << medianVal << endl;
    cout << "Min UMIs: " << minUMI << endl;
    
    return result;
}

OrdMagResult SimpleEmptyDropsStage::runCRSimpleFilter(
    const vector<uint32>& nUMIperCB,
    uint32 nCB,
    const OrdMagParams& params
) {
    OrdMagResult result;
    
    if (nCB == 0) {
        result.retainThreshold = 0;
        result.nCellsSimple = 0;
        result.minUMI = 0;
        result.medianVal = 0;
        result.candidateLastRank = 0;
        result.ambientRange = make_pair(0, 0);
        return result;
    }
    
    // Create index array sorted by UMI count (descending)
    typedef struct {uint32 index; uint32 count;} IndCount;
    vector<IndCount> indCount(nCB);
    for (uint32 ii = 0; ii < nCB; ii++) {
        indCount[ii].index = ii;
        indCount[ii].count = nUMIperCB[ii];
    }
    
    sort(indCount.begin(), indCount.end(), [](const IndCount& ic1, const IndCount& ic2) {
        return (ic1.count > ic2.count) || (ic1.count == ic2.count && ic1.index < ic2.index);
    });
    
    // Extract sorted totals
    vector<uint32> totalsSorted(nCB);
    for (uint32 ii = 0; ii < nCB; ii++) {
        totalsSorted[ii] = indCount[ii].count;
    }
    
    // Compute robust max index
    uint32 nExpectedCells = max(params.nExpectedCells, (uint32)1);
    uint32 maxInd = (uint32)round(nExpectedCells * (1.0 - params.maxPercentile));
    if (maxInd < 1) maxInd = 1;
    if (maxInd > nCB) maxInd = nCB;
    
    // Compute retain threshold
    uint32 nUMImax = totalsSorted[maxInd - 1];
    uint32 retain = max((uint32)round((double)nUMImax / params.maxMinRatio), (uint32)1);
    
    // Count cells passing simple filter
    uint32 ncellsSimple = 0;
    for (uint32 ii = 0; ii < nCB; ii++) {
        if (nUMIperCB[ii] >= retain) {
            ncellsSimple++;
        }
    }
    ncellsSimple = max(min(ncellsSimple, nCB), (uint32)0);
    
    // Fallback if no cells pass
    if (ncellsSimple == 0) {
        ncellsSimple = min(nExpectedCells, nCB);
        retain = totalsSorted[min(ncellsSimple, nCB) - 1];
        // Recompute ncellsSimple
        ncellsSimple = 0;
        for (uint32 ii = 0; ii < nCB; ii++) {
            if (nUMIperCB[ii] >= retain) {
                ncellsSimple++;
            }
        }
        ncellsSimple = max(min(ncellsSimple, nCB), (uint32)0);
    }
    
    // Compute median value
    uint32 medianVal;
    if (ncellsSimple == 0) {
        medianVal = totalsSorted[0];
    } else {
        uint32 medianRank = max((uint32)floor((double)ncellsSimple / 2.0), (uint32)1);
        medianRank = min(medianRank, nCB);
        medianVal = totalsSorted[medianRank - 1];
    }
    
    // Compute minUMI for candidates
    uint32 minUMIFromFrac = (uint32)round(params.umiMinFracMedian * medianVal);
    uint32 minUMI = max(params.umiMin, minUMIFromFrac);
    
    // Compute candidate limit
    uint32 candLimit = 0;
    for (uint32 ii = 0; ii < nCB; ii++) {
        if (nUMIperCB[ii] > minUMI) {
            candLimit++;
        }
    }
    candLimit = max(candLimit, ncellsSimple);
    uint32 iCandLast = min(ncellsSimple + params.candMaxN, min(candLimit, nCB));
    
    // Extract candidate indices (simple cells + tail) - matches Python behavior
    // Python: candidate_idx = np.concatenate([simple_idx, tail_idx])
    for (uint32 ii = 0; ii < iCandLast; ii++) {
        result.candidateIndices.push_back(indCount[ii].index);
    }
    
    // Extract passing indices (simple filter)
    for (uint32 ii = 0; ii < ncellsSimple; ii++) {
        result.passingIndices.push_back(indCount[ii].index);
    }
    
    // Extract ambient indices (1-based ranks ind_min .. ind_max inclusive)
    // Scale ambient range to tag matrix size: clamp to [min(nCells, indMin), min(nCells, indMax)]
    // This ensures we get a real ambient pool even for tag-specific matrices
    uint32 scaledIndMin = min(params.indMin, nCB);
    uint32 scaledIndMax = min(params.indMax, nCB);
    
    // Ensure we have at least some ambient cells (10% of matrix or at least 100 cells)
    uint32 minAmbientCells = (nCB >= 1000) ? (nCB / 10) : min((uint32)100, nCB);
    
    // Check if scaled range is valid and has enough cells
    // Note: scaledIndMin/scaledIndMax are 1-based inclusive, so window size = max - min + 1
    uint32 ambientWindowSize = (scaledIndMax >= scaledIndMin) ? (scaledIndMax - scaledIndMin + 1) : 0;
    
    // Always use fallback if range is invalid or too small (ensures we get enough ambient cells)
    if (scaledIndMax < scaledIndMin || ambientWindowSize < minAmbientCells || scaledIndMin == scaledIndMax) {
        // Fallback: use lowest-ranked cells (10% of matrix or at least 100 cells)
        uint32 fallbackSize = min(minAmbientCells, nCB);
        uint32 fallbackStart = (nCB >= fallbackSize) ? (nCB - fallbackSize) : 0;
        result.ambientIndices.clear();
        for (uint32 ii = fallbackStart; ii < nCB; ii++) {
            result.ambientIndices.push_back(indCount[ii].index);
        }
        result.ambientRange = make_pair(fallbackStart + 1, nCB);
    } else {
        // Use scaled range: convert 1-based inclusive to 0-based indices
        uint32 ambientStart = (scaledIndMin > 0) ? scaledIndMin - 1 : 0;
        uint32 ambientEnd = scaledIndMax;  // 1-based inclusive, so loop to ambientEnd (inclusive)
        // Clamp to valid range
        if (ambientEnd > nCB) ambientEnd = nCB;
        result.ambientIndices.clear();
        for (uint32 ii = ambientStart; ii < ambientEnd; ii++) {
            if (ii < nCB) {
                result.ambientIndices.push_back(indCount[ii].index);
            }
        }
        result.ambientRange = make_pair(scaledIndMin, ambientEnd);
    }
    
    // Store results
    result.retainThreshold = retain;
    result.nCellsSimple = ncellsSimple;
    result.minUMI = minUMI;
    result.medianVal = medianVal;
    result.candidateLastRank = iCandLast;
    
    return result;
}

void SimpleEmptyDropsStage::writeOutputs(
    const OrdMagResult& result,
    const vector<string>& barcodes,
    const string& outputDir,
    const OrdMagParams& params,
    bool writeFilteredMatrix
) {
    // Create output directory
    string cmd = "mkdir -p " + outputDir;
    system(cmd.c_str());
    
    // Write passing_barcodes.txt
    string passingFile = outputDir + "/passing_barcodes.txt";
    ofstream passingOut(passingFile);
    if (!passingOut.is_open()) {
        return; // Error - should log
    }
    
    for (auto idx : result.passingIndices) {
        if (idx < barcodes.size()) {
            passingOut << barcodes[idx] << "\n";
        }
    }
    passingOut.close();
    
    // Write filter_summary.json
    string summaryFile = outputDir + "/filter_summary.json";
    ofstream summaryOut(summaryFile);
    if (!summaryOut.is_open()) {
        return;
    }
    
    summaryOut << "{\n";
    summaryOut << "  \"retain_threshold\": " << result.retainThreshold << ",\n";
    summaryOut << "  \"ncells_simple\": " << result.nCellsSimple << ",\n";
    summaryOut << "  \"n_candidates\": " << result.candidateIndices.size() << ",\n";
    summaryOut << "  \"n_ambient\": " << result.ambientIndices.size() << ",\n";
    summaryOut << "  \"min_umi\": " << result.minUMI << ",\n";
    summaryOut << "  \"median_val\": " << result.medianVal << ",\n";
    summaryOut << "  \"candidate_last_rank\": " << result.candidateLastRank << ",\n";
    summaryOut << "  \"ambient_range\": [" << result.ambientRange.first << ", " << result.ambientRange.second << "],\n";
    summaryOut << "  \"parameters\": {\n";
    summaryOut << "    \"n_expected_cells\": " << params.nExpectedCells << ",\n";
    summaryOut << "    \"max_percentile\": " << fixed << setprecision(6) << params.maxPercentile << ",\n";
    summaryOut << "    \"max_min_ratio\": " << params.maxMinRatio << ",\n";
    summaryOut << "    \"umi_min\": " << params.umiMin << ",\n";
    summaryOut << "    \"umi_min_frac_median\": " << params.umiMinFracMedian << ",\n";
    summaryOut << "    \"cand_max_n\": " << params.candMaxN << ",\n";
    summaryOut << "    \"ind_min\": " << params.indMin << ",\n";
    summaryOut << "    \"ind_max\": " << params.indMax << "\n";
    summaryOut << "  }\n";
    summaryOut << "}\n";
    summaryOut.close();
}

OrdMagParams SimpleEmptyDropsStage::loadParamsFromSolo(const ParametersSolo& pSolo) {
    OrdMagParams params;
    // Use EmptyDrops_CR parameters (these match the Python defaults)
    params.nExpectedCells = (uint32)pSolo.cellFilter.knee.nExpectedCells;
    params.maxPercentile = pSolo.cellFilter.knee.maxPercentile;
    params.maxMinRatio = pSolo.cellFilter.knee.maxMinRatio;
    params.umiMin = pSolo.cellFilter.eDcr.umiMin;
    params.umiMinFracMedian = pSolo.cellFilter.eDcr.umiMinFracMedian;
    params.candMaxN = pSolo.cellFilter.eDcr.candMaxN;
    params.indMin = pSolo.cellFilter.eDcr.indMin;
    params.indMax = pSolo.cellFilter.eDcr.indMax;
    
    // Note: Flex-specific parameter overrides are now handled by FlexDriver directly
    // This function is kept for compatibility with Solo integration
    
    return params;
}
