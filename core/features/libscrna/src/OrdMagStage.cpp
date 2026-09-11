/**
 * @file OrdMagStage.cpp
 * @brief Implementation of Simple EmptyDrops (OrdMag) filtering
 */

#include "OrdMagStage.h"
#include "OrdMagRank.h"
#include "EmptyDropsMultinomial.h"
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
#include "ParallelTasks.h"
#include <vector>
#include <stdexcept>

using namespace std;

namespace {
// The estimator evaluates many baseline indices on the same bootstrap draw.
// Keep its rounding and inclusive threshold rule, but reuse one sorted copy.
uint32 findWithinSortedOrdmag(const vector<uint32>& sorted, uint32 baselineIdx) {
    if (sorted.empty()) return 0;
    const uint32 n = sorted.size();
    if (baselineIdx >= n) baselineIdx = n - 1;
    const uint32 baseline = sorted[n - baselineIdx - 1];
    const uint32 cutoff = max((uint32)1, (uint32)round(0.1 * baseline));
    return n - (lower_bound(sorted.begin(), sorted.end(), cutoff) - sorted.begin());
}
}

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
    
    return findWithinSortedOrdmag(sorted, baselineIdx);
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
    vector<uint32> sorted = nonzeroCounts;
    sort(sorted.begin(), sorted.end());
    
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
        uint32 filteredCells = findWithinSortedOrdmag(sorted, baselineIdx);
        
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
    return runCRSimpleFilterBootstrap(nUMIperCB, nCB, params,
                                     vector<uint32>(), vector<string>());
}

OrdMagResult SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
    const vector<uint32>& nUMIperCB,
    uint32 nCB,
    OrdMagParams& params,
    const vector<uint32>& detectedGenes,
    const vector<string>& barcodeIds,
    const vector<uint64_t>& nonMitoUMIs,
    OrdMagBootstrapTrace* trace
) {
    OrdMagResult result = {};
    if (trace) *trace = OrdMagBootstrapTrace();

    if (nUMIperCB.size() != nCB ||
        (!nonMitoUMIs.empty() && nonMitoUMIs.size() != nCB) ||
        (!detectedGenes.empty() && detectedGenes.size() != nCB) ||
        (!barcodeIds.empty() && barcodeIds.size() != nCB)) {
        throw std::invalid_argument("OrdMag tie metadata must match the barcode count");
    }
    
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

    // Make bootstrap samples depend on the count distribution, not the input
    // barcode order. The external caller already supplies descending counts.
    sort(nonzeroCounts.begin(), nonzeroCounts.end(), std::greater<uint32>());
    
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
        if (trace) trace->bootstrapThreads = nThreads;
        
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
        const uint32 workers = params.maxConcurrentThreads
            ? std::min(params.maxConcurrentThreads, nThreads) : nThreads;
        scrna::parallelFor(nThreads, workers, [&](size_t t, size_t) {
            uint32 startIdx = t * chunkSize;
            uint32 endIdx = min(startIdx + chunkSize, nBoot);
            if (startIdx < endIdx) {
                bootstrapWorker(startIdx, endIdx, baseSeed + t * 1000);
            }
        });
        
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
        if (trace) trace->bootstrapThreads = nThreads;
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
        const uint32 workers = params.maxConcurrentThreads
            ? std::min(params.maxConcurrentThreads, nThreads) : nThreads;
        scrna::parallelFor(nThreads, workers, [&](size_t t, size_t) {
            uint32 startIdx = t * chunkSize;
            uint32 endIdx = min(startIdx + chunkSize, nBoot);
            if (startIdx < endIdx) {
                topNWorker(startIdx, endIdx, baseSeed2 + t * 1000);
            }
        });
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
    if (trace) {
        trace->recoveredCells = recoveredCells;
        trace->meanRetained = meanTopN;
        trace->sdRetained = sdTopN;
    }
    
    // Round to get number of cells
    uint32 nCellsSimple = ordMagRetainCount(nNonzero, meanTopN);
    
    cout << "Bootstrap mean = " << meanTopN << ", sd = " << sdTopN << ", nCellsSimple = " << nCellsSimple << endl;
    
    // Step 5: Get actual filtered barcodes (top nCellsSimple by UMI)
    // Create index array sorted by UMI count (descending)
    typedef struct {uint32 index; uint32 count;} IndCount;
    vector<IndCount> indCount(nCB);
    for (uint32 ii = 0; ii < nCB; ii++) {
        indCount[ii].index = ii;
        indCount[ii].count = nUMIperCB[ii];
    }
    
    sort(indCount.begin(), indCount.end(), [&](const IndCount& ic1, const IndCount& ic2) {
        return ordMagRankBefore(ic1.index, ic2.index, nUMIperCB,
                                detectedGenes.empty() ? nullptr : &detectedGenes,
                                barcodeIds.empty() ? nullptr : &barcodeIds,
                                nonMitoUMIs.empty() ? nullptr : &nonMitoUMIs);
    });
    
    // The rank target does not override the configured primary floor. Apply it
    // before constructing the primary prefix, median and tail candidates:
    // primary cells receive automatic p=0 in EmptyDrops. Previously only
    // the legacy Flex wrapper removed these low-UMI passers afterward.
    const uint32 primaryFloor = ordMagPrimaryFloor(params);
    const uint32 estimatedRetainCount = nCellsSimple;
    while (nCellsSimple > 0 && indCount[nCellsSimple - 1].count < primaryFloor)
        --nCellsSimple;
    if (nCellsSimple != estimatedRetainCount)
        cout << "[OrdMag floor] primary_umi_min=" << primaryFloor
             << " removed=" << estimatedRetainCount - nCellsSimple
             << " nCellsSimple=" << nCellsSimple << endl;

    // Within the floor, keep the exact target and deterministic quality/barcode
    // ordering without expanding or dropping an entire boundary count group.
    
    // Compute retain threshold (UMI of last passing cell)
    uint32 retainThreshold = (nCellsSimple > 0 && nCellsSimple <= nCB) 
        ? indCount[nCellsSimple - 1].count : 0;

    if (nCellsSimple > 0) {
        const uint32 cutoff = indCount[nCellsSimple - 1].index;
        uint32 umiTie = 0, qualityTie = 0, selectedTie = 0;
        for (uint32 rank = 0; rank < nCB; ++rank) {
            const uint32 idx = indCount[rank].index;
            if (nUMIperCB[idx] != retainThreshold) continue;
            ++umiTie;
            if (rank < nCellsSimple) ++selectedTie;
            if ((detectedGenes.empty() || detectedGenes[idx] == detectedGenes[cutoff]) &&
                (nonMitoUMIs.empty() || nonMitoUMIs[idx] == nonMitoUMIs[cutoff])) ++qualityTie;
        }
        cout << "[OrdMag rank] total_UMI=" << retainThreshold
             << " detected_genes=" << (detectedGenes.empty() ? 0 : detectedGenes[cutoff])
             << " non_MT_UMI=" << (nonMitoUMIs.empty() ? retainThreshold : nonMitoUMIs[cutoff])
             << " umi_tie=" << umiTie << " quality_tie=" << qualityTie
             << " selected_umi_tie=" << selectedTie << endl;
    }
    
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
        if (meetsEmptyDropsCandidateFloor(nUMIperCB[ii], minUMI)) {
            candLimit++;
        }
    }
    candLimit = max(candLimit, nCellsSimple);
    uint32 iCandLast = min(nCellsSimple + params.candMaxN, min(candLimit, nCB));
    
    // Extract candidate indices
    for (uint32 ii = 0; ii < iCandLast; ii++) {
        result.candidateIndices.push_back(indCount[ii].index);
    }
    
    // Ambient ranks use UMIs and barcode identity only: tie quality must not
    // preferentially put low-quality cells into the ambient profile.
    sort(indCount.begin(), indCount.end(), [&](const IndCount& a, const IndCount& b) {
        return ordMagRankBefore(a.index, b.index, nUMIperCB, nullptr,
                                barcodeIds.empty() ? nullptr : &barcodeIds);
    });
    // Extract ambient indices (same window logic as non-bootstrap version)
    uint32 scaledIndMin = min(params.indMin, nCB);
    // STAR EmptyDrops_CR treats indMax as exclusive
    uint32 scaledIndMax = min(params.indMax, nCB);
    uint32 minAmbientCells = 0;
    if (nCB >= 1000) {
        // Guarded minimum: at least 2% or 5000 cells, but never more than 10%.
        uint32 minFrac = nCB / 50;  // 2%
        uint32 minAbs = 5000;
        uint32 minCandidate = max(minAbs, minFrac);
        minAmbientCells = min(nCB / 10, minCandidate);  // cap at 10%
    } else {
        minAmbientCells = min((uint32)100, nCB);
    }
    uint32 ambientWindowSize = (scaledIndMax > scaledIndMin) ? (scaledIndMax - scaledIndMin) : 0;
    
    if (scaledIndMax <= scaledIndMin || ambientWindowSize < minAmbientCells) {
        uint32 fallbackSize = min(minAmbientCells, nCB);
        uint32 fallbackStart = (nCB >= fallbackSize) ? (nCB - fallbackSize) : 0;
        result.ambientIndices.clear();
        for (uint32 ii = fallbackStart; ii < nCB; ii++) {
            result.ambientIndices.push_back(indCount[ii].index);
        }
        result.ambientRange = make_pair(fallbackStart + 1, nCB);
    } else {
        uint32 ambientStart = scaledIndMin;
        uint32 ambientEnd = scaledIndMax;
        result.ambientIndices.clear();
        for (uint32 ii = ambientStart; ii < ambientEnd; ii++) {
            if (ii < nCB) {
                result.ambientIndices.push_back(indCount[ii].index);
            }
        }
        result.ambientRange = make_pair(ambientStart, ambientEnd);
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
    OrdMagResult result = {};
    
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
    int maxInd = (int)llround(nExpectedCells * (1.0 - params.maxPercentile));
    if (maxInd < 0) maxInd = 0;
    if ((uint32)maxInd >= nCB) maxInd = (int)nCB - 1;
    
    // Compute retain threshold (match STAR EmptyDrops_CR indexing)
    uint32 nUMImax = totalsSorted[(uint32)maxInd];
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
    
    // Apply the same primary floor as the bootstrap path, after the fallback
    // so that it cannot reintroduce low-count cells (including zero counts).
    retain = max(retain, max(ordMagPrimaryFloor(params), (uint32)1));
    while (ncellsSimple > 0 && totalsSorted[ncellsSimple - 1] < retain)
        --ncellsSimple;

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
        if (meetsEmptyDropsCandidateFloor(nUMIperCB[ii], minUMI)) {
            candLimit++;
        }
    }
    candLimit = max(candLimit, ncellsSimple);
    uint32 iCandLast = min(ncellsSimple + params.candMaxN, min(candLimit, nCB));
    
    // Extract candidate indices (simple cells + tail) - matches Python behavior
    for (uint32 ii = 0; ii < iCandLast; ii++) {
        result.candidateIndices.push_back(indCount[ii].index);
    }
    
    // Extract passing indices (simple filter)
    for (uint32 ii = 0; ii < ncellsSimple; ii++) {
        result.passingIndices.push_back(indCount[ii].index);
    }
    
    // Extract ambient indices
    uint32 scaledIndMin = min(params.indMin, nCB);
    uint32 scaledIndMax = min(params.indMax, nCB);
    uint32 minAmbientCells = 0;
    if (nCB >= 1000) {
        // Guarded minimum: at least 2% or 5000 cells, but never more than 10%.
        uint32 minFrac = nCB / 50;  // 2%
        uint32 minAbs = 5000;
        uint32 minCandidate = max(minAbs, minFrac);
        minAmbientCells = min(nCB / 10, minCandidate);  // cap at 10%
    } else {
        minAmbientCells = min((uint32)100, nCB);
    }
    uint32 ambientWindowSize = (scaledIndMax >= scaledIndMin) ? (scaledIndMax - scaledIndMin) : 0;
    
    if (scaledIndMax <= scaledIndMin || ambientWindowSize < minAmbientCells) {
        uint32 fallbackSize = min(minAmbientCells, nCB);
        uint32 fallbackStart = (nCB >= fallbackSize) ? (nCB - fallbackSize) : 0;
        result.ambientIndices.clear();
        for (uint32 ii = fallbackStart; ii < nCB; ii++) {
            result.ambientIndices.push_back(indCount[ii].index);
        }
        result.ambientRange = make_pair(fallbackStart + 1, nCB);
    } else {
        uint32 ambientStart = scaledIndMin;
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
    summaryOut << "    \"primary_umi_min\": " << ordMagPrimaryFloor(params) << ",\n";
    summaryOut << "    \"umi_min_frac_median\": " << params.umiMinFracMedian << ",\n";
    summaryOut << "    \"cand_max_n\": " << params.candMaxN << ",\n";
    summaryOut << "    \"ind_min\": " << params.indMin << ",\n";
    summaryOut << "    \"ind_max\": " << params.indMax << "\n";
    summaryOut << "  }\n";
    summaryOut << "}\n";
    summaryOut.close();
}
