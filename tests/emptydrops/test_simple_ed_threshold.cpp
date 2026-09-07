/**
 * @file test_simple_ed_threshold.cpp
 * @brief Unit test: SimpleEmptyDropsStage thresholding
 * 
 * Target: SimpleEmptyDropsStage::runCRSimpleFilter
 * Test: Use tiny UMI list; verify retain threshold + simple passers match expected
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>

#include "OrdMagStage.h"
#include "AdaptiveAmbientWindow.h"
#include "EmptyDropsMultinomial.h"
#include "FlexTagGroup.h"

using namespace std;

int main() {
    cout << "Test: Simple EmptyDrops thresholding (SimpleEmptyDropsStage::runCRSimpleFilter)" << endl;
    
    // Create UMI counts simulating a typical experiment
    uint32_t nCells = 10000;
    vector<uint32_t> nUMIperCB(nCells);
    
    // Use fixed seed for reproducibility
    srand(42);
    
    // High-confidence cells (top 1000): 10k-50k UMI
    for (uint32_t i = 0; i < 1000; i++) {
        nUMIperCB[i] = 10000 + (rand() % 40000);
    }
    // Medium cells (2000): 1k-5k UMI
    for (uint32_t i = 1000; i < 3000; i++) {
        nUMIperCB[i] = 1000 + (rand() % 4000);
    }
    // Low cells (2000): 100-500 UMI
    for (uint32_t i = 3000; i < 5000; i++) {
        nUMIperCB[i] = 100 + (rand() % 400);
    }
    // Empty droplets (5000): 0-100 UMI
    for (uint32_t i = 5000; i < nCells; i++) {
        nUMIperCB[i] = rand() % 100;
    }
    
    // Set up parameters (matching CR defaults)
    SimpleEmptyDropsParams params;
    params.nExpectedCells = 3000;
    params.maxPercentile = 0.99;
    params.maxMinRatio = 10.0;
    params.umiMin = 500;
    params.umiMinFracMedian = 0.01;
    params.candMaxN = 20000;
    params.indMin = 45000;
    params.indMax = 90000;
    params.useBootstrap = false;
    
    // Run filter
    SimpleEmptyDropsResult result = SimpleEmptyDropsStage::runCRSimpleFilter(
        nUMIperCB, nCells, params);
    
    // Verify basic results
    if (result.nCellsSimple == 0) {
        cerr << "FAIL: No cells passed simple filter" << endl;
        return 1;
    }
    
    if (result.retainThreshold == 0) {
        cerr << "FAIL: Retain threshold is 0" << endl;
        return 1;
    }
    
    // Should detect at least some high-confidence cells (expect >= 500)
    if (result.nCellsSimple < 500) {
        cerr << "FAIL: Too few cells passed (" << result.nCellsSimple << ", expected >= 500)" << endl;
        return 1;
    }
    
    // Verify candidate indices are populated
    if (result.candidateIndices.empty()) {
        cerr << "FAIL: No candidate indices" << endl;
        return 1;
    }
    
    // Verify passing indices are subset of candidates
    if (result.passingIndices.size() > result.candidateIndices.size()) {
        cerr << "FAIL: More passing indices (" << result.passingIndices.size() 
             << ") than candidates (" << result.candidateIndices.size() << ")" << endl;
        return 1;
    }
    
    cout << "PASS: Simple filter: " << result.nCellsSimple << " cells, "
         << "threshold=" << result.retainThreshold << ", "
         << "candidates=" << result.candidateIndices.size() << endl;

    if (!meetsEmptyDropsCandidateFloor(500, 500) ||
        meetsEmptyDropsCandidateFloor(499, 500)) {
        cerr << "FAIL: the EmptyDrops candidate floor is not inclusive" << endl;
        return 1;
    }

    vector<pair<uint32_t, uint32_t>> ranked = {
        {100, 0}, {90, 1}, {40, 2}, {30, 3}, {20, 4}, {10, 5}
    };
    AdaptiveAmbientWindow fixed = selectAdaptiveAmbientWindow(ranked, 2, 4, 60);
    if (fixed.start != 2 || fixed.end != 4 || fixed.umiMass != 70) {
        cerr << "FAIL: established ambient window was not preserved" << endl;
        return 1;
    }
    AdaptiveAmbientWindow extended = selectAdaptiveAmbientWindow(ranked, 2, 4, 95);
    if (extended.start != 2 || extended.end != 6 || extended.umiMass != 100) {
        cerr << "FAIL: ambient endpoint did not extend to the requested mass" << endl;
        return 1;
    }
    AdaptiveAmbientWindow exhausted = selectAdaptiveAmbientWindow(ranked, 2, 4, 1000);
    if (exhausted.end != ranked.size() || exhausted.umiMass != 100) {
        cerr << "FAIL: ambient endpoint did not stop at the available droplets" << endl;
        return 1;
    }

    cout << "PASS: inclusive UMI floor and adaptive ambient window" << endl;

    const vector<string> fusedTags = {"ACTTTAGG", "AACGGGAA"};
    if (!barcodeHasAnyFlexTag("AAAAAAAAAAAAAAAAACTTTAGG", fusedTags) ||
        !barcodeHasAnyFlexTag("CCCCCCCCCCCCCCCCAACGGGAA-1", fusedTags) ||
        barcodeHasAnyFlexTag("GGGGGGGGGGGGGGGGAGTAGGCT", fusedTags)) {
        cerr << "FAIL: fused TAG8 membership did not preserve composite cell identity" << endl;
        return 1;
    }
    cout << "PASS: one- and two-tag composite barcode selection" << endl;
    return 0;
}
