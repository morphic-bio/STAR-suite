/**
 * Standalone test: verify bootstrap OrdMag estimates ~7000 recovered cells
 * on the UCSF full-sample raw matrix.
 *
 * Build:
 *   g++ -O2 -std=c++17 -I../core/features/libscrna/include \
 *       tests/test_ordmag_bootstrap.cpp \
 *       core/features/libscrna/src/OrdMagStage.cpp \
 *       -lpthread -o tests/test_ordmag_bootstrap
 *
 * Run:
 *   tests/test_ordmag_bootstrap <path-to-raw-matrix-dir>
 */

#include "OrdMagStage.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdlib>

using namespace std;

static vector<uint32_t> loadUMIcounts(const string& matrixDir) {
    string mtxPath = matrixDir + "/matrix.mtx";
    ifstream mtx(mtxPath);
    if (!mtx.is_open()) {
        cerr << "ERROR: Cannot open " << mtxPath << endl;
        exit(1);
    }

    // Skip comments
    string line;
    while (getline(mtx, line)) {
        if (line[0] != '%') break;
    }

    // Header: nFeatures nBarcodes nnz
    uint32_t nFeatures, nBarcodes, nnz;
    {
        istringstream iss(line);
        iss >> nFeatures >> nBarcodes >> nnz;
    }
    cerr << "Matrix: " << nFeatures << " features x " << nBarcodes << " barcodes, " << nnz << " entries" << endl;

    vector<uint32_t> umiCounts(nBarcodes, 0);
    uint32_t feat, bc, count;
    for (uint32_t i = 0; i < nnz; i++) {
        mtx >> feat >> bc >> count;
        if (bc >= 1 && bc <= nBarcodes) {
            umiCounts[bc - 1] += count;
        }
    }

    uint32_t nonzero = 0;
    for (auto c : umiCounts) if (c > 0) nonzero++;
    cerr << "Loaded " << nBarcodes << " barcodes, " << nonzero << " nonzero" << endl;

    return umiCounts;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <raw-matrix-dir>" << endl;
        cerr << "  e.g. " << argv[0] << " Solo.out/GeneFull/raw/" << endl;
        return 1;
    }

    string matrixDir = argv[1];
    vector<uint32_t> umiCounts = loadUMIcounts(matrixDir);
    uint32_t nCB = umiCounts.size();

    // Sort descending for reporting
    vector<uint32_t> sorted = umiCounts;
    sort(sorted.begin(), sorted.end(), greater<uint32_t>());

    cerr << "\nTop UMI counts: ";
    for (int i = 0; i < min(10, (int)nCB); i++) cerr << sorted[i] << " ";
    cerr << endl;

    // === Test 1: Old-style fixed nExpectedCells=3000 ===
    {
        cerr << "\n=== Test 1: Fixed nExpectedCells=3000 (current STAR behavior) ===" << endl;
        SimpleEmptyDropsParams params;
        params.nExpectedCells = 3000;
        params.maxPercentile = 0.99;
        params.maxMinRatio = 10.0;
        params.umiMin = 100;
        params.umiMinFracMedian = 0.01;
        params.candMaxN = 20000;
        params.indMin = 45000;
        params.indMax = 90000;

        auto t0 = chrono::steady_clock::now();
        SimpleEmptyDropsResult result = SimpleEmptyDropsStage::runCRSimpleFilter(umiCounts, nCB, params);
        auto t1 = chrono::steady_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();

        cerr << "  nCellsSimple = " << result.nCellsSimple << endl;
        cerr << "  retainThreshold = " << result.retainThreshold << endl;
        cerr << "  medianVal = " << result.medianVal << endl;
        cerr << "  minUMI = " << result.minUMI << endl;
        cerr << "  elapsed = " << ms << " ms" << endl;

        cout << "FIXED: nCellsSimple=" << result.nCellsSimple
             << " retainThreshold=" << result.retainThreshold << endl;
    }

    // === Test 2: Bootstrap estimation (CR9 style) ===
    {
        cerr << "\n=== Test 2: Bootstrap OrdMag (CR9 style) ===" << endl;
        SimpleEmptyDropsParams params;
        params.nExpectedCells = 0;  // trigger bootstrap estimation
        params.maxPercentile = 0.99;
        params.maxMinRatio = 10.0;
        params.umiMin = 100;
        params.umiMinFracMedian = 0.01;
        params.candMaxN = 20000;
        params.indMin = 45000;
        params.indMax = 90000;
        params.useBootstrap = true;
        params.nBootstrapSamples = 100;
        params.recoveredCellsQuantile = 0.99;
        params.minRecoveredCells = 50;
        params.maxExpectedCells = 0;  // auto: indMin/2

        auto t0 = chrono::steady_clock::now();
        SimpleEmptyDropsResult result = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(umiCounts, nCB, params);
        auto t1 = chrono::steady_clock::now();
        double ms = chrono::duration<double, milli>(t1 - t0).count();

        cerr << "  estimated nExpectedCells = " << params.nExpectedCells << endl;
        cerr << "  nCellsSimple = " << result.nCellsSimple << endl;
        cerr << "  retainThreshold = " << result.retainThreshold << endl;
        cerr << "  medianVal = " << result.medianVal << endl;
        cerr << "  minUMI = " << result.minUMI << endl;
        cerr << "  elapsed = " << ms << " ms" << endl;

        cout << "BOOTSTRAP: estimated_recovered=" << params.nExpectedCells
             << " nCellsSimple=" << result.nCellsSimple
             << " retainThreshold=" << result.retainThreshold << endl;
    }

    // === Test 3: Cross-check with CR expected values ===
    cerr << "\n=== Cross-check ===" << endl;
    cerr << "CR expected: ~7325 cells, cutoff ~78-80 UMIs" << endl;
    cerr << "UMI at rank 7325: " << (nCB > 7325 ? sorted[7324] : 0) << endl;
    cerr << "UMI at rank 7000: " << (nCB > 7000 ? sorted[6999] : 0) << endl;
    cerr << "UMI at rank 6601: " << (nCB > 6601 ? sorted[6600] : 0) << endl;

    return 0;
}
