/**
 * Unit test: the OrdMag primary floor is separate from the EmptyDrops
 * candidate floor. Synthetic counts, no data files.
 *
 * Build (after `make -C core/features/libscrna`):
 *   g++ -O2 -std=c++17 -Icore/features/libscrna/include \
 *       tests/test_ordmag_primary_floor.cpp \
 *       core/features/libscrna/lib/libscrna.a -lpthread \
 *       -o tests/test_ordmag_primary_floor
 *
 * Run:
 *   tests/test_ordmag_primary_floor
 */

#include "OrdMagStage.h"
#include <iostream>
#include <vector>

using namespace std;

static int failures = 0;

static void check(bool ok, const string& what) {
    cout << (ok ? "PASS " : "FAIL ") << what << endl;
    if (!ok) failures++;
}

static SimpleEmptyDropsParams baseParams(uint32 nExpected) {
    SimpleEmptyDropsParams p;
    p.nExpectedCells = nExpected;
    p.maxPercentile = 0.99;
    p.maxMinRatio = 10.0;
    p.umiMin = 500;
    p.umiMinFracMedian = 0.01;
    p.candMaxN = 20000;
    p.indMin = 45000;
    p.indMax = 90000;
    p.maxThreads = 1;
    return p;
}

// Low-depth library: 2,000 cells at 150-449 UMIs, 100,000 empties at 1-15.
static vector<uint32> lowDepth() {
    vector<uint32> v;
    for (uint32 i = 0; i < 2000; i++) v.push_back(150 + i % 300);
    for (uint32 i = 0; i < 100000; i++) v.push_back(1 + i % 15);
    return v;
}

// Deep library: 3,000 cells at 5,000-7,999, 500 barcodes at 200-499,
// 300 at 500-799, 100,000 empties at 1-15.
static vector<uint32> deep() {
    vector<uint32> v;
    for (uint32 i = 0; i < 3000; i++) v.push_back(5000 + i % 3000);
    for (uint32 i = 0; i < 500; i++) v.push_back(200 + i % 300);
    for (uint32 i = 0; i < 300; i++) v.push_back(500 + i);
    for (uint32 i = 0; i < 100000; i++) v.push_back(1 + i % 15);
    return v;
}

static uint32 tailBelow(const SimpleEmptyDropsResult& r, const vector<uint32>& umi, uint32 floor) {
    uint32 n = 0;
    for (size_t i = r.nCellsSimple; i < r.candidateIndices.size(); i++)
        if (umi[r.candidateIndices[i]] < floor) n++;
    return n;
}

int main() {
    const vector<uint32> low = lowDepth();
    const uint32 nLow = (uint32)low.size();

    // Fixed-knee path.
    {
        SimpleEmptyDropsParams p = baseParams(2000);
        SimpleEmptyDropsResult r = SimpleEmptyDropsStage::runCRSimpleFilter(low, nLow, p);
        check(r.nCellsSimple == 0, "fixed knee, primary floor unset: umiMin 500 trims all low-depth primaries (Flex)");

        p.primaryUmiMin = 1;
        r = SimpleEmptyDropsStage::runCRSimpleFilter(low, nLow, p);
        check(r.nCellsSimple == 2000, "fixed knee, primary floor 1: all 2,000 low-depth primaries kept");
        check(r.minUMI == 500, "fixed knee, primary floor 1: candidate floor stays 500");
    }

    // Bootstrap path.
    {
        SimpleEmptyDropsParams p = baseParams(0);
        p.useBootstrap = true;
        SimpleEmptyDropsResult r = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(low, nLow, p);
        check(r.nCellsSimple == 0, "bootstrap, primary floor unset: umiMin 500 trims all low-depth primaries (Flex)");

        p = baseParams(0);
        p.useBootstrap = true;
        p.primaryUmiMin = 1;
        r = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(low, nLow, p);
        check(r.nCellsSimple >= 1500, "bootstrap, primary floor 1: low-depth primaries kept (got " +
              to_string(r.nCellsSimple) + ")");
        check(r.minUMI == 500, "bootstrap, primary floor 1: candidate floor stays 500");
    }

    // Candidate floor on a deep library: no tail candidate below umiMin.
    {
        const vector<uint32> d = deep();
        const uint32 nDeep = (uint32)d.size();
        SimpleEmptyDropsParams p = baseParams(3000);
        p.primaryUmiMin = 1;
        SimpleEmptyDropsResult r500 = SimpleEmptyDropsStage::runCRSimpleFilter(d, nDeep, p);
        check(tailBelow(r500, d, 500) == 0 && r500.candidateIndices.size() > r500.nCellsSimple,
              "deep, umiMin 500: tail candidates present and none below 500");

        p.umiMin = 100;
        SimpleEmptyDropsResult r100 = SimpleEmptyDropsStage::runCRSimpleFilter(d, nDeep, p);
        check(r100.nCellsSimple == r500.nCellsSimple, "deep: primaries do not depend on umiMin when the primary floor is set");
        check(tailBelow(r100, d, 500) == 500, "deep, umiMin 100: the 500 barcodes at 200-499 become candidates");
    }

    cout << (failures ? "FAILED" : "ALL PASS") << " (" << failures << " failures)" << endl;
    return failures ? 1 : 0;
}
