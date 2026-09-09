#include "OrdMagStage.h"
#include "SimpleEDCaller.h"
#include <cassert>
#include <iostream>

static SimpleEmptyDropsParams parameters(uint32 floor) {
    SimpleEmptyDropsParams p;
    p.nExpectedCells = 3000;
    p.maxPercentile = .99;
    p.maxMinRatio = 10;
    p.umiMin = floor;
    p.umiMinFracMedian = .01;
    p.candMaxN = 20000;
    p.indMin = 45000;
    p.indMax = 90000;
    p.nBootstrapSamples = 10;
    p.maxThreads = 1;
    return p;
}

static void check(const vector<uint32>& counts, uint32 floor,
                  uint32 expected, bool bootstrap, bool forceFallback = false) {
    auto p = parameters(floor);
    if (forceFallback) p.maxMinRatio = .01;
    const auto result = bootstrap
        ? OrdMagStage::runCRSimpleFilterBootstrap(counts, counts.size(), p)
        : OrdMagStage::runCRSimpleFilter(counts, counts.size(), p);
    assert(result.nCellsSimple == expected);
    assert(result.passingIndices.size() == expected);
    assert(result.candidateIndices.size() == expected);
    assert(result.candidateLastRank == expected);
    for (size_t i = 0; i < expected; ++i) {
        assert(result.passingIndices[i] == result.candidateIndices[i]);
        assert(counts[result.passingIndices[i]] >= floor);
        assert(counts[result.passingIndices[i]] > 0);
    }
}

static void checkSharedCaller(bool bootstrap) {
    // Each unused tag can contain nonzero matrix entries, yet no eligible cell.
    // Exercise the actual shared path used by Flex, including EmptyDrops.
    vector<string> barcodes;
    vector<uint32> umi(3000), genes, counts, starts, nGenes(3000, 1);
    for (uint32 i = 0; i < umi.size(); ++i) {
        barcodes.push_back("unused_" + std::to_string(i));
        umi[i] = 1 + i % 3;
        starts.push_back(i);
        genes.push_back(i % 32);
        counts.push_back(umi[i]);
    }
    scrna_ed_config* config = scrna_ed_config_create();
    assert(config);
    config->use_bootstrap = bootstrap;
    config->n_expected_cells = 3000;
    config->umi_min = 500;
    config->sim_n = 128;
    SimpleEDOptions options;
    options.bootstrapThreads = 1;
    options.invariantChecks = true;
    scrna_ed_result result = {};
    SimpleEDRunInfo info;
    const int rc = runSimpleEDWithAmbient(barcodes, umi, genes, counts, starts,
        nGenes, 32, config, options, vector<uint8_t>(), &result, &info);
    if (rc) std::cerr << (result.error_message ? result.error_message : "caller failed") << '\n';
    assert(rc == 0);
    assert(result.n_simple_cells == 0);
    assert(result.n_candidates == 0);
    assert(result.n_barcodes == 0);
    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
}

int main() {
    vector<uint32> emptyTag(3000);
    for (uint32 i = 0; i < emptyTag.size(); ++i) emptyTag[i] = 1 + i % 3;
    vector<uint32> boundary(3000, 1);
    boundary[0] = 600;
    boundary[1] = 500;
    boundary[2] = 499;
    for (bool bootstrap : {false, true}) {
        check(emptyTag, 500, 0, bootstrap);
        check(boundary, 500, 2, bootstrap);
        check(emptyTag, 1, 3000, bootstrap); // Explicit low-count policy remains supported.
        check(vector<uint32>(100, 0), 500, 0, bootstrap);
        check(vector<uint32>(), 500, 0, bootstrap);
        checkSharedCaller(bootstrap);
    }
    check(emptyTag, 500, 0, false, true); // A fallback must not defeat the floor.
    std::cout << "OrdMag primary UMI floor regression tests PASSED\n";
}
