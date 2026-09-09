#include "OrdMagRank.h"
#include "OrdMagStage.h"
#include "scrna_api.h"
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>

static void require(bool ok, const char* message) {
    if (!ok) { std::cerr << "FAIL: " << message << '\n'; std::exit(1); }
}

static std::set<string> identities(const vector<uint32>& indices, const vector<string>& barcodes) {
    std::set<string> out;
    for (uint32 idx : indices) out.insert(barcodes[idx]);
    return out;
}

int main() {
    // Duplicate coordinates and explicit zeros do not inflate detected genes.
    // MT removal affects only scores; duplicate positive counts still sum.
    const vector<uint32> genes = {0, 1, 1, 2, 3};
    const vector<uint32> counts = {0, 200, 300, 500, 1000};
    vector<uint32> seen(4, UINT32_MAX);
    const uint8_t mt[] = {0, 0, 1, 0};
    const auto all = ordMagCellQuality(genes.data(), counts.data(), genes.size(), seen, 0);
    const auto masked = ordMagCellQuality(genes.data(), counts.data(), genes.size(), seen, 1, mt);
    require(all.detectedGenes == 3 && all.nonMitoUMIs == 2000, "positive distinct genes");
    require(masked.detectedGenes == 2 && masked.nonMitoUMIs == 1500, "MT-excluded scores");
    const uint32 oneGene = 2, oneCount = 1;
    const auto oneMT = ordMagCellQuality(&oneGene, &oneCount, 1, seen, 2, mt);
    require(oneMT.detectedGenes == 0 && oneMT.nonMitoUMIs == 0, "one MT UMI scores zero");

    // Higher total UMI remains primary; quality ranks equal totals. Preserve
    // TAG8 even when the first 16 barcode bases match.
    const vector<uint32> totals = {2000, 2000, 2000, 2000, 2001};
    const vector<uint32> diversity = {300, 100, 200, 200, 1};
    const vector<uint64_t> nonMT = {1500, 1900, 1900, 1900, 1};
    const vector<string> bcs = {"AAAAAAAAAAAAAAAAACAGTCTG", "AAAAAAAAAAAAAAAAAGTGAGTG",
        "CCCCCCCCCCCCCCCCAGTGAGTG", "CCCCCCCCCCCCCCCCACAGTCTG", "GGGGGGGGGGGGGGGGACAGTCTG"};
    vector<uint32> order = {0, 1, 2, 3, 4};
    std::sort(order.begin(), order.end(), [&](uint32 a, uint32 b) {
        return ordMagRankBefore(a, b, totals, &diversity, &bcs, &nonMT);
    });
    require(order == vector<uint32>({4, 3, 2, 1, 0}), "UMI, non-MT, genes, full barcode priority");
    std::sort(order.begin(), order.end(), [&](uint32 a, uint32 b) {
        return ordMagRankBefore(a, b, totals, &diversity, &bcs);
    });
    require(order == vector<uint32>({4, 0, 3, 2, 1}), "gene-only ranking");

    // User's midpoint counterexample, with realistic UMI values: 21,558 tied
    // cells at 2,000 UMIs must not jump together when the target changes by 0.2.
    require(ordMagRetainCount(34862, 24082.9) == 24083 &&
            ordMagRetainCount(34862, 24083.1) == 24083, "no midpoint group jump");
    uint32 previous = 0;
    for (uint32 tenth = 0; tenth <= 1200; ++tenth) {
        const uint32 now = ordMagRetainCount(100, tenth / 10.0);
        require(now >= previous && now <= previous + 1, "target moves at most one rank");
        previous = now;
    }
    require(ordMagRetainCount(0, 10) == 0 && ordMagRetainCount(100, -1) == 0, "target clamping");

    // Exercise the actual bootstrap caller at a 2,000-UMI boundary. Quality
    // metadata must not change its estimated count or ambient membership.
    vector<uint32> umi(100, 2000), ngenes(100);
    vector<uint64_t> nonmito(100);
    vector<string> barcodes(100);
    for (uint32 i = 0; i < 100; ++i) {
        umi[i] = i < 2 ? 100000 : i < 12 ? 10000 : 2000;
        ngenes[i] = 100 + (i * 17) % 41;
        nonmito[i] = umi[i] - (i % 5) * 100;
        barcodes[i] = "barcode" + std::to_string(1000 + i);
    }
    OrdMagParams params;
    params.nExpectedCells = 100;
    params.nBootstrapSamples = 40;
    params.maxThreads = 2;
    params.maxPercentile = .99;
    params.maxMinRatio = 10;
    params.umiMin = 500;
    params.umiMinFracMedian = .01;
    params.candMaxN = 100;
    params.indMin = 80;
    params.indMax = 100;
    auto a = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(umi, 100, params, ngenes, barcodes, nonmito);
    require(a.nCellsSimple > 12 && a.nCellsSimple < 100 && a.retainThreshold == 2000,
            "fixture actually splits realistic UMI tie");
    require(std::equal(a.passingIndices.begin(), a.passingIndices.end(), a.candidateIndices.begin()),
            "simple calls remain candidate prefix for EmptyDrops");
    const auto firstCalls = identities(a.passingIndices, barcodes);
    const auto firstAmbient = identities(a.ambientIndices, barcodes);
    std::reverse(umi.begin(), umi.end()); std::reverse(ngenes.begin(), ngenes.end());
    std::reverse(nonmito.begin(), nonmito.end()); std::reverse(barcodes.begin(), barcodes.end());
    auto b = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(umi, 100, params, ngenes, barcodes, nonmito);
    require(firstCalls == identities(b.passingIndices, barcodes), "permutation leaves calls unchanged");
    require(firstAmbient == identities(b.ambientIndices, barcodes), "permutation leaves ambient unchanged");
    auto c = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(umi, 100, params, ngenes, barcodes);
    require(c.nCellsSimple == a.nCellsSimple, "MT scores do not alter count estimation");
    require(firstAmbient == identities(c.ambientIndices, barcodes), "MT scores do not alter ambient ranks");
    bool invalid = false;
    try {
        SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(umi, 100, params, vector<uint32>(2), barcodes);
    } catch (const std::invalid_argument&) { invalid = true; }
    require(invalid, "reject misaligned metadata");

    // The new C ABI rejects an MT mask when counts needed to score it are absent.
    scrna_matrix_input input = {};
    scrna_ed_result result = {};
    scrna_ed_config* config = scrna_ed_config_create();
    require(scrna_emptydrops_run_with_rank_mask(&input, config, mt, &result) != 0 && result.error_message,
            "C API requires matrix evidence for MT scores");
    scrna_ed_result_free(&result);
    scrna_ed_config_destroy(config);
    std::cout << "PASS: quality ranking, exact targets, realistic ties, permutation, ambient isolation, C API guard\n";
}
