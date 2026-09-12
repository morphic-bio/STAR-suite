#include "libflex/FlexFilter.h"
#include "libflex/FlexFilterIO.h"
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <cmath>

static bool sameNumber(double a, double b) { return a == b || (std::isnan(a) && std::isnan(b)); }

static void require(bool value, const char* message) { if (!value) throw std::runtime_error(message); }
static std::string cb(uint32_t index) {
    std::string value(16, 'A');
    for (size_t i = 0; i < 16; ++i) { value[15-i] = "ACGT"[index & 3]; index >>= 2; }
    return value;
}

int main(int argc, char** argv) {
    require(argc == 2, "pass a fresh artifact directory");
    const std::string root = argv[1];
    require(FlexFilterIO::defaultCreateDirectory(root, 0775), "create fixture root");
    FlexFilter::MemoryInputs input;
    input.sampleLabels = {"single", "paired", "paired", "empty", "empty"};
    input.sampleTags = {"ACTTTAGG", "ACAGTCTG", "AGTGAGTG", "AGAGGCAA", "ACTACTCA"};
    auto& m = input.matrixData;
    m.nGenes = 33; m.countMatStride = 3;
    for (uint32_t g=0; g<31; ++g) m.features.push_back("G"+std::to_string(g));
    m.features.push_back("MT0");
    m.features.push_back("DEPRECATED_MODEL_ONLY");
    for (size_t t = 0; t < input.sampleTags.size(); ++t) {
        for (uint32_t i = 0; i < 45250; ++i) {
            m.barcodes.push_back(cb(i) + input.sampleTags[t]); // same CB16 in different tags is intentional
            m.countCellGeneUMIindex.push_back(m.countCellGeneUMI.size());
            const uint32_t total = t >= 3 ? 1 + i % 3 : i < 2 ? 100000 : i < 12 ? 10000 : i < 100 ? 2000 :
                i < 150 ? 600 : i < 175 ? 500 : i < 180 ? 499 : 1;
            uint32_t n = 0;
            auto add = [&](uint32_t gene, uint32_t count) {
                if (count) { m.countCellGeneUMI.insert(m.countCellGeneUMI.end(), {gene, count, 0}); ++n; }
            };
            if (total >= 500) {
                const bool ambientLike = i >= 100 && i % 2 == 0;
                const uint32_t mt = ambientLike ? 0 : (i % 3) * 10;
                // Half of the 500/600-UMI tail follows the ambient gene; the
                // other half expresses a gene absent from ambient droplets.
                if (ambientLike) {
                    uint32_t used=0;
                    for (uint32_t g=0; g<20; ++g) { uint32_t c=total*(g+1)/250; add(g,c); used+=c; }
                    add(20,total-used);
                } else {
                    const uint32_t modelOnly = i >= 150 && i < 175 ? 100 : 0;
                    add(30, total-mt-modelOnly); add(31, mt); add(32, modelOnly);
                }
            } else {
                // Twenty-one ambient genes with many count frequencies are
                // required for the existing Simple Good-Turing estimator.
                uint32_t position=i%250, gene=0;
                while (gene<20 && position>=gene+1) { position-=gene+1; ++gene; }
                add(t >= 3 ? 32 : gene,total);
            }
            m.nUMIperCB.push_back(total); m.nGenePerCB.push_back(n);
        }
    }
    m.nCells = m.barcodes.size(); m.countCellGeneUMIindex.push_back(m.countCellGeneUMI.size());
    std::ofstream features(root + "/features.tsv"), barcodes(root + "/barcodes.tsv"), matrix(root + "/matrix.mtx"), mt(root + "/mt.tsv");
    for (const auto& g : m.features) features << g << '\t' << g << "\tGene Expression\n";
    for (const auto& b : m.barcodes) barcodes << b << '\n';
    mt << "MT0\n";
    matrix << "%%MatrixMarket matrix coordinate integer general\n" << m.nGenes << " " << m.nCells << ' ' << m.countCellGeneUMI.size()/3 << '\n';
    for (uint32_t c = 0; c < m.nCells; ++c) for (uint32_t j = 0; j < m.nGenePerCB[c]; ++j) {
        const size_t pos = m.countCellGeneUMIindex[c] + j*3;
        matrix << m.countCellGeneUMI[pos]+1 << ' ' << c+1 << ' ' << m.countCellGeneUMI[pos+1] << '\n';
    }
    features.close();barcodes.close();matrix.close();mt.close();
    FlexFilter::Config config;
    config.tagAwareCaller = true;
    config.useThreadPermits = false;
    config.totalThreads = 1;
    config.emptydropsParams.simN = 128;
    config.emptydropsParams.FDR = .01;
    config.emptydropsParams.mcThreads = 2;
    config.simpleEmptyDropsParams.maxThreads = 2;
    config.mitochondrialGenesPath = root + "/mt.tsv";
    config.debugOutputDir = root + "/internal";
    config.enableInvariantChecks = true;
    FlexFilter filter; FlexFilter::Outputs result;
    require(filter.runFromMemory(input, &result, config) == 0, "internal grouped caller succeeds");
    require(result.tagResults.size() == 3, "duplicate labels produce one sample result");
    for (const auto& group : result.tagResults) {
        std::ofstream calls(root + "/internal/" + group.sampleLabel + "/calls.txt");
        std::set<std::string> seen;
        for (const auto& b : group.passingBarcodes) { require(b.size() == 24 && seen.insert(b).second, "full unique composite identities"); calls << b << '\n'; }
        require(group.occupancyRemoved == 0, "valid fixture calls survive observed-tag occupancy");
        if (group.sampleLabel == "empty") {
            require(group.nSimpleCells == 0 && group.nTailTested == 0 &&
                    group.edPasserBarcodes.empty() && group.passingBarcodes.empty(),
                    "unused low-UMI tags must never contribute cells to occupancy");
            continue;
        }
        require(group.nSimpleCells > 0, "OrdMag is primary");
        require(group.nTailTested > 0, "fixture exercises EmptyDrops tail");
    }
    // Vary group scheduling and MC workers while preserving bootstrap streams.
    const auto serial = result;
    for (uint32_t budget : {7u, 2u}) {
        if (budget == 2) {
            // Exercise a second physical layout against the same serial
            // scientific result, including noncontiguous paired-tag rows.
            std::vector<uint32_t> pairs;
            pairs.reserve(m.countCellGeneUMI.size() / 3 * 2);
            for (size_t i = 0; i < m.countCellGeneUMI.size(); i += 3) {
                pairs.push_back(m.countCellGeneUMI[i]);
                pairs.push_back(m.countCellGeneUMI[i + 1]);
            }
            for (auto& offset : m.countCellGeneUMIindex) offset = offset / 3 * 2;
            m.countCellGeneUMI.swap(pairs);
            m.countMatStride = 2;
        }
        config.totalThreads = budget;
        config.useThreadPermits = true;
        config.debugOutputDir = root + "/parallel_" + std::to_string(budget);
        require(filter.runFromMemory(input, &result, config) == 0, "parallel grouped caller succeeds");
        require(result.tagResults.size() == serial.tagResults.size(), "same ordered groups");
        for (size_t i = 0; i < result.tagResults.size(); ++i) {
            const auto& a = serial.tagResults[i]; const auto& b = result.tagResults[i];
            require(a.sampleLabel == b.sampleLabel && a.tag == b.tag && a.expectedCells == b.expectedCells &&
                a.nRetainWindow == b.nRetainWindow && a.nSimpleCells == b.nSimpleCells &&
                a.nTailTested == b.nTailTested && a.nTailPassers == b.nTailPassers &&
                a.occupancyRemoved == b.occupancyRemoved, "identical stage counts and occupancy");
            require(a.tagBarcodes == b.tagBarcodes && a.retainBarcodes == b.retainBarcodes &&
                a.passingBarcodes == b.passingBarcodes && a.filteredBarcodes == b.filteredBarcodes &&
                a.edPasserBarcodes == b.edPasserBarcodes, "identical ordered stage identities");
            require(a.emptydropsResults.size() == b.emptydropsResults.size(), "same candidate ledger size");
            for (size_t j = 0; j < a.emptydropsResults.size(); ++j) {
                const auto& x = a.emptydropsResults[j]; const auto& y = b.emptydropsResults[j];
                require(x.cellIndex == y.cellIndex && sameNumber(x.pValue, y.pValue) &&
                    sameNumber(x.pAdjusted, y.pAdjusted) && sameNumber(x.obsLogProb, y.obsLogProb) &&
                    x.passesRawP == y.passesRawP && x.passesFDR == y.passesFDR,
                    "identical candidate p-values, probabilities and decisions");
            }
        }
    }
    // A worker failure must clear partial results after joining all groups.
    const auto validGene = m.countCellGeneUMI[0];
    m.countCellGeneUMI[0] = m.nGenes;
    config.debugOutputDir.clear();
    require(filter.runFromMemory(input, &result, config) != 0 && result.tagResults.empty(),
        "parallel malformed-matrix error clears partial results");
    m.countCellGeneUMI[0] = validGene;
    // Duplicate tag assignments are rejected before starting a caller.
    input.sampleTags[2] = input.sampleTags[1];
    require(filter.runFromMemory(input, &result, config) != 0, "reject duplicate tags");
    std::cout << "PASS: internal sample grouping, composite identity, primary OrdMag, tail, empty-tag floor, duplicate-tag rejection\n";
}
