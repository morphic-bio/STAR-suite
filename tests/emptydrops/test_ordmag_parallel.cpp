#include "OrdMagStage.h"
#include "ParallelTasks.h"
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

static void sameResult(const OrdMagResult& result, const OrdMagResult& expected) {
    require(result.passingIndices == expected.passingIndices &&
        result.candidateIndices == expected.candidateIndices &&
        result.ambientIndices == expected.ambientIndices &&
        result.retainThreshold == expected.retainThreshold &&
        result.nCellsSimple == expected.nCellsSimple && result.minUMI == expected.minUMI &&
        result.medianVal == expected.medianVal &&
        result.candidateLastRank == expected.candidateLastRank &&
        result.ambientRange == expected.ambientRange, "worker-independent complete OrdMag result");
}

static void workerIndependentResamples(const vector<uint32>& umi,
                                      const vector<uint32>& genes,
                                      const vector<string>& barcodes) {
    const char* original = std::getenv("OMP_NUM_THREADS");
    const bool wasSet = original != nullptr;
    const std::string saved = original ? original : "";
    for (uint32 samples : {17u, 97u}) {
        for (uint32 expectedCells : {0u, 113u}) {
            OrdMagParams settings;
            settings.nBootstrapSamples = samples;
            settings.bootstrapSeed = samples == 17 ? 19 : 0;
            settings.nExpectedCells = expectedCells;
            settings.maxExpectedCells = 2000;
            settings.umiMin = 1;
            settings.umiMinFracMedian = .01;
            settings.candMaxN = 3000;
            settings.indMin = 3000;
            settings.indMax = 5000;
            settings.maxThreads = 1;
            ::unsetenv("OMP_NUM_THREADS");
            OrdMagParams params = settings;
            OrdMagBootstrapTrace reference;
            const auto expected = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
                umi, umi.size(), params, genes, barcodes, {}, &reference);
            for (const char* omp : {static_cast<const char*>(nullptr), "2", "32"}) {
                if (omp) ::setenv("OMP_NUM_THREADS", omp, 1);
                else ::unsetenv("OMP_NUM_THREADS");
                for (uint32 workers : {1u, 2u, 4u, 16u, 32u, 0u}) {
                    params = settings;
                    params.maxThreads = workers;
                    params.maxConcurrentThreads = 0;
                    OrdMagBootstrapTrace trace;
                    const auto result = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
                        umi, umi.size(), params, genes, barcodes, {}, &trace);
                    require(trace.recoveredCells == reference.recoveredCells &&
                        trace.meanRetained == reference.meanRetained &&
                        trace.sdRetained == reference.sdRetained,
                        "resample estimates ignore workers, hardware and OMP_NUM_THREADS");
                    sameResult(result, expected);
                }
            }
        }
    }
    if (wasSet) ::setenv("OMP_NUM_THREADS", saved.c_str(), 1);
    else ::unsetenv("OMP_NUM_THREADS");
}

int main() {
    // Uneven partitions and custom seeds must retain exactly the same samples.
    vector<uint32> umi(5000), genes(5000);
    vector<string> barcodes(5000);
    for (uint32 i = 0; i < umi.size(); ++i) {
        umi[i] = i < 7 ? 80000 + i * 731 : i < 53 ? 15000 + i * 53 :
            i < 280 ? 1000 + i % 17 : i < 490 ? 500 : 1 + i % 97;
        genes[i] = 1 + i % 137;
        barcodes[i] = "BC" + std::to_string(10000 + i);
    }
    for (uint32 samples : {17u, 97u}) {
        OrdMagParams params;
        params.maxThreads = 48; params.nBootstrapSamples = samples;
        params.bootstrapSeed = samples == 17 ? 19 : 0;
        params.maxExpectedCells = 2000;
        params.maxPercentile = .99; params.maxMinRatio = 10;
        params.umiMin = 500; params.umiMinFracMedian = .01;
        params.candMaxN = 3000; params.indMin = 3000; params.indMax = 5000;
        OrdMagBootstrapTrace reference;
        const auto expected = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
            umi, umi.size(), params, genes, barcodes, {}, &reference);
        for (uint32 workers : {1u, 6u, 13u}) {
            params.nExpectedCells = 0; // Re-run estimation, not only the top-N bootstrap.
            params.maxConcurrentThreads = workers;
            OrdMagBootstrapTrace trace;
            const auto result = SimpleEmptyDropsStage::runCRSimpleFilterBootstrap(
                umi, umi.size(), params, genes, barcodes, {}, &trace);
            require(trace.recoveredCells == reference.recoveredCells &&
                trace.meanRetained == reference.meanRetained && trace.sdRetained == reference.sdRetained &&
                trace.bootstrapThreads == reference.bootstrapThreads, "unchanged bootstrap streams and estimates");
            require(result.passingIndices == expected.passingIndices &&
                result.candidateIndices == expected.candidateIndices && result.ambientIndices == expected.ambientIndices &&
                result.retainThreshold == expected.retainThreshold, "unchanged ranks and stage membership");
        }
    }
    workerIndependentResamples(umi, genes, barcodes);
    std::atomic<unsigned> active(0), peak(0);
    vector<unsigned> visits(31, 0);
    scrna::parallelFor(visits.size(), 3, [&](size_t task, size_t worker) {
        require(worker < 3, "bounded worker IDs");
        const unsigned now = ++active;
        unsigned old = peak.load();
        while (old < now && !peak.compare_exchange_weak(old, now)) {}
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
        ++visits[task]; --active;
    });
    require(peak > 1 && peak <= 3 && active == 0, "parallel execution stays within budget");
    for (unsigned v : visits) require(v == 1, "each independent task runs exactly once");
    bool caught = false;
    try {
        scrna::parallelFor(31, 4, [&](size_t task, size_t) {
            if (task == 2) throw std::runtime_error("worker failure");
            ++active;
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
            --active;
        });
    } catch (const std::runtime_error&) { caught = true; }
    require(caught && active == 0, "join workers before propagating failure");
    std::cout << "PASS: bootstrap stream parity, bounded parallel tasks, joined exception propagation\n";
}
