#include "OrdMagStage.h"
#include "ParallelTasks.h"
#include <atomic>
#include <chrono>
#include <iostream>
#include <stdexcept>

static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
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
