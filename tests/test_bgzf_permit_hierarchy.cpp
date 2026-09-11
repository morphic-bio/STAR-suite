#include "BgzfRateController.h"
#include "ThreadControl.h"
#include "SaturationPermitController.h"
#include "WorkloadDrainEstimate.h"
#include <atomic>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using D = ThreadControl::PermitDomain;
using W = ThreadControl::PermitWork;

static void exerciseOuterPolicy() {
    using C = star::multiome::SaturationPermitController;
    for (unsigned mask : {1U, 2U, 3U}) {
        C::Config c;
        c.configuredPermits = 4; c.activeMask = mask; c.startFromFloors = true;
        c.initialFloors = {{mask == 1 ? 4 : (mask == 3 ? 2 : 0), mask == 2 ? 4 : (mask == 3 ? 2 : 0), 0}};
        C controller(c);
        assert(controller.initialDecision().phase == C::Phase::STEADY);
        C::Observation o;
        o.mapEtaSec = 100; o.featureEtaSec = 10;
        o.mapUnitsDelta = mask & 1U ? 10 : 0;
        o.featureUnitsDelta = mask & 2U ? 10 : 0;
        const auto d = controller.observe(o);
        assert(d.atacFloor == 0 && d.mapFloor + d.featureFloor <= 4);
        if (mask == 3) assert(d.mapFloor == 3 && d.featureFloor == 1);
    }
    C::Config c;
    c.configuredPermits = 4; c.activeMask = 3; c.startFromFloors = true;
    c.initialFloors = {{2, 2, 0}};
    C controller(c);
    C::Observation o;
    o.mapUnitsDelta = 10; o.mapEtaSec = 100; o.featureEtaSec = 10;
    // A delayed first library or gap between libraries must retain FEATURE.
    for (int i = 0; i < 5; ++i) assert(controller.observe(o).featureFloor == 2);
    o.featureEstimateComplete = true;
    controller.observe(o);
    assert(controller.observe(o).featureFloor == 0);
    WorkloadDrainEstimate estimate;
    auto first = estimate.observe(1000, 1000, 1, false);
    assert(first.estimate > 1000 && first.eta > 0 && std::isfinite(first.eta));
    auto final = estimate.observe(1000, 1000, 1, true);
    assert(final.eta == 0);
    WorkloadDrainEstimate unknown;
    assert(!std::isfinite(unknown.observe(100, 0, 1, false).eta));
}

static void waitFor(const std::function<bool()>& predicate) {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(5);
    while (!predicate()) {
        assert(std::chrono::steady_clock::now() < deadline);
        std::this_thread::yield();
    }
}

// Exercise the actual GEX/FEATURE sampler without verbose telemetry. Both
// workloads run long enough for an ETA transfer, then publish durable completion.
static void exerciseOuterRuntime(const std::string& log) {
    ThreadControl pool;
    pool.mapPermitConfigure(true, 4, 4, false, false);
    pool.mapPermitConfigureDomainFloors({2, 2, 0});
    pool.mapPermitStartHierarchy(log, true, true, 9000, 3600);
    std::atomic<int> remaining[2];
    remaining[0].store(3); remaining[1].store(2);
    std::vector<std::thread> actors;
    for (int i = 0; i < 2; ++i) {
        for (int worker = 0; worker < (i == 0 ? 3 : 2); ++worker) {
            actors.emplace_back([&, i] {
                const auto domain = i == 0 ? D::MAP : D::FEATURE;
                for (int n = 0; n < (i == 0 ? 3000 : 1800); ++n) {
                    const auto token = pool.mapPermitAcquireForDomain(domain);
                    std::this_thread::sleep_for(std::chrono::microseconds(2000));
                    pool.mapPermitReleaseForDomain(domain, token, 1, 1, 2000000);
                }
                if (remaining[i].fetch_sub(1) == 1) pool.mapPermitMarkDomainComplete(domain);
            });
        }
    }
    for (auto& actor : actors) actor.join();
    pool.mapPermitStopHierarchy();
    const auto snap = pool.mapPermitSnapshot();
    assert(snap.availablePermits == 4 && snap.inUsePermits == 0 && snap.currentWaiters == 0);
    assert(snap.mapDomain.completedReadPairs == 9000 && snap.featureDomain.completedReadPairs == 3600);
    assert(snap.mapDomain.complete && snap.featureDomain.complete);
    std::ifstream source(log + ".outer.tsv");
    const std::string trace((std::istreambuf_iterator<char>(source)), std::istreambuf_iterator<char>());
    assert(trace.find("map-eta-late") != std::string::npos);
}

// Reproduce long mapping leases queued ahead of needed decoding. A low input
// buffer must let the decoder use its split on the next release. A full buffer
// must let processing borrow that split, preserving ordinary FIFO order.
static void exerciseDecoderDemand(bool inputLow, const std::string& log) {
    ThreadControl pool;
    pool.mapPermitConfigure(true, 4, 4, true, false);
    pool.mapPermitConfigureDomainFloors({4, 0, 0});
    pool.mapPermitStartHierarchy(log, false);
    int reader = 0;
    pool.mapPermitObserveDecode(D::MAP, &reader, inputLow ? 0 : 8, 8, 8, 4, inputLow, 1);
    uint64_t leases[4];
    for (auto& lease : leases) lease = pool.mapPermitAcquireForDomain(D::MAP);
    std::atomic<int> first{0};
    std::atomic<bool> finish{false};
    std::thread processor([&] {
        auto token = pool.mapPermitAcquireForDomain(D::MAP);
        int empty = 0; first.compare_exchange_strong(empty, 1);
        waitFor([&] { return finish.load(); });
        pool.mapPermitReleaseForDomain(D::MAP, token, 1, 1, 1);
    });
    waitFor([&] { return pool.mapPermitSnapshot().mapDomain.currentWaiters == 1; });
    std::thread decoder([&] {
        auto token = pool.mapPermitAcquireForDomain(D::MAP, W::BGZF);
        int empty = 0; first.compare_exchange_strong(empty, 2);
        waitFor([&] { return finish.load(); });
        pool.mapPermitReleaseForDomain(D::MAP, token, 1, 1, 1, W::BGZF);
    });
    waitFor([&] { return pool.mapPermitSnapshot().mapDomain.currentWaiters == 2; });
    pool.mapPermitReleaseForDomain(D::MAP, leases[0], 1, 1, 1);
    waitFor([&] { return first.load() != 0; });
    assert(first.load() == (inputLow ? 2 : 1));
    finish.store(true);
    for (int i = 1; i < 4; ++i) pool.mapPermitReleaseForDomain(D::MAP, leases[i], 1, 1, 1);
    processor.join(); decoder.join();
    pool.mapPermitObserveDecode(D::MAP, &reader, 0, 0, 0, 0, 0, 0);
    pool.mapPermitMarkDomainComplete(D::MAP);
    pool.mapPermitStopHierarchy();
    assert(pool.mapPermitSnapshot().availablePermits == 4);
}

// A deterministic pipeline model checks the objective rather than reproducing
// controller branches: pipeline throughput is the slower of its two stages.
static int model(double decodeCost, double processCost) {
    BgzfRateController controller;
    BgzfRateController::Observation o;
    o.budget = 8; o.workers = 8; o.capacity = 16;
    o.decodePending = true; o.processorReady = true;
    int limit = controller.observe(o).limit;
    for (int i = 0; i < 160; ++i) {
        const double decodeRate = limit / decodeCost;
        const double processRate = (8 - limit) / processCost;
        o.completedPairs += static_cast<uint64_t>(std::min(decodeRate, processRate) * o.seconds);
        o.ready = decodeRate < processRate ? 0 : o.capacity;
        if (decodeRate < processRate) o.inputWaitNs += 100000000;
        limit = controller.observe(o).limit;
    }
    assert(limit >= 1 && limit < 8);
    o.budget = 1;
    assert(controller.observe(o).limit == 1); // shrink cannot strand the only permit
    return limit;
}

struct Pipeline {
    std::mutex mutex;
    std::condition_variable ready, space;
    std::deque<int> queue;
    std::atomic<int> next{0}, done{0};
    int producers = 2;
};

static void exercisePool(int budget, const std::string& log) {
    ThreadControl pool;
    pool.mapPermitConfigure(true, budget, budget, true, false);
    pool.mapPermitConfigureDomainFloors({budget / 2, budget - budget / 2, 0});
    pool.mapPermitStartHierarchy(log, true);
    Pipeline pipelines[2];
    std::vector<std::thread> actors;
    constexpr int reads = 96;
    for (int i = 0; i < 2; ++i) {
        auto domain = i ? D::FEATURE : D::MAP;
        auto& p = pipelines[i];
        for (int j = 0; j < 2; ++j) actors.emplace_back([&, domain, i] {
            auto& pipe = pipelines[i];
            for (;;) {
                int read = pipe.next.fetch_add(1);
                if (read >= reads) break;
                {
                    std::unique_lock<std::mutex> lock(pipe.mutex);
                    pipe.space.wait(lock, [&] { return pipe.queue.size() < 8; });
                }
                auto token = pool.mapPermitAcquireForDomain(domain, W::BGZF);
                std::this_thread::sleep_for(std::chrono::microseconds(100));
                pool.mapPermitReleaseForDomain(domain, token, 1, 64, 100000, W::BGZF);
                {
                    std::lock_guard<std::mutex> lock(pipe.mutex);
                    pipe.queue.push_back(read);
                    pool.mapPermitObserveDecode(domain, &pipe, pipe.queue.size(), pipe.queue.size(), 10, 2, 0, 1);
                }
                pipe.ready.notify_all();
            }
            {
                std::lock_guard<std::mutex> lock(pipe.mutex);
                --pipe.producers;
            }
            pipe.ready.notify_all();
        });
        for (int j = 0; j < 3; ++j) actors.emplace_back([&, domain, i] {
            auto& pipe = pipelines[i];
            for (;;) {
                {
                    std::unique_lock<std::mutex> lock(pipe.mutex);
                    pipe.ready.wait(lock, [&] { return !pipe.queue.empty() || !pipe.producers; });
                    if (pipe.queue.empty()) break;
                    pipe.queue.pop_front();
                    pool.mapPermitObserveDecode(domain, &pipe, pipe.queue.size(), pipe.queue.size(), 10, 2, pipe.queue.empty(), 1);
                }
                pipe.space.notify_all();
                auto token = pool.mapPermitAcquireForDomain(domain);
                std::this_thread::sleep_for(std::chrono::microseconds(400));
                const auto snap = pool.mapPermitSnapshot();
                assert(snap.inUsePermits <= budget && snap.availablePermits >= 0);
                assert(snap.mapDomain.decode.inUse <= snap.mapDomain.inUse);
                assert(snap.featureDomain.decode.inUse <= snap.featureDomain.inUse);
                pool.mapPermitReleaseForDomain(domain, token, 1, 64, 400000);
                pipe.done.fetch_add(1);
            }
        });
        (void)p;
    }
    for (auto& actor : actors) actor.join();
    for (int i = 0; i < 2; ++i) {
        assert(pipelines[i].done.load() == reads);
        auto domain = i ? D::FEATURE : D::MAP;
        pool.mapPermitObserveDecode(domain, &pipelines[i], 0, 0, 0, 0, 0, 0);
        pool.mapPermitMarkDomainComplete(domain);
    }
    pool.mapPermitStopHierarchy();
    const auto snap = pool.mapPermitSnapshot();
    assert(snap.availablePermits == budget && snap.inUsePermits == 0 && snap.currentWaiters == 0);
    assert(snap.mapDomain.completedReadPairs == reads && snap.featureDomain.completedReadPairs == reads);
    for (const auto* d : {&snap.mapDomain, &snap.featureDomain}) {
        assert(d->decode.blocks == reads && d->decode.inUse == 0 && d->decode.waiters == 0);
        assert(d->decode.acquireCalls == d->decode.releaseCalls);
        assert(d->floor == 0);
    }
}

int main(int argc, char** argv) {
    assert(argc == 2);
    exerciseOuterPolicy();
    assert(model(0.008, 0.002) == 6); // optimum with expensive decompression
    assert(model(0.001, 0.020) == 1); // optimum with expensive processing
    exercisePool(1, std::string(argv[1]) + ".one.tsv");
    exercisePool(4, std::string(argv[1]) + ".four.tsv");
    exerciseDecoderDemand(true, std::string(argv[1]) + ".starved.tsv");
    exerciseDecoderDemand(false, std::string(argv[1]) + ".full.tsv");
    exerciseOuterRuntime(std::string(argv[1]) + ".outer_runtime.tsv");
}
