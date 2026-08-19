#include "ThreadControl.h"
#include "SaturationPermitController.h"

#include <cassert>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>
#include <vector>

int main() {
    // Saturation-aware policy: observe the larger ATAC arm first, preserve
    // its sustained demand, then probe MAP with the remainder. When both
    // observed demands fit, ETA noise must not move the learned floors.
    {
        using Controller = star::multiome::SaturationPermitController;
        Controller controller(32, 0, 2);
        Controller::Decision initial = controller.initialDecision();
        assert(initial.phase == Controller::Phase::PROBE_ATAC);
        assert(initial.mapFloor == 1);
        assert(initial.atacFloor == 31);

        Controller::Observation probe{};
        probe.atacUnitsDelta = 1000;
        probe.atacOccupancy = 20.4;
        assert(!controller.observe(probe).floorsChanged);
        probe.atacOccupancy = 21.0;
        Controller::Decision atacLearned = controller.observe(probe);
        assert(atacLearned.floorsChanged);
        assert(atacLearned.phase == Controller::Phase::PROBE_MAP);
        assert(atacLearned.atacSaturation == 21);
        assert(atacLearned.atacSaturationKnown);
        assert(atacLearned.mapFloor == 11);

        probe = Controller::Observation{};
        probe.mapUnitsDelta = 1000;
        probe.mapOccupancy = 1.0;
        probe.mapInUse = 1;
        probe.mapWaiters = 10;
        probe.atacInUse = 31;
        assert(!controller.observe(probe).floorsChanged);

        probe = Controller::Observation{};
        probe.mapUnitsDelta = 1000;
        probe.mapOccupancy = 9.2;
        probe.mapInUse = 10;
        probe.mapWaiters = 1;
        probe.atacInUse = 21;
        assert(!controller.observe(probe).floorsChanged);
        probe.mapOccupancy = 9.6;
        Controller::Decision bothLearned = controller.observe(probe);
        assert(bothLearned.floorsChanged);
        assert(bothLearned.phase == Controller::Phase::STEADY);
        assert(bothLearned.mapSaturation == 10);
        assert(bothLearned.mapSaturationKnown);
        assert(!bothLearned.capacityLimited);
        assert(bothLearned.mapFloor == 10);
        assert(bothLearned.atacFloor == 21);

        probe.mapUnitsDelta = 1000;
        probe.atacUnitsDelta = 1000;
        probe.mapEtaSec = 1000.0;
        probe.atacEtaSec = 10.0;
        assert(!controller.observe(probe).floorsChanged);
    }

    // If both probes consume every permit made available to them, remaining
    // work becomes the coarse tie-breaker. It moves one borrowable floor per
    // window and retains at least one permit for each live arm.
    {
        using Controller = star::multiome::SaturationPermitController;
        Controller controller(32, 0, 2);
        Controller::Observation probe{};
        probe.atacUnitsDelta = 1000;
        probe.atacOccupancy = 31.0;
        controller.observe(probe);
        Controller::Decision atacLimited = controller.observe(probe);
        assert(!atacLimited.atacSaturationKnown);
        assert(atacLimited.mapFloor == 1);

        probe = Controller::Observation{};
        probe.mapUnitsDelta = 1000;
        probe.mapOccupancy = 1.0;
        controller.observe(probe);
        Controller::Decision capacityLimited = controller.observe(probe);
        assert(capacityLimited.capacityLimited);
        assert(capacityLimited.mapFloor == 1);
        assert(capacityLimited.atacFloor == 31);

        probe.atacUnitsDelta = 1000;
        probe.mapEtaSec = 300.0;
        probe.atacEtaSec = 100.0;
        Controller::Decision helpMap = controller.observe(probe);
        assert(helpMap.floorsChanged);
        assert(helpMap.reason == Controller::Reason::MAP_ETA_LATE);
        assert(helpMap.mapFloor == 2);
        assert(helpMap.atacFloor == 30);

        probe.mapEtaSec = 50.0;
        probe.atacEtaSec = 200.0;
        Controller::Decision helpAtac = controller.observe(probe);
        assert(helpAtac.floorsChanged);
        assert(helpAtac.reason == Controller::Reason::MAP_ETA_EARLY);
        assert(helpAtac.mapFloor == 1);
        assert(helpAtac.atacFloor == 31);
    }

    // Occupancy accounting is independent of floors/FIFO. The diagnostic
    // must still report who owns a permit when the legacy admission path is
    // active, otherwise a comparison run could hide the very imbalance being
    // investigated.
    {
        ThreadControl legacy;
        legacy.mapPermitConfigure(true, 4, 2, true, false);
        const uint64_t wait = legacy.mapPermitAcquireForDomain(
            ThreadControl::PermitDomain::ATAC);
        ThreadControl::MapPermitSnapshot held = legacy.mapPermitSnapshot();
        assert(!held.floorsActive);
        assert(!held.fifoEnabled);
        assert(held.atacDomain.inUse == 1);
        assert(held.atacDomain.maxInUse == 1);
        legacy.mapPermitReleaseForDomain(
            ThreadControl::PermitDomain::ATAC, wait, 1, 64, 1000);
        assert(legacy.mapPermitSnapshot().atacDomain.inUse == 0);
    }

    // Floors are minimum reservations, not a conserved allocation. With a
    // 22-permit pool and 1/1 MAP/ATAC floors, MAP may borrow all 22 while
    // ATAC has no waiter. A controller that changes only floors therefore
    // does not define an 11/11 target split.
    {
        ThreadControl shared;
        shared.mapPermitConfigure(true, 44, 22, true, false);
        shared.mapPermitConfigureDomainFloors(std::vector<int>{1, 0, 1});
        shared.mapPermitConfigureFifoWaiters(true);
        std::vector<uint64_t> waits;
        for (int i = 0; i < 22; ++i) {
            waits.push_back(shared.mapPermitAcquireForDomain(
                ThreadControl::PermitDomain::MAP));
        }
        ThreadControl::MapPermitSnapshot borrowed = shared.mapPermitSnapshot();
        assert(borrowed.mapDomain.inUse == 22);
        assert(borrowed.atacDomain.inUse == 0);
        assert(borrowed.mapDomain.floor == 1);
        assert(borrowed.atacDomain.floor == 1);
        assert(borrowed.floorChangeCalls == 1);
        for (uint64_t waitNs : waits) {
            shared.mapPermitReleaseForDomain(
                ThreadControl::PermitDomain::MAP, waitNs, 1, 64, 1000);
        }
        shared.mapPermitConfigureDomainFloors(std::vector<int>{0, 0, 6});
        assert(shared.mapPermitSnapshot().floorChangeCalls == 2);
    }

    ThreadControl permits;
    permits.mapPermitConfigure(true, 4, 2, true, false);
    permits.mapPermitConfigureDomainFloors(std::vector<int>{1, 0, 1});
    permits.mapPermitConfigureFifoWaiters(true);

    const uint64_t firstMapWait =
        permits.mapPermitAcquireForDomain(ThreadControl::PermitDomain::MAP);
    const uint64_t atacWait =
        permits.mapPermitAcquireForDomain(ThreadControl::PermitDomain::ATAC);

    ThreadControl::MapPermitSnapshot saturated = permits.mapPermitSnapshot();
    assert(saturated.configuredPermits == 2);
    assert(saturated.availablePermits == 0);
    assert(saturated.mapDomain.floor == 1);
    assert(saturated.atacDomain.floor == 1);
    assert(saturated.mapDomain.inUse == 1);
    assert(saturated.atacDomain.inUse == 1);
    assert(saturated.mapDomain.fastAcquireCalls == 1);
    assert(saturated.atacDomain.fastAcquireCalls == 1);

    std::mutex stateMutex;
    std::condition_variable stateCv;
    bool waiterStarted = false;
    bool waiterFinished = false;
    std::thread waiter([&]() {
        {
            std::lock_guard<std::mutex> lock(stateMutex);
            waiterStarted = true;
        }
        stateCv.notify_all();
        const uint64_t waitNs =
            permits.mapPermitAcquireForDomain(ThreadControl::PermitDomain::MAP);
        permits.mapPermitReleaseForDomain(ThreadControl::PermitDomain::MAP,
                                          waitNs, 1, 64, 1000);
        {
            std::lock_guard<std::mutex> lock(stateMutex);
            waiterFinished = true;
        }
        stateCv.notify_all();
    });

    {
        std::unique_lock<std::mutex> lock(stateMutex);
        stateCv.wait(lock, [&]() { return waiterStarted; });
    }
    ThreadControl::MapPermitSnapshot queued{};
    for (int attempt = 0; attempt < 100; ++attempt) {
        queued = permits.mapPermitSnapshot();
        if (queued.mapDomain.currentWaiters == 1) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    assert(queued.mapDomain.currentWaiters == 1);
    assert(queued.fifoQueueDepth == 1);
    assert(queued.mapDomain.blockedAcquireCalls == 1);

    permits.mapPermitReleaseForDomain(ThreadControl::PermitDomain::ATAC,
                                      atacWait, 1, 64, 1000);
    {
        std::unique_lock<std::mutex> lock(stateMutex);
        stateCv.wait(lock, [&]() { return waiterFinished; });
    }
    waiter.join();
    permits.mapPermitReleaseForDomain(ThreadControl::PermitDomain::MAP,
                                      firstMapWait, 1, 64, 1000);

    const ThreadControl::MapPermitSnapshot done = permits.mapPermitSnapshot();
    assert(done.availablePermits == 2);
    assert(done.inUsePermits == 0);
    assert(done.mapDomain.inUse == 0);
    assert(done.atacDomain.inUse == 0);
    assert(done.mapDomain.currentWaiters == 0);
    assert(done.fifoQueueDepth == 0);
    assert(done.mapDomain.acquireCalls == 2);
    assert(done.mapDomain.fastAcquireCalls == 1);
    assert(done.mapDomain.queuedGrantCalls == 1);
    assert(done.mapDomain.releaseCalls == 2);
    assert(done.atacDomain.releaseCalls == 1);
    assert(done.mapDomain.maxInUse == 2);
    assert(done.mapDomain.maxWaiters == 1);
    assert(done.mapDomain.inUsePermitNs > 0);
    assert(done.mapDomain.waiterNs > 0);
    assert(done.noAdmissibleGrantEvents == 0);
    return 0;
}
