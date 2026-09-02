#include "ThreadControl.h"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>

namespace {
inline void atomicStoreMax(std::atomic<uint64_t> &target, uint64_t value) {
    uint64_t current = target.load(std::memory_order_relaxed);
    while (current < value && !target.compare_exchange_weak(current, value, std::memory_order_relaxed)) {
    }
}

inline int normalizePermitCount(int configuredPermits, int totalThreads) {
    const int normalizedTotal = std::max(1, totalThreads);
    int normalizedPermits = configuredPermits;
    if (normalizedPermits <= 0) {
        normalizedPermits = normalizedTotal;
    }
    return std::max(1, std::min(normalizedPermits, normalizedTotal));
}

inline size_t permitDomainIndex(ThreadControl::PermitDomain domain) {
    switch (domain) {
        case ThreadControl::PermitDomain::FEATURE: return 1;
        case ThreadControl::PermitDomain::ATAC:    return 2;
        case ThreadControl::PermitDomain::MAP:
        default:                                   return 0;
    }
}

inline const char* permitDomainName(ThreadControl::PermitDomain domain) {
    switch (domain) {
        case ThreadControl::PermitDomain::FEATURE: return "feature";
        case ThreadControl::PermitDomain::ATAC:    return "atac";
        case ThreadControl::PermitDomain::MAP:
        default:                                   return "map";
    }
}

inline uint64_t steadyNowNs() {
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now().time_since_epoch()
        ).count()
    );
}

constexpr uint64_t kPermitWaitPollNs = 250000000ULL;      // 250ms
constexpr uint64_t kPermitStallWarnEveryNs = 5000000000ULL; // 5s

bool readProcStatCpuTotals(uint64_t &idleAll, uint64_t &totalAll) {
    std::ifstream in("/proc/stat");
    if (!in.is_open()) {
        return false;
    }

    std::string label;
    uint64_t user = 0;
    uint64_t nice = 0;
    uint64_t system = 0;
    uint64_t idle = 0;
    uint64_t iowait = 0;
    uint64_t irq = 0;
    uint64_t softirq = 0;
    uint64_t steal = 0;
    uint64_t guest = 0;
    uint64_t guestNice = 0;
    in >> label >> user >> nice >> system >> idle >> iowait >> irq >> softirq >> steal >> guest >> guestNice;
    if (!in.good() || label != "cpu") {
        return false;
    }

    idleAll = idle + iowait;
    totalAll = user + nice + system + idle + iowait + irq + softirq + steal + guest + guestNice;
    return totalAll >= idleAll;
}
}

ThreadControl::ThreadControl() {
    chunkInN=0;
    chunkOutN=0;
    // Initialize auto-trim synchronization state
    autoTrimPending.store(false);
    autoTrimFileIndex.store(0);
    autoTrimThreadsDone.store(0);
    autoTrimBarrierInitialized = false;
//     chunkOutBAMposition=new uint [MAX_chunkOutBAMposition];
};

void ThreadControl::mapPermitConfigure(
    bool enabled,
    int totalThreads,
    int configuredPermits,
    bool telemetryEnabled,
    bool variableThreads
) {
    std::lock_guard<std::mutex> lock(mapPermitMutex);

    const uint64_t telemetryStartNs = steadyNowNs();

    mapPermitEnabledFlag = enabled;
    mapPermitTelemetryEnabledFlag = telemetryEnabled;
    mapPermitTotalThreads = std::max(1, totalThreads);
    mapPermitVariableThreadsEnabledFlag = mapPermitEnabledFlag && variableThreads;
    mapPermitRetuneEveryAcquires = 0;
    mapPermitRetuneSequence.clear();
    mapPermitRetuneTraceEnabled = mapPermitTelemetryEnabledFlag && mapPermitVariableThreadsEnabledFlag;
    mapPermitRetuneTraceTargets.clear();
    mapPermitRetuneTraceDropped = 0;
    mapPermitFloorsActive = false;
    mapPermitFifoEnabled = false;
    mapPermitWaitQueue.clear();

    if (!mapPermitEnabledFlag) {
        mapPermitConfigured = mapPermitTotalThreads;
        mapPermitTargetPermits = mapPermitConfigured;
        mapPermitAvailable = mapPermitConfigured;
    } else {
        const int normalizedPermits = normalizePermitCount(configuredPermits, mapPermitTotalThreads);
        mapPermitTargetPermits = normalizedPermits;
        mapPermitConfigured = normalizedPermits;
        mapPermitAvailable = normalizedPermits;
    }

    mapPermitAcquireCalls.store(0, std::memory_order_relaxed);
    mapPermitAcquireOrdinal.store(0, std::memory_order_relaxed);
    mapPermitRetuneStep.store(0, std::memory_order_relaxed);
    mapPermitRetuneCalls.store(0, std::memory_order_relaxed);
    mapPermitBlockedAcquireCalls.store(0, std::memory_order_relaxed);
    mapPermitWaitTimeoutEvents.store(0, std::memory_order_relaxed);
    mapPermitStallWarnEvents.store(0, std::memory_order_relaxed);
    mapPermitCurrentWaiters.store(0, std::memory_order_relaxed);
    mapPermitMaxWaiters.store(0, std::memory_order_relaxed);
    mapPermitLastReleaseNs.store(steadyNowNs(), std::memory_order_relaxed);
    mapPermitWaitNsTotal.store(0, std::memory_order_relaxed);
    mapPermitWaitNsMax.store(0, std::memory_order_relaxed);
    mapPermitWorkUnitsTotal.store(0, std::memory_order_relaxed);
    mapPermitWorkBytesTotal.store(0, std::memory_order_relaxed);
    mapPermitWorkNsTotal.store(0, std::memory_order_relaxed);
    mapPermitWorkNsMax.store(0, std::memory_order_relaxed);
    mapPermitTelemetryStartNs = telemetryStartNs;
    mapPermitTelemetryLastStateNs = telemetryStartNs;
    mapPermitAvailablePermitNs = 0;
    mapPermitContendedIdlePermitNs = 0;
    mapPermitNoAdmissibleGrantEvents = 0;
    mapPermitFloorChangeCalls = 0;
    for (size_t domain = 0; domain < mapPermitDomainCount; ++domain) {
        mapPermitDomainFloor[domain] = 0;
        mapPermitDomainInUse[domain] = 0;
        mapPermitDomainWaiters[domain] = 0;
        mapPermitDomainComplete[domain] = false;
        mapPermitAcquireCallsByDomain[domain].store(0, std::memory_order_relaxed);
        mapPermitWaitNsTotalByDomain[domain].store(0, std::memory_order_relaxed);
        mapPermitWaitNsMaxByDomain[domain].store(0, std::memory_order_relaxed);
        mapPermitWorkUnitsTotalByDomain[domain].store(0, std::memory_order_relaxed);
        mapPermitWorkBytesTotalByDomain[domain].store(0, std::memory_order_relaxed);
        mapPermitWorkNsTotalByDomain[domain].store(0, std::memory_order_relaxed);
        mapPermitWorkNsMaxByDomain[domain].store(0, std::memory_order_relaxed);
        mapPermitDomainInUsePermitNs[domain] = 0;
        mapPermitDomainWaiterNs[domain] = 0;
        mapPermitDomainMaxInUse[domain] = 0;
        mapPermitDomainMaxWaiters[domain] = 0;
        mapPermitDomainBlockedAcquireCalls[domain] = 0;
        mapPermitDomainFastAcquireCalls[domain] = 0;
        mapPermitDomainQueuedGrantCalls[domain] = 0;
        mapPermitDomainReleaseCalls[domain] = 0;
    }
}

void ThreadControl::mapPermitAccountStateLocked(uint64_t nowNs) const {
    if (mapPermitTelemetryLastStateNs == 0 ||
        nowNs < mapPermitTelemetryLastStateNs) {
        mapPermitTelemetryLastStateNs = nowNs;
        return;
    }
    const uint64_t elapsedNs = nowNs - mapPermitTelemetryLastStateNs;
    mapPermitTelemetryLastStateNs = nowNs;
    if (!mapPermitTelemetryEnabledFlag || elapsedNs == 0) {
        return;
    }

    const uint64_t available = static_cast<uint64_t>(
        std::max(0, mapPermitAvailable));
    mapPermitAvailablePermitNs += elapsedNs * available;
    int totalWaiters = 0;
    for (size_t domain = 0; domain < mapPermitDomainCount; ++domain) {
        const uint64_t inUse = static_cast<uint64_t>(
            std::max(0, mapPermitDomainInUse[domain]));
        const uint64_t waiters = static_cast<uint64_t>(
            std::max(0, mapPermitDomainWaiters[domain]));
        mapPermitDomainInUsePermitNs[domain] += elapsedNs * inUse;
        mapPermitDomainWaiterNs[domain] += elapsedNs * waiters;
        totalWaiters += mapPermitDomainWaiters[domain];
    }
    if (available > 0 && totalWaiters > 0) {
        // Permit-slot nanoseconds that were simultaneously free and wanted.
        // Any material value here is admission/controller loss, not work-size
        // variation inside a permit holder.
        mapPermitContendedIdlePermitNs += elapsedNs * available;
    }
}

void ThreadControl::mapPermitConfigureCpuAware(bool enabled, int sampleIntervalMs, double emaAlpha) {
    std::lock_guard<std::mutex> lock(mapPermitMutex);
    mapPermitCpuAwareEnabledFlag = enabled;
    mapPermitCpuInitializedFlag = false;
    mapPermitCpuSampleIntervalMs = std::max(1, sampleIntervalMs);
    mapPermitCpuEmaAlpha = std::max(0.0, std::min(1.0, emaAlpha));
    mapPermitCpuLastIdleAll = 0;
    mapPermitCpuLastTotalAll = 0;
    mapPermitCpuLastSampleNs = 0;
    mapPermitCpuSampleCount = 0;
    mapPermitCpuBusyInstant = 0.0;
    mapPermitCpuBusyEma = 0.0;
}

bool ThreadControl::mapPermitCpuMaybeSample() {
    std::lock_guard<std::mutex> lock(mapPermitMutex);
    if (!mapPermitEnabledFlag || !mapPermitCpuAwareEnabledFlag) {
        return false;
    }

    const uint64_t nowNs = steadyNowNs();
    const uint64_t sampleIntervalNs =
        static_cast<uint64_t>(std::max(1, mapPermitCpuSampleIntervalMs)) * 1000000ULL;
    if (mapPermitCpuInitializedFlag &&
        nowNs >= mapPermitCpuLastSampleNs &&
        (nowNs - mapPermitCpuLastSampleNs) < sampleIntervalNs) {
        return false;
    }

    uint64_t idleAll = 0;
    uint64_t totalAll = 0;
    if (!readProcStatCpuTotals(idleAll, totalAll) || totalAll == 0) {
        return false;
    }

    if (!mapPermitCpuInitializedFlag) {
        mapPermitCpuLastIdleAll = idleAll;
        mapPermitCpuLastTotalAll = totalAll;
        mapPermitCpuLastSampleNs = nowNs;
        mapPermitCpuSampleCount = 1;
        mapPermitCpuBusyInstant = 0.0;
        mapPermitCpuBusyEma = 0.0;
        mapPermitCpuInitializedFlag = true;
        return true;
    }

    if (totalAll <= mapPermitCpuLastTotalAll || idleAll < mapPermitCpuLastIdleAll) {
        mapPermitCpuLastIdleAll = idleAll;
        mapPermitCpuLastTotalAll = totalAll;
        mapPermitCpuLastSampleNs = nowNs;
        return false;
    }

    const uint64_t idleDelta = idleAll - mapPermitCpuLastIdleAll;
    const uint64_t totalDelta = totalAll - mapPermitCpuLastTotalAll;
    if (totalDelta == 0) {
        mapPermitCpuLastSampleNs = nowNs;
        return false;
    }

    double busyInstant = 1.0 - (static_cast<double>(idleDelta) / static_cast<double>(totalDelta));
    if (busyInstant < 0.0) {
        busyInstant = 0.0;
    } else if (busyInstant > 1.0) {
        busyInstant = 1.0;
    }

    mapPermitCpuBusyInstant = busyInstant;
    if (mapPermitCpuSampleCount <= 1) {
        mapPermitCpuBusyEma = busyInstant;
    } else {
        mapPermitCpuBusyEma =
            mapPermitCpuEmaAlpha * busyInstant + (1.0 - mapPermitCpuEmaAlpha) * mapPermitCpuBusyEma;
    }
    mapPermitCpuLastIdleAll = idleAll;
    mapPermitCpuLastTotalAll = totalAll;
    mapPermitCpuLastSampleNs = nowNs;
    ++mapPermitCpuSampleCount;
    return true;
}

void ThreadControl::mapPermitConfigureRetunePlan(const std::vector<int> &permitSequence, int retuneEveryAcquires) {
    std::lock_guard<std::mutex> lock(mapPermitMutex);

    mapPermitRetuneSequence.clear();
    mapPermitRetuneEveryAcquires = 0;
    mapPermitRetuneTraceTargets.clear();
    mapPermitRetuneTraceDropped = 0;

    if (!mapPermitEnabledFlag || !mapPermitVariableThreadsEnabledFlag) {
        mapPermitAcquireOrdinal.store(0, std::memory_order_relaxed);
        mapPermitRetuneStep.store(0, std::memory_order_relaxed);
        return;
    }

    if (retuneEveryAcquires <= 0 || permitSequence.empty()) {
        mapPermitAcquireOrdinal.store(0, std::memory_order_relaxed);
        mapPermitRetuneStep.store(0, std::memory_order_relaxed);
        return;
    }

    mapPermitRetuneEveryAcquires = retuneEveryAcquires;
    mapPermitRetuneSequence.reserve(permitSequence.size());
    for (size_t i = 0; i < permitSequence.size(); ++i) {
        mapPermitRetuneSequence.push_back(normalizePermitCount(permitSequence[i], mapPermitTotalThreads));
    }

    mapPermitAcquireOrdinal.store(0, std::memory_order_relaxed);
    mapPermitRetuneStep.store(0, std::memory_order_relaxed);
}

void ThreadControl::mapPermitConfigureDomainFloors(const std::vector<int> &floorsByDomainIndex) {
    std::vector<PermitWaiter*> toWake;
    {
        std::lock_guard<std::mutex> lock(mapPermitMutex);
        mapPermitAccountStateLocked(steadyNowNs());
        int oldFloors[mapPermitDomainCount]{};
        for (size_t i = 0; i < mapPermitDomainCount; ++i) {
            oldFloors[i] = mapPermitDomainFloor[i];
        }
        const bool oldFloorsActive = mapPermitFloorsActive;
        bool anyFloor = false;
        int sumFloors = 0;
        for (size_t i = 0; i < mapPermitDomainCount; ++i) {
            const int requested = (i < floorsByDomainIndex.size()) ? floorsByDomainIndex[i] : 0;
            const int normalized = mapPermitDomainComplete[i]
                ? 0
                : std::max(0, std::min(requested, mapPermitConfigured));
            mapPermitDomainFloor[i] = normalized;
            if (normalized > 0) {
                anyFloor = true;
            }
            sumFloors += normalized;
        }
        // Sanity-clip: combined floors must not exceed the pool itself; if they
        // do, scale them down proportionally so each is at most a fair share.
        if (sumFloors > mapPermitConfigured && mapPermitConfigured > 0) {
            const double scale =
                static_cast<double>(mapPermitConfigured) / static_cast<double>(sumFloors);
            for (size_t i = 0; i < mapPermitDomainCount; ++i) {
                mapPermitDomainFloor[i] =
                    std::max(0, static_cast<int>(mapPermitDomainFloor[i] * scale));
            }
            // Recompute anyFloor in case scaling drove a floor to 0.
            anyFloor = false;
            for (size_t i = 0; i < mapPermitDomainCount; ++i) {
                if (mapPermitDomainFloor[i] > 0) { anyFloor = true; break; }
            }
        }
        mapPermitFloorsActive = mapPermitEnabledFlag && anyFloor;
        bool changed = oldFloorsActive != mapPermitFloorsActive;
        for (size_t i = 0; i < mapPermitDomainCount; ++i) {
            changed = changed || oldFloors[i] != mapPermitDomainFloor[i];
        }
        if (changed) {
            ++mapPermitFloorChangeCalls;
            if (mapPermitFifoEnabled) {
                grantFifoWaitersLocked(toWake);
            }
        }
    }
    for (PermitWaiter *waiter : toWake) {
        waiter->cv.notify_one();
    }
    mapPermitCv.notify_all();
    for (size_t i = 0; i < mapPermitDomainCount; ++i) {
        mapPermitDomainCv[i].notify_all();
    }
}

void ThreadControl::mapPermitConfigureFifoWaiters(bool enabled) {
    std::lock_guard<std::mutex> lock(mapPermitMutex);
    mapPermitFifoEnabled = enabled && mapPermitEnabledFlag;
}

void ThreadControl::mapPermitMarkDomainComplete(PermitDomain domain) {
    std::vector<PermitWaiter*> toWake;
    {
        std::lock_guard<std::mutex> lock(mapPermitMutex);
        mapPermitAccountStateLocked(steadyNowNs());
        const size_t index = permitDomainIndex(domain);
        if (mapPermitDomainComplete[index]) {
            return;
        }
        mapPermitDomainComplete[index] = true;
        if (mapPermitDomainFloor[index] != 0) {
            mapPermitDomainFloor[index] = 0;
            ++mapPermitFloorChangeCalls;
        }
        mapPermitFloorsActive = false;
        for (size_t i = 0; i < mapPermitDomainCount; ++i) {
            mapPermitFloorsActive = mapPermitFloorsActive ||
                mapPermitDomainFloor[i] > 0;
        }
        if (mapPermitFifoEnabled) {
            grantFifoWaitersLocked(toWake);
        }
    }
    for (PermitWaiter *waiter : toWake) {
        waiter->cv.notify_one();
    }
    mapPermitCv.notify_all();
    for (size_t i = 0; i < mapPermitDomainCount; ++i) {
        mapPermitDomainCv[i].notify_all();
    }
}

void ThreadControl::grantFifoWaitersLocked(std::vector<PermitWaiter*> &toWake) {
    // When floors are inactive: strict FIFO from the head.
    // When floors are active: scan the queue and grant the first waiter
    // whose domain can be admitted without starving another domain that
    // has waiters AND is below its floor. Within a single domain this
    // preserves arrival order; across domains it lets a tail waiter on a
    // below-floor domain bypass a head waiter that's already at-floor.
    while (mapPermitAvailable > 0 && !mapPermitWaitQueue.empty()) {
        auto pickIt = mapPermitWaitQueue.end();
        if (!mapPermitFloorsActive) {
            pickIt = mapPermitWaitQueue.begin();
        } else {
            for (auto it = mapPermitWaitQueue.begin();
                 it != mapPermitWaitQueue.end(); ++it) {
                const size_t di = permitDomainIndex((*it)->domain);
                if (mapPermitDomainInUse[di] < mapPermitDomainFloor[di]) {
                    pickIt = it; break;
                }
                bool blockedByFloor = false;
                for (size_t j = 0; j < mapPermitDomainCount; ++j) {
                    if (j == di) continue;
                    if (mapPermitDomainWaiters[j] > 0
                            && mapPermitDomainInUse[j]
                                    < mapPermitDomainFloor[j]) {
                        blockedByFloor = true; break;
                    }
                }
                if (!blockedByFloor) {
                    pickIt = it; break;
                }
            }
            if (pickIt == mapPermitWaitQueue.end()) {
                // No queued waiter can be admitted under current floors;
                // hold the freed permit until a release/retune changes
                // either the floors or the queue composition.
                if (mapPermitTelemetryEnabledFlag) {
                    ++mapPermitNoAdmissibleGrantEvents;
                }
                break;
            }
        }
        PermitWaiter *waiter = *pickIt;
        const size_t di = permitDomainIndex(waiter->domain);
        mapPermitWaitQueue.erase(pickIt);
        --mapPermitDomainWaiters[di];
        ++mapPermitDomainInUse[di];
        if (mapPermitTelemetryEnabledFlag) {
            ++mapPermitDomainQueuedGrantCalls[di];
            mapPermitDomainMaxInUse[di] = std::max<uint64_t>(
                mapPermitDomainMaxInUse[di], mapPermitDomainInUse[di]);
        }
        waiter->granted = true;
        --mapPermitAvailable;
        toWake.push_back(waiter);
    }
}

void ThreadControl::mapPermitSetTargetPermits(int targetPermits) {
    if (!mapPermitEnabledFlag || !mapPermitVariableThreadsEnabledFlag) {
        return;
    }

    std::vector<PermitWaiter*> toWake;
    {
        std::lock_guard<std::mutex> lock(mapPermitMutex);
        mapPermitAccountStateLocked(steadyNowNs());
        const int normalizedTarget = normalizePermitCount(targetPermits, mapPermitTotalThreads);
        mapPermitTargetPermits = normalizedTarget;

        if (normalizedTarget == mapPermitConfigured) {
            return;
        }

        int inUse = mapPermitConfigured - mapPermitAvailable;
        if (inUse < 0) {
            inUse = 0;
        }

        mapPermitConfigured = normalizedTarget;
        int newAvailable = normalizedTarget - inUse;
        if (newAvailable < 0) {
            newAvailable = 0;
        }
        if (newAvailable > normalizedTarget) {
            newAvailable = normalizedTarget;
        }
        mapPermitAvailable = newAvailable;

        // If FIFO is on and the retune raised available permits, hand them
        // to queued waiters under the same lock. On shrink (newAvailable <
        // previous), the helper is a no-op because mapPermitAvailable is
        // already clamped down. The legacy mapPermitCv only wakes the
        // floor-aware/legacy paths; FIFO waiters need their per-waiter CV.
        if (mapPermitFifoEnabled) {
            grantFifoWaitersLocked(toWake);
        }

        if (mapPermitRetuneTraceEnabled) {
            if (mapPermitRetuneTraceTargets.size() < mapPermitRetuneTraceLimit) {
                mapPermitRetuneTraceTargets.push_back(normalizedTarget);
            } else {
                ++mapPermitRetuneTraceDropped;
            }
        }
    }

    for (PermitWaiter *w : toWake) {
        w->cv.notify_one();
    }
    mapPermitRetuneCalls.fetch_add(1, std::memory_order_relaxed);
    mapPermitCv.notify_all();
}

bool ThreadControl::mapPermitEnabled() const {
    std::lock_guard<std::mutex> lock(mapPermitMutex);
    return mapPermitEnabledFlag;
}

uint64_t ThreadControl::mapPermitAcquire() {
    return mapPermitAcquireForDomain(PermitDomain::MAP);
}

uint64_t ThreadControl::mapPermitAcquireForDomain(PermitDomain domain) {
    if (!mapPermitEnabledFlag) {
        return 0;
    }

    if (mapPermitVariableThreadsEnabledFlag && mapPermitRetuneEveryAcquires > 0 && !mapPermitRetuneSequence.empty()) {
        const uint64_t ordinal = mapPermitAcquireOrdinal.fetch_add(1, std::memory_order_relaxed) + 1;
        if (ordinal % static_cast<uint64_t>(mapPermitRetuneEveryAcquires) == 0) {
            const uint64_t step = mapPermitRetuneStep.fetch_add(1, std::memory_order_relaxed);
            const int nextTarget = mapPermitRetuneSequence[step % mapPermitRetuneSequence.size()];
            mapPermitSetTargetPermits(nextTarget);
        }
    }

    const auto waitStart = std::chrono::steady_clock::now();
    bool blocked = false;
    uint64_t nextWarnNs = kPermitStallWarnEveryNs;
    const size_t domainIndex = permitDomainIndex(domain);

    // FIFO waiter-queue path. New arrivals can only fast-path when the
    // pool has a permit AND the queue is empty AND (when floors are
    // active) admitting this domain wouldn't violate another domain's
    // floor. Otherwise they enqueue at the tail and wait for an
    // under-the-lock grant from release/retune, eliminating the
    // notify_one() bypass race.
    if (mapPermitFifoEnabled) {
        std::unique_lock<std::mutex> lock(mapPermitMutex);
        // Domain-aware fast path: only take a free permit if the queue is
        // empty AND admitting this domain wouldn't violate another domain's
        // floor. When floors are inactive this collapses to the simple check.
        auto canAdmitFifo = [&]() -> bool {
            if (mapPermitAvailable <= 0) return false;
            if (!mapPermitFloorsActive) return true;
            if (mapPermitDomainInUse[domainIndex]
                    < mapPermitDomainFloor[domainIndex]) {
                return true;
            }
            for (size_t i = 0; i < mapPermitDomainCount; ++i) {
                if (i == domainIndex) continue;
                if (mapPermitDomainWaiters[i] > 0
                        && mapPermitDomainInUse[i]
                                < mapPermitDomainFloor[i]) {
                    return false;
                }
            }
            return true;
        };
        if (mapPermitWaitQueue.empty() && canAdmitFifo()) {
            mapPermitAccountStateLocked(steadyNowNs());
            --mapPermitAvailable;
            ++mapPermitDomainInUse[domainIndex];
            if (mapPermitTelemetryEnabledFlag) {
                ++mapPermitDomainFastAcquireCalls[domainIndex];
                mapPermitDomainMaxInUse[domainIndex] = std::max<uint64_t>(
                    mapPermitDomainMaxInUse[domainIndex],
                    mapPermitDomainInUse[domainIndex]);
            }
        } else {
            PermitWaiter waiter;
            waiter.domain = domain;
            mapPermitAccountStateLocked(steadyNowNs());
            mapPermitWaitQueue.push_back(&waiter);
            mapPermitBlockedAcquireCalls.fetch_add(1, std::memory_order_relaxed);
            ++mapPermitDomainWaiters[domainIndex];
            if (mapPermitTelemetryEnabledFlag) {
                ++mapPermitDomainBlockedAcquireCalls[domainIndex];
                mapPermitDomainMaxWaiters[domainIndex] = std::max<uint64_t>(
                    mapPermitDomainMaxWaiters[domainIndex],
                    mapPermitDomainWaiters[domainIndex]);
            }
            const uint64_t waiters =
                mapPermitCurrentWaiters.fetch_add(1, std::memory_order_relaxed) + 1;
            atomicStoreMax(mapPermitMaxWaiters, waiters);
            while (!waiter.granted) {
                if (waiter.cv.wait_for(
                        lock, std::chrono::nanoseconds(kPermitWaitPollNs)) ==
                    std::cv_status::timeout) {
                    if (mapPermitTelemetryEnabledFlag) {
                        const uint64_t elapsedNs = static_cast<uint64_t>(
                            std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now() - waitStart)
                                .count());
                        if (elapsedNs >= nextWarnNs) {
                            const int inUse =
                                mapPermitConfigured - mapPermitAvailable;
                            const uint64_t sinceReleaseNs =
                                steadyNowNs() -
                                mapPermitLastReleaseNs.load(
                                    std::memory_order_relaxed);
                            mapPermitStallWarnEvents.fetch_add(
                                1, std::memory_order_relaxed);
                            std::fprintf(
                                stderr,
                                "WARNING: dynamic-permit stall domain=%s "
                                "waitMs=%.3f configured=%d target=%d "
                                "available=%d inUse=%d waiters=%llu retunes=%llu "
                                "sinceReleaseMs=%.3f fifo=on\n",
                                permitDomainName(domain), elapsedNs / 1.0e6,
                                mapPermitConfigured, mapPermitTargetPermits,
                                mapPermitAvailable, inUse,
                                static_cast<unsigned long long>(
                                    mapPermitCurrentWaiters.load(
                                        std::memory_order_relaxed)),
                                static_cast<unsigned long long>(
                                    mapPermitRetuneCalls.load(
                                        std::memory_order_relaxed)),
                                sinceReleaseNs / 1.0e6);
                            nextWarnNs += kPermitStallWarnEveryNs;
                        }
                    }
                }
            }
            // Release/grantFifoWaitersLocked already popped us from the
            // queue, pre-decremented mapPermitAvailable, decremented our
            // domain waiter count, and incremented our domain in-use count
            // — all under the same lock.
            mapPermitCurrentWaiters.fetch_sub(1, std::memory_order_relaxed);
        }
        // Common acquire-end accounting (telemetry).
        const auto waitEnd = std::chrono::steady_clock::now();
        const uint64_t waitNs = static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(waitEnd -
                                                                  waitStart)
                .count());
        if (mapPermitTelemetryEnabledFlag) {
            mapPermitAcquireCalls.fetch_add(1, std::memory_order_relaxed);
            mapPermitAcquireCallsByDomain[domainIndex].fetch_add(
                1, std::memory_order_relaxed);
        }
        return waitNs;
    }

    {
        std::unique_lock<std::mutex> lock(mapPermitMutex);
        // Floor-aware admission: a domain that has reached its floor cannot
        // acquire more permits while another domain with waiters is still
        // below its floor. Otherwise (no floors set, or all floors satisfied),
        // standard FCFS-ish admission via the legacy global cv.
        auto canAdmit = [&]() -> bool {
            if (mapPermitAvailable <= 0) {
                return false;
            }
            if (!mapPermitFloorsActive) {
                return true;
            }
            // Always admit up to D's own floor.
            if (mapPermitDomainInUse[domainIndex] < mapPermitDomainFloor[domainIndex]) {
                return true;
            }
            // D is at or above its floor. Yield if any other domain has
            // waiters AND is below its floor.
            for (size_t i = 0; i < mapPermitDomainCount; ++i) {
                if (i == domainIndex) continue;
                if (mapPermitDomainWaiters[i] > 0 &&
                    mapPermitDomainInUse[i] < mapPermitDomainFloor[i]) {
                    return false;
                }
            }
            return true;
        };
        while (!canAdmit()) {
            if (!blocked) {
                blocked = true;
                mapPermitAccountStateLocked(steadyNowNs());
                mapPermitBlockedAcquireCalls.fetch_add(1, std::memory_order_relaxed);
                const uint64_t waiters = mapPermitCurrentWaiters.fetch_add(1, std::memory_order_relaxed) + 1;
                atomicStoreMax(mapPermitMaxWaiters, waiters);
                ++mapPermitDomainWaiters[domainIndex];
                if (mapPermitTelemetryEnabledFlag) {
                    ++mapPermitDomainBlockedAcquireCalls[domainIndex];
                    mapPermitDomainMaxWaiters[domainIndex] = std::max<uint64_t>(
                        mapPermitDomainMaxWaiters[domainIndex],
                        mapPermitDomainWaiters[domainIndex]);
                }
            }

            if (mapPermitFloorsActive) {
                // Wait on this domain's CV — release wakes the specific
                // domain whose floor is unmet, eliminating the notify_one()
                // wakeup-fairness pathology.
                mapPermitDomainCv[domainIndex].wait_for(
                    lock, std::chrono::nanoseconds(kPermitWaitPollNs));
            } else {
                mapPermitCv.wait_for(lock, std::chrono::nanoseconds(kPermitWaitPollNs));
            }
            if (canAdmit()) {
                break;
            }

            mapPermitWaitTimeoutEvents.fetch_add(1, std::memory_order_relaxed);

            if (mapPermitTelemetryEnabledFlag) {
                const uint64_t elapsedNs = static_cast<uint64_t>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - waitStart
                    ).count()
                );
                if (elapsedNs >= nextWarnNs) {
                    const int inUse = mapPermitConfigured - mapPermitAvailable;
                    const uint64_t sinceReleaseNs = steadyNowNs() - mapPermitLastReleaseNs.load(std::memory_order_relaxed);
                    mapPermitStallWarnEvents.fetch_add(1, std::memory_order_relaxed);
                    std::fprintf(
                        stderr,
                        "WARNING: dynamic-permit stall domain=%s waitMs=%.3f configured=%d target=%d available=%d inUse=%d waiters=%llu retunes=%llu sinceReleaseMs=%.3f\n",
                        permitDomainName(domain),
                        elapsedNs / 1.0e6,
                        mapPermitConfigured,
                        mapPermitTargetPermits,
                        mapPermitAvailable,
                        inUse,
                        static_cast<unsigned long long>(mapPermitCurrentWaiters.load(std::memory_order_relaxed)),
                        static_cast<unsigned long long>(mapPermitRetuneCalls.load(std::memory_order_relaxed)),
                        sinceReleaseNs / 1.0e6
                    );
                    nextWarnNs += kPermitStallWarnEveryNs;
                }
            }
        }
        mapPermitAccountStateLocked(steadyNowNs());
        if (blocked) {
            mapPermitCurrentWaiters.fetch_sub(1, std::memory_order_relaxed);
            --mapPermitDomainWaiters[domainIndex];
            if (mapPermitTelemetryEnabledFlag) {
                ++mapPermitDomainQueuedGrantCalls[domainIndex];
            }
        } else if (mapPermitTelemetryEnabledFlag) {
            ++mapPermitDomainFastAcquireCalls[domainIndex];
        }
        --mapPermitAvailable;
        ++mapPermitDomainInUse[domainIndex];
        if (mapPermitTelemetryEnabledFlag) {
            mapPermitDomainMaxInUse[domainIndex] = std::max<uint64_t>(
                mapPermitDomainMaxInUse[domainIndex],
                mapPermitDomainInUse[domainIndex]);
        }
    }
    const auto waitEnd = std::chrono::steady_clock::now();
    const uint64_t waitNs = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(waitEnd - waitStart).count());

    if (mapPermitTelemetryEnabledFlag) {
        mapPermitAcquireCalls.fetch_add(1, std::memory_order_relaxed);
        mapPermitAcquireCallsByDomain[permitDomainIndex(domain)].fetch_add(1, std::memory_order_relaxed);
    }

    return waitNs;
}

void ThreadControl::mapPermitRelease(uint64_t waitNs, uint64_t workUnits, uint64_t workBytes, uint64_t workNs) {
    mapPermitReleaseForDomain(PermitDomain::MAP, waitNs, workUnits, workBytes, workNs);
}

void ThreadControl::mapPermitReleaseForDomain(
    PermitDomain domain,
    uint64_t waitNs,
    uint64_t workUnits,
    uint64_t workBytes,
    uint64_t workNs
) {
    if (!mapPermitEnabledFlag) {
        return;
    }

    // FIFO release: grant the next permit directly to the head waiter,
    // under the lock, before any new arrival can race in.
    if (mapPermitFifoEnabled) {
        std::vector<PermitWaiter*> toWake;
        {
            std::lock_guard<std::mutex> lock(mapPermitMutex);
            mapPermitAccountStateLocked(steadyNowNs());
            ++mapPermitAvailable;
            if (mapPermitAvailable > mapPermitConfigured) {
                mapPermitAvailable = mapPermitConfigured;
            }
            const size_t releasingIdx = permitDomainIndex(domain);
            if (mapPermitDomainInUse[releasingIdx] > 0) {
                --mapPermitDomainInUse[releasingIdx];
            }
            if (mapPermitTelemetryEnabledFlag) {
                ++mapPermitDomainReleaseCalls[releasingIdx];
            }
            grantFifoWaitersLocked(toWake);
        }
        mapPermitLastReleaseNs.store(steadyNowNs(), std::memory_order_relaxed);
        for (PermitWaiter *w : toWake) {
            w->cv.notify_one();
        }

        if (mapPermitTelemetryEnabledFlag) {
            const size_t domainIndex = permitDomainIndex(domain);
            mapPermitWaitNsTotal.fetch_add(waitNs, std::memory_order_relaxed);
            atomicStoreMax(mapPermitWaitNsMax, waitNs);
            mapPermitWaitNsTotalByDomain[domainIndex].fetch_add(waitNs, std::memory_order_relaxed);
            atomicStoreMax(mapPermitWaitNsMaxByDomain[domainIndex], waitNs);
            mapPermitWorkUnitsTotal.fetch_add(workUnits, std::memory_order_relaxed);
            mapPermitWorkUnitsTotalByDomain[domainIndex].fetch_add(workUnits, std::memory_order_relaxed);
            mapPermitWorkBytesTotal.fetch_add(workBytes, std::memory_order_relaxed);
            mapPermitWorkBytesTotalByDomain[domainIndex].fetch_add(workBytes, std::memory_order_relaxed);
            mapPermitWorkNsTotal.fetch_add(workNs, std::memory_order_relaxed);
            mapPermitWorkNsTotalByDomain[domainIndex].fetch_add(workNs, std::memory_order_relaxed);
            atomicStoreMax(mapPermitWorkNsMax, workNs);
            atomicStoreMax(mapPermitWorkNsMaxByDomain[domainIndex], workNs);
        }
        return;
    }

    int notifyDomain = -1; // -1 = use legacy global cv
    {
        std::lock_guard<std::mutex> lock(mapPermitMutex);
        mapPermitAccountStateLocked(steadyNowNs());
        ++mapPermitAvailable;
        if (mapPermitAvailable > mapPermitConfigured) {
            mapPermitAvailable = mapPermitConfigured;
        }
        const size_t releasingIdx = permitDomainIndex(domain);
        if (mapPermitDomainInUse[releasingIdx] > 0) {
            --mapPermitDomainInUse[releasingIdx];
        }
        if (mapPermitTelemetryEnabledFlag) {
            ++mapPermitDomainReleaseCalls[releasingIdx];
        }
        if (mapPermitFloorsActive) {
            // Pick the wake target: priority is any domain with waiters AND
            // inUse < floor (under-floor with waiters). Else any domain with
            // waiters (most-waiters first as a simple fairness tiebreak).
            int underFloorIdx = -1;
            int mostWaitersIdx = -1;
            int mostWaiters = 0;
            for (size_t i = 0; i < mapPermitDomainCount; ++i) {
                if (mapPermitDomainWaiters[i] <= 0) continue;
                if (mapPermitDomainInUse[i] < mapPermitDomainFloor[i] &&
                    underFloorIdx < 0) {
                    underFloorIdx = static_cast<int>(i);
                }
                if (mapPermitDomainWaiters[i] > mostWaiters) {
                    mostWaiters = mapPermitDomainWaiters[i];
                    mostWaitersIdx = static_cast<int>(i);
                }
            }
            notifyDomain = (underFloorIdx >= 0) ? underFloorIdx : mostWaitersIdx;
        }
    }
    mapPermitLastReleaseNs.store(steadyNowNs(), std::memory_order_relaxed);
    if (notifyDomain >= 0) {
        mapPermitDomainCv[notifyDomain].notify_one();
    } else {
        mapPermitCv.notify_one();
    }

    if (mapPermitTelemetryEnabledFlag) {
        const size_t domainIndex = permitDomainIndex(domain);
        mapPermitWaitNsTotal.fetch_add(waitNs, std::memory_order_relaxed);
        atomicStoreMax(mapPermitWaitNsMax, waitNs);
        mapPermitWaitNsTotalByDomain[domainIndex].fetch_add(waitNs, std::memory_order_relaxed);
        atomicStoreMax(mapPermitWaitNsMaxByDomain[domainIndex], waitNs);
        mapPermitWorkUnitsTotal.fetch_add(workUnits, std::memory_order_relaxed);
        mapPermitWorkUnitsTotalByDomain[domainIndex].fetch_add(workUnits, std::memory_order_relaxed);
        mapPermitWorkBytesTotal.fetch_add(workBytes, std::memory_order_relaxed);
        mapPermitWorkBytesTotalByDomain[domainIndex].fetch_add(workBytes, std::memory_order_relaxed);
        mapPermitWorkNsTotal.fetch_add(workNs, std::memory_order_relaxed);
        mapPermitWorkNsTotalByDomain[domainIndex].fetch_add(workNs, std::memory_order_relaxed);
        atomicStoreMax(mapPermitWorkNsMax, workNs);
        atomicStoreMax(mapPermitWorkNsMaxByDomain[domainIndex], workNs);
    }
}

ThreadControl::MapPermitSnapshot ThreadControl::mapPermitSnapshot() const {
    MapPermitSnapshot snapshot{};
    snapshot.enabled = mapPermitEnabledFlag;
    snapshot.telemetryEnabled = mapPermitTelemetryEnabledFlag;
    snapshot.variableThreadsEnabled = mapPermitVariableThreadsEnabledFlag;
    snapshot.cpuAwareEnabled = mapPermitCpuAwareEnabledFlag;
    uint64_t cpuLastSampleNs = 0;

    {
        std::lock_guard<std::mutex> lock(mapPermitMutex);
        const uint64_t snapshotNowNs = steadyNowNs();
        mapPermitAccountStateLocked(snapshotNowNs);
        snapshot.retuneEveryAcquires = mapPermitRetuneEveryAcquires;
        snapshot.sequenceLength = static_cast<int>(mapPermitRetuneSequence.size());
        snapshot.targetPermits = mapPermitTargetPermits;
        snapshot.configuredPermits = mapPermitConfigured;
        snapshot.availablePermits = mapPermitAvailable;
        snapshot.inUsePermits = mapPermitConfigured - mapPermitAvailable;
        snapshot.floorsActive = mapPermitFloorsActive;
        snapshot.fifoEnabled = mapPermitFifoEnabled;
        snapshot.fifoQueueDepth = mapPermitWaitQueue.size();
        snapshot.telemetryElapsedNs =
            snapshotNowNs >= mapPermitTelemetryStartNs
                ? snapshotNowNs - mapPermitTelemetryStartNs
                : 0;
        snapshot.availablePermitNs = mapPermitAvailablePermitNs;
        snapshot.contendedIdlePermitNs = mapPermitContendedIdlePermitNs;
        snapshot.noAdmissibleGrantEvents = mapPermitNoAdmissibleGrantEvents;
        snapshot.floorChangeCalls = mapPermitFloorChangeCalls;
        snapshot.retuneTraceTargets = mapPermitRetuneTraceTargets;
        snapshot.retuneTraceDropped = mapPermitRetuneTraceDropped;
        snapshot.cpuInitialized = mapPermitCpuInitializedFlag;
        snapshot.cpuSampleIntervalMs = mapPermitCpuSampleIntervalMs;
        snapshot.cpuSampleCount = mapPermitCpuSampleCount;
        snapshot.cpuBusyInstant = mapPermitCpuBusyInstant;
        snapshot.cpuBusyEma = mapPermitCpuBusyEma;
        cpuLastSampleNs = mapPermitCpuLastSampleNs;
        PermitDomainSnapshot *domainSnapshots[mapPermitDomainCount] = {
            &snapshot.mapDomain, &snapshot.featureDomain, &snapshot.atacDomain};
        for (size_t domain = 0; domain < mapPermitDomainCount; ++domain) {
            domainSnapshots[domain]->complete = mapPermitDomainComplete[domain];
            domainSnapshots[domain]->floor = mapPermitDomainFloor[domain];
            domainSnapshots[domain]->inUse = mapPermitDomainInUse[domain];
            domainSnapshots[domain]->currentWaiters =
                mapPermitDomainWaiters[domain];
            domainSnapshots[domain]->maxInUse =
                mapPermitDomainMaxInUse[domain];
            domainSnapshots[domain]->maxWaiters =
                mapPermitDomainMaxWaiters[domain];
            domainSnapshots[domain]->blockedAcquireCalls =
                mapPermitDomainBlockedAcquireCalls[domain];
            domainSnapshots[domain]->fastAcquireCalls =
                mapPermitDomainFastAcquireCalls[domain];
            domainSnapshots[domain]->queuedGrantCalls =
                mapPermitDomainQueuedGrantCalls[domain];
            domainSnapshots[domain]->releaseCalls =
                mapPermitDomainReleaseCalls[domain];
            domainSnapshots[domain]->inUsePermitNs =
                mapPermitDomainInUsePermitNs[domain];
            domainSnapshots[domain]->waiterNs =
                mapPermitDomainWaiterNs[domain];
        }
    }

    snapshot.acquireCalls = mapPermitAcquireCalls.load(std::memory_order_relaxed);
    snapshot.retuneCalls = mapPermitRetuneCalls.load(std::memory_order_relaxed);
    snapshot.blockedAcquireCalls = mapPermitBlockedAcquireCalls.load(std::memory_order_relaxed);
    snapshot.waitTimeoutEvents = mapPermitWaitTimeoutEvents.load(std::memory_order_relaxed);
    snapshot.stallWarnEvents = mapPermitStallWarnEvents.load(std::memory_order_relaxed);
    snapshot.currentWaiters = mapPermitCurrentWaiters.load(std::memory_order_relaxed);
    snapshot.maxWaiters = mapPermitMaxWaiters.load(std::memory_order_relaxed);
    const uint64_t lastReleaseNs = mapPermitLastReleaseNs.load(std::memory_order_relaxed);
    const uint64_t nowNs = steadyNowNs();
    snapshot.lastReleaseAgoNs = (nowNs >= lastReleaseNs) ? (nowNs - lastReleaseNs) : 0;
    snapshot.cpuLastSampleAgoNs = (snapshot.cpuInitialized && nowNs >= cpuLastSampleNs)
        ? (nowNs - cpuLastSampleNs)
        : 0;
    snapshot.cpuIdleEma = std::max(0.0, 1.0 - snapshot.cpuBusyEma);
    snapshot.waitNsTotal = mapPermitWaitNsTotal.load(std::memory_order_relaxed);
    snapshot.waitNsMax = mapPermitWaitNsMax.load(std::memory_order_relaxed);
    snapshot.workUnitsTotal = mapPermitWorkUnitsTotal.load(std::memory_order_relaxed);
    snapshot.workBytesTotal = mapPermitWorkBytesTotal.load(std::memory_order_relaxed);
    snapshot.workNsTotal = mapPermitWorkNsTotal.load(std::memory_order_relaxed);
    snapshot.workNsMax = mapPermitWorkNsMax.load(std::memory_order_relaxed);
    snapshot.mapDomain.acquireCalls = mapPermitAcquireCallsByDomain[permitDomainIndex(PermitDomain::MAP)].load(std::memory_order_relaxed);
    snapshot.mapDomain.waitNsTotal = mapPermitWaitNsTotalByDomain[permitDomainIndex(PermitDomain::MAP)].load(std::memory_order_relaxed);
    snapshot.mapDomain.waitNsMax = mapPermitWaitNsMaxByDomain[permitDomainIndex(PermitDomain::MAP)].load(std::memory_order_relaxed);
    snapshot.mapDomain.workUnitsTotal = mapPermitWorkUnitsTotalByDomain[permitDomainIndex(PermitDomain::MAP)].load(std::memory_order_relaxed);
    snapshot.mapDomain.workBytesTotal = mapPermitWorkBytesTotalByDomain[permitDomainIndex(PermitDomain::MAP)].load(std::memory_order_relaxed);
    snapshot.mapDomain.workNsTotal = mapPermitWorkNsTotalByDomain[permitDomainIndex(PermitDomain::MAP)].load(std::memory_order_relaxed);
    snapshot.mapDomain.workNsMax = mapPermitWorkNsMaxByDomain[permitDomainIndex(PermitDomain::MAP)].load(std::memory_order_relaxed);
    snapshot.featureDomain.acquireCalls = mapPermitAcquireCallsByDomain[permitDomainIndex(PermitDomain::FEATURE)].load(std::memory_order_relaxed);
    snapshot.featureDomain.waitNsTotal = mapPermitWaitNsTotalByDomain[permitDomainIndex(PermitDomain::FEATURE)].load(std::memory_order_relaxed);
    snapshot.featureDomain.waitNsMax = mapPermitWaitNsMaxByDomain[permitDomainIndex(PermitDomain::FEATURE)].load(std::memory_order_relaxed);
    snapshot.featureDomain.workUnitsTotal = mapPermitWorkUnitsTotalByDomain[permitDomainIndex(PermitDomain::FEATURE)].load(std::memory_order_relaxed);
    snapshot.featureDomain.workBytesTotal = mapPermitWorkBytesTotalByDomain[permitDomainIndex(PermitDomain::FEATURE)].load(std::memory_order_relaxed);
    snapshot.featureDomain.workNsTotal = mapPermitWorkNsTotalByDomain[permitDomainIndex(PermitDomain::FEATURE)].load(std::memory_order_relaxed);
    snapshot.featureDomain.workNsMax = mapPermitWorkNsMaxByDomain[permitDomainIndex(PermitDomain::FEATURE)].load(std::memory_order_relaxed);
    snapshot.atacDomain.acquireCalls = mapPermitAcquireCallsByDomain[permitDomainIndex(PermitDomain::ATAC)].load(std::memory_order_relaxed);
    snapshot.atacDomain.waitNsTotal = mapPermitWaitNsTotalByDomain[permitDomainIndex(PermitDomain::ATAC)].load(std::memory_order_relaxed);
    snapshot.atacDomain.waitNsMax = mapPermitWaitNsMaxByDomain[permitDomainIndex(PermitDomain::ATAC)].load(std::memory_order_relaxed);
    snapshot.atacDomain.workUnitsTotal = mapPermitWorkUnitsTotalByDomain[permitDomainIndex(PermitDomain::ATAC)].load(std::memory_order_relaxed);
    snapshot.atacDomain.workBytesTotal = mapPermitWorkBytesTotalByDomain[permitDomainIndex(PermitDomain::ATAC)].load(std::memory_order_relaxed);
    snapshot.atacDomain.workNsTotal = mapPermitWorkNsTotalByDomain[permitDomainIndex(PermitDomain::ATAC)].load(std::memory_order_relaxed);
    snapshot.atacDomain.workNsMax = mapPermitWorkNsMaxByDomain[permitDomainIndex(PermitDomain::ATAC)].load(std::memory_order_relaxed);
    return snapshot;
}
