#ifndef THREAD_CONTROL_DEF
#define THREAD_CONTROL_DEF

#include "ReadAlignChunk.h"
#include <pthread.h>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <condition_variable>
#include <mutex>
#include <vector>

#define MAX_chunkOutBAMposition 100000

class ThreadControl {
public:
    enum class PermitDomain : uint8_t {
        MAP = 0,
        FEATURE = 1,
        ATAC = 2
    };

    struct PermitHookContext {
        PermitDomain domain;
    };

    struct PermitDomainSnapshot {
        uint64_t acquireCalls;
        uint64_t waitNsTotal;
        uint64_t waitNsMax;
        uint64_t workUnitsTotal;
        uint64_t workBytesTotal;
        uint64_t workNsTotal;
        uint64_t workNsMax;
    };

    bool threadBool;

    pthread_t *threadArray;
    pthread_mutex_t mutexInRead, mutexOutSAM, mutexOutBAM1, mutexOutChimSAM, mutexOutChimJunction, mutexOutUnmappedFastx, mutexOutFilterBySJout;
    pthread_mutex_t mutexStats, mutexLogMain, mutexBAMsortBins, mutexError;
    pthread_mutex_t mutexOutYFastq[MAX_N_MATES], mutexOutNoYFastq[MAX_N_MATES];  // Y/noY FASTQ output mutexes per mate
    
    // Auto-trim synchronization
    pthread_barrier_t autoTrimBarrier;      // Barrier for file boundary synchronization
    pthread_mutex_t mutexAutoTrim;          // Mutex for auto-trim computation
    std::atomic<bool> autoTrimPending;      // Flag indicating trim computation needed
    std::atomic<uint32_t> autoTrimFileIndex; // File index that triggered trim
    std::atomic<int> autoTrimThreadsDone;   // Threads that completed current file
    bool autoTrimBarrierInitialized;        // Whether barrier was initialized

    uint chunkInN,chunkOutN;

    ThreadControl();
    
    struct MapPermitSnapshot {
        bool enabled;
        bool telemetryEnabled;
        bool variableThreadsEnabled;
        bool cpuAwareEnabled;
        bool cpuInitialized;
        int retuneEveryAcquires;
        int sequenceLength;
        int targetPermits;
        int configuredPermits;
        int availablePermits;
        int inUsePermits;
        int cpuSampleIntervalMs;
        std::vector<int> retuneTraceTargets;
        uint64_t retuneTraceDropped;
        uint64_t acquireCalls;
        uint64_t retuneCalls;
        uint64_t blockedAcquireCalls;
        uint64_t waitTimeoutEvents;
        uint64_t stallWarnEvents;
        uint64_t currentWaiters;
        uint64_t maxWaiters;
        uint64_t lastReleaseAgoNs;
        uint64_t waitNsTotal;
        uint64_t waitNsMax;
        uint64_t workUnitsTotal;
        uint64_t workBytesTotal;
        uint64_t workNsTotal;
        uint64_t workNsMax;
        uint64_t cpuSampleCount;
        uint64_t cpuLastSampleAgoNs;
        double cpuBusyInstant;
        double cpuBusyEma;
        double cpuIdleEma;
        PermitDomainSnapshot mapDomain;
        PermitDomainSnapshot featureDomain;
        PermitDomainSnapshot atacDomain;
    };

    void mapPermitConfigure(bool enabled, int totalThreads, int configuredPermits, bool telemetryEnabled, bool variableThreads);
    void mapPermitConfigureCpuAware(bool enabled, int sampleIntervalMs, double emaAlpha);
    void mapPermitConfigureRetunePlan(const std::vector<int> &permitSequence, int retuneEveryAcquires);
    // Per-domain borrowable floors (Step 5a in
    // multiomic-atac-scrna plans/2026-04-27-atac-permits-controller-followups.md
    // v6). Each domain reserves at LEAST `floor` concurrent permits whenever
    // it has waiters; surplus permits are fully shared. floors must be
    // length mapPermitDomainCount, indexed by permitDomainIndex(domain).
    // Pass all-zero floors (or skip configuring) to disable the floor-aware
    // path and fall through to the legacy global-cv behavior.
    void mapPermitConfigureDomainFloors(const std::vector<int> &floorsByDomainIndex);
    void mapPermitSetTargetPermits(int targetPermits);
    bool mapPermitEnabled() const;
    bool mapPermitCpuMaybeSample();
    uint64_t mapPermitAcquireForDomain(PermitDomain domain);
    uint64_t mapPermitAcquire();
    void mapPermitReleaseForDomain(PermitDomain domain, uint64_t waitNs, uint64_t workUnits, uint64_t workBytes, uint64_t workNs);
    void mapPermitRelease(uint64_t waitNs, uint64_t workUnits, uint64_t workBytes, uint64_t workNs);
    MapPermitSnapshot mapPermitSnapshot() const;

    static void* threadRAprocessChunks(void *RAchunk) {
        ( (ReadAlignChunk*) RAchunk )->processChunks();
        pthread_exit(0);
        return NULL;
    };

private:
    static constexpr size_t mapPermitDomainCount = 3;
    bool mapPermitEnabledFlag = false;
    bool mapPermitTelemetryEnabledFlag = false;
    bool mapPermitVariableThreadsEnabledFlag = false;
    int mapPermitRetuneEveryAcquires = 0;
    std::vector<int> mapPermitRetuneSequence;
    bool mapPermitRetuneTraceEnabled = false;
    size_t mapPermitRetuneTraceLimit = 64;
    std::vector<int> mapPermitRetuneTraceTargets;
    uint64_t mapPermitRetuneTraceDropped = 0;
    bool mapPermitCpuAwareEnabledFlag = false;
    bool mapPermitCpuInitializedFlag = false;
    int mapPermitCpuSampleIntervalMs = 500;
    double mapPermitCpuEmaAlpha = 0.2;
    uint64_t mapPermitCpuLastIdleAll = 0;
    uint64_t mapPermitCpuLastTotalAll = 0;
    uint64_t mapPermitCpuLastSampleNs = 0;
    uint64_t mapPermitCpuSampleCount = 0;
    double mapPermitCpuBusyInstant = 0.0;
    double mapPermitCpuBusyEma = 0.0;
    int mapPermitTotalThreads = 1;
    int mapPermitTargetPermits = 1;
    int mapPermitConfigured = 0;
    int mapPermitAvailable = 0;
    mutable std::mutex mapPermitMutex;
    std::condition_variable mapPermitCv;
    // Per-domain borrowable-floor state. Active iff any element of
    // mapPermitDomainFloor is > 0 (mapPermitFloorsActive). The floor-aware
    // acquire/release path uses one CV per domain so release can wake the
    // specific domain whose floor is unmet, eliminating the notify_one()
    // wakeup-fairness pathology under contention.
    bool mapPermitFloorsActive = false;
    int mapPermitDomainFloor[mapPermitDomainCount]{};
    int mapPermitDomainInUse[mapPermitDomainCount]{};
    int mapPermitDomainWaiters[mapPermitDomainCount]{};
    std::condition_variable mapPermitDomainCv[mapPermitDomainCount];

    std::atomic<uint64_t> mapPermitAcquireCalls{0};
    std::atomic<uint64_t> mapPermitAcquireOrdinal{0};
    std::atomic<uint64_t> mapPermitRetuneStep{0};
    std::atomic<uint64_t> mapPermitRetuneCalls{0};
    std::atomic<uint64_t> mapPermitBlockedAcquireCalls{0};
    std::atomic<uint64_t> mapPermitWaitTimeoutEvents{0};
    std::atomic<uint64_t> mapPermitStallWarnEvents{0};
    std::atomic<uint64_t> mapPermitCurrentWaiters{0};
    std::atomic<uint64_t> mapPermitMaxWaiters{0};
    std::atomic<uint64_t> mapPermitLastReleaseNs{0};
    std::atomic<uint64_t> mapPermitWaitNsTotal{0};
    std::atomic<uint64_t> mapPermitWaitNsMax{0};
    std::atomic<uint64_t> mapPermitWorkUnitsTotal{0};
    std::atomic<uint64_t> mapPermitWorkBytesTotal{0};
    std::atomic<uint64_t> mapPermitWorkNsTotal{0};
    std::atomic<uint64_t> mapPermitWorkNsMax{0};
    std::atomic<uint64_t> mapPermitAcquireCallsByDomain[mapPermitDomainCount]{};
    std::atomic<uint64_t> mapPermitWaitNsTotalByDomain[mapPermitDomainCount]{};
    std::atomic<uint64_t> mapPermitWaitNsMaxByDomain[mapPermitDomainCount]{};
    std::atomic<uint64_t> mapPermitWorkUnitsTotalByDomain[mapPermitDomainCount]{};
    std::atomic<uint64_t> mapPermitWorkBytesTotalByDomain[mapPermitDomainCount]{};
    std::atomic<uint64_t> mapPermitWorkNsTotalByDomain[mapPermitDomainCount]{};
    std::atomic<uint64_t> mapPermitWorkNsMaxByDomain[mapPermitDomainCount]{};
};

#endif
