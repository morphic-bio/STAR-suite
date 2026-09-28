#ifndef STAR_HOST_PERMIT_TYPES_H
#define STAR_HOST_PERMIT_TYPES_H

// Plain snapshot types of STAR's shared permit pool (ThreadControl).
// Public part of the STAR Suite host interface (docs/HOST_API.md); this
// header has no STAR dependencies so a host program can include it alone.

#include <cstdint>
#include <vector>

namespace star {
namespace host {

// The three permit domains of STAR's pool. MAP and FEATURE are STAR's own
// mapping and feature-assignment work. EXTERNAL is lent to a host program
// for work it runs beside STAR in the same process; STAR never uses it.
enum class Domain : uint8_t { Map = 0, Feature = 1, External = 2 };

struct DecodeSnapshot {
    int inUse = 0, waiters = 0, limit = 1;
    uint64_t acquireCalls = 0, releaseCalls = 0, maxInUse = 0;
    uint64_t blocks = 0, bytes = 0, workNs = 0, waitNs = 0;
    uint64_t ready = 0, outstanding = 0, capacity = 0, workers = 0;
    uint64_t waitingReaders = 0, inputWaitNs = 0;
};

struct PermitDomainSnapshot {
    DecodeSnapshot decode;
    uint64_t completedReadPairs = 0;
    bool complete;
    int floor;
    int inUse;
    int currentWaiters;
    uint64_t maxInUse;
    uint64_t maxWaiters;
    uint64_t blockedAcquireCalls;
    uint64_t fastAcquireCalls;
    uint64_t queuedGrantCalls;
    uint64_t releaseCalls;
    uint64_t inUsePermitNs;
    uint64_t waiterNs;
    uint64_t acquireCalls;
    uint64_t waitNsTotal;
    uint64_t waitNsMax;
    uint64_t workUnitsTotal;
    uint64_t workBytesTotal;
    uint64_t workNsTotal;
    uint64_t workNsMax;
};

struct PermitSnapshot {
    bool enabled;
    bool telemetryEnabled;
    bool variableThreadsEnabled;
    bool cpuAwareEnabled;
    bool cpuInitialized;
    bool floorsActive;
    bool fifoEnabled;
    int retuneEveryAcquires;
    int sequenceLength;
    int targetPermits;
    int configuredPermits;
    int availablePermits;
    int inUsePermits;
    uint64_t fifoQueueDepth;
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
    uint64_t telemetryElapsedNs;
    uint64_t availablePermitNs;
    uint64_t contendedIdlePermitNs;
    uint64_t noAdmissibleGrantEvents;
    uint64_t floorChangeCalls;
    uint64_t cpuSampleCount;
    uint64_t cpuLastSampleAgoNs;
    double cpuBusyInstant;
    double cpuBusyEma;
    double cpuIdleEma;
    PermitDomainSnapshot mapDomain;
    PermitDomainSnapshot featureDomain;
    PermitDomainSnapshot externalDomain;
};

}  // namespace host
}  // namespace star

#endif  // STAR_HOST_PERMIT_TYPES_H
