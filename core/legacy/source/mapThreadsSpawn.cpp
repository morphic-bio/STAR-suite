#include "mapThreadsSpawn.h"
#include "FlexPipeline.h"
#include "FlexHashScreen.h"
#include "SoloReadFeature.h"
#include "SoloFeatureTypes.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include <algorithm>
#include <sstream>
#include <zlib.h>

namespace {
std::string serializeRetuneTrace(const std::vector<int> &traceTargets) {
    if (traceTargets.empty()) {
        return "-";
    }

    std::ostringstream traceStream;
    for (size_t i = 0; i < traceTargets.size(); ++i) {
        if (i > 0) {
            traceStream << "|";
        }
        traceStream << traceTargets[i];
    }
    return traceStream.str();
}
}

static bool flexPipelineActivationGuard(Parameters &P) {
    const auto &ps = P.pSolo;
    if (ps.flexPipelineStr == "no") {
        P.inOut->logMain << "Flex pipeline: disabled by --flexPipeline no\n" << std::flush;
        return false;
    }
    if (!ps.flexMode) {
        P.inOut->logMain << "Flex pipeline: not active (flexMode=false)\n" << std::flush;
        return false;
    }
    if (!ps.hashScreenEnabled) {
        P.inOut->logMain << "Flex pipeline: not active (hash screen not enabled)\n" << std::flush;
        return false;
    }
    if (P.outSAMtype.at(0) != "None") {
        P.inOut->logMain << "Flex pipeline: not active (outSAMtype is not None)\n" << std::flush;
        return false;
    }
    int nLanes = static_cast<int>(P.readFilesN);
    int nTriage = ps.flexPipelineNTriage;
    int nSolo = ps.flexPipelineNSolo;
    bool fullyFused = (nTriage == 0 && nSolo == 0);
    int minThreads = fullyFused ? 1 : (nLanes + nTriage + nSolo + 1);
    if (P.runThreadN < minThreads) {
        P.inOut->logMain << "Flex pipeline: not active (runThreadN=" << P.runThreadN
                         << " < minimum " << minThreads << " for " << nLanes << " lanes + "
                         << nTriage << " triage + "
                         << nSolo << " solo + 1 worker)\n" << std::flush;
        return false;
    }
    return true;
}

static void mapThreadsSpawnFlexPipeline(Parameters &P, ReadAlignChunk** RAchunk) {
    // Kill legacy FIFO decompression children (zcat/gzip) — the pipeline
    // opens its own gzFile handles and does not use the FIFO streams.
    P.closeReadsFiles();

    const int nLanes = static_cast<int>(P.readFilesN);
    const int nTriage = P.pSolo.flexPipelineNTriage;
    const int nSolo = P.pSolo.flexPipelineNSolo;

    // Mode selection: nTriage=0 + nSolo=0 → fully fused (read+route+solo inline)
    //                 nTriage=0 + nSolo>0 → fused reader+router, separate solo consumers
    //                 nTriage>0            → separate reader, triage, solo threads
    const bool fullyFused = (nTriage == 0 && nSolo == 0);
    const bool fusedReaderRouter = (nTriage == 0 && nSolo > 0);
    const int actualNTriage = (nTriage == 0) ? 0 : nTriage;
    const int actualNSolo = fullyFused ? 0 : nSolo;
    const int nWorkers = fullyFused ? 0 : (P.runThreadN - nLanes - actualNTriage - actualNSolo);
    const int nFusedThreads = fullyFused ? P.runThreadN : 0;

    const bool noAlign = (P.pSolo.flexNoAlign != 0);
    P.inOut->logMain << "Flex pipeline: runThreadN=" << P.runThreadN
                     << ", nLanes=" << nLanes
                     << ", triage=" << actualNTriage
                     << (fullyFused ? " (fully fused + lane steal + role switch)" :
                         fusedReaderRouter ? " (fused reader+router)" : "")
                     << ", soloConsumers=" << actualNSolo
                     << ", fusedThreads=" << nFusedThreads
                     << ", dedicatedWorkers=" << nWorkers
                     << (noAlign ? ", noAlign=ON (alignment skipped)" : "")
                     << "\n" << std::flush;

    FlexPipelineState state;
    state.init(nLanes, actualNSolo, actualNTriage);

    // Store lane file paths for dynamic claiming (fully-fused mode)
    if (fullyFused) {
        state.laneFiles.resize(nLanes);
        for (int lane = 0; lane < nLanes; ++lane) {
            state.laneFiles[lane].r2path = P.readFilesNames[0][lane];
            state.laneFiles[lane].r1path = P.readFilesNames[1][lane];
        }
        state.nFusedThreads = nFusedThreads;
    }

    // Per-thread SoloReadFeature + Stats for fully-fused mode
    std::vector<SoloReadFeature *> fusedFeats(nFusedThreads);
    std::vector<Stats> fusedStats(nFusedThreads);

    if (fullyFused) {
        // All runThreadN threads are fused: they steal lanes, then switch to alignment
        std::vector<pthread_t> fusedThreads(nFusedThreads);
        std::vector<FlexLaneReaderArgs> fusedArgs(nFusedThreads);

        for (int i = 0; i < nFusedThreads; ++i) {
            fusedFeats[i] = new SoloReadFeature(SoloFeatureTypes::Gene, P, -(200 + i));
            fusedStats[i] = Stats();
            fusedArgs[i].state = &state;
            fusedArgs[i].P = &P;
            fusedArgs[i].gzR2 = nullptr;
            fusedArgs[i].gzR1 = nullptr;
            fusedArgs[i].laneId = -1;
            fusedArgs[i].readFeat = fusedFeats[i];
            fusedArgs[i].stats = &fusedStats[i];
            fusedArgs[i].RA = RAchunk[i]->RA;
            fusedArgs[i].threadId = i;
            pthread_create(&fusedThreads[i], nullptr, flexLaneReaderFullThread, &fusedArgs[i]);
        }
        P.inOut->logMain << "  " << nFusedThreads << " fused threads started (lane-steal + role-switch)\n" << std::flush;

        // Stats reporter
        pthread_t reporterThread;
        FlexStatsReporterArgs reporterArgs;
        reporterArgs.state = &state;
        reporterArgs.P = &P;
        pthread_create(&reporterThread, nullptr, flexStatsReporterThread, &reporterArgs);

        // Join all fused threads (they self-transition from reader to aligner)
        for (int i = 0; i < nFusedThreads; ++i) pthread_join(fusedThreads[i], nullptr);
        P.inOut->logMain << "  All fused threads joined\n" << std::flush;

        state.pipelineDone.store(true, std::memory_order_relaxed);
        pthread_join(reporterThread, nullptr);
        P.inOut->logMain << "  Stats reporter joined\n" << std::flush;

        // Merge stats from all fused threads
        for (int i = 0; i < nFusedThreads; ++i) {
            g_statsAll.addStats(fusedStats[i]);
            if (!noAlign)
                g_statsAll.addStats(RAchunk[i]->RA->statsRA);
        }
        g_statsAll.readN = state.counters.readsTotal.load();

        // Merge Solo hashes
        if (P.pSolo.inlineHashMode && P.pSolo.featureYes[SoloFeatureTypes::Gene]) {
            SoloReadFeature *workerGeneFeat = RAchunk[0]->RA->soloRead->readFeat[P.pSolo.featureInd[SoloFeatureTypes::Gene]];
            for (int i = 0; i < nFusedThreads; ++i) {
                workerGeneFeat->mergeInlineHash(*fusedFeats[i]);
            }
        }

        for (int i = 0; i < nFusedThreads; ++i) delete fusedFeats[i];

        P.inOut->logMain << "Flex pipeline complete: total=" << state.counters.readsTotal.load()
                         << ", triageKeep=" << state.counters.triageKeep.load()
                         << ", triageDeny=" << state.counters.triageDeny.load()
                         << ", triageMiss=" << state.counters.triageMiss.load() << "\n";
        for (int i = 0; i < nLanes; ++i) {
            P.inOut->logMain << "  Lane " << i << ": " << state.counters.perLaneReads[i] << " reads\n";
        }
        P.inOut->logMain << std::flush;
        return;
    }

    // --- Non-fused paths (legacy pipeline modes) ---

    // Per-lane SoloReadFeature + Stats (not used in non-fused)
    std::vector<SoloReadFeature *> laneFeats;
    std::vector<Stats> laneStats;

    std::vector<pthread_t> readerThreads(nLanes);
    std::vector<FlexLaneReaderArgs> readerArgs(nLanes);
    for (int lane = 0; lane < nLanes; ++lane) {
        const std::string &r2path = P.readFilesNames[0][lane];
        const std::string &r1path = P.readFilesNames[1][lane];
        gzFile gzR2 = gzopen(r2path.c_str(), "rb");
        gzFile gzR1 = gzopen(r1path.c_str(), "rb");
        if (!gzR2 || !gzR1) {
            std::ostringstream errOut;
            errOut << "EXITING: Flex pipeline cannot open lane " << lane
                   << " files: " << r2path << " / " << r1path << "\n";
            exitWithError(errOut.str(), std::cerr, P.inOut->logMain, 1, P);
        }
        gzbuffer(gzR2, 1 << 20);
        gzbuffer(gzR1, 1 << 20);

        readerArgs[lane].state = &state;
        readerArgs[lane].P = &P;
        readerArgs[lane].gzR2 = gzR2;
        readerArgs[lane].gzR1 = gzR1;
        readerArgs[lane].laneId = lane;
        readerArgs[lane].readFeat = nullptr;
        readerArgs[lane].stats = nullptr;
        readerArgs[lane].RA = nullptr;
        readerArgs[lane].threadId = lane;

        if (fusedReaderRouter) {
            pthread_create(&readerThreads[lane], nullptr, flexLaneReaderRouterThread, &readerArgs[lane]);
            P.inOut->logMain << "  Flex reader+router lane " << lane << " started\n" << std::flush;
        } else {
            pthread_create(&readerThreads[lane], nullptr, flexLaneReaderThread, &readerArgs[lane]);
            P.inOut->logMain << "  Flex reader lane " << lane << " started\n" << std::flush;
        }
    }

    // --- Stage 2: Triage (skipped when any fused mode) ---
    std::vector<pthread_t> triageThreads(actualNTriage);
    std::vector<FlexTriageArgs> triageArgs(actualNTriage);
    for (int i = 0; i < actualNTriage; ++i) {
        triageArgs[i].state = &state;
        triageArgs[i].P = &P;
        pthread_create(&triageThreads[i], nullptr, flexTriageThread, &triageArgs[i]);
        P.inOut->logMain << "  Flex triage " << i << " started\n" << std::flush;
    }

    // --- Stage 3a: Solo consumers ---
    std::vector<pthread_t> soloThreads(actualNSolo);
    std::vector<FlexSoloConsumerArgs> soloArgs(actualNSolo);
    std::vector<SoloReadFeature *> soloFeats(actualNSolo);
    std::vector<Stats> soloStats(actualNSolo);

    for (int i = 0; i < actualNSolo; ++i) {
        soloFeats[i] = new SoloReadFeature(SoloFeatureTypes::Gene, P, -(100 + i));
        soloArgs[i].state = &state;
        soloArgs[i].P = &P;
        soloArgs[i].consumerId = i;
        soloArgs[i].readFeat = soloFeats[i];
        soloArgs[i].stats = &soloStats[i];
        pthread_create(&soloThreads[i], nullptr, flexSoloConsumerThread, &soloArgs[i]);
        P.inOut->logMain << "  Flex solo consumer " << i << " started\n" << std::flush;
    }

    // --- Stage 3b: Alignment workers ---
    std::vector<pthread_t> workerThreads(nWorkers);
    std::vector<FlexAlignWorkerArgs> workerArgs(nWorkers);
    for (int i = 0; i < nWorkers; ++i) {
        workerArgs[i].state = &state;
        workerArgs[i].RA = RAchunk[i]->RA;
        pthread_create(&workerThreads[i], nullptr, flexAlignWorkerThread, &workerArgs[i]);
        P.inOut->logMain << "  Flex align worker " << i << " started\n" << std::flush;
    }

    // --- Stats reporter (does not count against runThreadN) ---
    pthread_t reporterThread;
    FlexStatsReporterArgs reporterArgs;
    reporterArgs.state = &state;
    reporterArgs.P = &P;
    pthread_create(&reporterThread, nullptr, flexStatsReporterThread, &reporterArgs);
    P.inOut->logMain << "  Flex stats reporter started\n" << std::flush;

    // --- Join all threads ---
    for (int i = 0; i < nLanes; ++i) pthread_join(readerThreads[i], nullptr);
    P.inOut->logMain << "  All readers joined\n" << std::flush;

    for (int i = 0; i < actualNTriage; ++i) pthread_join(triageThreads[i], nullptr);
    if (actualNTriage > 0)
        P.inOut->logMain << "  All triage threads joined\n" << std::flush;

    for (int i = 0; i < actualNSolo; ++i) pthread_join(soloThreads[i], nullptr);
    if (actualNSolo > 0)
        P.inOut->logMain << "  All solo consumers joined\n" << std::flush;

    for (int i = 0; i < nWorkers; ++i) pthread_join(workerThreads[i], nullptr);
    P.inOut->logMain << "  All alignment workers joined\n" << std::flush;

    // Stop and join stats reporter
    state.pipelineDone.store(true, std::memory_order_relaxed);
    pthread_join(reporterThread, nullptr);
    P.inOut->logMain << "  Stats reporter joined\n" << std::flush;

    // --- Merge stats ---
    for (int i = 0; i < actualNSolo; ++i) {
        g_statsAll.addStats(soloStats[i]);
    }
    for (int i = 0; i < nWorkers; ++i) {
        g_statsAll.addStats(RAchunk[i]->RA->statsRA);
    }

    g_statsAll.readN = state.counters.readsTotal.load();

    // --- Merge Solo hashes into per-feature sum ---
    if (P.pSolo.inlineHashMode && P.pSolo.featureYes[SoloFeatureTypes::Gene]) {
        SoloReadFeature *workerGeneFeat = RAchunk[0]->RA->soloRead->readFeat[P.pSolo.featureInd[SoloFeatureTypes::Gene]];
        for (int i = 0; i < actualNSolo; ++i) {
            workerGeneFeat->mergeInlineHash(*soloFeats[i]);
        }
    }

    // Cleanup
    for (int i = 0; i < actualNSolo; ++i) delete soloFeats[i];

    // Log pipeline counters
    P.inOut->logMain << "Flex pipeline complete: total=" << state.counters.readsTotal.load()
                     << ", triageKeep=" << state.counters.triageKeep.load()
                     << ", triageDeny=" << state.counters.triageDeny.load()
                     << ", triageMiss=" << state.counters.triageMiss.load() << "\n";
    for (int i = 0; i < nLanes; ++i) {
        P.inOut->logMain << "  Lane " << i << ": " << state.counters.perLaneReads[i] << " reads\n";
    }
    P.inOut->logMain << std::flush;
}

void mapThreadsSpawn (Parameters &P, ReadAlignChunk** RAchunk) {
    // Check activation guard for Flex pipeline mode
    if (flexPipelineActivationGuard(P)) {
        mapThreadsSpawnFlexPipeline(P, RAchunk);
        return;
    }

    const bool interfaceEnabled = (P.dynamicThreadInterface == 1);
    const bool telemetryEnabled = (P.dynamicThreadTelemetry == 1);
    const bool variableThreadsEnabled = (P.variableThreads == 1);
    // Permit-pool budget. With chromapAtac concurrent the pool spans STAR's
    // GEX MAP/FEATURE workers AND chromap's ATAC workers, so the budget is
    // runThreadN + chromapAtac.threads (a separate thread budget than runThreadN
    // alone). When chromapAtac is off, pool stays at runThreadN. The user can
    // override via --dynamicThreadConstMapPermits, which is now honored as-is
    // (no clamp to runThreadN); pass 0 for the auto-sized default.
    const int permitTotalThreads = (P.chromapAtac.enabled == 1)
        ? (P.runThreadN + std::max(0, P.chromapAtac.threads))
        : P.runThreadN;
    const int configuredPermits = (P.dynamicThreadConstMapPermits > 0)
        ? P.dynamicThreadConstMapPermits
        : permitTotalThreads;
    g_threadChunks.mapPermitConfigure(interfaceEnabled, permitTotalThreads, configuredPermits, telemetryEnabled, variableThreadsEnabled);
    g_threadChunks.mapPermitConfigureCpuAware(
        interfaceEnabled && P.dynamicThreadPfControllerCpuAware == 1,
        P.dynamicThreadPfControllerCpuSampleMs,
        P.dynamicThreadPfControllerCpuEmaAlpha
    );
    g_threadChunks.mapPermitConfigureRetunePlan(P.variableThreadsPermitSequence, P.variableThreadsRetuneEveryAcquires);

    // Per-domain borrowable floors (Step 5a). Index order must match
    // ThreadControl::permitDomainIndex(): MAP=0, FEATURE=1, ATAC=2.
    {
        std::vector<int> domainFloors(3, 0);
        domainFloors[0] = std::max(0, P.dynamicThreadMapFloor);
        domainFloors[1] = std::max(0, P.dynamicThreadFeatureFloor);
        domainFloors[2] = std::max(0, P.dynamicThreadAtacFloor);
        g_threadChunks.mapPermitConfigureDomainFloors(domainFloors);
    }
    // FIFO waiter queue (Step 7). When enabled, ThreadControl serves
    // queued waiters in strict arrival order; new arrivals cannot
    // fast-path past existing waiters.
    g_threadChunks.mapPermitConfigureFifoWaiters(P.dynamicThreadFifoWaiters == 1);

    if (interfaceEnabled) {
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Dynamic thread interface enabled: map permits=" << configuredPermits
                         << " (runThreadN=" << P.runThreadN
                         << ", telemetry=" << (telemetryEnabled ? "on" : "off")
                         << ", variableThreads=" << (variableThreadsEnabled ? "on" : "off")
                         << ", cpuAware=" << ((P.dynamicThreadPfControllerCpuAware == 1) ? "on" : "off")
                         << ", cpuSampleMs=" << P.dynamicThreadPfControllerCpuSampleMs
                         << ", cpuEmaAlpha=" << P.dynamicThreadPfControllerCpuEmaAlpha
                         << ", retuneEveryAcquires=" << P.variableThreadsRetuneEveryAcquires
                         << ", retuneSequenceLength=" << P.variableThreadsPermitSequence.size() << ")\n" << flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    }

    for (int ithread=1;ithread<P.runThreadN;ithread++) {//spawn threads
        int threadStatus=pthread_create(&g_threadChunks.threadArray[ithread], NULL, &g_threadChunks.threadRAprocessChunks, (void *) RAchunk[ithread]);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while creating thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Created thread # " <<ithread <<"\n"<<flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    };

    RAchunk[0]->processChunks(); //start main thread

    for (int ithread=1;ithread<P.runThreadN;ithread++) {//wait for all threads to complete
        int threadStatus = pthread_join(g_threadChunks.threadArray[ithread], NULL);
        if (threadStatus>0) {//something went wrong with one of threads
                ostringstream errOut;
                errOut << "EXITING because of FATAL ERROR: phtread error while joining thread # " << ithread <<", error code: "<<threadStatus ;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, 1, P);
        };
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Joined thread # " <<ithread <<"\n"<<flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    };

    if (interfaceEnabled && telemetryEnabled) {
        ThreadControl::MapPermitSnapshot snapshot = g_threadChunks.mapPermitSnapshot();
        const double waitMsTotal = snapshot.waitNsTotal / 1.0e6;
        const double waitMsMax = snapshot.waitNsMax / 1.0e6;
        const double workMsTotal = snapshot.workNsTotal / 1.0e6;
        const double workMsMax = snapshot.workNsMax / 1.0e6;
        const double lastReleaseAgoMs = snapshot.lastReleaseAgoNs / 1.0e6;
        const double avgWaitMs = snapshot.acquireCalls > 0 ? waitMsTotal / static_cast<double>(snapshot.acquireCalls) : 0.0;
        const double avgWorkMs = snapshot.acquireCalls > 0 ? workMsTotal / static_cast<double>(snapshot.acquireCalls) : 0.0;

        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Dynamic thread telemetry: acquires=" << snapshot.acquireCalls
                         << ", retunes=" << snapshot.retuneCalls
                         << ", blockedAcquires=" << snapshot.blockedAcquireCalls
                         << ", waitTimeouts=" << snapshot.waitTimeoutEvents
                         << ", stallWarnings=" << snapshot.stallWarnEvents
                         << ", waiters(current/max)=" << snapshot.currentWaiters << "/" << snapshot.maxWaiters
                         << ", retuneEveryAcquires=" << snapshot.retuneEveryAcquires
                         << ", retuneSequenceLength=" << snapshot.sequenceLength
                         << ", targetPermits=" << snapshot.targetPermits
                         << ", configuredPermits=" << snapshot.configuredPermits
                         << ", inUsePermits=" << snapshot.inUsePermits
                         << ", cpuAware=" << (snapshot.cpuAwareEnabled ? "on" : "off")
                         << ", cpuInitialized=" << (snapshot.cpuInitialized ? "yes" : "no")
                         << ", cpuSampleCount=" << snapshot.cpuSampleCount
                         << ", cpuBusyInstant=" << snapshot.cpuBusyInstant
                         << ", cpuBusyEma=" << snapshot.cpuBusyEma
                         << ", cpuIdleEma=" << snapshot.cpuIdleEma
                         << ", retuneTrace=" << serializeRetuneTrace(snapshot.retuneTraceTargets)
                         << ", retuneTraceDropped=" << snapshot.retuneTraceDropped
                         << ", lastReleaseAgoMs=" << lastReleaseAgoMs
                         << ", workUnits=" << snapshot.workUnitsTotal
                         << ", workBytes=" << snapshot.workBytesTotal
                         << ", waitMs(total/avg/max)=" << waitMsTotal << "/" << avgWaitMs << "/" << waitMsMax
                         << ", workMs(total/avg/max)=" << workMsTotal << "/" << avgWorkMs << "/" << workMsMax
                         << "\n" << flush;
        pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
    }
};
