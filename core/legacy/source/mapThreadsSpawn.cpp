#include "mapThreadsSpawn.h"
#include "ThreadControl.h"
#include "GlobalVariables.h"
#include "ErrorWarning.h"
#include <algorithm>
#include <sstream>

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

void mapThreadsSpawn (Parameters &P, ReadAlignChunk** RAchunk) {
    const bool interfaceEnabled = (P.dynamicThreadInterface == 1);
    const bool telemetryEnabled = (P.dynamicThreadTelemetry == 1);
    const bool variableThreadsEnabled = (P.variableThreads == 1);
    const int configuredPermits = (P.dynamicThreadConstMapPermits > 0)
        ? std::min(P.dynamicThreadConstMapPermits, P.runThreadN)
        : P.runThreadN;
    g_threadChunks.mapPermitConfigure(interfaceEnabled, P.runThreadN, configuredPermits, telemetryEnabled, variableThreadsEnabled);
    g_threadChunks.mapPermitConfigureRetunePlan(P.variableThreadsPermitSequence, P.variableThreadsRetuneEveryAcquires);

    if (interfaceEnabled) {
        pthread_mutex_lock(&g_threadChunks.mutexLogMain);
        P.inOut->logMain << "Dynamic thread interface enabled: map permits=" << configuredPermits
                         << " (runThreadN=" << P.runThreadN
                         << ", telemetry=" << (telemetryEnabled ? "on" : "off")
                         << ", variableThreads=" << (variableThreadsEnabled ? "on" : "off")
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
