// STAR Suite host interface: permit facade, log access and the run view.
// runMain itself is defined in STAR.cpp. See host/StarHost.h and
// docs/HOST_API.md.

#include "host/StarHost.h"
#include "host/StarHostInternal.h"

#include "GlobalVariables.h"
#include "Parameters.h"
#include "ThreadControl.h"
#include "TimeFunctions.h"

#include <atomic>
#include <cctype>
#include <pthread.h>

namespace star {
namespace host {

namespace {

std::atomic<bool> g_runMainEntered{false};
std::atomic<std::ofstream*> g_hostLogMain{nullptr};
std::string g_externalLabel = "external";

ThreadControl::PermitDomain toPermitDomain(Domain domain) {
    switch (domain) {
        case Domain::Feature: return ThreadControl::PermitDomain::FEATURE;
        case Domain::External: return ThreadControl::PermitDomain::EXTERNAL;
        case Domain::Map:
        default: return ThreadControl::PermitDomain::MAP;
    }
}

}  // namespace

namespace detail {

bool enterRunMain() {
    return !g_runMainEntered.exchange(true);
}

const char* externalLabel(const Hooks* hooks) {
    if (hooks != nullptr && hooks->externalLabel != nullptr && hooks->externalLabel[0] != '\0') {
        g_externalLabel = hooks->externalLabel;
    }
    return g_externalLabel.c_str();
}

void attachLog(std::ofstream* logMain) {
    g_hostLogMain.store(logMain);
}

RunView makeRunView(const Parameters& P) {
    RunView view;
    view.runThreadN = P.runThreadN;
    view.batchMode = P.batchMode;
    view.outFileNamePrefix = P.outFileNamePrefix;
    view.outFileTmp = P.outFileTmp;
    view.commandLine = P.commandLine;
    view.commandLineFull = P.commandLineFull;
    view.permitInterface = (P.dynamicThreadInterface == 1);
    view.permitTelemetry = (P.dynamicThreadTelemetry == 1);
    view.fifoWaiters = (P.dynamicThreadFifoWaiters == 1);
    view.bgzfHierarchy = (P.dynamicThreadBgzfHierarchy == 1);
    view.variableThreads = (P.variableThreads == 1);
    view.constMapPermits = P.dynamicThreadConstMapPermits;
    view.pfControllerMode = P.dynamicThreadPfControllerMode;
    view.mapFloor = P.dynamicThreadMapFloor;
    view.featureFloor = P.dynamicThreadFeatureFloor;
    view.mapWorkEstimate = P.dynamicThreadMapWorkEstimate;
    view.featureWorkEstimate = P.dynamicThreadFeatureWorkEstimate;
    return view;
}

std::string hostErrorText(const std::string& error, const char* fallback) {
    return error.empty() ? std::string(fallback) : error;
}

std::string upperCopy(const std::string& text) {
    std::string out(text);
    for (char& c : out) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return out;
}

}  // namespace detail

bool permitsEnabled() {
    return g_threadChunks.mapPermitEnabled();
}

uint64_t permitAcquire(Domain domain) {
    return g_threadChunks.mapPermitAcquireForDomain(toPermitDomain(domain));
}

void permitRelease(Domain domain, uint64_t waitNs, uint64_t units, uint64_t bytes, uint64_t workNs) {
    g_threadChunks.mapPermitReleaseForDomain(toPermitDomain(domain), waitNs, units, bytes, workNs);
}

void permitMarkComplete(Domain domain) {
    g_threadChunks.mapPermitMarkDomainComplete(toPermitDomain(domain));
}

void permitSetFloors(const int floors[3]) {
    g_threadChunks.mapPermitConfigureDomainFloors(std::vector<int>(floors, floors + 3));
}

PermitSnapshot permitSnapshot() {
    return g_threadChunks.mapPermitSnapshot();
}

void logMain(const std::string& text) {
    std::ofstream* log = g_hostLogMain.load();
    if (log == nullptr) {
        return;
    }
    pthread_mutex_lock(&g_threadChunks.mutexLogMain);
    *log << text;
    log->flush();
    pthread_mutex_unlock(&g_threadChunks.mutexLogMain);
}

std::string timestamp() {
    return timeMonthDayTime();
}

}  // namespace host
}  // namespace star
