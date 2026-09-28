// Test host for the STAR Suite host interface (core/legacy/source/host).
// Links libstar_suite.a and runs STAR through star::host::runMain. The mode
// comes from STAR_HOST_TEST_MODE so that argv is passed to STAR unchanged:
//   null      runMain(argc, argv, nullptr)
//   empty     a Hooks object with no callbacks (default)
//   external  a dummy external domain: two threads take and return permits
//             in the External domain while STAR maps, joined at finish
//   params    parameter pass-through for names starting with "hostTest"
//   fail-preflight | fail-start | fail-finish   the callback fails
// Events are written to <outFileNamePrefix>host_events.tsv at finish (and at
// preflight for the failure modes, which never reach finish).
#include "StarHost.h"
#include "SaturationPermitController.h"

static_assert(star::host::kApiVersion == 1, "host test expects host API version 1");

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace {

struct TestHost {
    std::string mode;
    std::mutex eventsMutex;
    std::vector<std::string> events;
    std::string outPrefix;
    std::vector<std::thread> workers;
    std::atomic<bool> stop{false};
    std::atomic<uint64_t> acquired{0};
    int iterations = 200;
    bool permitInterface = false;

    void event(const std::string& text) {
        std::lock_guard<std::mutex> lock(eventsMutex);
        events.push_back(text);
    }

    void writeEvents() {
        std::lock_guard<std::mutex> lock(eventsMutex);
        std::ofstream out((outPrefix + "host_events.tsv").c_str());
        for (const auto& e : events) out << e << "\n";
    }
};

std::string join(const std::vector<std::string>& values) {
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) out << (i ? "|" : "") << values[i];
    return out.str();
}

bool onParameter(void* ctx, const std::string& name, const std::vector<std::string>& values,
                 int inputLevel, std::string* error) {
    TestHost* host = static_cast<TestHost*>(ctx);
    if (name.compare(0, 8, "hostTest") != 0) {
        return false;  // not ours: STAR reports an unrecognized parameter
    }
    if (name == "hostTestBad") {
        *error = "hostTestBad is always rejected by the test host";
        return false;
    }
    host->event("parameter\t" + name + "\t" + join(values) + "\t" + std::to_string(inputLevel));
    return true;
}

bool onPreflight(void* ctx, const star::host::RunView& run, std::string* error) {
    TestHost* host = static_cast<TestHost*>(ctx);
    host->outPrefix = run.outFileNamePrefix;
    host->permitInterface = run.permitInterface;
    std::ostringstream line;
    line << "preflight\trunThreadN=" << run.runThreadN << "\tpermitInterface=" << run.permitInterface
         << "\ttelemetry=" << run.permitTelemetry << "\tbatchMode=" << run.batchMode
         << "\tcommandLineFullHasHost=" << (run.commandLineFull.find("--hostTest") != std::string::npos);
    host->event(line.str());
    star::host::logMain("[host test] preflight at " + star::host::timestamp() + "\n");
    if (host->mode == "fail-preflight") {
        host->writeEvents();
        *error = "host test preflight failure\nSOLUTION: this failure is intentional\n";
        return false;
    }
    if (host->mode == "external" && !run.permitInterface) {
        *error = "external mode needs --dynamicThreadInterface 1";
        return false;
    }
    return true;
}

void externalWorker(TestHost* host, int id) {
    // Never acquire before the pool is enabled (see StarHost.h).
    while (!star::host::permitsEnabled() && !host->stop.load()) {
        std::this_thread::sleep_for(std::chrono::microseconds(200));
    }
    for (int i = 0; i < host->iterations && star::host::permitsEnabled(); ++i) {
        const uint64_t waitNs = star::host::permitAcquire(star::host::Domain::External);
        const auto start = std::chrono::steady_clock::now();
        volatile uint64_t sink = 0;
        for (int k = 0; k < 2000; ++k) sink += static_cast<uint64_t>(k * (id + 1));
        const uint64_t workNs = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - start).count());
        star::host::permitRelease(star::host::Domain::External, waitNs, 1, 64, workNs);
        host->acquired.fetch_add(1);
    }
}

bool onStart(void* ctx, const star::host::RunView& run, std::string* error) {
    TestHost* host = static_cast<TestHost*>(ctx);
    host->event("start\toutFileNamePrefix=" + run.outFileNamePrefix);
    if (host->mode == "fail-start") {
        host->writeEvents();
        *error = "host test start failure\n";
        return false;
    }
    if (host->mode == "external") {
        for (int id = 0; id < 2; ++id) host->workers.emplace_back(externalWorker, host, id);
    }
    return true;
}

int onExtraPermitThreads(void* ctx) {
    TestHost* host = static_cast<TestHost*>(ctx);
    host->event("extraPermitThreads\t2");
    return 2;
}

bool onInitialFloors(void* ctx, int configuredPermits, int floors[3], std::string*) {
    TestHost* host = static_cast<TestHost*>(ctx);
    std::ostringstream line;
    line << "initialFloors\tconfigured=" << configuredPermits << "\tin=" << floors[0] << "/" << floors[1]
         << "/" << floors[2];
    floors[2] = 1;
    line << "\tout=" << floors[0] << "/" << floors[1] << "/" << floors[2];
    host->event(line.str());
    return true;
}

bool onExternalActive(void* ctx) {
    return static_cast<TestHost*>(ctx)->mode == "external";
}

bool onRequiresFullPoolAtExit(void* ctx) {
    return static_cast<TestHost*>(ctx)->mode == "external";
}

bool onFinish(void* ctx, const star::host::RunView&, std::string* error) {
    TestHost* host = static_cast<TestHost*>(ctx);
    host->stop.store(true);
    for (auto& worker : host->workers) worker.join();
    host->workers.clear();
    if (host->mode == "external" && star::host::permitsEnabled()) {
        star::host::permitMarkComplete(star::host::Domain::External);
    }
    const star::host::PermitSnapshot snap = star::host::permitSnapshot();
    std::ostringstream line;
    line << "finish\tacquired=" << host->acquired.load() << "\texternalAcquireCalls="
         << snap.externalDomain.acquireCalls << "\texternalInUse=" << snap.externalDomain.inUse
         << "\tconfigured=" << snap.configuredPermits << "\tavailable=" << snap.availablePermits
         << "\tenabled=" << snap.enabled;
    host->event(line.str());
    host->event(std::string("controllerLabels\t") +
                star::permits::SaturationPermitController::phaseName(
                    star::permits::SaturationPermitController::Phase::PROBE_EXTERNAL, "hosttest") + "\t" +
                star::permits::SaturationPermitController::domainName(
                    star::permits::SaturationPermitController::Domain::EXTERNAL));
    star::host::logMain("[host test] finish\n");
    host->writeEvents();
    if (host->mode == "fail-finish") {
        *error = "host test finish failure\n";
        return false;
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    const char* modeEnv = std::getenv("STAR_HOST_TEST_MODE");
    TestHost host;
    host.mode = modeEnv != nullptr ? modeEnv : "empty";
    if (host.mode == "null") {
        return star::host::runMain(argc, argv, nullptr);
    }
    star::host::Hooks hooks;
    hooks.ctx = &host;
    if (host.mode != "empty") {
        hooks.externalLabel = "hosttest";
        hooks.parameter = onParameter;
        hooks.preflight = onPreflight;
        hooks.start = onStart;
        hooks.finish = onFinish;
        hooks.externalActive = onExternalActive;
        hooks.requiresFullPoolAtExit = onRequiresFullPoolAtExit;
        if (host.mode == "external") {
            hooks.extraPermitThreads = onExtraPermitThreads;
            hooks.initialFloors = onInitialFloors;
        }
    }
    const int rc = star::host::runMain(argc, argv, &hooks);
    host.stop.store(true);
    for (auto& worker : host.workers) worker.join();
    return rc;
}
