#ifndef STAR_HOST_STARHOST_H
#define STAR_HOST_STARHOST_H

// STAR Suite host interface, API version 1 (STAR Suite 1.10).
//
// A host program links libstar_suite.a, calls runMain() instead of STAR's
// main(), and may run its own work beside STAR in the same process on STAR's
// shared permit pool (the External domain). The interface is generic: STAR
// knows nothing about what the host runs. See docs/HOST_API.md.
//
// With hooks == nullptr, runMain() is exactly the standalone STAR executable.

#include <cstdint>
#include <string>
#include <vector>

#include "PermitTypes.h"

namespace star {
namespace host {

constexpr int kApiVersion = 1;

// Read-only view of STAR's run configuration. STAR fills a fresh copy before
// each lifecycle callback; the host may keep the values it needs.
struct RunView {
    int runThreadN = 1;
    bool batchMode = false;
    std::string outFileNamePrefix;
    std::string outFileTmp;
    std::string commandLine;       // raw command line as STAR records it
    std::string commandLineFull;   // effective command line, host parameters included
    // Permit pool settings (STAR's --dynamicThread* and --variableThreads).
    bool permitInterface = false;  // --dynamicThreadInterface 1
    bool permitTelemetry = false;  // --dynamicThreadTelemetry 1
    bool fifoWaiters = false;      // --dynamicThreadFifoWaiters 1
    bool bgzfHierarchy = false;    // --dynamicThreadBgzfHierarchy 1
    bool variableThreads = false;  // --variableThreads 1
    int constMapPermits = 0;       // --dynamicThreadConstMapPermits (0 = pool size)
    std::string pfControllerMode;  // --dynamicThreadPfControllerMode
    int mapFloor = 0;              // --dynamicThreadMapFloor
    int featureFloor = 0;          // --dynamicThreadFeatureFloor
    uint64_t mapWorkEstimate = 0;      // --dynamicThreadMapWorkEstimate
    uint64_t featureWorkEstimate = 0;  // --dynamicThreadFeatureWorkEstimate
};

// Host callbacks. Every member is optional; an absent callback leaves STAR's
// standalone behaviour at that point. Callbacks returning bool report failure
// with false and, where an error string is offered, a message; STAR then
// prints its usual "EXITING because of ... ERROR: " prefix followed by the
// message (verbatim, so a trailing "\nSOLUTION: ...\n" is kept) and exits with
// the exit code given below. All callbacks run on STAR's main thread.
struct Hooks {
    void* ctx = nullptr;
    // Name of the External domain in STAR's permit log lines, e.g. "atac".
    const char* externalLabel = "external";

    // Parameter pass-through. Called for each parameter name STAR does not
    // know, from the command line (inputLevel 2) or a --parametersFiles file
    // (inputLevel >= 5), in input order; a name given in a file and again on
    // the command line is delivered twice and the later value wins. values
    // are the whitespace-separated tokens after the name (a double-quoted
    // token is one value). Return true to accept: STAR records the parameter
    // in Log.out and in its effective command line (BAM @PG, Log.out) and the
    // raw command line (BAM @CO) already contains it. Return false with an
    // empty error to let STAR report an unrecognized parameter as usual, or
    // with a message to report an invalid value. Parameters never reach the
    // host from STAR's defaults or from genomeParameters.txt.
    bool (*parameter)(void* ctx, const std::string& name,
                      const std::vector<std::string>& values, int inputLevel,
                      std::string* error) = nullptr;

    // Lifecycle. preflight and start run once, back to back, after the
    // genome, transcriptome and run-time inputs are loaded and before read
    // mapping starts (exit code 102 on preflight failure, 103 on start
    // failure). finish runs once after all of STAR's own work (mapping,
    // Solo, feature assignment, BAM sorting, wiggle output) and before the
    // permit exit invariant and Log.final.out (exit code 103 on failure).
    bool (*preflight)(void* ctx, const RunView& run, std::string* error) = nullptr;
    bool (*start)(void* ctx, const RunView& run, std::string* error) = nullptr;
    bool (*finish)(void* ctx, const RunView& run, std::string* error) = nullptr;

    // Permit plan, consulted when STAR sizes its pool as mapping starts.
    // extraPermitThreads: threads the host runs in the External domain; the
    //   pool is runThreadN + extra (unless --dynamicThreadConstMapPermits).
    // initialFloors: floors[0..2] (Map, Feature, External) arrive filled
    //   with STAR's configured floors (External 0); the host may change them.
    //   Called once when STAR configures the pool, also with the permit
    //   interface off; not called when STAR's BGZF hierarchy owns the floors.
    //   Failure exits with code 1, as STAR's other thread-setup failures.
    // externalActive: whether the host has External work in this run,
    //   answered from its parameters (it can be called before preflight);
    //   STAR then refuses routes that cannot share the pool (fused Flex).
    // requiresFullPoolAtExit: when true, the exit invariant also requires
    //   every configured permit to be available again.
    int (*extraPermitThreads)(void* ctx) = nullptr;
    bool (*initialFloors)(void* ctx, int configuredPermits, int floors[3],
                          std::string* error) = nullptr;
    bool (*externalActive)(void* ctx) = nullptr;
    bool (*requiresFullPoolAtExit)(void* ctx) = nullptr;
};

// Runs STAR with argv exactly as the STAR executable would. Returns STAR's
// exit status; fatal errors exit the process, as in STAR. One run per process.
int runMain(int argc, char** argv, const Hooks* hooks);

// Permit facade over STAR's pool, callable from any thread. The pool is
// configured as mapping starts, and only when --dynamicThreadInterface 1 or
// BGZF input enables it. Until permitsEnabled() is true, acquire returns 0
// at once and release is ignored, so a permit "acquired" then must not be
// released after the pool is enabled: wait for permitsEnabled() before the
// first acquire, or skip permits when the interface is off.
bool permitsEnabled();
// Blocks until a permit in `domain` is granted; returns the wait in ns.
uint64_t permitAcquire(Domain domain);
// Returns one permit with the work it covered (units, bytes, ns).
void permitRelease(Domain domain, uint64_t waitNs, uint64_t units,
                   uint64_t bytes, uint64_t workNs);
// Declares durable completion of a domain and releases its floor.
void permitMarkComplete(Domain domain);
// Replaces the borrowable floors, indexed Map, Feature, External.
void permitSetFloors(const int floors[3]);
PermitSnapshot permitSnapshot();

// Appends text to STAR's Log.out under STAR's log mutex and flushes.
// No-op before Log.out is open.
void logMain(const std::string& text);
// STAR's log timestamp, e.g. "Sep 28 12:00:00".
std::string timestamp();

}  // namespace host
}  // namespace star

#endif  // STAR_HOST_STARHOST_H
