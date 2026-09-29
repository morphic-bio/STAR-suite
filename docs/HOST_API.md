# STAR Suite Host Interface

STAR Suite 1.10 can be linked into another program as a library. The host
program calls `star::host::runMain()` in place of STAR's `main()` and can run
its own work beside STAR in the same process, sharing STAR's thread-permit
pool. The interface is generic: STAR knows nothing about what the host runs.
Multiomics Suite uses it to build the joint RNA + ATAC (multiome) binary.

Host API version: **1** (`star::host::kApiVersion`).

## Build and link

```bash
make star-host-lib                   # STAR's bundled HTSlib
make star-host-lib HTSLIB=external   # installed HTSlib (pkg-config htslib)
```

This writes, under `core/legacy/source/`:

- `libstar_suite.a`: every STAR object except `STARmain.o`, the STAR
  executable's `main()`.
- `libstar_suite.link`: one line with what must follow the archive on the
  host's link line: STAR's dependent archives (`libflex`, `libtrim`, `libem`,
  `libprocess_features_integrated`, `libscrna`), HTSlib and system libraries.

Public headers are in `core/legacy/source/host/` and need no other STAR
header: `StarHost.h` (with `PermitTypes.h`) and `SaturationPermitController.h`.

```bash
S=/path/to/STAR-suite/core/legacy/source
g++ -std=c++11 -I"$S/host" host.cpp "$S/libstar_suite.a" $(cat "$S/libstar_suite.link") -o host
```

Use `HTSLIB=external` when the host also links code built against an installed
HTSlib, so the executable has one HTSlib ABI. The `STAR` executable itself is
built by `make core` exactly as before; `STARmain.cpp` is
`return star::host::runMain(argc, argv, nullptr);`.

## Entry point

```cpp
int star::host::runMain(int argc, char** argv, const star::host::Hooks* hooks);
```

`argv` is STAR's command line. With `hooks == nullptr`, or with a `Hooks`
object whose callbacks are all null, the run is the standalone STAR run
(checked byte for byte by `tests/host_api/`). STAR still exits the process on
fatal errors, and `runMain` may be called once per process because STAR keeps
global state. `--version`, `--help` and the other informational options behave
as in STAR, so a host handles its own informational options before calling
`runMain`.

## Callbacks

All members of `Hooks` are optional and run on STAR's main thread. `ctx` is
passed back unchanged.

| Callback | When STAR calls it | Standalone equivalent |
|---|---|---|
| `parameter(name, values, inputLevel, error)` | For each parameter name STAR does not know, from the command line or a `--parametersFiles` file | STAR rejects the name |
| `preflight(run, error)` | Once, after the genome, transcriptome and run-time inputs are loaded, before read mapping | none |
| `start(run, error)` | Right after `preflight` | none |
| `extraPermitThreads()` | When the permit pool is sized, as mapping starts | 0 |
| `initialFloors(configuredPermits, floors[3], error)` | Once when the pool is configured (not when the BGZF hierarchy owns the floors) | floors unchanged |
| `externalActive()` | When STAR decides whether the fused Flex route may run (also before `preflight`) | false |
| `requiresFullPoolAtExit()` | In the permit exit invariant | false |
| `finish(run, error)` | After all STAR work (mapping, Solo, feature assignment, BAM sorting, wiggle), before the permit exit invariant and `Log.final.out` | none |

`RunView` is a read-only copy of the run configuration the host needs:
`runThreadN`, `batchMode`, output prefix and temporary directory, the raw and
effective command lines, and STAR's own permit settings (interface, telemetry,
FIFO, BGZF hierarchy, variable threads, constant permits, PF controller mode,
MAP and FEATURE floors and work estimates).

A callback that returns false fails the run. STAR prints its usual prefix and
the host's message verbatim, so a message ending in `"\nSOLUTION: ...\n"` keeps
that form:

| Callback | Message prefix | Exit code |
|---|---|---|
| `preflight` | `EXITING because of fatal ERROR: ` | 102 (`EXIT_CODE_PARAMETER`) |
| `start`, `finish` | `EXITING because of fatal ERROR: ` | 103 (`EXIT_CODE_RUNTIME`) |
| `initialFloors` | `EXITING because of FATAL ERROR: ` | 1, as STAR's other thread-setup failures |

## Parameters

A host parameter uses STAR's syntax (`--name value ...`) on the same command
line or in a `--parametersFiles` file. STAR delivers `values` as the
whitespace-separated tokens after the name, a double-quoted token being one
value; a scalar parameter uses `values[0]` and ignores the rest, as STAR does.
`inputLevel` is STAR's input level: 2 for the command line, 5 and above for
parameter files. A name given in a file and again on the command line is
delivered twice, file first, and the later value wins; a name repeated within
one input is a duplicate error, as for STAR's own parameters. A parameter the
host never receives keeps the host's default, so `inputLevel > 0` marks a value
the user set.

Return true to accept. STAR then logs the parameter in `Log.out` like its own
and appends it, in first-seen order, to the effective command line (`Log.out`
"Final effective command line" and BAM `@PG CL`). The raw command line (BAM
`@CO`, `Log.out`) contains it anyway. Return false with an empty error for
"not mine" (STAR reports an unrecognized parameter as usual) or with a message
for an invalid value (exit 102). Parameters never reach the host from STAR's
defaults or from `genomeParameters.txt`.

## Shared permits

STAR's pool has three domains: `Map` and `Feature` (STAR's mapping and feature
assignment) and `External`, which STAR never uses and lends to the host.

- The pool is `runThreadN + extraPermitThreads()` permits, unless
  `--dynamicThreadConstMapPermits` fixes it.
- `permitAcquire(Domain::External)` blocks until a permit is granted and
  returns the wait in ns; `permitRelease(Domain::External, waitNs, units,
  bytes, workNs)` returns it with the work it covered. Hold at most one permit
  per worker thread and release it between batches.
- The pool exists only when `--dynamicThreadInterface 1` (or BGZF input)
  enables it, from the start of mapping. Until `permitsEnabled()` is true,
  acquire returns 0 at once and release is ignored, so a permit taken then must
  not be released later: wait for `permitsEnabled()` before the first acquire,
  or run without permits when `RunView::permitInterface` is false.
- `permitMarkComplete(Domain::External)` declares the host's work finished and
  releases its floor. `permitSetFloors(floors)` replaces the borrowable floors
  (Map, Feature, External); `permitSnapshot()` returns counters for all three
  domains (`externalDomain` for the host's).
- At exit STAR requires every domain's in-use count and the waiter count to be
  zero; with `requiresFullPoolAtExit()` true it also requires all configured
  permits to be available.

`SaturationPermitController.h` (namespace `star::permits`) is the saturation
policy STAR uses for MAP/FEATURE balancing; a host may drive all three domains
with it. Its `EXTERNAL` domain, phases and reasons print as "external"
(`probe-external`, `external-complete`); the overloads taking a label print the
host's name instead.

## Logging

- `externalLabel` names the External domain in STAR's permit log lines
  (`floors(map/feature/<label>)`, `<label>State(...)`, `<label>InUse=` in the
  final invariant, stall warnings). The default is `external`; a multiome host
  passes `atac`.
- `logMain(text)` appends to `Log.out` under STAR's log mutex and flushes; it
  is a no-op before `Log.out` is open. `timestamp()` is STAR's log timestamp.

## Environment

- `genomeGenerate` looks for its `compute_expected_gc` helper on `PATH`, then
  relative to the executable's directory in a STAR-suite tree. A host
  executable elsewhere should put the helper on `PATH`, or the TranscriptVB
  expected-GC sidecar is skipped with a warning.
- Output files are those of the same STAR run. Only the lines that record the
  command line (the header of `genomeParameters.txt`, BAM `@PG`/`@CO`, the
  `Chimeric.out.junction` header) and the logs show the host's `argv[0]` and
  parameters.

## Tests

`make host-api-tests` builds STAR and the library and runs
`tests/host_api/run_host_api_tests.sh`: a no-callback host reproduces STAR on
the synthetic scRNA fixture (eight UMI modes and the genome index); a host
with a dummy External domain (two threads, 400 permits) leaves Solo outputs
identical and the exit invariant clean; unknown parameters are rejected
without a host and recorded with one; failing callbacks exit with the codes
above; and the permit unit tests pass against the library. CI runs it as the
`host-api` partial build.

## From the 1.9.5 Chromap orchestration

The Chromap integration left STAR in 1.10.0 (see
`docs/HANDOVER_MULTIOMICS_1.10.md`). Its STAR dependencies map one to one:

| 1.9.5 orchestration | Host interface |
|---|---|
| `preflightStarChromapAtacIfEnabled`, `startStarChromapAtacIfEnabled` (`STAR.cpp`) | `preflight`, `start` |
| `runStarChromapAtacIfEnabled` (`STAR.cpp`) | `finish` |
| `P.chromapAtac.threads` in the pool size | `extraPermitThreads` |
| mode-2 initial floors, `dynamicThreadAtacFloor` (`mapThreadsSpawn.cpp`) | `initialFloors` |
| Flex fused-route refusal when ATAC is on | `externalActive` |
| `dynamicThreadAtacController == 2` full-pool exit check | `requiresFullPoolAtExit` |
| `g_threadChunks.mapPermitAcquireForDomain(ATAC)` / `…ReleaseForDomain` | `permitAcquire` / `permitRelease` (`Domain::External`) |
| `…MarkDomainComplete`, `…ConfigureDomainFloors`, `…Enabled`, `…Snapshot` | `permitMarkComplete`, `permitSetFloors`, `permitsEnabled`, `permitSnapshot` |
| `snapshot.atacDomain` | `snapshot.externalDomain` |
| `P.inOut->logMain` under `mutexLogMain`; `timeMonthDayTime()` | `logMain`; `timestamp` |
| `P.runThreadN`, `P.outFileNamePrefix`, `P.dynamicThread*` (STAR's own) | `RunView` |
| the 57 `chromapAtac*`, `multiomeAtac*` and ATAC-only permit parameters; `parameterInputLevel()` | the host's own options through `parameter`, with `inputLevel` |
| `star::multiome::SaturationPermitController` (`ATAC`, `atac*`) | `star::permits::SaturationPermitController` (`EXTERNAL`, `external*`) |
| Parameter checks in `Parameters.cpp` (ATAC floor and controller, BGZF hierarchy with ATAC, mode 2 with an applying PF controller) | host `preflight`, reading `RunView` |
