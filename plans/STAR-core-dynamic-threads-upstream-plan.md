# STAR Core Dynamic Threads Upstream Plan

Date: 2026-02-17  
Branch target: `core-dynamic-threads`  
Status: staged design for upstream-friendly implementation

Iteration runbook: `plans/dynamic_scheduler_runbook.md`

## CRISPR Calling Threshold Policy (2026-02-18)
1. Default `crMinUmi` is reset to `3` for general CRISPR feature calling.
2. A375 parity runs are pinned to `--crMinUmi 10` via:
   `tests/run_a375_gex_features_cr_parity.sh` and
   `tests/run_a375_gex_features_cr_parity_genefull.sh`
   (`A375_CR_MIN_UMI`, default `10`).
3. GMM library now short-circuits trivial EM cases:
   no counts above min-UMI floor or single non-zero cell.

## PF Preload Phase (2026-02-18)
1. Added asynchronous pf-multi preload API:
   `startPfMultiConfigPreload(const Parameters&)`.
2. STAR now starts preload during initialization (after parameter parse/filter).
3. `processPfMultiConfig(...)` consumes preloaded prep artifacts with safe
   synchronous fallback on preload error.
4. Preloaded artifacts:
   parsed multi config, chemistry resolution, CR assign output dirs,
   filtered feature refs, and anchored map/feature read estimates.
5. Current gate:
   preload runs when `--pfMultiConfig` is set and `--dynamicThreadInterface 1`.
6. Smoke validation:
   `tests/test_cr_compat_crispr_calling.sh` with
   `--readMapNumber 100 --dynamicThreadInterface 1`.
   Artifact: `/tmp/pf_preload_smoke_20260218_190530_696476/Log.out`
   (`pf-preload: consumed async preparation (sec=6.92337, featureLibraries=1)`).

## Objective
Build a dynamic threading controller in STAR core that can be reused by other
modules now, then evolve into a producer/consumer execution model later.

## Upstream Constraints
1. Default behavior must remain byte-for-byte compatible unless an explicit new
   flag is enabled.
2. Existing STAR thread creation model stays intact in initial phases.
3. New logic must be isolated and minimally invasive in `core/legacy`.
4. Metrics and control logic must be optional and low overhead when disabled.

## Scope
1. Near term: dynamic active-thread permits over existing worker pools.
2. Mid term: module-facing permit API so feature modules can participate.
3. Long term: staged migration to producer/consumer execution.

## Non-Goals (Initial Upstream Slice)
1. No full chunk-manager rewrite in the first merged phases.
2. No required behavior change for users not enabling dynamic mode.
3. No mandatory dependency on process_features internals.

## Core Architecture
1. Introduce a central `ThreadPermitController` in STAR core.
2. Manage a global permit budget:
   `total_permits = --runThreadN` by default.
3. Expose workload-scoped permits:
   `map`, `feature`, and reserved buckets for future modules.
4. Keep workers fixed at launch; only gate active work acquisition.
5. Add explicit telemetry for wait time, active time, and queue pressure.

## API Direction (Integration-Ready)
1. Internal C++ API for STAR core:
   `acquire(workload, timeout)`, `release(workload)`, `snapshot()`.
2. Thin C ABI bridge for C modules (process_features) in a later phase.
3. No module may directly mutate permit counts; only controller updates split.
4. Controller updates must be bounded by min/max floors and cooldown windows.

## Interface Stability for Producer/Consumer Migration
1. Keep a stable module-facing contract in all phases:
   `acquire`, `release`, `snapshot`, and optional `submit`/`complete`.
2. Treat scheduler model as internal:
   fixed-thread loop now, producer/consumer queue later.
3. Keep CLI behavior and flags consistent across both scheduler modes.
4. Add an adapter layer in core so existing paths can call the same contract.
5. Preserve ordering semantics as explicit policy, not implicit thread behavior.

## Expected Interface Delta if STAR Goes Full Producer/Consumer
1. Low change for module callers if Stage 2 contract is adopted first.
2. Internal STAR core change is high (loop and queue ownership rewiring).
3. Unavoidable additions likely include:
   work-item identity, completion state, and cancellation/flush hooks.
4. Permit controller API should stay compatible if queue operations are wrapped.

## Staged Implementation Plan

## Stage 0: Baseline Telemetry and Guardrails
1. Add timing and queue-depth telemetry in existing paths only.
2. Add no-op controller scaffolding compiled but disabled by default.
3. Validate default-mode zero diff on representative datasets.
4. Exit criteria:
   telemetry visible, no functional output changes, no measurable slowdown.

## Stage 1: STAR-Core Permit Gating (Static Split)
1. Add flags:
   `--dynamicThreadPermits`, `--dynamicThreadPermitsMapMin`,
   `--dynamicThreadPermitsFeatureMin`, `--dynamicThreadPermitsMode`,
   `--variableThreads`.
2. Gate STAR chunk processing in
   `core/legacy/source/ReadAlignChunk_processChunks.cpp`.
3. Ensure permit release on all success/failure exits.
4. Keep split static at first (`map=runThreadN`, `feature=0` unless overridden).
5. Exit criteria:
   determinism holds, no deadlocks, no permit leaks.

## Stage 1a: Interface Smoke + Mock Consumer (Repro Harness)
1. Add tiny fixture smoke harness:
   `tests/run_dynamic_threads_tiny_fixture.sh`.
2. Add a mock consumer report generator that parses telemetry from `Log.out`:
   `tests/dynamic_threads/mock_consumer_report.py`.
3. Add parser unit test:
   `tests/dynamic_threads/test_mock_consumer_report.py`.
4. Add reproducibility runbook with fixed fixture paths and commands:
   `docs/DYNAMIC_THREADS_TINY_FIXTURE_RUNBOOK_20260217.md`.
5. Exit criteria:
   smoke passes on tiny fixture and emits JSON/TXT telemetry reports.
6. Variable sequence smoke targets:
   `3->2->4` and `1->2->1` permit transitions.
7. Retune trace requirement:
   mock consumer validates exact cyclic `retuneTrace` behavior, not only
   aggregate `retunes` count.

## Stage 2: Module Integration Surface
1. Add stable controller access points in `core/legacy` for non-core modules.
2. Provide optional C ABI wrappers for `process_features` integration.
3. Add integration tests for controller lifecycle and ABI safety.
4. Exit criteria:
   module can read metrics and acquire/release permits without race regressions.
5. Deliverable:
   freeze and document Interface v1 so later producer/consumer work is
   implementation-internal for module consumers.

## Stage 3: Feature Consumer Permit Gating
1. Integrate permit gating into
   `core/features/process_features/src/assignBarcodes.c`.
2. Preserve existing queue semantics and ordering.
3. Add conservative fallback split if no telemetry is available.
4. Exit criteria:
   parity in strict mode, expected throughput gain in mixed workloads.

## Stage 4: Dynamic Rebalancer
1. Add periodic controller loop that rebalances map/feature permits.
2. Inputs:
   service rates, recent permit wait time, and estimated remaining work.
3. Work sizing policy:
   anchored estimate: exact read count from the smallest primary input file
   (`R1` for map, barcode file for PF), then scale by total compressed bytes
   across the domain; monotonic inflation if remaining estimate goes negative.
4. Stability controls:
   hysteresis, min hold time, bounded step size, starvation cap.
5. Retune policy:
   move one permit per interval toward the side with larger ETA when ETA gap
   exceeds 10%.
6. Exit criteria:
   no oscillatory behavior and no sustained starvation in stress tests.

## Stage 5: Producer/Consumer Readiness Refactor
1. Introduce internal work-item abstractions while preserving current semantics.
2. Move chunk lifecycle metadata behind explicit interfaces.
3. Keep old execution path selectable for regression safety.
4. Exit criteria:
   refactor-only parity with no scheduler behavior change required.

## Stage 6: Producer/Consumer Pilot (Flagged)
1. Add opt-in producer/consumer execution mode for STAR core.
2. Reuse permit controller as backpressure and fairness mechanism.
3. Benchmark against Stage 4 mode for throughput and determinism.
4. Exit criteria:
   pilot mode stable enough for extended validation, not default.

## PR Slices for Upstream Review
1. PR1: Stage 0 telemetry + controller scaffolding + default-off flags.
2. PR2: Stage 1 STAR-only permit gating (static).
3. PR3: Stage 2 module API and ABI bridge.
4. PR4: Stage 3 feature-side permit gating.
5. PR5: Stage 4 dynamic rebalancer.
6. PR6: Stage 5 readiness refactor.
7. PR7: Stage 6 producer/consumer pilot (experimental).

## Validation Matrix
1. Functional invariants:
   `PairedKeepInputOrder`, `readFilesIndex` transitions, chunk counters, and
   output determinism.
2. Regression suites:
   `tests/run_cbub_regression_test.sh`,
   `tests/test_cr_compat_crispr_calling.sh`,
   perturb compatibility datasets (A375, UCSF contexts).
3. Performance metrics:
   wallclock, reads/min, CPU util, iowait, queue occupancy, permit history.
4. Stress/failure cases:
   end-of-input partial chunks, FILE marker transitions, low/high feature load,
   starvation windows, error-path release guarantees.

## Risk Register
1. Hidden ordering regressions under contention.
2. Permit leaks on uncommon error paths.
3. Rebalancer oscillation under bursty queue behavior.
4. ABI drift between core controller and C modules.

## TODO Backlog (Cross-Cutting)
1. Consolidate CB correction to a single implementation across Solo/PF-facing
   paths; remove dual-system behavior over time.
2. Candidate direction: keep exact-match hash fast path, move 1MM and N rescue
   to on-demand evaluation instead of global precomputed 1MM variant map.
3. Add migration flagging:
   new unified path default-on after validation, temporary legacy fallback
   flag for rollback.
4. Add acceptance gates before removing legacy path:
   output parity on perturb fixtures, deterministic rescue behavior, startup
   time reduction, and no regression in steady-state throughput.
5. Track this work as independent from dynamic-permit scheduling so it can land
   in upstream-friendly PR slices without blocking permit integration.

## CB Init Benchmark Snapshot (2026-02-18)
Artifact root:
`/tmp/cb_init_bench_20260218_174640`

Setup:
1. Dataset: UCSF GEX lane (`/storage/ucsf-2M/GEX/iPSC2_1_AALG2/...L001...`).
2. Whitelist: 3M (`3M-february-2018_TRU.txt`).
3. Probe run: `--readMapNumber 1`, `--runThreadN 1` (startup-focused).

Results:
1. `default`:
   - `elapsed_sec=83.98`
   - `maxrss_gb=36.11`
2. `STAR_BUILD_LEGACY_AMBIG_HASH=1`:
   - `elapsed_sec=115.90` (`+31.92s`, `+38.0%` vs default)
   - `maxrss_gb=40.10` (`+3.99 GB`)
3. `STAR_CBCORRECTOR_MAX_HAMMING=0`:
   - `elapsed_sec=9.96` (`-74.02s`, `-88.1%` vs default)
   - `maxrss_gb=28.44` (`-7.67 GB`)
4. `STAR_SKIP_CBCORRECTOR_INIT=1`:
   - `elapsed_sec=9.30` (`-74.68s`, `-88.9%` vs default)
   - `maxrss_gb=28.08` (`-8.03 GB`)

Note:
1. These env toggles were temporary benchmark instrumentation and were removed
   from the runtime code path after measurement.

Interpretation:
1. Dominant startup cost is global 1MM precompute in `CbCorrector`.
2. Legacy `ambiguousCbByKey` precompute adds further startup+memory cost.
3. Exact-only init (`maxHamming=0`) is close to no-`CbCorrector` init, which
   supports migration to on-demand/lazy 1MM correction with cache.

## Shared CB Correction Library Direction (Decision)
Decision:
1. Build a shared CB correction library and use it from both Solo and PF.
2. Keep one correction semantics surface and one optimization path.
3. Decouple this from dynamic-thread scheduling rollout so work can land in
   parallel.

Library scope (v1):
1. Exact whitelist lookup.
2. On-demand 1MM candidate generation and resolution (no global 1MM prebuild).
3. N-base expansion with configurable cap.
4. Optional bounded variant cache (LRU or shard-local cache).
5. Deterministic ambiguity handling hooks (including frozen-count mode).

Integration model:
1. Solo adapter: replace current `CbCorrector` hot path with shared library
   calls while preserving existing CLI behavior.
2. PF adapter: replace local barcode-correction routines with shared library
   calls while preserving assignBarcodes/pf_api interfaces.
3. Keep legacy fallbacks behind temporary flags until parity gates are met.

Phased delivery:
1. Phase A:
   introduce shared library with Solo-only integration (startup reduction win).
2. Phase B:
   integrate PF adapter and align ambiguity/posterior behavior.
3. Phase C:
   remove duplicated correction implementations after parity + perf gates.

Acceptance gates:
1. Output parity on A375 and UCSF fixtures in default modes.
2. Deterministic rescue behavior under multithreaded schedules.
3. Startup and memory improvements on 3M whitelist runs.
4. No steady-state throughput regression in Solo or PF.

## Success Criteria
1. Default mode remains unchanged.
2. Dynamic mode yields reproducible throughput improvements on perturb-seq
   workloads.
3. Architecture is upstream-acceptable and extensible to producer/consumer
   mode without another deep rewrite.
