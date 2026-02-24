# STAR Core Dynamic Threads Upstream Plan

Date: 2026-02-17  
Branch target: `core-dynamic-threads`  
Status: staged design for upstream-friendly implementation

Iteration runbook: `plans/dynamic_scheduler_runbook.md`

## C/C++ LOC Delta Snapshot (2026-02-23)
1. Counter script: `core/legacy/scripts/count_cpp_lines.sh`.
2. Current STAR core scope (`core/legacy/source` + `core/features`,
   excluding only `htslib/` and `opal/`): `444` files, `108,813` lines.
3. Included integrated library example:
   `core/features/libscrna` = `24` files, `7,236` lines.
4. Upstream STAR `2.7.11b` (`source/` with the same exclusions):
   `250` files, `28,228` lines.
5. Net delta vs upstream: `+194` files, `+80,585` lines.

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

## UCSF Timing Update (2026-02-18, CR9 32-core rerun)
Dataset/scope:
1. UCSF iPSC2_1_AALG2 comparison slice used in the existing STAR dynamic
   benchmark (`/storage/ucsf-2M/...` inputs).
2. CR9 reference rerun kept the same sample/feature/transcriptome inputs and
   changed scheduler resources to `--localcores 32`.
3. Consolidated report:
   `comparisons/ucsf_dynamic_threads_vs_cr9_20260218/RESULTS.md`

Timing summary:

| Tool/Mode | Threads | Wall time | Peak RSS | Artifact |
|---|---:|---:|---:|---|
| STAR `dynamic_on` | 32 | `1:56.30` | `42,077,388 KB` | `/storage/ucsf-2M/star_runs/compare_pfopt_dynamic_vs_baseline_20260218_192106/PARITY_REPORT.txt` |
| STAR `dynamic_off` | 32 | `1:56.30` | `42,076,712 KB` | `/storage/ucsf-2M/star_runs/compare_pfopt_dynamic_vs_baseline_20260218_192106/PARITY_REPORT.txt` |
| Cell Ranger 9 (`crstar_sameidx`, prior) | 16 | `real 118.07s` | `13.616 GB` | `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813.log` |
| Cell Ranger 9 (`crstar_sameidx`, rerun) | 32 | `1:57.83` (`117.83s`) | `13.618 GB` | `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260218_32c.log` |

CR9 16-core vs 32-core reproducibility check:
1. CRISPR call outputs are byte-identical:
   `protospacer_calls_per_cell.csv`, `protospacer_calls_summary.csv`,
   `protospacer_umi_thresholds.csv`, `protospacer_umi_thresholds.json`.
2. Filtered MEX gzip files have different compressed MD5 values but identical
   decompressed contents (gzip header/metadata variance only).

32-core CR9 run provenance:
1. Run id:
   `cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260218_32c`
2. Resources:
   `--localcores 32 --localmem 128`
3. Martian fork walltime:
   `116.712882166s` from
   `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260218_32c/_perf`
4. Repro command:
   ```bash
   /usr/bin/time -v cellranger count \
     --id cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260218_32c \
     --transcriptome /storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar \
     --libraries /storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260218_32c_libraries.csv \
     --feature-ref /mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv \
     --include-introns true --min-crispr-umi 3 --create-bam false \
     --localcores 32 --localmem 128 --disable-ui
   ```

## UCSF 2M Fixture Refresh (2026-02-24)
1. Active UCSF 2M fixture baseline is now the 2026-02-24 rerun set:
   - sequential: `/storage/ucsf-2M/star_runs/fixture_ucsf2m_current_sequential`
   - dynamic: `/storage/ucsf-2M/star_runs/fixture_ucsf2m_current_dynamic`
2. New-run consistency checks show no regressions across current code paths:
   - sequential vs dynamic: `17/17` key outputs identical
   - sequential (`search=1`) vs sequential (`consumer=4/search=4`): `17/17` identical
3. Old 2026-02-18 fixture remains a historical benchmark reference but is no
   longer the canonical fixture for current validation.
4. Evidence artifacts:
   - `/tmp/ucsf2m_seq_vs_fixture_20260224_090754/COMPARE_REPORT.txt`
   - `/tmp/ucsf2m_seq_fixture_check_20260224_091148/SUMMARY.txt`

## UCSF Full Benchmark Refresh (2026-02-24)
1. Dynamic full-sample run (provision both map and PF to 32):
   - `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_full_dynamic_32x32_20260224_092512`
2. Effective dynamic profile used:
   - `--runThreadN 32`
   - `--dynamicThreadInterface 1`
   - `--dynamicThreadConstMapPermits 32`
   - `--dynamicThreadTelemetry 1`
   - `--crAssignConsumerThreads 32`
   - `--crAssignSearchThreads 1`
3. Wall-time comparison:
   - STAR dynamic `32x32`: `20:26.93`
   - STAR full fixture baseline: `30:20`
   - Cell Ranger 9 full (`32 cores`): `33:57.18`
4. Performance summary:
   - vs CR full: `~1.66x` throughput, `~39.8%` less wall time
   - vs STAR full fixture baseline: `~1.48x` throughput, `~32.6%` less wall time
5. Parity summary:
   - New full dynamic run matches STAR full fixture on checked core GEX/feature
     matrices and `crispr_analysis` outputs.

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

## Stage 4a: Phase-Aware Permit Handoff + Expandable PF Pools
Goal:
1. Increase available worker capacity, then throttle with permits so map and PF
   can overlap and finish at similar times.
2. After alignment completes, shift almost all permits to PF for tail collapse.

Implementation slices:
1. PF worker model update:
   split PF from fixed-size workers to max-size workers + permit-gated active
   work. Keep one global permit budget across map and PF.
2. Programmatic phase machine:
   `INIT_PRELOAD`, `ALIGN_ACTIVE`, `SOLO_POST_ALIGN`, `PF_ONLY`.
3. Deterministic transition hooks:
   transition on explicit STAR markers (`started mapping`, `finished mapping`,
   Solo finalize boundaries), not only periodic heuristics.
4. Target allocation policy:
   - `INIT_PRELOAD`: map low floor, PF gets remainder for preload/read-side prep.
   - `ALIGN_ACTIVE`: ETA-based balancing with bounded retune step.
   - `SOLO_POST_ALIGN`: keep small map reserve (2-4), hand remainder to PF.
   - `PF_ONLY`: map 0 (or reserve 1 if required), PF gets all remaining permits.
5. Stability controls:
   hysteresis threshold (10-15%), cooldown window, min floors, and starvation
   guards.

Required guards:
1. `runThreadN == 1`:
   disable dynamic mode and execute legacy sequential path.
2. Zero-work domains:
   if PF inputs absent, route all permits to map; if map done, route all permits
   to PF.
3. Error paths:
   guarantee permit release on all exits (success, early return, exception).

Telemetry requirements:
1. Emit phase transitions with timestamp and permit split.
2. Emit per-domain ETA snapshots used for retune decisions.
3. Emit final summary:
   time spent per phase, permit occupancy, retune count, and blocked-acquire
   metrics.

Acceptance criteria:
1. No deadlock/livelock under 100-read and 500-read smoke cases.
2. Deterministic outputs vs dynamic-off baseline (GEX + feature outputs).
3. `dynamic_on` wall time improves vs `dynamic_off` on UCSF perturb runs.
4. Post-alignment tail is materially reduced by PF permit expansion.

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

## Planned Full-Sample Benchmark (Not Executed Yet)
Goal:
1. Repeat the timing/parity comparison on full UCSF sample inputs to validate
   scaling behavior beyond the current comparison slice.

Execution plan:
1. Freeze binaries and provenance:
   `STAR` commit/build tag, `cellranger-9.0.1`, reference paths, feature ref.
2. Stage full-sample GEX/guide FASTQs to SSD-backed `/storage` before running
   (copy from source roots under `/mnt/pikachu/ucsf-perturb-seq/...`).
3. Build one shared input manifest after staging (absolute path, file size,
   checksum) and use this same manifest to generate both:
   STAR `--readFilesIn` lists and CR `--libraries` CSV roots.
4. Lock parameters to the downsampled comparison values and do not tune per run:
   STAR: `--runThreadN 32`, `--dynamicThreadConstMapPermits 24`,
   `--dynamicThreadInterface {0|1}`, PF assign flags
   (`crAssignMaxHamming=1`, `crAssignFeatureOffset=0`,
   `crAssignLimitSearch=-1`, `crAssignMinCounts=0`,
   `crAssignMaxBarcodeMismatches=5`, `crAssignFeatureN=1`,
   `crAssignBarcodeN=2`, `crAssignConsumerThreads=4`,
   `crAssignSearchThreads=4`), and same Solo/GEX flags as downsampled run.
   CR: `--include-introns true --min-crispr-umi 3 --create-bam false` with
   `--localcores 32`.
5. Run three jobs on the same host under minimal background load:
   STAR `dynamic_off` (32), STAR `dynamic_on` (32, map permits 24 start),
   Cell Ranger 9 (32 cores, same reference family).
6. Collect timing/memory:
   STAR `/usr/bin/time -v` and `Log.out` phase markers;
   CR9 wrapper `time -v`, `_log`, and `_perf`.
7. Run parity checks on shared outputs:
   GEX shared matrices, CRISPR calls, threshold files, and feature matrices.
8. Record semantic differences with explicit `L1`/`max_abs` summaries if byte
   parity fails.

Acceptance criteria:
1. No crashes/deadlocks.
2. `dynamic_on` remains at least neutral to `dynamic_off` on wall time.
3. CRISPR outputs remain deterministic across repeated CR9 runs at fixed inputs.
4. STAR and CR are confirmed to use the same staged FASTQ set via shared
   manifest checks.
5. Any STAR dynamic-vs-baseline differences are bounded and documented with
   per-feature/per-barcode localization.

## Success Criteria
1. Default mode remains unchanged.
2. Dynamic mode yields reproducible throughput improvements on perturb-seq
   workloads.
3. Architecture is upstream-acceptable and extensible to producer/consumer
   mode without another deep rewrite.
