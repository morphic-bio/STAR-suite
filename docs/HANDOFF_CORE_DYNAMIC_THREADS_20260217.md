# Handoff: STAR Core Dynamic Threads for Upstreaming

Date: 2026-02-17 (updated 2026-02-18)  
Current working branch: `core-alignment-threads-integration`  
Original dynamic branch: `core-dynamic-threads`  
Previous context branch: `noBCfix`

## Resume Snapshot (2026-02-18)
1. `master` has been fast-forwarded to include compatibility fixes and docs:
   - `5a76763` (`docs: record 2026-02-18 compatibility gate and status`)
2. Dynamic-thread work from `core-dynamic-threads` was ported onto the current branch:
   - `474f3d0` (dynamic-thread handoff/plan docs)
   - `a882fe7` (interface stability docs)
   - `40c4bc2` (core dynamic-thread retune trace + harness)
3. Remote branch for continued work:
   - `origin/core-alignment-threads-integration` @ `40c4bc2`
4. Recommendation:
   - Continue all new implementation on `core-alignment-threads-integration`
   - Treat `core-dynamic-threads` as historical source context only.

## Why This Exists
This handoff reframes the perturb parallelization work as an upstream-ready STAR
core initiative:
1. First deliver dynamic active-thread control in core.
2. Expose integration hooks for feature and future modules.
3. Then transition toward producer/consumer architecture in staged phases.

## Canonical Documents
1. Upstream plan:
   `plans/STAR-core-dynamic-threads-upstream-plan.md`
2. Original perturb-focused plan (source context):
   `plans/STAR-perturb-parallel.plans.md`
3. Original handoff (source context):
   `docs/HANDOFF_STAR_PERTURB_PARALLEL_20260217.md`
4. UCSF parity context:
   `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/RESULTS.md`

## Architecture Direction
1. Use fixed worker pools with dynamic permits first.
2. Centralize permit logic in STAR core (`ThreadPermitController` direction).
3. Keep module integration optional through a stable interface.
4. Keep all behavior default-off behind explicit flags.

## Interface Migration Expectation
1. If Interface v1 is frozen in Stage 2, module-facing changes for a later full
   producer/consumer STAR are expected to be small.
2. Most churn will be internal to STAR core execution and queue ownership.
3. Likely additive API elements for producer/consumer mode:
   work-item identity, completion status, cancellation/flush controls.
4. CLI contract and default behavior should remain unchanged across models.

## Immediate Implementation Targets
1. STAR core gating path:
   `core/legacy/source/ReadAlignChunk_processChunks.cpp`
2. Control state and synchronization:
   `core/legacy/source/ThreadControl.h`,
   `core/legacy/source/ThreadControl.cpp`
3. CLI parameter wiring:
   `core/legacy/source/Parameters.h`,
   `core/legacy/source/Parameters.cpp`,
   `core/legacy/source/parametersDefault`
4. Module integration path (later stage):
   `core/features/process_features/src/assignBarcodes.c`

## Non-Negotiable Invariants
1. Default mode unchanged when dynamic flags are off.
2. Ordering/determinism preserved, including `PairedKeepInputOrder`.
3. `FILE` marker and `readFilesIndex` semantics preserved.
4. Global chunk/counter behavior preserved (`P.iReadAll`, `chunkInN/chunkOutN`).
5. No permit leaks on error or early-return paths.

## Staged Delivery Checklist

## Stage 0 Checklist
1. Add telemetry-only metrics and disabled controller scaffolding.
2. Confirm no output diffs with flags off.
3. Capture baseline timing and queue metrics.

## Stage 1 Checklist
1. Add dynamic permit flags (default off).
2. Gate STAR chunk processing with static split.
3. Validate repeated-run determinism and deadlock safety.

## Stage 2 Checklist
1. Add module-facing permit API surface.
2. Add C ABI bridge for C modules when enabled.
3. Add unit/integration checks for controller lifecycle.
4. Status (2026-02-18): `pf-api` permit hook surface and C bridge are now wired
   for perturb feature assignment (`PfMultiAssign` -> `pf_api` -> `assignBarcodes`
   consumer loop), with runtime gating preserved behind dynamic-thread flags.

## Stage 3 Checklist
1. Integrate feature consumer permit gating.
2. Keep existing queue behavior unchanged.
3. Validate parity and throughput on perturb workloads.

## Stage 4 Checklist
1. Add queue-aware dynamic rebalancer with bounds.
2. Add anti-oscillation controls (cooldown/hysteresis).
3. Stress test low/high feature load and starvation limits.

## Stage 5-6 Checklist
1. Refactor for producer/consumer readiness behind flags.
2. Add producer/consumer pilot mode as experimental.
3. Benchmark pilot against Stage 4 dynamic mode.

## Validation Gates
1. Core regressions:
   `tests/run_cbub_regression_test.sh`
2. CR-compat gate:
   `tests/test_cr_compat_crispr_calling.sh`
3. Perturb parity/performance:
   A375 and UCSF contexts with dynamic mode off/on comparisons.

## Artifact Hygiene
1. Keep generated outputs in `/tmp` or `tests/*_output*/`.
2. Add new stable output paths to `tests/ARTIFACTS.md`.
3. Do not commit datasets, binaries, or large generated artifacts.

## First PR Recommendation
1. Stage 0 + Stage 1 only:
   flags, scaffolding, STAR-only gating, telemetry.
2. Explicitly include validation proof that flag-off mode is unchanged.
3. Defer module API and rebalancer to follow-up PRs for cleaner upstream review.

## Immediate Next Actions for Next Agent
1. Rebuild from clean state before evaluating behavior:
   - `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`
2. Re-run dynamic tiny-fixture validation on this branch:
   - `tests/run_dynamic_threads_tiny_fixture.sh`
   - `tests/run_dynamic_threads_variable_sequences_smoke.sh`
3. Confirm no regressions with dynamic flags off (parity baseline).
4. Start Stage 2 implementation slices only after Stage 1 rerun is green:
   - module-facing permit API scaffolding
   - C ABI bridge shape for feature consumers
5. Keep default behavior unchanged unless new flags are explicitly set.
