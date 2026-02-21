# Handoff: STAR Perturb Parallel Integration (Dynamic Thread Permits)

Date: 2026-02-17 (updated 2026-02-18)  
Current branch for continuation: `core-alignment-threads-integration`  
Legacy source branch: `core-dynamic-threads`  
Current tip: `40c4bc2`

## Branch Update (2026-02-18)
This handoff was originally written against `core-dynamic-threads`. The active
continuation branch is now `core-alignment-threads-integration`, which already
contains the dynamic-thread commits and sits on top of the latest compatibility
baseline merged into `master` (`5a76763`).

## Scope
This handoff is for implementing the STAR perturb parallelization sketch for STAR mapping + feature assignment coordination.

Primary objective:
- Add dynamic thread-permit scheduling across STAR mapping and feature assignment without rewriting STAR core into a full producer/consumer engine.

Non-goals for first implementation:
- No full STAR chunk manager rewrite.
- No default behavior changes.

## Canonical References
- Runbook: `plans/STAR-perturb-parallel.plans.md`
- Upstream staged plan: `plans/STAR-core-dynamic-threads-upstream-plan.md`
- Tiny fixture reproducibility runbook: `docs/DYNAMIC_THREADS_TINY_FIXTURE_RUNBOOK_20260217.md`
- Todo anchor: `docs/todos` (dynamic permit item)
- UCSF call-only parity context: `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/RESULTS.md`

## Key Decisions Already Made
1. Use a hybrid model with fixed worker pools and dynamic active permits.
2. Keep STAR thread creation unchanged; dynamic behavior means gating active work, not spawning/killing threads mid-run.
3. Roll out in phases behind a feature flag (default off).
4. Validate STAR-only permit gating first, then cross-pipeline rebalancing.

## Implemented on `core-dynamic-threads` (current staged state)
1. STAR-core dynamic-thread interface and telemetry are wired behind:
   - `--dynamicThreadInterface`
   - `--dynamicThreadConstMapPermits`
   - `--dynamicThreadTelemetry`
2. Variable permit mode is wired behind:
   - `--variableThreads`
   - `--variableThreadsRetuneEveryAcquires`
   - `--variableThreadsPermitSequence`
3. Deterministic permit retuning is implemented in `ThreadControl` and logged via
   interface/telemetry lines in `Log.out` (including bounded `retuneTrace`).
4. Repro harness and mock consumer are added:
   - `tests/run_dynamic_threads_tiny_fixture.sh`
   - `tests/dynamic_threads/mock_consumer_report.py`
   - `tests/dynamic_threads/test_mock_consumer_report.py`
   - `tests/run_dynamic_threads_variable_sequences_smoke.sh`
5. Requested smoke transition scenarios pass:
   - `3->2->4` (with `runThreadN=4`, initial permits `3`)
   - `1->2->1` (with `runThreadN=2`, initial permits `1`)
6. Tiny harness now supports optional sorted-BAM parity checks and canonical
   `Log.final.out` parity diffing for dynamic-vs-baseline validation.
7. Phase 2 (`pf-api` first) integration now wires permit hooks through:
   - `core/features/process_features/include/pf_api.h`
   - `core/features/process_features/src/pf_api.c`
   - `core/features/process_features/src/assignBarcodes.c`
   - `core/legacy/source/PfMultiAssign.cpp`
   - `core/legacy/source/PfMultiProcess.cpp`
   This remains gated by the dynamic-thread interface (`--dynamicThreadInterface 1`);
   default behavior is unchanged when disabled.
8. Phase 3a Slice A/B groundwork is now in place:
   - domain-aware permit telemetry in `ThreadControl` (map vs feature)
   - typed feature hook context in `PfMultiAssign` callback bridge
   - 100K perturb smoke harness: `tests/run_pf_dynamic_permit_100k_smoke.sh`
   - latest validation (2026-02-18): dynamic hook-on run reports non-zero
     `dynamicPermitDelta.feature.*`, and baseline-vs-dynamic raw/filtered MEX
     parity checks pass on tier-100000 A375 fixture.
9. Phase 3a Slice C test coverage now includes:
   - variable permit retune sequences in pf stage (`3->2->4`, `1->2->1`)
   - acquire/work parity assertions (`acquires == workUnits`) for aggregate and
     feature-domain counters
   - map-domain zero-delta assertion in pf-only stage
   - forced early-exit timeout probe followed by successful recovery run.
10. Phase 3a Slice D external-controller scaffold is now wired:
   - new STAR params:
     - `--dynamicThreadPfControllerMode off|shadow|active`
     - `--dynamicThreadPfControllerIntervalMs`
     - `--dynamicThreadPfControllerSequence`
     - `--dynamicThreadPfControllerMaxUpdates`
   - `PfMultiProcess` can run a pf-local controller thread in `shadow` (propose-only)
     or `active` (applies `mapPermitSetTargetPermits`) mode.
   - 100K pf smoke now includes controller `shadow` and `active` cases behind
     `PF_DYNAMIC_100K_RUN_PF_CONTROLLER_CASES`.

## Validation Notes (2026-02-18)
- CLI parameter contract for the new controller flags is validated:
  - invalid mode token fails with explicit error
  - `active` mode requires `--variableThreads 1`
  - non-`off` mode requires non-empty controller sequence + interval
  - valid `shadow`/`active` combos pass parameter checks and continue to normal
    file validation.
- Full 100K fixture smoke is currently blocked by a mapping-stage crash on this
  branch before pf stage is reached (`std::length_error` / segfault observed in
  baseline and dynamic runs), so controller runtime assertions were not
  re-collected in this pass.

## Current Runtime/Architecture Observations
- Feature pipeline already has producer/consumer queues:
  - `core/features/process_features/src/assignBarcodes.c`
- STAR mapping is thread-parallel with shared input lock, not central queue:
  - `core/legacy/source/mapThreadsSpawn.cpp`
  - `core/legacy/source/ReadAlignChunk_processChunks.cpp`

## Implementation Entry Points

### STAR-side permit gating (Phase 1)
- `core/legacy/source/ReadAlignChunk_processChunks.cpp`
  - Gate chunk processing before `mapChunk()` with map permits.
  - Guarantee permit release on all exits.
- `core/legacy/source/ThreadControl.h` / `core/legacy/source/ThreadControl.cpp`
  - Add shared permit counters + synchronization.
- `core/legacy/source/Parameters.h` + `core/legacy/source/Parameters.cpp` + `core/legacy/source/parametersDefault`
  - Add/validate dynamic-thread and variable-thread flags (default off).

### Feature-side permit gating (Phase 2)
- `core/features/process_features/src/assignBarcodes.c`
  - Gate consumer work loop with feature permits.
  - Keep existing queue semantics unchanged.

### Rebalancer (Phase 3)
- Add a controller loop to adjust permit split from queue pressure and service rates.
- Guardrails:
  - min permits for both sides
  - bounded step size
  - hysteresis/cooldown to avoid oscillation

## Phase 3a Hook-Up Plan (Tracked Mirror)
This mirrors the expanded Phase 3a details in `plans/STAR-perturb-parallel.plans.md`
so the handoff stays self-contained in tracked docs.

### Interface v1 freeze
1. Keep one shared scheduler contract for map + feature clients:
   - `acquire(domain)`
   - `release(domain, wait_ns, work_units, work_bytes, work_ns)`
   - `set_target(domain, permits)` (controller-owned updates)
   - `snapshot()`
2. Do not fork a separate feature scheduler API; adapt per domain on top of one controller.
3. Keep `pf_api` callback ABI additive/backward-compatible.

### Adapter ownership
1. Map adapter remains STAR map-thread path (`ReadAlignChunk_processChunks.cpp`).
2. Feature adapter remains `PfMultiAssign` callback bridge into `assignBarcodes`.
3. Feature bridge should pass typed hook context (`domain=feature`) rather than null context.
4. Controller target updates stay STAR-core-owned; pf consumers are permit clients.

### Producer/consumer integration details
1. Keep `consume_reads` permit acquire after dequeue/copy and before match/assign work.
2. Keep permit release after one consumed work unit is fully processed.
3. Keep one consumed block as one permit work unit.
4. Add starvation guard assertions in smoke tests (`acquire == release`, bounded wait).

### Telemetry and validation gates
1. Add per-domain snapshot counters (map + feature) while preserving aggregate counters.
2. Keep current `dynamicPermitDelta.*` fields and add per-domain deltas in `assignBarcodes.api_run.txt`.
3. Validate with:
   - tiny callback smoke fixture
   - 100K perturb integration smoke (`--dynamicThreadInterface 1 --dynamicThreadTelemetry 1`)
   - hook-disabled parity check (`--dynamicThreadInterface 0`)

## Invariants That Must Not Regress
1. Determinism and ordering behavior, especially `PairedKeepInputOrder`.
2. File-boundary semantics (`FILE` markers, `readFilesIndex`, per-file stop logic).
3. Global counters/chunk numbering (`P.iReadAll`, `chunkInN/chunkOutN`).
4. Default mode behavior/performance with feature flag off.

## Validation Sequence (Required)

### Step A: STAR-only gating
- Enable flag with fixed map permits.
- Verify output parity against baseline (flag off).
- Verify no deadlocks/permit leaks.

### Step B: A375 regression gate
- Run CR-compatible A375 parity checks before enabling any dynamic rebalance by default.

### Step C: Feature gating + rebalance
- Add feature permit consumption.
- Add dynamic controller.
- Re-run parity + throughput on UCSF sample/full set.

## Benchmarks to Capture
- End-to-end wallclock
- Reads/minute (mapping + feature)
- CPU utilization
- I/O wait
- Queue occupancy/permit history

Use `/tmp` or `tests/*_output*/` for outputs. Add stable artifact paths to `tests/ARTIFACTS.md` as needed.

## Current Comparison Context (Do Not Re-solve First)
From UCSF call-only parity checks:
- With `--min-umi 3`, common-cell set-equivalent `feature_call` concordance is 100%.
- Residual row-count differences are mostly zero-call row emission policy (`None`) and one STAR-only non-`None` edge call.

This is useful downstream context, but not blocking for permit scheduler implementation.

## Recommended First PR Slice
1. Add feature flag + permit structs.
2. STAR-only permit gating with static map permit count.
3. Add telemetry for permit wait/acquire/release timings.
4. Provide validation note showing no output diffs with feature flag off.

Then follow with feature-side gating and the rebalancer.

## UCSF Parity Debug Status (2026-02-18, unresolved)
Context:
- Requested target: reproduce prior UCSF call parity where calls are identical
  except one outlier (`--min-umi 3`: `5766/5767` set-equivalent).
- This target is documented in:
  `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/RESULTS.md`.
- Additional GEX parity goal for recovery automation:
  Pearson `>= 0.93`, Spearman `>= 0.94` on normalized/common barcodes.

What was run (ad hoc path, not accepted as gate):
1. `scripts/run_gex_feature_parity_checks.sh` against:
   - STAR:
     `/storage/ucsf-2M/star_runs/star_baseline_iPSC2_1_AALG2_1M_nxt_20260217_160217`
   - CR:
     `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813`
2. Translation mode tested with CR-side mapping
   (`--translate-side cr --translation-direction left-to-right`).
3. Artifacts:
   - `/tmp/ucsf2m_oldscripts_try_20260218_222834/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
   - `/tmp/ucsf2m_oldscripts_try_20260218_222834/FILTERED_BARCODE_SET_OVERLAP.txt`
4. Same run repeated against CR run `...064209`:
   - `/tmp/ucsf2m_oldscripts_try_064209_20260218_222907/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`

Observed outcome from ad hoc path:
- Feature-call parity on common rows remained `831/937` (`88.6873%`), not the
  expected one-outlier behavior.
- This confirms the current mixed GEX+feature parity harness is not the right
  gate for the requested call-only parity target.

Action taken:
1. Added dedicated recovery runbook:
   `docs/UCSF_PARITY_RECOVERY_RUNBOOK_20260218.md`
2. Runbook defines:
   - locked inputs
   - clean-build requirement
   - script plan for a canonical parity wrapper
   - explicit pass/fail thresholds aligned to prior known-good results
   - explicit GEX parity goals (`>=0.93` Pearson, `>=0.94` Spearman)
3. Next agent should implement:
   `scripts/run_ucsf_call_parity_recovery.sh`
   by wrapping (not replacing):
   `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/run_callonly_gmm_parity.sh`

## Related Follow-Up (2026-02-19)
For UCSF/A375 parity continuation specifically focused on Cell Ranger multimap
compatibility policy, use:

- `docs/HANDOFF_CR_COMPAT_MULTIMAP_POLICY_20260219.md`

That handoff tracks the updated analysis scripts, matrix-selection controls,
and the proposed compatibility-mode implementation plan (default OFF).
