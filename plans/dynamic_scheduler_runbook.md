# Dynamic Scheduler Runbook (Alignment + PF)

Date: 2026-02-18  
Branch: `core-alignment-threads-integration`  
Latest commit with current wiring: `595e280`

## Purpose
Capture the working design and guardrails for moving from dynamic permit
throttling to a true scheduler-backed producer/consumer model, while keeping
current behavior safe and testable at each step.

## Current State
1. Alignment and PF can share dynamic permit telemetry/hooks.
2. PF read cap is wired from STAR `--readMapNumber` into pf-api `maxReads`.
3. PF controller (`off`/`shadow`/`active`) is instrumented into
   `assignBarcodes.api_run.txt`.
4. Dynamic permit stall diagnostics are available:
   `blockedAcquires`, `waitTimeoutEvents`, `stallWarnEvents`,
   `waiters`, `lastReleaseAgoNs`.

## Confirmed Gaps
1. No central scheduler queue for cross-module load balancing yet.
2. Feature worker topology is still static at startup.
3. Merge/finalize remain serial tail work (not scheduler tasks).
4. Map side has no native work queue; pressure must be inferred from chunk
   telemetry and permit wait metrics.

## Design Decisions (Locked for Next Iteration)
1. Queue only normal work:
   `MAP_CHUNK`, `PF_ASSIGN_BLOCK`.
2. Do not queue merge/finalize initially.
3. Treat merge/finalize as post-processing exclusive phases:
   pause normal dispatch, run immediately, resume.
4. Single-thread mode bypasses dynamic scheduler entirely:
   legacy `align -> assign -> finalize`.
5. Zero-work fast path exits immediately (or runs post-only work directly).
6. Rebalancing is done by retargeting permits, not by moving OS threads.
7. Work-size sizing uses an anchored estimate:
   exact count on the smallest input file, then byte scaling to total input.
8. If estimated remaining work goes negative, inflate estimate by 10% until
   remaining is non-negative.

## Task Model
1. `MAP_CHUNK`
   one STAR mapping chunk.
2. `PF_ASSIGN_BLOCK`
   one PF consumer block (bounded batch size).
3. Post-processing (not queued in MVP):
   merge + finalize phases as exclusive execution windows.

## Scheduler Policy (MVP)
1. Start fixed worker pools at launch; do not spawn/teardown for retunes.
2. Use one global permit budget `N=runThreadN`.
3. Track per-domain targets `mapTarget` and `pfTarget` with
   `mapTarget + pfTarget = N`.
4. `acquire(domain)` requires:
   global permit availability and domain `inUse < target`.
5. Retune targets at coarse intervals (2-5s), moving at most one permit per
   tick, with cooldown/hysteresis.
6. If one side is done, immediately hand all permits to the active side.
7. Post-processing still preempts normal dispatch.

## Work Estimation Policy (Locked)
1. Build an anchored estimate per domain (`map` and `pf`):
   pick the smallest file from the domain's primary file set
   (`R1` for map, barcode read for PF), then compute exact reads with
   `zcat | wc -l` divided by 4.
2. Scale anchor reads by compressed byte ratio:
   `estimate = anchorReads * (totalBytesAcrossDomain / anchorBytes)`.
3. Prefer `--readMapNumber` as exact total when provided.
4. Maintain monotonic estimate:
   `estimate = max(estimate, done * 1.05)`.
5. If `remaining = estimate - done < 0`, repeatedly apply `estimate *= 1.10`
   until `remaining >= 0`.

## Edge-Case Handling
1. `runThreadN <= 1`:
   force legacy sequential path.
2. No tasks available:
   if producers done and no post pending, exit cleanly.
3. If both sides have pending work, enforce floors (`mapMin>=1`, `pfMin>=1`).
4. If one side has no pending work, floors can be relaxed for the active side.
5. Exit predicate:
   map+pf work complete and no post pending.

## Invariants
1. No permit leaks on any success/failure path.
2. No new work dispatch during exclusive post-processing.
3. Deterministic outputs in default mode and single-thread fallback.
4. Dynamic mode never deadlocks on empty queues or zero budgets.

## Implementation Phases
1. Phase A: scheduler skeleton + single-thread bypass + zero-work fast path.
2. Phase B: domain permit caps with static split (`mapTarget`, `pfTarget`).
3. Phase C: ETA-based retune loop with estimate-only work sizing.
4. Phase D: optional producer/consumer refinement (if needed later).

## Tests Required Per Phase
1. `runThreadN=1`, no work.
2. `runThreadN=1`, map-only work.
3. `runThreadN=1`, PF-only work.
4. `runThreadN=2+`, empty queues after startup.
5. `runThreadN=2+`, mixed map+PF, post phase triggered.
6. Dynamic controller `shadow` and `active` parity checks.
7. Full 100K smoke uncapped.
8. Estimate overshoot path:
   verify estimate auto-inflates and controller remains stable.

## Runtime Commands
1. Fast capped probe (100):
```bash
PF_DYNAMIC_100K_READ_MAP_NUMBER=100 \
PF_DYNAMIC_100K_RUN_BASELINE=0 \
PF_DYNAMIC_100K_RUN_VARIABLE_CASES=0 \
PF_DYNAMIC_100K_RUN_FORCED_EXIT_PROBE=0 \
PF_DYNAMIC_100K_RUN_PF_CONTROLLER_CASES=1 \
tests/run_pf_dynamic_permit_100k_smoke.sh
```
2. Full uncapped 100K:
```bash
PF_DYNAMIC_100K_RUN_BASELINE=0 \
PF_DYNAMIC_100K_RUN_VARIABLE_CASES=0 \
PF_DYNAMIC_100K_RUN_FORCED_EXIT_PROBE=0 \
PF_DYNAMIC_100K_RUN_PF_CONTROLLER_CASES=0 \
tests/run_pf_dynamic_permit_100k_smoke.sh
```

## Artifact Keys to Track
1. `dynamicPermitDelta.acquires`
2. `dynamicPermitDelta.feature.acquires`
3. `dynamicPermitDelta.blockedAcquires`
4. `dynamicPermitDelta.waitTimeoutEvents`
5. `dynamicPermitDelta.stallWarnEvents`
6. `dynamicPermitAfter.waiters.max`
7. `pfController.mode`, `pfController.ticks`, `pfController.applied`
8. `maxReads`

## Known Caveat
1. PF `maxReads` currently applies per producer/lane-set in practice.
   On 2 lanes with `maxReads=100`, observed PF work is 200 reads total.
   Keep this behavior explicit until we decide between:
   per-lane cap vs global cap.

## Next Actionable Step
1. Implement Phase A guardrails first:
   single-thread hard bypass + zero-work fast path + explicit scheduler state
   predicates.
