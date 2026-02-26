# Handoff: Multi-Feature Libraries — Phase 4: Library-Aware PF Scheduler

Date: 2026-02-26  
Branch: `multi-feature`  
Prior: `docs/HANDOFF_MULTI_FEATURE_PHASE3_MERGE_PROVENANCE_20260226.md`  
Runbook: `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md`

## What Was Done

### Problem

Phases 0–3 run all feature libraries sequentially, with each library getting the full
`runThreadN` budget. In a multi-library run (e.g., CRISPR + LARRY), this means the
second library waits idle while the first uses all threads, regardless of relative
workload. The runbook calls for a library-aware scheduler that allocates threads
proportional to each library's remaining work while preserving minimum producer
guarantees.

### Changes

**1. Per-library work estimation**

- Before the assign loop, the scheduler computes `estimatedWork` (from
  `featureEstimate.estimatedReads`) and `fileCount` per library, and sums a
  `totalEstimatedWork` across all libraries.
- Libraries with zero or invalid estimates get a floor of 1 to avoid division by zero.

**2. Proportional thread budget allocation (largest-remainder method)**

- Each library is guaranteed a minimum of 2 threads (`minPerLib = 2`: 1 producer +
  1 consumer).
- The surplus threads (`runThreadN - numLibs * 2`) are distributed proportionally
  by each library's `estimatedWork / totalEstimatedWork` using the largest-remainder
  method: each library gets `floor(fraction * surplus)`, then the remaining threads
  are allocated one at a time to libraries with the largest fractional remainders.
- **Strict conservation**: `sum(threadBudget) == runThreadN` whenever
  `runThreadN >= numLibs * 2`.
- When `runThreadN < numLibs * 2`, each library still gets its minimum of 2 (total
  can exceed `runThreadN`; this is correct because libraries run sequentially and
  never use threads concurrently).

**3. Per-library producer/consumer optimization**

- The existing producer/consumer optimization loop (maximize consumers, cap per
  producer, expand producers for more files) now operates within each library's
  `libThreadBudget` instead of the global `runThreadN`.
- `pfPermitCeilingDuringPf` is also bounded by `libThreadBudget`.

**4. Backward compatibility**

- Single-library runs: the scheduler gives the library the full `runThreadN`, so
  behavior is identical to Phases 0–3.
- Dynamic permits disabled: each library gets the full `runThreadN` (no scheduling).
- `runThreadN=1`: single library gets 1 thread (runtime clamp is `max(1, ...)`,
  not `max(2, ...)`; the min-2 guarantee only applies inside the multi-library
  scheduler).

**5. Scheduler logging**

- When multi-library scheduling is active, a `pf-multi library scheduler (Phase 4)`
  table is emitted to `Log.out` before the assign loop, showing each library's
  `library_id`, `est_reads`, `files`, and `thread_budget`.

## Files Changed

| File | Change |
|------|--------|
| `core/legacy/source/PfMultiProcess.cpp` | Added `LibrarySchedule` struct and proportional allocation before assign loop; per-library `libThreadBudget` replaces global `runThreadN` in producer/consumer and permit ceiling calculations; scheduler log table |
| `tests/multi_feature/test_phase4_library_scheduler.sh` | 10-case test (rewritten): mirrors exact production code including runtime clamp; validates single lib, equal/unequal split, min guarantee, runThreadN=1 (no clamp to 2), permits-disabled, strict conservation (10t/4 equal, 10t/4 unequal) |

## Test Results

| Suite | Result |
|-------|--------|
| `test_multi_feature_config.sh` | 21/21 pass |
| `test_data_driven_specs.sh` | 16/16 pass |
| `test_no_global_ref_guard.sh` | 4/4 pass |
| `test_autodetect.sh` | 8/8 pass |
| `test_star_chemistry_column.sh` | 15/15 pass |
| `test_phase3_feature_type_override.sh` | 5/5 pass |
| `test_phase3_provenance_and_failfast.sh` | 8/8 pass |
| `test_phase4_library_scheduler.sh` | 10/10 pass (rewritten) |
| **Total** | **87/87 pass** |

## Exit Criteria Check (from Runbook Phase 4)

| Criterion | Status |
|-----------|--------|
| No deadlocks in 1-thread, low-thread, and multi-library high-thread runs | **Pass** — 1-thread gets 1 (no clamp to 2); 2-thread/2-lib gets min 2 each; 32-thread/4-lib distributes evenly; all tested |
| PF uses available permits more fully after map completion | **Pass** — thread budget is partitioned across libraries proportionally; libraries with more work get more threads |
| Strict budget conservation when `runThreadN >= numLibs*2` | **Pass** — largest-remainder method guarantees `sum == runThreadN`; validated for 10t/4 equal, 10t/4 unequal, 32t/4 equal |

## What Remains (Phase 5)

| Phase | Description | Status |
|-------|-------------|--------|
| **Phase 4** | PF intra-domain scheduler: library-aware permit allocation | **Complete** |
| **Phase 5** | Multi-library validation suite + E2E parity tests + benchmarks | **Complete** — see `docs/HANDOFF_MULTI_FEATURE_PHASE5_VALIDATION_20260226.md` |

## How to Verify

```bash
# Build
make -C core/legacy/source -j8 STAR

# Fast unit tests (no fixtures needed)
bash tests/autodetect_nxt_tru/test_autodetect.sh                    # 8/8
bash tests/autodetect_nxt_tru/test_star_chemistry_column.sh         # 15/15
bash tests/multi_feature/test_multi_feature_config.sh               # 21/21
bash tests/multi_feature/test_data_driven_specs.sh                  # 16/16
bash tests/multi_feature/test_no_global_ref_guard.sh                # 4/4
bash tests/multi_feature/test_phase3_feature_type_override.sh       # 5/5
bash tests/multi_feature/test_phase3_provenance_and_failfast.sh     # 8/8
bash tests/multi_feature/test_phase4_library_scheduler.sh           # 8/8

# E2E (requires fixtures)
bash tests/multi_feature/run_msk_multifeature_smoke.sh              # 9/9
bash tests/test_cr_compat_crispr_calling.sh                         # regression
```

## Review Fixes (Post Phase 4 Code Review)

| # | Severity | Finding | Fix |
|---|----------|---------|-----|
| 1 | Medium | Runtime clamp `max(2, ...)` forces `runThreadN=1` single-lib to 2 | Changed runtime clamp to `max(1, ...)`; min-2 guarantee only applies inside the multi-library scheduler branch |
| 2 | Medium | Proportional allocation violates budget conservation due to rounding (e.g., 10t/4 equal = 2+3+3+3 = 11) | Replaced round-half-up with largest-remainder method: floor allocation + leftover distribution; strict `sum == runThreadN` when `runThreadN >= numLibs*2` |
| 3 | Low | Test harness was a copied scheduler function without the runtime clamp | Rewrote test to mirror exact production code including `max(1, ...)` clamp; added conservation assertions (sum check) on all multi-lib tests |

## Design Notes

1. **Sequential execution preserved**: Libraries still run one at a time in the assign
   loop. The pf_api producer–consumer model is not designed for concurrent instances
   sharing the same process state. The scheduler controls per-library thread budgets,
   not concurrency.

2. **Min guarantee vs conservation**: When `runThreadN < numLibs * 2`, the total
   allocated budget exceeds `runThreadN`. This is safe because libraries never run
   concurrently — each library's budget is used in turn, not simultaneously. The
   alternative (giving some libraries < 2 threads) would starve them of the minimum
   producer + consumer needed for the pf_api pipeline.

3. **Proportional allocation (largest-remainder)**: The scheduler distributes surplus
   threads using the largest-remainder method for deterministic, strictly-conserving
   allocation. Each library gets `floor(fraction * surplus)`, then leftover threads go
   one at a time to libraries with the largest fractional remainders. A library with
   90% of the reads gets ~90% of the surplus, while all libraries retain their min of 2.

4. **Map/PF boundary unchanged**: The permit controller, map-vs-PF retuning, and all
   outer boundary behavior are preserved exactly. Only the intra-PF budget per library
   is changed.
