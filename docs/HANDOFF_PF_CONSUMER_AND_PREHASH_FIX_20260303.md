# PF Consumer Cap & d2 Prehash Budget Fix

Date: 2026-03-03
Branch: `CR-Larry-perturb`
Parent commit: `2e0b93c`

## Problem Statement

Two issues were identified from the MSK 30-KO two-feature benchmark runs:

1. **Capped consumers**: With 32 threads and 2 feature libraries, each library
   only got 3 consumer threads (vs ~7 for a single-library run). The root
   cause was twofold: the thread budget was split proportionally between
   libraries (even though they run sequentially), and each consumer's budget
   was divided by a "search thread quantum" of 4.

2. **d2 prehash never built when requested**: When a user sets
   `--crAssignMaxHamming 2`, the d2 prehash was blocked by a 50M-entry budget
   cap. For large feature sets like LARRY (246K barcodes, ~435M estimated d2
   entries), the budget check prevented the d2 hash from being built, silently
   falling back to d1-only matching.

## Changes

### 1. `core/features/process_features/src/io.c` — Remove d2 entry budget cap

**Before**: If the estimated d2 entry count exceeded `feature_prehash_max_entries`
(default 50M), the d2 hash was skipped entirely and only d1 was built.

**After**: The d2 hash is always built when `feature_prehash_max_hamming >= 2`.
A warning is emitted if the estimated entry count exceeds 500M.

The d1 entry budget check is retained (if d1 itself would exceed the budget,
prehash is disabled entirely — this guards against truly pathological inputs).

### 2. `core/features/process_features/src/pf_api.c` — No change to default

The `pf_config` default `max_hamming_distance` remains **1**. This is the
correct default: most runs use d1, and d2 should be an explicit opt-in via
`--crAssignMaxHamming 2`.

The existing pre-build clamp in `pf_load_feature_ref` ensures the prehash
level matches the per-library assignment `max_hamming_distance`:

- Default (`maxHamming=1`): prehash built to d1 only, d1 assignment.
- Explicit `--crAssignMaxHamming 2`: prehash built to d1+d2, d2 assignment.

The clamp prevents building d2 prehash that would never be used at runtime.

### 3. `core/legacy/source/PfMultiProcess.cpp` — Uncap consumer threads

**Before** (two problems):

a) Multi-library budget split: When `numFeatureLibs > 1`, the scheduler
   divided `P.runThreadN` proportionally between libraries using the
   largest-remainder method. With 32 threads and 2 equal-work libraries, each
   got 16 threads. But libraries run sequentially (pf_api holds a mutex), so
   half the machine sat idle during each library's processing.

b) Quantum divisor: The auto consumer calculation divided the remaining budget
   by `pfSearchThreadQuantum` (4), assuming each consumer permanently occupies
   4 OMP threads for feature matching. With prehash active, matching is O(1)
   hash lookup — the OMP loop rarely fires.

c) Per-producer cap: `pfConsumerCapPerProducer = 8` limited consumers to
   8 per producer slot, further restricting parallelism.

**After**:

- Each library gets the full `P.runThreadN` budget regardless of library count.
- Auto mode: `pfConsumerThreadsForRun = libThreadBudget - 1` (reserve 1 for
  the producer, rest are consumers). No quantum divisor, no per-producer cap.
- Explicit mode (`--crAssignConsumerThreads N`): honored directly.
- The `pfConsumerCapPerProducer` constant and associated cap-lifting logic are
  removed.

**Effect with 32 threads, 2 feature libraries:**

| Metric                      | Before | After |
|-----------------------------|--------|-------|
| Per-library thread budget   | 16     | 32    |
| Consumer threads per library| 3      | 31    |

## Behavioral Summary

| Scenario | d1 hash | d2 hash | Assignment | Consumers (32 threads) |
|----------|---------|---------|------------|----------------------|
| Default (maxHamming=1), 1 lib | built | skipped | d1 | 31 |
| Default (maxHamming=1), 2 libs | built | skipped | d1 | 31 each |
| `--crAssignMaxHamming 2`, small ref (29 guides) | built | built | d2 | 31 |
| `--crAssignMaxHamming 2`, large ref (246K barcodes) | built | built (~435M entries, warning emitted) | d2 | 31 |

## Files Modified

| File | Change |
|------|--------|
| `core/features/process_features/src/io.c` | Replace d2 budget gate with 500M warning |
| `core/legacy/source/PfMultiProcess.cpp` | Remove budget split, remove consumer caps and quantum divisor |

## Validation Plan

1. **Clean rebuild** (mandatory per AGENTS.md):
   ```
   make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
   make feature-barcodes-tools
   ```

2. **Smoke test — single feature library** (gRNA, 29 guides):
   Run with default settings. Verify logs show:
   - `Feature prehash clamp: requested max_hamming=2, assignment max_hamming=1, effective=1`
   - `Feature prehash enabled: <=1 only`
   - `Will use 31 threads to process the reads` (with 32 total threads)

3. **Smoke test — two feature libraries** (gRNA + LARRY):
   Run with `--crAssignMaxHamming 2`. Verify logs show:
   - No prehash clamp message (prehash=2 matches assignment=2)
   - `Feature prehash enabled: <=1+<=2` for both libraries
   - For LARRY: `WARNING: d2 prehash is large (N estimated entries)`
   - `Will use 31 threads to process the reads` for each library
   - `policy=full_budget_sequential` in library scheduler log

4. **MSK 30-KO benchmark rerun**: Compare wall time and parity against CR
   baselines per the plan in `docs/HANDOFF_MSK_2FEATURE_BENCHMARK_20260303.md`.

## Risk Assessment

- **Memory**: The d2 prehash for LARRY (246K barcodes, 20bp) builds ~435M
  entries. At ~12 bytes/entry (key + value + overhead in khash), this is
  roughly 5-6 GB. Machines running these workloads typically have 110+ GB RAM,
  so this is acceptable.

- **Consumer oversubscription**: 31 consumer threads each capable of spawning
  up to 4 OMP threads (for the rare reads that miss the prehash) could
  momentarily reach 124 threads on a 32-core machine. Since the OMP regions
  are sub-microsecond bursts and the OS scheduler handles this efficiently,
  this is not a practical concern.
