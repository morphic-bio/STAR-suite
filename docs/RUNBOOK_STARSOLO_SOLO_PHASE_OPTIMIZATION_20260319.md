# Runbook: STARsolo Solo-Phase Optimization

**Date**: 2026-03-19  
**Branch context**: `benchmark-flex`  
**Scope**: standard STARsolo `Gene` / `GeneFull` post-map path  
**Out of scope**: Flex inline-hash fast path, FlexFilter occupancy tuning, EmptyDrops semantics

## Goal

Optimize the standard STARsolo post-map "Solo phase" after mapping is complete.

The target is the classic STARsolo counting path used for scRNA-style
`Gene/GeneFull` processing:

1. `sumThreads()`
2. `countCBgeneUMI()`
3. `outputResults()`
4. `cellFiltering()`

The objective is to reduce wall time and memory churn without changing MEX
contents, filtered-cell outputs, or CLI behavior.

## Current Architecture

The main orchestration is in:

- `core/legacy/source/SoloFeature_processRecords.cpp`

For the standard non-inline path, the relevant flow is:

1. `SoloFeature::sumThreads()`
2. `SoloFeature::countCBgeneUMI()`
3. `SoloFeature::outputResults(false, ...raw/)`
4. `SoloFeature::cellFiltering()`

The important implementation files are:

- `core/legacy/source/SoloFeature_processRecords.cpp`
- `core/legacy/source/SoloFeature_countCBgeneUMI.cpp`
- `flex/source/SoloReadInfoLoader.cpp`
- `flex/source/SoloReadInfoLoader.h`
- `flex/source/SoloReadInfoSink.cpp`
- `flex/source/SoloReadInfoSink.h`
- `core/legacy/source/SoloFeature_collapseUMIall.cpp`
- `core/legacy/source/SoloFeature_outputResults.cpp`
- `core/legacy/source/MexWriter.cpp`
- `core/legacy/source/SoloFeature_cellFiltering.cpp`

## Source-Level Findings

These are the main optimization opportunities identified from code review.

### 1. `collapseUMIperCB()` is doing avoidable full-copy sorting

In `core/legacy/source/SoloFeature_collapseUMIall.cpp`, the old in-place
`qsort()` path is commented out and replaced with:

- allocate `std::vector<Triplet>`
- copy each `(gene, umi, read)` tuple out of `rGU`
- `std::sort()` the vector
- write everything back into `rGU`

This is in the hot collapse path and is the strongest obvious optimization
target.

### 2. `CountingSink` double-materializes the read records

The current standard path in `countCBgeneUMI()` does:

- parse packed records through `SoloReadInfoLoader`
- push `ReadInfoRecord` objects into `CountingSink::perWL`
- later allocate legacy `rGeneUMI` / `rCBp`
- copy all buffered records again into `rGeneUMI`
- then run `collapseUMIall()`

This is structurally similar to the extra materialization that was expensive on
the Flex side.

### 3. Per-record sink dispatch is more generic than the hot path needs

`SoloReadInfoLoader` uses:

- `using RecordSink = std::function<void(const ReadInfoRecord&)>;`

and calls it for every record in the main loader loop.

That abstraction is convenient, but it is expensive in a hot loop.

### 4. `CountingSink` does an unconditional conflict-guard hash lookup

For every record with a real `readIndex`, `CountingSink::onRecord()` checks
`readToCb` and enforces a conflict guard.

That looks useful for debugging and parity validation, but it is extra hot-path
hash work that may not need to run unconditionally in production.

### 5. MEX output still has streaming inefficiencies

There are two distinct output-side issues:

- `SoloFeature_outputResults.cpp` builds the Unique+Mult matrix body into a
  full `std::ostringstream` before writing it to disk.
- `core/legacy/source/MexWriter.cpp` writes one `fprintf()` per matrix entry,
  barcode, and feature row.

These are real cleanup opportunities, but they are probably secondary to the
collapse/materialization path.

## Strategic Decision

Do **not** start by parallelizing EmptyDrops or cell filtering.

For the standard STARsolo path, the most plausible wins are:

1. collapse hot-path cleanup
2. removing extra record materialization
3. reducing per-record sink overhead
4. then output-path cleanup

Only after these are measured and cleaned up should deeper parallelization of
post-map finalize be considered.

## Blast Radius Guide

Before making changes, classify them by how much they can affect Flex.

### Solo-only or effectively Solo-only

These are the safest starting points if the goal is to optimize standard
STARsolo without disturbing the current Flex inline-hash path:

- `core/legacy/source/SoloFeature_collapseUMIall.cpp`
- `core/legacy/source/SoloFeature_outputResults.cpp`
- `core/legacy/source/MexWriter.cpp`

Why:

- the standard non-inline post-map path uses these directly
- the Flex inline-hash fast path returns early before `outputResults()` in
  `SoloFeature_processRecords.cpp`
- `collapseUMIall()` is the legacy packed-array collapse path, distinct from
  `collapseUMIall_fromHash()`

Practical interpretation:

- `P1` is effectively Solo-only
- `P4` is effectively Solo-only

### Shared with Flex-adjacent plumbing

These files are used by the standard Solo path but also participate in shared
or Flex-adjacent flows:

- `flex/source/SoloReadInfoLoader.cpp`
- `flex/source/SoloReadInfoLoader.h`
- `flex/source/SoloReadInfoSink.cpp`
- `flex/source/SoloReadInfoSink.h`
- `core/legacy/source/SoloFeature_countCBgeneUMI.cpp`

Why:

- the standard packed-record counting path uses them
- Flex helper paths such as `prepareReadInfoOnly()` also use the loader/sink
  plumbing
- `countCBgeneUMI()` contains both the standard non-inline path and the
  `inlineHashMode` branch

Practical interpretation:

- `P2` is shared code and must be validated on both Solo and Flex-adjacent
  paths

### Not part of the Solo-only optimization track

These are explicitly out of scope for the current runbook:

- `collapseUMIall_fromHash()`
- Flex inline-hash record/aggregation paths
- hash-screen logic
- FlexFilter / occupancy / EmptyDrops tuning

If a proposed patch touches these surfaces, it is no longer a "Solo-only"
optimization patch and should be treated as a mixed-scope change.

## Optimization Priority

### P0: Add direct Solo-phase timers and counters

**Target**: make each optimization measurable on the real path before and after
patches.

Add explicit timers around these stages:

- `SoloFeature::sumThreads()`
- `SoloFeature::countCBgeneUMI()`
- loader parse loop only
- `CountingSink::finalize()`
- `runCliqueCorrection()` if enabled in this path
- `collapseUMIall()`
- `outputResults(false, raw/)`
- `cellFiltering()`

Inside `countCBgeneUMI()`, split timing further into:

- loader time
- sink ingest time
- sink finalize/materialize time
- collapse time

Also record structural counters:

- total records emitted by loader
- number of detected CBs
- total `rGeneUMI` triplets
- total per-CB collapse groups
- nnz written to raw MEX

Implementation notes:

- Prefer `Stats` or `logMain`-style timers already used elsewhere in STAR.
- Keep the default logging concise.
- If needed, gate verbose sub-phase timers behind an env var such as
  `STAR_DEBUG_SOLO_TIMERS=1`.

Acceptance:

- timing appears in `Log.out` or `Log.final.out`
- timers are deterministic enough to compare repeated 100K runs
- counters line up with existing raw MEX / summary outputs

### P1: Fix the collapse hot path first

**Blast radius**: `Solo-only / low Flex risk`

**Target**: remove the `Triplet` copy/sort/writeback path in
`collapseUMIperCB()`.

Current problem:

- per-CB extra heap allocation
- full copy out of `rGU`
- full copy back into `rGU`
- extra cache traffic in the hottest part of collapse

Recommended patch order:

1. Restore a true in-place sort for the existing `rGU` layout.
2. If `qsort()` is undesirable, replace it with a zero-copy typed view over the
   underlying buffer, not a copied vector.
3. Preserve the current effective ordering semantics for downstream code.

Do **not** combine this with unrelated algorithmic changes in the first patch.

Acceptance:

- exact raw MEX parity on a representative 100K fixture
- exact filtered MEX parity
- measurable reduction in `collapseUMIall()` time from P0 timers

### P2: Remove double materialization from the counting path

**Blast radius**: `Shared with Flex-adjacent plumbing`

**Target**: reduce memory traffic and copy volume in the standard loader path.

Current path:

1. loader parses records
2. `CountingSink::perWL[cbIdx].push_back(rec)`
3. `finalize()` allocates `rGeneUMI`
4. `finalize()` copies per-WL vectors into `rGeneUMI`
5. `collapseUMIall()` consumes `rGeneUMI`

Recommended implementation path:

#### P2a: cheap cleanup first

- Gate the `readToCb` conflict guard behind a debug/parity flag.
- Keep it available for validation runs.
- Remove it from the default hot path.

#### P2b: replace `std::function` sink dispatch

- Replace `RecordSink = std::function<...>` with a templated or specialized
  sink path.
- Keep the public behavior the same.
- Do not over-abstract the hot loop.

#### P2c: flatten the buffering model

Preferred direction:

- first pass: count records per detected CB
- allocate contiguous target storage once
- second pass: fill final packed `rGeneUMI` directly

Alternative direction:

- keep one flat vector of packed records and an index over CB boundaries
- avoid nested `vector<vector<ReadInfoRecord>>`

The important point is to stop buffering every record as a heavy object and
then copying it into the actual collapse layout.

Acceptance:

- exact raw/filtered MEX parity
- lower peak RSS or lower loader/finalize time
- no changes to CB/UB or multimapper semantics
- no Flex regression on any path that uses `SoloReadInfoLoader` /
  `CountingSink` / `prepareReadInfoOnly`

### P3: Flatten expensive nested containers in collapse-only paths

**Blast radius**: `Mostly Solo-only, but depends on correction/multimap usage`

**Target**: reduce allocator and hash-map churn inside `collapseUMIall()`.

Candidate surfaces:

- `umiGeneMapCount`
- `umiGeneMapCount0`
- `umiCorrected`

These appear in the multi-gene / correction subpaths and likely contribute
substantial allocator overhead on large runs.

Recommended approach:

- instrument these paths first with P0 timers and counters
- if hot, replace nested `unordered_map<umi, unordered_map<gene, count>>`
  patterns with flatter packed-key structures
- keep semantics unchanged for MultiGeneUMI and correction modes

This is a second-wave optimization, not the first patch.

Acceptance:

- exact MEX parity under the affected multimapper/correction settings
- measurable improvement in the specific timed sub-phase

### P4: Optimize MEX output after collapse/materialization cleanup

**Blast radius**: `Solo-only / low Flex risk`

**Target**: reduce output-side overhead once the main compute path is cleaned up.

Two concrete patches:

#### P4a: remove full-file `ostringstream` staging

In `SoloFeature_outputResults.cpp`, the Unique+Mult path currently builds the
entire matrix body in memory before writing it.

Replace this with:

- direct streaming to the output file, or
- buffered chunk writes using a reusable byte buffer

Do not change output ordering.

#### P4b: buffer the shared MEX writer

In `core/legacy/source/MexWriter.cpp`, replace per-line `fprintf()` loops with
buffered batched writes.

This is low-risk and transferable to multiple paths, but it is probably not the
main Solo win by itself.

Acceptance:

- exact file parity for raw and filtered MEX outputs
- lower `outputResults()` time from P0 timers

### P5: Only then consider parallelizing post-map finalize

**Blast radius**: `Mixed / depends on implementation`

**Target**: improve scaling after obvious serial inefficiencies are removed.

Possible directions:

- parallelize collapse by CB ranges after `sumThreads()`
- shard output matrix assembly by CB ranges
- parallelize parts of `sumThreads()`

Do **not** start here. First remove the single-thread waste in P1-P4.

## Validation Matrix

Every optimization patch should use the same validation ladder.

### Level 1: build + smoke

Clean rebuild first:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
```

Smoke:

```bash
tests/run_solo_smoke.sh
```

### Level 2: raw MEX parity on a representative fixture

Use an existing 100K-scale scRNA fixture when available.

Recommended starting point:

- `tests/run_cr_parity_100k.sh`

Primary checks:

- raw MEX exact diff or high-precision compare
- filtered MEX exact diff or high-precision compare
- barcode count stability
- gene-level Pearson / barcode Jaccard if comparing against CR-format outputs

For shared-plumbing patches (`P2`), also run a Flex-adjacent validation surface
before merging:

- a Flex smoke or existing internal Flex regression harness
- any path that exercises `SoloReadInfoLoader` / `CountingSink`
- confirm no CB/UB tag regressions if sorted BAM tag injection is relevant

### Level 3: runtime measurement

Run before/after on the same binary configuration and same fixture:

- same thread count
- same output settings
- no concurrent benchmark jobs on the machine

Capture:

- total wall time
- P0 sub-phase timers
- peak RSS if available

### Level 4: optional larger-scale confirmation

After a patch is clean on 100K:

- run a larger fixture
- confirm that the same phase stays improved
- confirm no raw/filtered output drift

## Recommended Patch Sequence

This is the recommended order for a coding agent.

1. Add P0 timers and counters.
2. Patch P1 only: remove `Triplet` copy/sort/writeback.
3. Validate parity and benchmark.
4. Patch P2a and P2b: guard conflict map, remove `std::function` sink overhead.
5. Validate parity and benchmark.
6. Patch P2c: flatten buffering / direct materialization.
7. Validate parity and benchmark.
8. Patch P4a/P4b for output cleanup.
9. Only then consider P3/P5.

If the goal is to avoid Flex risk entirely, stop after:

1. `P0`
2. `P1`
3. `P4`

and defer `P2` until a broader validation window is available.

## What Not To Change In Early Patches

Do not mix these into the first optimization patch:

- EmptyDrops logic
- cell-calling thresholds
- multimapper resolution semantics
- CB/UB validity rules
- CR-compat output formatting
- Flex inline-hash code paths

The first Solo-phase optimization patches should be pure performance work.

## Success Criteria

This optimization branch is successful if it achieves all of the following:

- exact raw MEX parity
- exact filtered MEX parity
- no filtered-cell regressions
- measurable reduction in Solo post-map wall time
- lower memory churn in the counting/collapse path
- no CLI or output-format regression

## Short Version For The Coding Agent

If handing this to another agent, the short instruction is:

1. Instrument `sumThreads`, loader ingest, `CountingSink::finalize`,
   `collapseUMIall`, `outputResults`, and `cellFiltering`.
2. Fix the obvious regression in `collapseUMIperCB()` first.
3. Then remove the `CountingSink` double-buffering and hot-path generic sink
   overhead.
4. Only after that, optimize MEX output.
5. Validate every step on `run_solo_smoke.sh` plus a 100K parity fixture.
