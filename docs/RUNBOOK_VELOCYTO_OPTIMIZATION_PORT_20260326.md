# Runbook: Porting Velocyto-Side Optimizations to New Paths

Date: 2026-03-26  
Branch intent: `feature/velocyto-optimizations-20260326`

Related docs:
- `docs/RUNBOOK_VELOCYTO_COUNT_RESOLUTION_20260326.md`
- `docs/RUNBOOK_GZIP_PIPELINE_ITERATION_20260326.md`
- `tests/ARTIFACTS.md` section `Flex Full-Scale Benchmark Artifacts`

## Goal

Capture the main optimization ideas that were validated on the Flex / no-align /
Velocyto-adjacent path so they can be ported deliberately into other paths
without re-learning the same lessons.

This runbook is about **how to port the optimization ideas**, not about exact
Velocyto counting semantics. Exact Velocyto semantics are documented separately
in `RUNBOOK_VELOCYTO_COUNT_RESOLUTION_20260326.md`.

## Core Idea

The important result from the recent work is not one isolated micro-change. It
is a sequence of structural simplifications:

1. stop materializing and regrouping the same information repeatedly
2. convert global hash / map traversals into linear CSR-style scans
3. emit final matrix structures directly instead of building temporary copies
4. fuse counting work into earlier passes when the information is already hot
5. keep replay / snapshot harnesses so each step can be validated in isolation

These ideas were validated on the full Flex `--flexNoAlign 1` surface and are
good candidates to reuse elsewhere.

## Validated Optimization Sequence

From `tests/ARTIFACTS.md`, the successful sequence was:

### 1. Replace regroup-heavy collapse with CB-first CSR layout

Main idea:

- count hash entries per CB first
- build one flat CSR-like layout
- do per-CB local work on contiguous spans

Why it helped:

- removed repeated unordered-map churn
- converted scattered work into linear memory access
- made the expensive phase proportional to contiguous spans rather than global
  hash regrouping

Observed effect on the Flex no-align snapshot harness:

- collapse total dropped from the older `extract/sort/group` shape to about
  `5.5 s`

### 2. Direct emit into the final inline matrix bundle

Main idea:

- do not rebuild composite barcode keys and triplets in a second pass
- emit `InlineMatrixBundle` contents directly during the CSR scan

Why it helped:

- removed dead post-regroup code
- eliminated extra sorts and extra temporary vectors
- reduced duplicated work between collapse and output preparation

Observed effect:

- `buildInlineMatrixFromHash` became removable dead code
- `collapseUMIall_fromHash` dropped further, from about `52.6 s` to `49.6 s`

### 3. Fuse pre-dedup counting into the first pass

Main idea:

- when Phase 1 already walks all entries, accumulate the counts needed later
  there instead of scanning the entire hash again

Why it helped:

- removed one full extra pass over the large inline-hash state
- improved cache locality by consuming the data while it is already hot

Observed effect:

- `collapseUMIall_fromHash` dropped again, to about `47.6 s`

### 4. Output-path cleanup

Main idea:

- use larger stdio buffers for MEX writing
- pass through already-built barcode strings / views instead of copying them

Why it helped:

- output code is usually not the dominant problem, but once the main collapse
  work is reduced, this cleanup becomes measurable

Observed effect:

- `processRecords` dropped modestly after the larger structural wins landed

### 5. Keep a snapshot / replay harness the whole time

Main idea:

- preserve a reproducible boundary where mapping can be skipped
- compare count-surface or full-output parity after each optimization

Why it mattered:

- made it possible to optimize collapse and output paths without paying full
  mapping cost every iteration
- prevented “faster but wrong” regressions

## The Main Porting Principles

When moving these ideas into another path, keep these rules:

### Prefer structural wins over micro-tuning

Good:

- fewer passes
- fewer regroups
- CSR / flat spans
- direct emit
- fused counting

Lower value unless the structure is already clean:

- tiny branch tweaks
- container-reserve fiddling without pass reduction
- output buffering before the main collapse path is fixed

### Avoid rebuilding the same representation twice

If a pass already knows:

- CB
- gene or feature id
- dedup outcome
- final barcode/sample slot

then the next question should be:

- can the final bundle be emitted here directly?

### Make per-CB work contiguous

The successful pattern was:

- global count per CB
- prefix sum
- scatter into flat storage
- local sort/dedup within each CB span

This is the pattern to try first in other large-cell-count paths.

### Fuse only after invariants are understood

The order matters:

1. first reproduce the old behavior cleanly
2. then replace regrouping with CSR/flat spans
3. then direct-emit
4. only then fuse earlier counting or metadata collection

If fusion happens too early, correctness becomes hard to reason about.

## What To Reuse First

For the new branch, the most reusable ideas are:

1. **CB-first CSR collapse**
2. **direct emit into final matrix/output bundles**
3. **fused pre-dedup counting**
4. **zero-copy / low-copy output plumbing**
5. **snapshot harnesses for count-surface replay**

These are more important than path-specific details like the exact Flex output
layout.

## Expected Port Order

Use this order unless a path-specific blocker forces a change.

### Phase 0. Baseline and instrumentation

- identify the current hot phase
- record wall, RSS, and per-phase timing
- add one clear replay boundary if none exists

### Phase 1. CSR / flat-layout conversion

- replace global regrouping with count -> prefix sum -> scatter
- keep behavior otherwise unchanged

### Phase 2. Per-span local dedup

- sort and dedup only within the per-CB span
- avoid a second global sort if possible

### Phase 3. Direct emit

- emit the final matrix/output bundle from the local scan
- remove dead intermediate bundle builders

### Phase 4. Fused counting

- pull later “just counting” passes back into an earlier full walk
- remove redundant full-state scans

### Phase 5. Output cleanup

- larger write buffers
- barcode / feature passthrough
- remove avoidable string rebuilding

## What Not To Mix In

Do not combine this work with:

- gzip-path tuning
- alignment changes
- EmptyDrops / cell-calling policy changes
- multimapper policy changes
- new output semantics

Those each need their own harness and attribution.

## Validation Standard

Each optimization step should preserve one of:

1. **count-surface replay parity**
2. **full-output diff parity**
3. both, when practical

For performance validation, always record:

- wall time
- peak RSS
- the specific phase timing that was supposed to improve

## Current Working Hypothesis

For the new branch, the likely highest-yield reuse is:

- take the **CSR collapse**
- then **direct emit**
- then **fused counting**

before doing any smaller cleanups.

That is the shortest path to reproducing the largest validated wins from the
Velocyto-side / Flex no-align optimization thread.

## First Concrete Step

Before editing code:

1. identify the target path to optimize
2. locate the equivalent of:
   - hash/materialized record store
   - regroup phase
   - final bundle construction
3. write down the exact pass structure
4. mark which passes can become:
   - CSR build
   - direct emit
   - fused count collection

Only after that should implementation begin.
