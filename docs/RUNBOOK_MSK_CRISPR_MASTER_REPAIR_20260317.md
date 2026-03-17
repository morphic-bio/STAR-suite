# Runbook: MSK CRISPR Master Repair

**Date**: 2026-03-17  
**Status**: Proposed plan for review  
**Primary source**: `docs/HANDOFF_MSK_CRISPR_REGRESSION_DEBUG_20260317.md`

## Goal

Repair the MSK 30polyKO mixed-chemistry CRISPR path on the current `master`
checkout by restoring the missing `CR-Larry-perturb` behavior that is still
required for correctness, while avoiding a destructive overwrite of newer local
work already present in this tree.

## What We Know

- The paper benchmark regression is real and is not explained by GEX/index/poly-G
  differences. The bad run used the `master` binary and lost major chunks of the
  CR-Larry feature-assignment path.
- The current checkout is **not** plain `master`. It already contains local WIP
  in the same files, plus post-squash fixes on `master`.
- Because the worktree is dirty and the overlap is non-trivial, the repair
  should be a **manual backport/file-by-file overlay**, not a blind merge or
  cherry-pick of the 15 missing commits.
- GEX expression looks unaffected. The immediate problem is feature assignment
  and downstream CRISPR calling on the MSK mixed TRU/NXT/TRU run.

## Current Tree vs `CR-Larry-perturb`

Already present or partially restored in the current checkout:

- Tiered hash search entry points exist in `assignBarcodes.c`.
- Bootstrap offset finalization exists.
- A brute-force bootstrap fallback already exists in `assignBarcodes.c`.
- `PfMultiProcess.cpp` already sets bootstrap/anchor options and has env
  overrides for bootstrap debugging.
- Async feature assignment launch is already present in `STAR.cpp` /
  `PfMultiProcess.cpp`.
- Consumer hot-path fixes on `master` are already newer than the missing branch
  in some places.

Still clearly missing relative to `CR-Larry-perturb`:

- Bootstrap replay buffer + some tiered-path bookkeeping/diagnostics in
  `core/features/process_features/src/assignBarcodes.c`
- `use_hot_hash` / `skip_heatmaps` config plumbing across `pf_api`,
  globals, and `PfMultiAssign`
- The later `PfMultiMerge` / `PfMultiProcess` MEX optimization stack
  (`pruneZeroCountFeatures`, vector remap, streaming gzip, in-memory GEX/MEX,
  Phase A/B split result object)

## Merge Strategy

- Do **not** run a raw `git merge` into the current dirty tree.
- Do **not** replace entire files with branch versions.
- Backport by topic, reconciling each change against:
  - current local WIP
  - newer `master` consumer-path fixes
  - current CLI/API options such as `probeOnly` and `skipQcOutputs`

## Proposed Changes

### P0: Must Restore for MSK Correctness

- [ ] Backport the remaining bootstrap/tiered assignment logic from
  `CR-Larry-perturb` into
  `core/features/process_features/src/assignBarcodes.c`.
  Scope:
  - bootstrap replay buffer for reads that miss during bootstrap
  - `feature_mode_search_offsets_reset()`
  - guarded bootstrap histogram/finalize behavior from the branch
  - branch-side tiered-path safety helpers (`pf_bloom_maybe_*`,
    replay-drain path, tiered diagnostics where useful)
  Constraint:
  - preserve the newer consumer hot-path fixes already on `master`
  - preserve the current local brute-force bootstrap fallback unless testing
    shows the branch replay path supersedes it cleanly

- [ ] Backport the missing `process_features` config plumbing for
  `use_hot_hash` and `skip_heatmaps`.
  Files:
  - `core/features/process_features/include/globals.h`
  - `core/features/process_features/src/globals.c`
  - `core/features/process_features/include/pf_api.h`
  - `core/features/process_features/src/pf_api.c`
  - `core/features/process_features/src/main.c`
  - `core/features/process_features/src/mex_writer.c`
  - `core/features/process_features/src/assignBarcodes.c`
  - `core/legacy/source/PfMultiAssign.h`
  - `core/legacy/source/PfMultiAssign.cpp`
  Constraint:
  - keep current `probeOnly` and `skipQcOutputs` behavior on `master`
  - treat this as an overlay, not a rollback

- [ ] Reconcile `PfMultiAssign` option wiring with the current tree.
  Concrete changes:
  - wire `useHotHash` through to `pf_config_set_use_hot_hash()`
  - add `skipHeatmaps` and wire it through to `pf_config_set_skip_heatmaps()`
  - keep current `skipQcOutputs` support
  - keep current api-run summary fields and add branch fields only if they
    still match current permit-controller data

- [ ] Audit the mixed-chemistry namespace path for the MSK layout and patch only
  the deltas that are still missing on the current tree.
  Focus areas:
  - NXT PolyIII assignment against TRU GEX whitelist
  - `translateNxt` behavior through assign output and MEX merge
  - filtered-barcode normalization for mixed namespace runs
  - per-library `star_max_hamming`
  Rationale:
  - MSK is `GEX=TRU`, `gRNA=NXT`, `LARRY=TRU`
  - even if the big count regression is primarily assignment-side, we should
    re-verify this path before calling the repair complete

### P1: Likely Useful but Not Required for First Correctness Pass

- [ ] Bring over the low-risk CR-compat behavior from `CR-Larry-perturb` in
  `core/legacy/source/PfMultiProcess.cpp`.
  Candidates:
  - `assignOpts.skipHeatmaps = true`
  - the `kPfConsumersPerProducer` producer/consumer apportioning logic
  - any branch-only logging that helps validate mixed-chem assignment
  Note:
  - this is more about runtime/QC noise than the primary MSK correctness bug

### P2: Defer Unless Post-Fix Runtime Is Still Bad

- [ ] Defer the large MEX optimization stack until after correctness is back.
  Files:
  - `core/legacy/source/PfMultiMerge.cpp`
  - `core/legacy/source/PfMultiMerge.h`
  - `core/legacy/source/PfMultiProcess.cpp`
  - `core/legacy/source/PfMultiProcess.h`
  - `core/legacy/source/STAR.cpp`
  - `core/legacy/source/SoloFeature.h`
  - `core/legacy/source/SoloFeature_outputResults.cpp`
  - `core/legacy/source/SoloFeature_processRecords.cpp`
  Features:
  - vector-remap merge
  - streaming gzip
  - pre-prune feature MEX
  - in-memory GEX/MEX caching
  - fuller Phase A/B split result object
  Rationale:
  - these are valuable, but they are not the first thing blocking MSK CRISPR
    correctness

## Verification Gates

### Gate 0: Clean Rebuild First

Required before any regression validation:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

### Gate 1: Small/Targeted Validation

- Run a mixed-chemistry smoke before the full paper benchmark.
  Suggested surfaces:
  - `tests/multi_feature/run_msk_multifeature_smoke.sh`
  - `tests/multi_feature/run_mixed_chemistry_filtered_smoke.sh`
- Validate using:
  - feature `stats.txt`
  - merged filtered barcodes overlap
  - final CRISPR/LARRY outputs
- Do **not** use only `resolve_resolved` counters as the success criterion.
  The bootstrap handoff already showed those counters are not the right parity
  surface.

### Gate 2: Full MSK Benchmark

Run:

```bash
scripts/paper/run_msk_30polyko_benchmark.sh
```

Expected direction of change versus the bad `master` run:

- total feature counts move back toward the reference scale
  (`~100M`, not `~34M`)
- total deduped feature counts move back toward `~4.7M`
- gRNA-assigned cells recover dramatically
- zero-molecule filtered cells collapse from the catastrophic `28,550` level

### Gate 3: Regression Guard

After MSK repair:

- confirm GEX outputs are unchanged within expected tolerance
- confirm UCSF CRISPR path still behaves
- confirm mixed-chemistry filtered-barcode normalization still works

## Open Questions

- Is the replay buffer strictly required once the current local brute-force
  bootstrap fallback is present, or is it just the cleanest way to recover the
  original branch behavior?
- Do we want `useHotHash` enabled by default, opt-in only, or not at all for
  the first repair pass?
- Do we want to pull the `PfMultiMerge` optimization stack into the first repair,
  or keep the first pass tightly scoped to assignment correctness?

## Recommended Execution Order

1. Finish the `process_features` assignment-path backport.
2. Finish `PfMultiAssign`/`PfMultiProcess` option plumbing reconciliation.
3. Clean rebuild.
4. Run targeted mixed-chem validation.
5. Run the full MSK paper benchmark.
6. Only then decide whether to port the larger MEX optimization stack.
