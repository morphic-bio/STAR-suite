# Runbook: Flex Integration with Modular Components

## Status: COMPLETED (2026-01-25)

This runbook has been executed. libflex now uses libscrna for EmptyDrops, OrdMag, and Occupancy stages.

### Summary of Changes
- **Stage 0**: Baseline created at `tests/baselines/flex_pre_modular_20260125_060551/`
- **Stage 1**: Added libscrna include paths to STAR and libflex Makefiles
- **Stage 2**: Verified API compatibility (headers match except for unused STAR-specific declarations)
- **Stage 3**: Removed OrdMagStage.o, EmptyDropsMultinomial.o, EmptyDropsCRSampler.o, OccupancyGuard.o from libflex.a (now provided by libscrna.a)
- **Stage 4**: All regression tests pass (flex_inline_test, nonflex_mex_smoke, flex_smoke)

### Current libflex.a Contents
```
FlexFilter.o
FlexFilterAdapters.o
FlexFilterIO.o
```

### Current libscrna.a Contents
```
EmptyDropsCRSampler.o
EmptyDropsMultinomial.o
OrdMagStage.o
OccupancyGuard.o
scrna_api.o
```

### Header Separation

Flex now includes libscrna headers explicitly via a `libscrna/` prefix (wrapper headers live in
`core/features/libscrna/include/libscrna/`). This avoids ambiguous includes from `flex/source/libflex`.

---

## Goal
Move Flex to the modular single-cell components (libscrna + process_features) while preserving Flex outputs and defaults.

## Scope
- Flex EmptyDrops/OrdMag/Occupancy should be sourced from `core/features/libscrna`.
- Flex feature calling stays behaviorally identical (seed defaults, thresholds, occupancy behavior).
- Non-Flex STARsolo behavior must not regress.

## References
- `plans/scrna_flex_plan.md` (staged libscrna integration)
- `plans/emptydrops_refactor_plan.md` (design/guardrails for EmptyDrops + Simple ED)
- `plans/compat_flex_handoff.md` (compat state + Flex modular goal)
- `docs/TODO_crispr_feature_calling.md` (crMinUmi regression matrix)
- `docs/todos` (Flex-related runtime checks)

## Preconditions
- Clean tree (no untracked build artifacts).
- Flex smoke fixtures available (see `tests/flex_smoke/run_flex_smoke.sh`).
- Capture and track artifacts in `tests/ARTIFACTS.md`.

## Stage 0: Baseline + Inventory
1. Freeze Flex outputs:
   - `flex/tools/flexfilter/test_smoke.sh`
   - `tests/flex_smoke/run_flex_smoke.sh`
   - Optional: `tests/run_flex_inline_test.sh`, `tests/run_flex_with_tags_test.sh`, `tests/run_flex_multisample_test.sh`
2. Store outputs in a dated baseline folder:
   - `tests/baselines/flex_pre_modular_<stamp>/`
3. Inventory Flex sources to be replaced by libscrna:
   - `flex/source/libflex/EmptyDropsCRSampler.*`
   - `flex/source/libflex/EmptyDropsMultinomial.*`
   - `flex/source/libflex/OrdMagStage.*`
   - `flex/source/libflex/OccupancyGuard.*`
   - `flex/source/libflex/CRLogProb.h`
   - `flex/source/libflex/SampleMatrixData.h`
   - `flex/source/libflex/pcg_*.hpp`, `flex/source/libflex/SimpleGoodTuring/`
4. Confirm Flex link flags (avoid `--whole-archive` conflicts when adding `libscrna.a`).

### 100K Flex Baselines (Stage 0 add-on)
- STAR-suite scripts (100K downsampled SC2300771):
  - `tests/flex_smoke/run_flex_smoke.sh` (uses `reference/tests/100K` by default)
  - `tests/run_flex_inline_test.sh`
  - `tests/run_flex_cbub_validation_test.sh`
  - `tests/run_flex_with_tags_test.sh`
  - `tests/run_flex_multisample_test.sh`
  - Fixture env defaults: `tests/external_fixtures_env.sh` and `tests/setup_parms_tests.sh`
- STAR-Flex scripts/docs (source-of-truth baselines):
  - `tests/run_nonflex_solo_regression_test.sh`
  - `reference/runSTAR.sh` (100K Flex run)
  - `reference/run100K_test.sh`
  - `flex/docs/TESTING_flex.md`
  - `reference/README_MEX_INTEGRATION.md`, `reference/FLEXFILTER_COMPLETE_SUMMARY.md`
  - Baseline outputs referenced under `/storage/100K/SC2300771/` and `/storage/downsampled_100K/SC2300771/`

## Stage 1: Wire libscrna into Flex build (no behavior change)
1. Build `core/features/libscrna/libscrna.a` standalone.
2. Update `flex/source/Makefile`:
   - Add `core/features/libscrna/include` to includes.
   - Add `libscrna.a` to link line (if not using `--whole-archive`).
3. Rebuild Flex and rerun baseline smoke tests.
4. Validate outputs identical to Stage 0 (hash/diff).

## Stage 2: Switch Flex to libscrna implementations
1. Add thin adapter headers in `flex/source/libflex/` that `#include` libscrna headers.
2. Update Flex orchestration files to include adapters:
   - `flex/source/libflex/FlexFilter.cpp`
   - `flex/source/libflex/FlexFilterIO.cpp`
   - `flex/source/libflex/FlexFilterAdapters.cpp`
3. Ensure config passthrough remains identical:
   - seed default = 1
   - Simple ED is fallback-only (unless explicitly enabled)
   - occupancy filtering remains enabled for Flex
4. Rerun Flex smoke tests; compare to Stage 0 baseline.

## Stage 3: Remove duplicate libflex sources
1. Remove duplicated Flex sources from build:
   - `EmptyDropsCRSampler.*`, `EmptyDropsMultinomial.*`, `OrdMagStage.*`, `OccupancyGuard.*`
   - `CRLogProb.*`, `SampleMatrixData.h`, `pcg_*.hpp`, `SimpleGoodTuring/`
2. Keep orchestration files only:
   - `FlexFilter.{h,cpp}`, `FlexFilterIO.{h,cpp}`, `FlexFilterAdapters.{h,cpp}`
3. Confirm Flex links cleanly against `libscrna.a` only (no duplicate symbols).

## Stage 4: Flex-specific regression checks
1. Validate `--soloAddTagsToUnsorted` behavior on Flex runs.
2. Verify CB/UB strict coupling:
   - Keep strict policy for Flex.
   - Relax for non-Flex (if required by `docs/todos`).
3. Run CRISPR feature calling regression matrix for `--crMinUmi`:
   - Use `docs/TODO_crispr_feature_calling.md` matrix.
   - Update assay recommendations in `tests/crispr_feature_calling_comparison_report.md`
     and `docs/feature_barcodes.md` after validation.
4. Re-run 100K Flex tests and compare to STAR-Flex baselines:
   - Compare output hashes and key logs against STAR-Flex outputs when available.
   - Record results in `tests/ARTIFACTS.md`.

## Stage 5: Documentation + Cleanup
1. Update `tests/ARTIFACTS.md` with new baseline paths.
2. Update `docs/todos` for completed items and remaining follow-ups.
3. Add a short summary doc (e.g., `docs/FLEX_MODULAR_INTEGRATION_SUMMARY.md`).

## Acceptance Criteria
- Flex builds and links against `libscrna.a` without duplicate symbol warnings.
- Flex smoke tests match Stage 0 baseline outputs.
- Flex occupancy behavior unchanged.
- No regressions in non-Flex STARsolo tests.

## Rollback Plan
- Revert to Stage 0 baseline (restore libflex sources in build).
- Remove libscrna link line and adapters.
- Re-run smoke tests to confirm original outputs restored.
