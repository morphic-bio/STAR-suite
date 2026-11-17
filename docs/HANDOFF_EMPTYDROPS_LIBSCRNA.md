# EmptyDrops (libscrna) Migration Handoff

## Summary

EmptyDrops_CR now uses **libscrna** directly in STAR core (scRNA + perturb),
with a legacy escape hatch for the old STAR implementation.

**Default behavior**
- `--soloCellFilter EmptyDrops_CR` uses libscrna backend.
- Legacy path is available via:
  - `--soloEmptyDropsLegacy yes`

## Implementation Notes

- Core integration: `core/legacy/source/SoloFeature_emptyDrops_libscrna.cpp`
- Switch logic: `core/legacy/source/SoloFeature_cellFiltering.cpp`
- Legacy toggle: `--soloEmptyDropsLegacy` (default `no`)
- OrdMag/ambient indexing aligned to STAR in:
  - `core/features/libscrna/src/OrdMagStage.cpp`

## Parity & Regression Tests

### A375 EmptyDrops_CR Parity (legacy vs libscrna)

**Script**: `tests/run_emptydrops_parity.sh`  
**Default dataset**: `/storage/A375/star_gex_features_cr_parity_20260129_174602/Solo.out/Gene/raw`  
**Outputs**: `tests/emptydrops_parity_output_YYYYMMDD_HHMMSS/`

Latest run (2026‑01‑31):
- legacy count = 1168
- libscrna count = 1168
- intersection = 1168
- diff = 0
- Simple filter parity: `nUMImax=25422`, `nUMImin=2542`, `nCellsSimple=1132`

Note added 2026-11-17. In the outputs of that run the log of the "legacy" arm
names the libscrna backend, so the run compared the libscrna backend with
itself. Both arms are re-calls with `--runMode soloCellFiltering`, which at the
time left the last barcode of the matrix out of the call. The script now asks
both arms for the same knee (`--soloEmptyDropsLegacyKnee yes`) and fails
unless each arm's `Log.out` names its backend. Result on the same matrix
after the last-barcode fix:
- legacy count = 1167
- libscrna count = 1168
- one barcode differs (GGGGTTAAGGTCTAGC, 948 UMIs: adjusted p 0.0214 in the
  legacy backend, 0.0085 in libscrna), within the tolerance of one
- Simple filter parity: `nUMImax=25422`, `nUMImin=2542`, `nCellsSimple=1132`

The libscrna backend is the default and the reference; the script's header
describes how the two backends adjust p-values.

### Flex EmptyDrops Parity (large MEX)

**Script**: `tests/run_flexfilter_parity.sh`  
**Dataset**: `/mnt/pikachu/prod-12-28/SC2300771/Solo.out/Gene/raw`  
**Gold summary**: `/mnt/pikachu/prod-12-28/SC2300771/per_sample/flexfilter_summary.tsv`  
**Outputs**: `tests/flexfilter_parity_output_YYYYMMDD_HHMMSS/`

Latest run (2026‑01‑31):
- 0 diffs vs gold summary across all key columns

## Notes

- If any parity issues arise, first check `summary_diff.tsv` in the test output.
- For debugging tail/ambient behavior, compare `EmptyDrops/emptydrops_results.tsv`
  between legacy vs libscrna outputs.

