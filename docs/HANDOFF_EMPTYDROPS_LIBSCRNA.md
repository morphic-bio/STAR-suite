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

