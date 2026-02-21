# CR-Compat Parity Results: UCSF iPSC2_1_AALG2 Full Sample

**Date**: 2026-02-21
**Branch**: `core-alignment-threads-integration`
**Dataset**: UCSF iPSC2_1_AALG2 full sample (8 GEX lanes + 2 guide lanes)

## Summary

Three fixes bring STAR Solo CR-compat mode to near-identical cell calling
and expression quantification versus Cell Ranger 9:

| Fix | Component | Impact |
|-----|-----------|--------|
| Poly-G trimming | `ClipCR4::polyTail3p` | Gene Pearson 0.30 → 0.997 |
| EmptyDrops BH + umiMin | `libscrna` | 8 additional rescued cells |
| Bootstrap OrdMag | `libscrna` | Cell gap 716 → 39 (Jaccard 0.99) |

## Cell Calling Comparison

| Metric | STAR (before) | STAR (after) | CR9 |
|--------|--------------|--------------|-----|
| Filtered cells | 6,609 | **7,286** | 7,325 |
| Cell count gap | 716 (9.8%) | **39 (0.5%)** | — |
| OrdMag cutoff (UMIs) | 111 | **88** | ~80-85 |
| Bootstrap recovered_cells | — | **7,263** | ~7,325 |

## Barcode Overlap

| | Count |
|---|---:|
| STAR filtered | 7,286 |
| CR filtered | 7,325 |
| Both | 7,268 |
| Only STAR | 18 |
| Only CR | 57 |
| **Jaccard** | **0.9898** |

All 6,609 previously called cells are retained (proper superset).

## Expression Correlation (7,268 common barcodes)

| Metric | Value |
|--------|------:|
| Gene Pearson (all genes) | 0.9969 |
| Gene Spearman (all genes) | 0.9731 |
| Gene Pearson (filtered, 21,521 genes) | 0.9969 |
| Gene Spearman (filtered) | 0.9789 |
| Cell Pearson | 0.9987 |
| Cell Spearman | 0.9983 |

## Run Details

### STAR
- Binary: `/mnt/pikachu/STAR-suite/core/legacy/source/STAR`
- Output: `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_full_iPSC2_1_AALG2_forward_rescue_guides_bootstrap_20260221_010635/`
- Key params: `--clipAdapterType CellRanger4 --clip3pPolyG auto --soloCellFilter EmptyDrops_CR --soloCrMultimapRescue yes --soloStrand Forward --soloFeatures GeneFull`

### Cell Ranger 9
- Output: `/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804/`

### Reference genome
- Index: `/storage/autoindex_110_44/bulk_index`

## Technical Details

### Bootstrap OrdMag Algorithm
Matches CR9's `filter_cellular_barcodes_ordmag()`:
1. Generate log2-spaced range of candidate `recovered_cells` (1 to `maxExpectedCells`)
2. For each candidate, compute baseline UMI at rank `round(recovered * 0.01)`
3. Count cells within order of magnitude (UMI >= baseline/10)
4. Minimize loss `(filtered - recovered)^2 / recovered`
5. Bootstrap 100 samples for robustness; take mean of estimates
6. Second bootstrap pass: sample top-N with variance for final cell count

### Flex Isolation
Flex path (`FlexFilter.cpp`) is completely untouched:
- Fixed `nExpectedCells=3000`, `umiMin=500`
- Raw p-value gating (no BH)
- 10K MC simulations
- No bootstrap OrdMag

### Files Changed
- `core/features/libscrna/include/scrna_api.h` — `use_bootstrap` config field
- `core/features/libscrna/src/scrna_api.cpp` — bootstrap dispatch
- `core/legacy/source/SoloFeature_emptyDrops_libscrna.cpp` — non-Flex defaults
- `core/features/libscrna/include/EmptyDropsMultinomial.h` — BH correction param
- `core/features/libscrna/src/EmptyDropsMultinomial.cpp` — BH implementation
- `core/legacy/source/ClipCR4.cpp` — poly-G trimming
- `core/legacy/source/ClipCR4.h` — `clipPolyG` flag
- `core/legacy/source/ParametersClip.h` — poly-G parameter
- `core/legacy/source/ParametersClip_initialize.cpp` — poly-G resolution
- `core/legacy/source/ClipMate.h` — `cr4` nullptr init (segfault fix)
- `core/legacy/source/ClipMate_initialize.cpp` — allocate cr4 for type==11
- `core/legacy/source/Parameters.cpp` — register `clip3pPolyG`
- `core/legacy/source/parametersDefault` — `clip3pPolyG auto` default
