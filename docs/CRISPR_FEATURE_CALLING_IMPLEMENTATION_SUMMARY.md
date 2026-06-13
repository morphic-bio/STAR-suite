# CRISPR Feature Calling Implementation Summary

**Date:** January 24, 2026  
**Status:** Complete (pending regression testing for non-CRISPR assays)

## Overview

Implemented automatic CRISPR feature calling in STAR's CR-compat mode (`--pfMultiConfig`). When processing multi-configs with CRISPR Guide Capture libraries, STAR now automatically generates `crispr_analysis/*` outputs matching Cell Ranger's format. An optional ambient-FDR guide caller is also available as a secondary STAR-suite output; GMM remains the default compatibility caller.

## Key Results

### Validation (A375 Dataset)

| Metric | Value |
|--------|-------|
| STAR cells | 1,167 |
| CR cells | 1,083 |
| Common barcodes | 1,083 |
| **Exact matches** | **1,083 (100.0%)** |

**UMI Thresholds:**
- RAB1A-2: STAR=225, CR=217
- Non_Target-1_MS: STAR=33, CR=32

## Changes Made

### 1. Code Changes

#### `core/legacy/source/star_feature_call.cpp`
- Kept default `--min-umi` at 3 for general use
- A375 parity workflows explicitly pin `--min-umi 10`

#### `core/legacy/source/PfMultiProcess.cpp`
- Added `runCrisprFeatureCalling()` function
- Filters merged MEX to CRISPR Guide Capture features only
- Runs GMM calling with configurable min_umi threshold
- Writes to `outs/crispr_analysis/`
- Optional `--crGuideCaller ambient-fdr|both` path writes ambient-FDR outputs under `outs/crispr_analysis/ambient_fdr/`

#### `core/legacy/source/Parameters.h`
- Added `crMinUmi` field to `pfMulti` struct
- Added ambient-FDR guide-caller fields (`crGuideCaller`, `crGuideFdr`, `crGuideFdrMinUmi`, `crGuideFdrEmitQvalues`)

#### `core/legacy/source/Parameters.cpp`
- Added `--crMinUmi` parameter (default: 3)
- Added optional ambient-FDR guide-caller parameters; current CR-compatible defaults remain `--crGuideCaller gmm`

#### `core/legacy/source/Makefile`
- Added `libprocess_features.a` and `libscrna.a` as STAR dependencies
- Added `-I$(PROCESS_FEATURES_INC)` to include path

### 2. New Parameter

**`--crMinUmi N`** - Minimum UMI threshold for CRISPR feature calling

- Default: 3 (general default for perturb/lineage contexts)
- A375 CR-parity workflows: pin to 10
- For lineage barcodes/stable features: use 2-3
- For FLEX probes: assay-specific validation still required

### 3. Documentation

| File | Description |
|------|-------------|
| `docs/todos` | Top-level TODO with regression testing notes |
| `docs/TODO_crispr_feature_calling.md` | Detailed regression test plan |
| `docs/feature_barcodes.md` | User-facing documentation with assay recommendations |
| `tests/crispr_feature_calling_comparison_report.md` | Validation report and implementation details |

### 4. Test Scripts

| File | Description |
|------|-------------|
| `tests/test_cr_compat_crispr_calling.sh` | Integration test for CR-compat CRISPR calling |

## Output Files (CR-compat mode)

When CRISPR Guide Capture features are present:

```
<outPrefix>/outs/
├── raw_feature_bc_matrix/          # Combined GEX + features
├── filtered_feature_bc_matrix/     # Filtered by GEX EmptyDrops
└── crispr_analysis/                # NEW: CRISPR feature calls
    ├── protospacer_calls_per_cell.csv
    ├── protospacer_calls_summary.csv
    ├── protospacer_umi_thresholds.csv
    ├── protospacer_umi_thresholds.json
    └── ambient_fdr/                # optional when --crGuideCaller ambient-fdr|both
        ├── guide_fdr_calls_per_cell.csv
        ├── guide_fdr_summary.json
        ├── guide_ambient_rates.tsv
        ├── guide_qvalues.mtx
        ├── guide_qvalues_barcodes.tsv
        └── guide_qvalues_features.tsv
```

Ambient-FDR estimates guide ambient rates from raw non-cell barcodes and writes
sparse q-values only for observed guide counts in filtered cells. Missing
cell-guide entries imply q-value 1 and are not materialized. BH correction uses
`n_filtered_cells * n_guides` tests.

`guide_fdr_calls_per_cell.csv` carries both the statistical call and the UMI
support needed for downstream manual tuning: `num_umis`, `min_called_umi`,
`max_called_umi`, `min_qvalue`, `call_status`, `default_fdr`, and `caller`.
`scripts/integrate_feature_library.py --ambient-fdr-calls-csv` imports these
fields into h5ad `obs` as `guide_fdr_*` columns, and
`scripts/build_multiome_mudata.py` propagates those columns to top-level
MuData `obs`.

## Assay-Specific Recommendations

| Assay Type | Recommended `--crMinUmi` | Status |
|------------|--------------------------|--------|
| CRISPR Guide Capture | 3 (default) | ✓ Validated baseline |
| A375 CR-parity fixture | 10 (override) | ✓ Validated |
| Lineage Barcodes | 2-3 | **Needs Testing** |
| FLEX Probes | 10 (TBD) | **Needs Testing** |
| Antibody Capture | TBD | **Needs Testing** |

**Important:** Default `3` is the general STAR-suite setting. Use higher values only for assay-specific parity targets (for example, A375 CR guide parity at `10`).

## Pending Work

See `docs/todos` for the full list:

1. **Regression Testing Required**
   - Test `--crMinUmi` with FLEX probe sets
   - Test with lineage barcode datasets at thresholds 2, 3, 5, 10
   - Update assay recommendations table after testing

2. **libflex Migration**
   - Change libflex to use modular libscrna
   - Include minUMI regression testing in migration

## How to Test

```bash
# Run the integration test
bash tests/test_cr_compat_crispr_calling.sh

# Or manually with custom threshold
STAR ... --pfMultiConfig config.csv --crMinUmi 3
```

## Files Modified (Summary)

```
core/legacy/source/
├── PfMultiProcess.cpp      # CRISPR calling integration
├── Parameters.cpp          # --crMinUmi parameter
├── Parameters.h            # crMinUmi field
├── Makefile               # libprocess_features linkage
└── star_feature_call.cpp  # Default min-umi change

docs/
├── todos                              # Top-level TODOs
├── TODO_crispr_feature_calling.md     # Regression test plan
├── feature_barcodes.md                # User documentation
└── CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md  # This file

tests/
├── crispr_feature_calling_comparison_report.md  # Validation report
└── test_cr_compat_crispr_calling.sh             # Integration test
```

## Review Checklist

- [x] CRISPR-only MEX filtering implemented
- [x] GMM calling with configurable min_umi
- [x] CR-compat output format (`crispr_analysis/*`)
- [x] 100% concordance with Cell Ranger on A375
- [x] `--crMinUmi` parameter exposed
- [x] Documentation updated
- [x] Integration test created
- [x] TODO items for regression testing added
- [ ] Regression testing for lineage barcodes (pending)
- [ ] Regression testing for FLEX (pending)

## Related NB-EM / SCEPTRE-Parity Results

See `tests/nbem_sceptre_parity_results.md` for NB‑EM / Poisson‑EM / SCEPTRE‑parity results on A375 and vignette datasets, including test set locations and outputs.
