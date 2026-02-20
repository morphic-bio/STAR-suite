# GEX + Feature Parity Runbook

This runbook captures the reusable parity scripts for STAR vs Cell Ranger MEX
comparisons on:
1. Raw matrices
2. Raw matrices restricted to Cell Ranger filtered barcodes
3. Filtered-vs-filtered matrices

## Scripts

1. `scripts/run_gex_feature_parity_checks.sh`
2. `scripts/compare_barcode_sets.py`
3. `scripts/report_additional_parity_metrics.py`
4. Existing dependency: `tests/compare_feature_mex.py`

## Default Usage

```bash
bash scripts/run_gex_feature_parity_checks.sh
```

Defaults are set to the current UCSF full-sample benchmark paths:
1. STAR: `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_dynamic_off_full_20260218_203459`
2. CR: `/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804`

## Custom Usage

```bash
bash scripts/run_gex_feature_parity_checks.sh \
  --star-run /path/to/star/run \
  --cr-run /path/to/cellranger/run \
  --out-dir /path/to/output_dir
```

## Matrix Selection (Unique vs Rescue Policy Probes)

The runner supports selecting matrix basenames independently for CR and STAR.
This is useful when evaluating STAR multimapper policies (for example
`UniqueAndMult-Rescue.mtx`) against Cell Ranger `matrix.mtx`.

```bash
--cr-raw-matrix-basename matrix.mtx
--cr-filtered-matrix-basename matrix.mtx
--star-raw-matrix-basename UniqueAndMult-Rescue.mtx
--star-filtered-matrix-basename matrix.mtx
```

Notes:
1. `UniqueAndMult-*` files are typically available in raw output directories.
2. Filtered directories usually provide `matrix.mtx`, so Rescue probes should
   be interpreted primarily from raw and raw-restricted sections.

## Barcode Normalization

By default, the runner normalizes barcode namespaces using the CR9 translation
table and compares in a shared namespace (left-to-right mapping, i.e. col1->col2).
This is required for UCSF-style parity where Cell Ranger outputs are TRU-formatted
(column 2 namespace).

Useful flags:

```bash
--barcode-translation /home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz
--translation-direction left-to-right
--translate-side star
```

When CR barcodes are already TRU-formatted and STAR is NXT-formatted, prefer
`--translate-side star` (NXT -> TRU). Avoid remapping CR-side barcodes in this
case.

Disable translation-based normalization only for debugging:

```bash
--no-barcode-normalization
```

If CR has per-sample outputs, the runner auto-detects
`outs/per_sample_outs/*/count/sample_filtered_barcodes.csv` and uses it as the
filtered-cell barcode list. Otherwise it falls back to:
`outs/filtered_feature_bc_matrix/barcodes.tsv.gz`.

## Outputs

The runner writes:
1. `PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
   - six matrix parity sections (GEX + CRISPR features; raw, CR-filtered, filtered)
   - additional metrics section:
     - GEX Pearson/Spearman on common barcodes
     - feature-call parity from `protospacer_calls_per_cell.csv`
2. `FILTERED_BARCODE_SET_OVERLAP.txt`
   - normalized filtered barcode-set overlaps

Default output root:
`/tmp/gex_feature_parity_<timestamp>/`
