# UCSF STAR 1-Hamming vs Cell Ranger Analysis Scripts

This folder contains reproducible scripts used to analyze UCSF CRISPR feature-call
parity differences between STAR (`process_features`/`assignBarcodes`) and Cell Ranger.

## Scripts

1. `build_star_m1_delta_vs_cr.py`
- Builds the STAR 1-Hamming rescue proxy table by taking `STAR(m=1) - STAR(m=0)`
  per `(feature, barcode)` pair.
- Maps STAR NXT barcodes to TRU with the 2-column translation whitelist.
- Joins each delta pair to CR raw MEX counts and CR call-list membership.
- Outputs:
  - `STAR_M1_DELTA_VS_CR.tsv`
  - `STAR_M1_DELTA_VS_CR_SUMMARY.txt`

2. `classify_star_m1_cr_misses.py`
- Takes `STAR_M1_DELTA_VS_CR.tsv`, filters rows with `cr_raw_count == 0`, and
  classifies missing rows into:
  - `barcode_absent_from_cr_matrix`
  - `shift_to_partner_feature` (`_A` <-> `_B`)
  - `no_partner_signal_in_cr`
- Outputs:
  - `STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.tsv`
  - `STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.summary.txt`

3. `inspect_barcode_feature_totals.py`
- For selected TRU barcodes, reports STAR and CR raw CRISPR feature totals side
  by side (for discrepancy deep dives such as FOXD3).
- Outputs:
  - `BARCODE_INSPECTION.tsv`
  - `BARCODE_INSPECTION_SUMMARY.txt`

4. `build_star_exact_vs_cr.py`
- Compares STAR exact matches (`m=0`) directly against CR raw MEX and CR call list.
- Outputs:
  - `STAR_EXACT_VS_CR.tsv`
  - `STAR_EXACT_VS_CR_SUMMARY.txt`
  - `STAR_EXACT_VS_CR_MISSING_TOP.tsv`

5. `export_paper_miss_tables.py`
- Combines H0 and H1 misses into publication-ready tables.
- Inputs:
  - `STAR_EXACT_VS_CR.tsv`
  - `STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.tsv`
- Outputs:
  - `UCSF_STAR_vs_CR_MISSES_H0_H1.tsv`
  - `UCSF_STAR_vs_CR_MISSES_H0_H1_SUMMARY.tsv`
  - `UCSF_STAR_vs_CR_MISSES_H0_H1_PREVIEW.md`

6. `audit_cr_star_zero_cells.py`
- Audits `CR>0, STAR=0` barcodes on a shared barcode set (typically CR filtered
  barcodes) for GEX MEX comparisons.
- Supports optional side translation (`--translate-side`) to reproduce or
  diagnose namespace effects.
- Supports optional FASTQ R1 probing to confirm exact CB presence for STAR-zero
  rows (`--fastq-r1` + `--barcode-translation`).
- Outputs:
  - `CR_FILTERED_GEX_UMI_DELTA_RANKED.tsv`
  - `CR_FILTERED_GEX_UMI_DELTA_TOP{25,50,100}.tsv`
  - `SUMMARY.txt`
  - `STAR_ZERO_FASTQ_PROBE.tsv` (when FASTQ probe is enabled)


## Using Cell Ranger config.csv as the UCSF input source

Use `scripts/run_ucsf_full_compat_forward_rescue_guides.sh --cr-config <config.csv>` to
derive the exact GEX/guides FASTQ inputs and feature reference from an existing
Cell Ranger `config.csv`. This preserves the local pinned STAR `genomeDir` policy
while making the CR config the source of truth for library inputs.

The script writes a generated `pf_multi_config.from_cr.csv` into the run output
directory and records the original `cr_config`, `cr_gene_expression_reference`,
and `cr_gene_expression_chemistry` in `RUN_MANIFEST.txt`.

The same CR-config input mode is also available for the UCSF perturb Y-removal
runner:
- `scripts/run_ucsf_perturb_yremove_batch.sh --cr-config <config.csv>`

Smoke coverage for this path:
- `tests/test_ucsf_cr_config_input.sh`
- `tests/test_ucsf_batch_cr_config_input.sh`
- `tests/test_ucsf_batch_cr_config_tiny_smoke.sh`
- `tests/run_ucsf_cr_config_1m_smoke.sh`

## End-to-end wrapper

Use `scripts/run_ucsf_star1h_cr_analysis.sh` to run all three scripts with
canonical UCSF defaults.

## Notes

- Inputs are matrix-market directories produced by STAR and Cell Ranger.
- CR barcodes are normalized by stripping `-1` before joins.
- The STAR `m=1 - m=0` delta is a practical proxy for 1-Hamming rescued feature
  assignments, not a per-read alignment trace.

## Multimap Policy Evaluation (GEX)

For GEX parity, the main driver is:
- `scripts/run_gex_feature_parity_checks.sh`

It now supports explicit matrix basename selection for each side:
- `--cr-raw-matrix-basename`
- `--cr-filtered-matrix-basename`
- `--star-raw-matrix-basename`
- `--star-filtered-matrix-basename`

This is required when comparing STAR `Unique` vs `UniqueAndMult-*` outputs
against Cell Ranger behavior without rewriting analysis logic.

Low-level MEX compare helper:
- `tests/compare_feature_mex.py`
- Supports:
  - `--matrix-basename` (single basename for both sides)
  - `--matrix-basename-a` / `--matrix-basename-b` (side-specific)
  - Real-valued MatrixMarket entries (needed by `UniqueAndMult-*` matrices)

Example (raw GEX parity on CR filtered barcodes, STAR Rescue matrix):

```bash
bash scripts/run_gex_feature_parity_checks.sh \
  --star-run /path/to/star_run \
  --cr-run /path/to/cr_run \
  --out-dir /tmp/ucsf_multimap_rescue_check \
  --star-raw-matrix-basename UniqueAndMult-Rescue.mtx \
  --star-filtered-matrix-basename matrix.mtx \
  --cr-raw-matrix-basename matrix.mtx \
  --cr-filtered-matrix-basename matrix.mtx \
  --translate-side star \
  --translation-direction left-to-right
```

Interpretation note:
- `UniqueAndMult-*` matrices are emitted for unfiltered raw outputs.
- Filtered STAR matrices are typically `matrix.mtx` unless a custom export is
  produced, so filtered-vs-filtered sections should be interpreted accordingly.
