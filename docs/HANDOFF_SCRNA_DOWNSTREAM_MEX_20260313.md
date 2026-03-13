# Handoff: scRNA Downstream Analysis from MEX Output

Date: 2026-03-13
Branch: `scrna-downstream`
Base: `master` at `4b7dc99`

## Objective

Start downstream single-cell analysis from the finished STAR-suite MEX outputs,
without rerunning alignment or feature calling.

Primary dataset for the next agent:

- UCSF perturb full run:
  `/mnt/pikachu/ucsf-perturb-yremove_20260311_224631`

This run is complete end to end and is the best current starting point for
downstream work on real scRNA/perturb outputs.

## Current Status

The upstream UCSF run is finished.

Final state:

- samples completed: `20 / 20`
- failed samples: `0`
- incomplete samples: `0`
- Y/noY transfer cleanup completed for all samples

Per-sample downstream outputs were kept locally on `/mnt/pikachu`.
Only BAM/YFASTQ payloads were transferred/cleaned.

## Primary Output Layout

For each sample:

```text
/mnt/pikachu/ucsf-perturb-yremove_20260311_224631/samples/<sample>/run/outs/
  raw_feature_bc_matrix/
    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz
  filtered_feature_bc_matrix/
    barcodes.tsv.gz
    features.tsv.gz
    matrix.mtx.gz
  crispr_analysis/
    protospacer_calls_per_cell.csv
    protospacer_calls_summary.csv
    protospacer_umi_thresholds.csv
    protospacer_umi_thresholds.json
```

The 20 sample IDs are:

- `EBs1_1_AALG1`
- `EBs1_1_AALG2`
- `EBs1_2_AALG1`
- `EBs1_2_AALG2`
- `EBs1_3_AALG1`
- `EBs1_3_AALG2`
- `EBs1_4_AALG1`
- `EBs1_4_AALG2`
- `EBs1_5_AALG1`
- `EBs1_5_AALG2`
- `EBs2_1_AALG1`
- `EBs2_1_AALG2`
- `EBs2_2_AALG1`
- `EBs2_2_AALG2`
- `EBs2_3_AALG1`
- `EBs2_3_AALG2`
- `EBs2_4_AALG1`
- `EBs2_4_AALG2`
- `EBs2_5_AALG1`
- `EBs2_5_AALG2`

## Matrix Schema

These are standard gzipped 10x-style MEX outputs.

Important details:

- barcodes are `TRU` namespace
- barcodes keep the `-1` suffix
- `features.tsv.gz` has three columns:
  - feature ID
  - feature name
  - feature type

For the UCSF perturb run, the filtered matrix contains:

- `38606` `Gene Expression` rows
- `548` `CRISPR Guide Capture` rows

Representative examples:

- gene rows at the top of `features.tsv.gz`
- guide rows at the bottom of `features.tsv.gz`

So downstream code should split on the third column rather than assume the
matrix is gene-only.

## Recommended Starting Surface

Use:

- `filtered_feature_bc_matrix/` for ordinary downstream cell-level analysis
- `raw_feature_bc_matrix/` only if you need to revisit cell filtering, ambient
  behavior, or custom thresholding

Use `crispr_analysis/protospacer_calls_per_cell.csv` as the guide-call table
for downstream perturb analyses. It has one row per filtered barcode.

## Quick Size Summary

Filtered barcode counts per sample:

```text
EBs1_1_AALG1  4519
EBs1_1_AALG2  7752
EBs1_2_AALG1  4492
EBs1_2_AALG2  6724
EBs1_3_AALG1  4111
EBs1_3_AALG2  7877
EBs1_4_AALG1  4324
EBs1_4_AALG2  7814
EBs1_5_AALG1  4528
EBs1_5_AALG2  4782
EBs2_1_AALG1  4426
EBs2_1_AALG2  8093
EBs2_2_AALG1  4229
EBs2_2_AALG2  8531
EBs2_3_AALG1  4406
EBs2_3_AALG2  7293
EBs2_4_AALG1  4301
EBs2_4_AALG2  5305
EBs2_5_AALG1  4453
EBs2_5_AALG2  7433
```

`protospacer_calls_per_cell.csv` has the same row count as the filtered barcode
file for each sample.

## Existing Repo Material To Reuse

For reading/comparing MEX:

- `scripts/run_ucsf_star1h_cr_analysis.sh`
- `scripts/ucsf_parity/build_star_exact_vs_cr.py`
- `scripts/ucsf_parity/build_star_m1_delta_vs_cr.py`
- `scripts/ucsf_parity/classify_star_m1_cr_misses.py`
- `scripts/ucsf_parity/inspect_barcode_feature_totals.py`

For UCSF parity context:

- `docs/UCSF_PARITY_RECOVERY_RUNBOOK_20260218.md`
- `docs/CR_CONFIG_COMPATIBILITY_MATRIX_20260312.md`
- `docs/CR_COMPAT_PARITY_RESULTS_UCSF_iPSC2_20260221.md`

For existing downstream-style MEX artifacts and examples:

- `tests/ARTIFACTS.md`
- `tests/sceptre_example_*_mex`
- `tests/crispat_vignette_mex/`
- `tests/nbem_crispat_vignette_output/`

These are useful if the next agent wants to compare STAR-suite MEX outputs to
existing downstream consumers or reuse prior MEX-loading patterns.

## Recommended First Steps For The Next Agent

1. Pick one UCSF sample first, not all 20.
   Start with `EBs1_1_AALG1` or `EBs2_5_AALG2`.

2. Load `filtered_feature_bc_matrix/` in the downstream framework of choice
   and confirm:
   - barcode import is correct
   - feature type split is correct
   - gene rows and guide rows are separated cleanly

3. Join `crispr_analysis/protospacer_calls_per_cell.csv` to the filtered
   barcode table using `cell_barcode`.

4. Only after that, decide how to merge samples:
   - keep sample ID as metadata
   - do not strip `-1`
   - do not assume barcodes are globally unique across samples without a sample
     prefix

5. If a CR comparison is needed, use the UCSF parity scripts above rather than
   writing new MEX readers from scratch.

## Important Constraints

- Do not rerun the UCSF upstream job just to start downstream analysis.
  The outputs are already complete.
- Treat the MEX outputs as the source of truth for downstream work.
- For perturb analyses, do not infer guide calls from the matrix alone before
  checking whether `crispr_analysis/` already gives the needed call surface.

## Non-Goals For The Next Agent

Unless the new task explicitly asks for it, the next agent should not start by:

- rerunning STAR
- rerunning Y removal
- rebuilding the UCSF reference
- re-deriving Cell Ranger compatibility inputs

The immediate work is downstream use of the finished MEX outputs.

## Follow-Up Findings

Additional downstream and velocyto-specific findings from the follow-up review
are documented here:

- `docs/HANDOFF_SCRNA_DOWNSTREAM_MEX_VELOCYTO_FINDINGS_20260313.md`
