# UCSF 2M Output-Namespace A/B Baseline (2026-02-20)

## Purpose

Establish a strict A/B baseline for the CR-compat output-namespace behavior on
UCSF 2M and confirm no feature/GEX regressions when switching output namespace.

The only changed argument between runs was:

- A: `--crOutputChemistry NXT`
- B: `--crOutputChemistry TRU`

All other STAR arguments, inputs, and references were identical.

## Runs

- NXT run:
  - `/storage/ucsf-2M/star_runs/star_abtest_iPSC2_1_AALG2_1M_nxt_outputchem_NXT_20260220_214005`
- TRU run:
  - `/storage/ucsf-2M/star_runs/star_abtest_iPSC2_1_AALG2_1M_nxt_outputchem_TRU_20260220_214846`

Both runs completed successfully.

## Direct Output Equivalence (STAR artifacts)

The following were identical between A/B:

- `outs/raw_feature_bc_matrix/matrix.mtx.gz`
  - dims/nnz/sum: `(39154, 16096, 71448, 1586203)`
- `outs/filtered_feature_bc_matrix/matrix.mtx.gz`
  - dims/nnz/sum: `(39154, 6069, 56772, 1352165)`
- CRISPR calls table cardinality:
  - `outs/crispr_analysis/protospacer_calls_per_cell.csv`
  - rows: `6069`, non-`None`: `5085`

Expected namespace artifact difference:

- NXT run: no namespace-conversion sidecars emitted
- TRU run: sidecars emitted (conversion active)
  - `outs/raw_feature_bc_matrix/barcodes.native.tsv.gz`
  - `outs/raw_feature_bc_matrix/barcodes.namespace_map.tsv.gz`
  - `outs/filtered_feature_bc_matrix/barcodes.native.tsv.gz`
  - `outs/filtered_feature_bc_matrix/barcodes.namespace_map.tsv.gz`

## Parity Reports vs Cell Ranger Baseline

CR baseline:

- `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813`

Parity runs (`--translate-side star`, historical UCSF parity mode):

- NXT parity:
  - `/tmp/ucsf2m_abtest_parity_nxt_20260220_215212/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
- TRU parity:
  - `/tmp/ucsf2m_abtest_parity_tru_20260220_215212/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`

Metrics most relevant to call/count consistency:

- `rows_star`: `6069` (both)
- `rows_star_non_none`: `5085` (both)
- `total_assigned_umis_star`: `1304633` (both)
- `total_assigned_umis_delta_star_minus_cr`: `-142751` (both)

Small differences in some parity fields were observed when translation is
applied despite TRU output already matching CR namespace (translation-mode
artifact), not from STAR output differences.

Sanity check without translation on TRU run:

- `/tmp/ucsf2m_abtest_parity_tru_notrans_20260220_215435/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`

## Conclusion (New Baseline)

UCSF 2M A/B confirms that switching CR-compat output namespace from NXT to TRU
introduces no regression in core STAR GEX/feature matrices or CRISPR call
cardinality under fixed inputs and parameters.

This A/B pair is the new baseline reference for output-namespace behavior.
