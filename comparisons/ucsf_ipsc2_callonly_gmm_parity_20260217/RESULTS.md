# UCSF iPSC2_1_AALG2 Call-Only GMM Parity (CR vs STAR)

Date: 2026-02-17

## Run 1: Original UCSF Full-Sample Check (Historical)

Branch: `noBCfix`  
Tool: `core/legacy/source/star_feature_call --call-only --compat-perturb`

Inputs:
- Cell Ranger filtered MEX:
  `/mnt/pikachu/ucsf-perturb-seq/iPSC2_1_AALG2_autoindex11044_run3_20260216/outs/per_sample_outs/iPSC2_1_AALG2_autoindex11044_run3_20260216/count/sample_filtered_feature_bc_matrix`
- Cell Ranger CRISPR calls:
  `/mnt/pikachu/ucsf-perturb-seq/iPSC2_1_AALG2_autoindex11044_run3_20260216/outs/per_sample_outs/iPSC2_1_AALG2_autoindex11044_run3_20260216/count/crispr_analysis`

CRISPR-only MEX extracted from CR filtered MEX:
- Features: `548`
- Matrix nonzeros: `45,248`

Results:
- `--min-umi 10`:
  - CR rows `6,119`, STAR rows `7,325`, common `6,119`
  - Set-equivalent calls on common rows: `6,118 / 6,119` (`99.9837%`)
- `--min-umi 3`:
  - CR rows `6,119`, STAR rows `7,325`, common `6,119`
  - Set-equivalent calls on common rows: `6,119 / 6,119` (`100.0000%`)

## Run 2: UCSF 1M Downsample, Clean-Master Feature-Only Check

Code/binary provenance:
- Source tree: `/mnt/pikachu/STAR-suite-master-pristine` (detached `origin/master`)
- Commit: `a470c8c6202cb202973d53961227c1bb82679f62`
- Binary: `/mnt/pikachu/STAR-suite-master-pristine/core/legacy/source/star_feature_call`

Inputs:
- CR filtered MEX:
  `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_20260217_064209/outs/filtered_feature_bc_matrix`
- CR CRISPR calls:
  `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_20260217_064209/outs/crispr_analysis`

Run command:

```bash
OUT_BASE="/storage/ucsf-2M/callonly_parity_masterclean_20260217_073728"
STAR_FEATURE_CALL_BIN="/mnt/pikachu/STAR-suite-master-pristine/core/legacy/source/star_feature_call" \
CR_MEX_DIR="/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_20260217_064209/outs/filtered_feature_bc_matrix" \
CR_CALLS_DIR="/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_20260217_064209/outs/crispr_analysis" \
OUT_BASE="$OUT_BASE" \
MIN_UMI_LIST="10 3" \
bash comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/run_callonly_gmm_parity.sh
```

CRISPR-only MEX extracted from CR filtered MEX:
- Features: `548`
- Barcodes: `7,346`
- Matrix nonzeros: `16,850`

Artifacts:
- `parity_minumi10.tsv`, `parity_minumi3.tsv`
- `star_feature_call_minumi10.log`, `star_feature_call_minumi3.log`
- `call_only_minumi10/crispr_analysis/*`
- `call_only_minumi3/crispr_analysis/*`

Row-policy normalization artifacts:
- STAR filtered to CR barcode set:
  - `/storage/ucsf-2M/callonly_parity_masterclean_20260217_073728/call_only_minumi10/crispr_analysis/protospacer_calls_per_cell.cr_row_filtered.csv`
  - `/storage/ucsf-2M/callonly_parity_masterclean_20260217_073728/call_only_minumi3/crispr_analysis/protospacer_calls_per_cell.cr_row_filtered.csv`

### Guide-Call Parity

`--min-umi 10`:
- CR rows: `5,767`
- STAR rows: `7,346`
- Common rows: `5,767`
- Set-equivalent calls on common rows: `5,457 / 5,767` (`94.6246%`)
- Real mismatches on common rows: `310`
- STAR-only rows: `1,579` (all `None`)

`--min-umi 3`:
- CR rows: `5,767`
- STAR rows: `7,346`
- Common rows: `5,767`
- Set-equivalent calls on common rows: `5,766 / 5,767` (`99.9827%`)
- Real mismatches on common rows: `1`
- STAR-only rows: `1,579` (`1,577` `None`, `2` non-`None`)

### Overall UMI Count Parity

Totals on CR barcode set (`n=5,767`, using `num_umis` total per row):
- CR total assigned UMIs: `1,447,384`
- STAR `--min-umi 10`: `1,445,234` (delta `-2,150`, `-0.148544%`)
- STAR `--min-umi 3`: `1,447,391` (delta `+7`, `+0.000484%`)

Totals on full STAR rows:
- STAR `--min-umi 10`: `1,445,234` (delta vs CR total `-2,150`)
- STAR `--min-umi 3`: `1,447,402` (delta vs CR total `+18`)

### Correlation Metrics (Per-Cell Total Assigned UMIs)

Correlation summary file:
- `/storage/ucsf-2M/callonly_parity_masterclean_20260217_073728/umi_correlation_summary.tsv`

Values:
- `--min-umi 10`, all CR rows (`n=5,767`): Pearson `0.999963`, Spearman `0.999940`
- `--min-umi 10`, STAR-called rows (`n=5,713`): Pearson `0.999971`, Spearman `0.999939`
- `--min-umi 10`, set-equivalent rows (`n=5,457`): Pearson `1.000000`, Spearman `1.000000`
- `--min-umi 3`, all CR rows (`n=5,767`): Pearson `1.000000`, Spearman `1.000000`
- `--min-umi 3`, STAR-called rows (`n=5,767`): Pearson `1.000000`, Spearman `1.000000`
- `--min-umi 3`, set-equivalent rows (`n=5,766`): Pearson `1.000000`, Spearman `1.000000`

### Runtime (Feature-Only Call Stage)

From `star_feature_call_minumi*.log`:
- `--min-umi 10`: `real 5.25s`
- `--min-umi 3`: `real 5.27s`

## Reproducibility

Primary script:
- `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/run_callonly_gmm_parity.sh`

The script writes outputs under `OUT_BASE` and emits:
- `parity_minumi*.tsv`
- `call_only_minumi*/crispr_analysis/*`
- `star_feature_call_minumi*.log`
