# UCSF Parity Recovery Runbook (2026-02-18)

Date: 2026-02-18
Owner handoff: next agent

## Goal
Recover the previously achieved UCSF call-level parity baseline and stop using
ad hoc mixed-metric comparisons as pass/fail.

Default policy reminder:
- Global STAR-suite default for CRISPR feature calling is `min-umi/crMinUmi = 3`.
- A375 parity workflows are the explicit exception and pin `10`.

Primary target is **call-only parity** (not mixed GEX+feature matrix overlap).

Secondary target is **GEX parity correlation** on normalized/common barcodes.

## Step 0: Barcode Namespace / Whitelist Preflight (Do First)
UCSF parity must normalize barcode namespaces before any GEX/feature overlap
judgement.

Verified local whitelist sources:
- 2-column translation whitelist:
  `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz`
- 1-column TRU whitelist:
  `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt`

Verified facts (2026-02-18):
- Translation file has `3,686,400` rows and exactly 2 columns on every row.
- TRU whitelist has `3,686,400` rows.
- Translation column 2 matches TRU whitelist namespace (all rows).
- Translation direction for normalization:
  - `left-to-right` = NXT -> TRU
  - `right-to-left` = TRU -> NXT
- Example verified pair in translation whitelist:
  - NXT: `AAACCCATCCAGAAAG`
  - TRU: `AAACCCAAGCAGAAAG`

Observed namespace behavior on UCSF artifacts:
- CR9 filtered barcodes
  (`/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/filtered_feature_bc_matrix/barcodes.tsv.gz`)
  are TRU-formatted (`7346/7346`, 100% hit on translation column 2 after
  stripping `-1`).
- STAR NXT feature outputs
  (`/storage/ucsf-2M/star_runs/star_baseline_iPSC2_1_AALG2_1M_nxt_20260217_160217/Solo.out/GeneFull/filtered/barcodes.tsv`)
  are primarily NXT-formatted (`6677/6816`, 97.9607% hit on translation
  column 1).
- Raw barcode outputs can be mixed namespace, so always detect namespace from
  the actual barcode files used in that run.

Preflight command (run per input/output barcode file before comparisons):

```bash
python3 - <<'PY' /path/to/barcodes.tsv[.gz] \
  /home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz
import gzip, sys
bc_path, tr_path = sys.argv[1], sys.argv[2]
op = gzip.open if bc_path.endswith(".gz") else open
left, right = set(), set()
with gzip.open(tr_path, "rt") as f:
    for line in f:
        p = line.split()
        if len(p) >= 2:
            left.add(p[0]); right.add(p[1])
bcs = set()
with op(bc_path, "rt") as f:
    for line in f:
        bc = line.strip()
        if not bc:
            continue
        if bc.endswith("-1"):
            bc = bc[:-2]
        bcs.add(bc)
lh = sum(1 for b in bcs if b in left)
rh = sum(1 for b in bcs if b in right)
print(f"unique={len(bcs)} left_hits={lh} ({100*lh/len(bcs):.4f}%) right_hits={rh} ({100*rh/len(bcs):.4f}%)")
PY
```

Normalization rule for this runbook:
- If one side is mostly NXT and the other is mostly TRU, normalize to TRU using
  translation `left-to-right`.
- Use `--translate-side star` when STAR is NXT and CR is already TRU.
- Do not use `--translate-side cr` when CR is TRU: `577,920` barcodes are
  present in both translation columns, so CR can be remapped unintentionally.
- If both sides already match the same namespace, do not translate.

## Known-Good Reference (Must Match)
Source of truth:
- `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/RESULTS.md`

Expected from compatibility-gate rerun:
- `--min-umi 3`: set-equivalent `5766/5767` (`99.9827%`), one real mismatch.
- `--min-umi 3`: total UMI delta on common rows `+7` vs CR.
- `--min-umi 3`: Pearson/Spearman ~`1.0`.

This is the acceptance baseline for script automation.

GEX parity goal (from prior recollection for UCSF parity runs):
- Pearson around `0.93`
- Spearman around `0.94`

For automation, treat these as target floors:
- `pearson >= 0.93`
- `spearman >= 0.94`

## Problem Statement
Recent ad hoc runs used `scripts/run_gex_feature_parity_checks.sh` as the
gating comparison and produced lower feature-call parity (`831/937 = 88.6873%`).
That output is not equivalent to the call-only parity workflow and should not be
used as the release gate for this question.

## Locked Inputs
Use these exact paths first:
- CR filtered MEX:
  `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/filtered_feature_bc_matrix`
- CR calls:
  `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/crispr_analysis`
- Barcode translation whitelist (2-column, NXT/TRU):
  `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz`
- Driver script:
  `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/run_callonly_gmm_parity.sh`

## Build Hygiene (Required)
Before parity judgement, clean rebuild:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR star_feature_call
```

## Script Plan (Call-Only Wrapper)
Implemented:
- `scripts/run_ucsf_call_parity_recovery.sh`

Script behavior:
1. Validate required inputs and binaries.
2. Run the existing canonical driver (no re-implementation):
   - `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/run_callonly_gmm_parity.sh`
   - `MIN_UMI_LIST="10 3"`.
3. Parse outputs:
   - `parity_minumi10.tsv`
   - `parity_minumi3.tsv`
   - `umi_correlation_summary.tsv`
4. Emit one normalized summary:
   - `PARITY_RECOVERY_SUMMARY.md`
   - `PARITY_RECOVERY_SUMMARY.tsv`
5. Enforce acceptance checks for `--min-umi 3`:
   - `real_mismatch_count <= 1`
   - `set_equivalent_pct >= 99.98`
   - `abs(delta_umi_common_rows) <= 10`
   - `pearson_all_cr_rows >= 0.99999`
   - `spearman_all_cr_rows >= 0.99999`
6. If failed, print actionable diffs:
   - `real_mismatch_barcodes`
   - `only_star_non_none_detail`
   - threshold mismatch detail rows

Current scope note:
- The implemented wrapper currently gates call-only parity first (hard gate).
- GEX parity remains diagnostic/secondary and can be added as a next extension.
- Barcode namespace preflight in Step 0 remains required before mixed GEX/feature
  comparisons.

## Diagnostic Section (Non-Gating)
If needed, run `scripts/run_gex_feature_parity_checks.sh` only as a diagnostic
for matrix overlap trends. Do not use it as pass/fail for the one-outlier
call-parity expectation.

Exception for this recovery runbook:
- GEX correlation goals (`>=0.93` Pearson, `>=0.94` Spearman) are tracked as
  explicit goals and reported in the recovery summary.

Latest diagnostic rerun (2026-02-18, sequential STAR GEX+PF):
- Command:
  `bash scripts/run_gex_feature_parity_checks.sh --star-run /storage/ucsf-2M/star_runs/star_baseline_iPSC2_1_AALG2_1M_nxt_20260217_160217 --cr-run /storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813 --translation-direction left-to-right --translate-side star --out-dir /tmp/ucsf2m_seq_gex_pf_minumi3_20260218_232125`
- Output report:
  `/tmp/ucsf2m_seq_gex_pf_minumi3_20260218_232125/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
- Result summary:
  - Feature-call parity on common rows: `5338/5379` (`99.2378%`)
  - GEX correlations (filtered vs filtered): Pearson `0.791041`, Spearman `0.701676`
  - Filtered barcode Jaccard: `0.631048`
- Conclusion:
  - Sequential GEX+PF parity does **not** recover the known call-only baseline.
  - Keep call-only parity as hard gate until this gap is explained.

## STAR vs CR Feature-Selection Discrepancy Method (Publication-Preserved)
Date: 2026-02-19

This section captures the reproducible methodology used to analyze why STAR and
Cell Ranger differ on a subset of CRISPR feature assignments.

### Scripts Saved In-Repo
- End-to-end wrapper:
  `scripts/run_ucsf_star1h_cr_analysis.sh`
- Analysis modules:
  - `scripts/ucsf_parity/build_star_m1_delta_vs_cr.py`
  - `scripts/ucsf_parity/classify_star_m1_cr_misses.py`
  - `scripts/ucsf_parity/inspect_barcode_feature_totals.py`
  - `scripts/ucsf_parity/build_star_exact_vs_cr.py`
  - `scripts/ucsf_parity/export_paper_miss_tables.py`
- Script notes:
  `scripts/ucsf_parity/README.md`

### Inputs Used For Reproduced Run
- STAR m=1 PF output:
  `/tmp/ucsf_pf_trace_foxd3_nxtbc_20260219_004502/iPSC2_1_AALG2`
- STAR m=0 PF output:
  `/tmp/ucsf_pf_trace_foxd3_nxtbc_m0_20260219_010611/iPSC2_1_AALG2`
- CR raw MEX:
  `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/raw_feature_bc_matrix`
- CR calls CSV:
  `/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/crispr_analysis/protospacer_calls_per_cell.csv`
- NXT/TRU translation whitelist:
  `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz`

### Reproduction Command
```bash
OUT_DIR=/tmp/ucsf_star1h_vs_cr_repro_$(date +%Y%m%d_%H%M%S)
scripts/run_ucsf_star1h_cr_analysis.sh --out-dir "$OUT_DIR"
```

### Analysis Design
1. Compute STAR 1-Hamming rescue proxy:
   - For each `(feature, barcode)`, compute `STAR(m=1) - STAR(m=0)`.
   - Keep rows with positive delta (`>=1`) as rescued assignments.
2. Normalize barcode namespace:
   - STAR barcodes are NXT; CR barcodes are TRU.
   - Translate NXT -> TRU before joining STAR rows to CR raw MEX/calls.
3. Join each rescued pair to CR:
   - `cr_raw_count` from CR raw matrix at same `(feature, barcode_tru)`.
   - `in_cr_call_list` from CR `protospacer_calls_per_cell.csv`.
4. Classify rescued rows where `cr_raw_count == 0`:
   - `barcode_absent_from_cr_matrix`
   - `shift_to_partner_feature` (`_A` <-> `_B`)
   - `no_partner_signal_in_cr`
5. Independently test exact-match misses:
   - Compare STAR m=0 exact counts vs CR raw for same `(feature, barcode_tru)`.
   - Quantify exact STAR rows that are absent in CR raw.

### Key Findings (Reproduced)
From `STAR_M1_DELTA_VS_CR_SUMMARY.txt`:
- `pairs_total=11579`
- `delta_sum=88171`
- `pairs_cr_raw_gt0=11470` (`99.06%`)
- `pairs_cr_raw_eq0=109` (`0.94%`)
- `delta_sum_cr_raw_eq0=999` (`1.13%` of rescued UMI)

From `STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.summary.txt`:
- `missing_pairs=109`
- `missing_delta_sum=999`
- `barcode_absent_from_cr_matrix_pairs=12` (`delta_sum=12`)
- `shift_to_partner_feature_pairs=15` (`delta_sum=480`)
- `no_partner_signal_in_cr_pairs=82` (`delta_sum=507`)

From `STAR_EXACT_VS_CR_SUMMARY.txt` (STAR m=0 exact rows):
- `pairs_total=34814`
- `star_exact_sum=1498143`
- `pairs_cr_raw_eq0=352` (`1.0111%`)
- `star_sum_cr_raw_eq0=2043` (`0.1364%` of STAR exact UMI)
- Missing exact rows split:
  - `missing_shift_to_partner_pairs=36` (`star_sum=896`)
  - `missing_no_partner_signal_pairs=316` (`star_sum=1147`)

### Example Confirming Non-Trivial Miss
Barcode `AGGTAGGAGAAGATCT` (TRU):
- STAR raw includes `FOXD3_P1_B=101` and `KLF2_P1P2_A=62`
- CR raw includes `KLF2_P1P2_A=62` but `FOXD3_P1_B=0`

This pattern indicates selective drop/reassignment behavior in CR for a subset
of rows, not a blanket inability to recover 1-Hamming events.

### Output Files For Manuscript/Audit Trail
- `STAR_M1_DELTA_VS_CR.tsv`
- `STAR_M1_DELTA_VS_CR_SUMMARY.txt`
- `STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.tsv`
- `STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.summary.txt`
- `STAR_EXACT_VS_CR.tsv`
- `STAR_EXACT_VS_CR_SUMMARY.txt`
- `STAR_EXACT_VS_CR_MISSING_TOP.tsv`
- `BARCODE_INSPECTION.tsv`
- `BARCODE_INSPECTION_SUMMARY.txt`

### Paper Table Export
After running the main analysis wrapper, export combined H0/H1 miss tables:

```bash
python3 scripts/ucsf_parity/export_paper_miss_tables.py \
  --analysis-dir "$OUT_DIR" \
  --out-dir comparisons/ucsf_star_cr_miss_tables_20260219 \
  --top-n 25
```

Produced in this repository:
- `comparisons/ucsf_star_cr_miss_tables_20260219/UCSF_STAR_vs_CR_MISSES_H0_H1.tsv`
- `comparisons/ucsf_star_cr_miss_tables_20260219/UCSF_STAR_vs_CR_MISSES_H0_H1_SUMMARY.tsv`
- `comparisons/ucsf_star_cr_miss_tables_20260219/UCSF_STAR_vs_CR_MISSES_H0_H1_PREVIEW.md`

## CR Multimap Compatibility Workstream (2026-02-19)

Objective:
- Mirror Cell Ranger multimap resolution behavior in a compatibility mode
  without changing default STAR behavior.

Key observations to preserve:
- Some large STAR-vs-CR GEX discrepancies are not barcode-namespace issues and
  not feature-leakage regressions.
- A subset of discordant genes is driven by multimapping classification policy
  differences (including equal-best alignments where one candidate is
  gene-annotated and another is non-gene in the filtered reference context).
- `UniqueAndMult-*` matrices are required to inspect policy effects in STAR; the
  default `matrix.mtx` alone is insufficient for this comparison.

Script updates now available:
- `tests/compare_feature_mex.py` accepts side-specific matrix basenames:
  - `--matrix-basename-a`
  - `--matrix-basename-b`
- `scripts/run_gex_feature_parity_checks.sh` accepts:
  - `--cr-raw-matrix-basename`
  - `--cr-filtered-matrix-basename`
  - `--star-raw-matrix-basename`
  - `--star-filtered-matrix-basename`
- `scripts/report_additional_parity_metrics.py` accepts the same four basenames
  and can consume real-valued matrices (for `UniqueAndMult-*`).

Recommended multimap policy comparison commands:

1. STAR Unique baseline (current compatibility baseline)
```bash
bash scripts/run_gex_feature_parity_checks.sh \
  --star-run /path/to/star_run \
  --cr-run /path/to/cr_run \
  --out-dir /tmp/ucsf_multimap_unique_$(date +%Y%m%d_%H%M%S) \
  --cr-raw-matrix-basename matrix.mtx \
  --cr-filtered-matrix-basename matrix.mtx \
  --star-raw-matrix-basename matrix.mtx \
  --star-filtered-matrix-basename matrix.mtx \
  --translate-side star \
  --translation-direction left-to-right
```

2. STAR Rescue-policy probe (raw GEX only)
```bash
bash scripts/run_gex_feature_parity_checks.sh \
  --star-run /path/to/star_run \
  --cr-run /path/to/cr_run \
  --out-dir /tmp/ucsf_multimap_rescue_$(date +%Y%m%d_%H%M%S) \
  --cr-raw-matrix-basename matrix.mtx \
  --cr-filtered-matrix-basename matrix.mtx \
  --star-raw-matrix-basename UniqueAndMult-Rescue.mtx \
  --star-filtered-matrix-basename matrix.mtx \
  --translate-side star \
  --translation-direction left-to-right
```

Interpretation guardrail:
- Filtered STAR outputs generally provide `matrix.mtx`; `UniqueAndMult-*` are
  emitted in raw outputs. Use raw-restricted-to-CR-filtered sections as the
  primary apples-to-apples check for Rescue-policy probes.

## Suggested Command Skeleton

```bash
OUT_BASE=/tmp/ucsf_call_parity_recovery_$(date +%Y%m%d_%H%M%S)
STAR_FEATURE_CALL_BIN=/mnt/pikachu/STAR-suite/core/legacy/source/star_feature_call \
CR_MEX_DIR=/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/filtered_feature_bc_matrix \
CR_CALLS_DIR=/storage/ucsf-2M/cellranger_runs/cr_baseline_iPSC2_1_AALG2_1M_crstar_sameidx_20260217_200813/outs/crispr_analysis \
OUT_BASE="$OUT_BASE" MIN_UMI_LIST="10 3" \
bash comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/run_callonly_gmm_parity.sh
```

## Output Contract for Next Agent
The new script must produce:
- `PARITY_RECOVERY_SUMMARY.md`
- `PARITY_RECOVERY_SUMMARY.tsv`
- explicit PASS/FAIL line with threshold reasons
- pointer to raw artifacts directory
- explicit GEX parity block with:
  - Pearson
  - Spearman
  - common barcode count
  - pass/fail against `0.93/0.94` thresholds

## Notes for Handoff
- Keep call-only parity as the hard gate.
- Treat mixed GEX+feature parity script as secondary diagnostics.
- Do not change STAR/PF behavior until this parity gate is reproducibly green.
