# SLAM-seq Blank Model Stage 1+2 Handoff (2026-02-01)

## Summary

Stage 1 (blank transition model) and Stage 2 (p_e regression) have been completed.
A critical finding emerged: **the 0h sample is NOT a true blank** - it shows ~2% T>C
conversion rate, characteristic of labeled RNA. The no4su (QC blank) shows the
expected error rate of ~0.08%.

## Execution

```bash
ARTIFACT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201 \
scripts/run_slam_blank_model.sh
```

- **Duration**: ~40 minutes (processing ~150M reads total)
- **Exit code**: 0 (success)

## Output Location

`/mnt/pikachu/slam_blank_artifacts_20260201/blank_model/`

### Files Generated

| File | Description |
|------|-------------|
| `blank_transition_model.tsv` | 0h transition rates + regression coefficients |
| `blank_transition_model.global.tsv` | no4su (QC blank) transition rates |
| `pe_estimates.tsv` | p_e estimates for all samples |
| `blank_model_summary.json` | Machine-readable summary |

## Key Results

### Transition Rates Comparison

| Sample | T>C Rate | T>A Rate | T>G Rate | Interpretation |
|--------|----------|----------|----------|----------------|
| **no4su (QC blank)** | 0.000835 (0.08%) | 0.000667 (0.07%) | 0.000169 (0.02%) | True error rate |
| **0h** | 0.020418 (2.04%) | 0.000763 (0.08%) | 0.000230 (0.02%) | NOT a blank! |
| **6h** | 0.019958 (2.00%) | 0.000707 (0.07%) | 0.000222 (0.02%) | Labeled |
| **24h** | 0.023351 (2.34%) | 0.000649 (0.06%) | 0.000239 (0.02%) | Labeled |

### p_e Estimates (Stage 2)

| Sample | p_e_blank | p_e_regression | p_e_direct | Comment |
|--------|-----------|----------------|------------|---------|
| 0h | 0.02042 | 0.00099 | 0.02042 | Regression better matches QC blank |
| 6h | 0.02042 | 0.00093 | 0.01996 | Signal = ~1.9% above regression p_e |
| 24h | 0.02042 | 0.00089 | 0.02335 | Signal = ~2.2% above regression p_e |

## Critical Finding

**The 0h sample is NOT suitable as a blank.** Evidence:

1. **T>C rate is 24x higher** than the QC blank (no4su): 2.04% vs 0.08%
2. **T>A and T>G rates are similar** across all samples (~0.07% and ~0.02%), confirming
   these represent sequencing errors
3. **The regression model (p_e_regression)** predicts p_e of ~0.09-0.10%, which matches
   the QC blank rate of 0.08%

This suggests the 0h timepoint already has significant label incorporation. The no4su
sample should be used as the primary blank for this dataset.

## Validation

| Metric | Expected | Observed | Status |
|--------|----------|----------|--------|
| p_e_regression close to QC blank | ~0.0008 | 0.00089-0.00099 | PASS |
| T>A/T>G similar across samples | Yes | Yes (0.06-0.08%) | PASS |
| QC blank T>C << 0h T>C | Yes | 0.08% << 2.04% | **WARNING** |

## Regression Model

```
T>C_predicted = 0.0 + 1.0 * T>A + 1.0 * T>G
```

The regression coefficients (1.0 for both T>A and T>G) are the default Grand-SLAM
approximation. For production use, these should be calibrated from multiple samples.

## Recommendations

1. **Use no4su as primary blank** for this dataset, not 0h
2. **Investigate 0h sample**: Why does it show 2% T>C? Possible causes:
   - Accidental labeling at t=0
   - Batch effect or contamination
   - Different cell preparation protocol
3. **For production**, fit regression coefficients from multiple unlabeled samples
4. **p_e_regression provides better estimates** than using 0h as blank

## Implementation Notes

### Tools Created

1. `slam/tools/blank_model/build_blank_model.py` - Main Python script
   - Counts all 12 transitions from BAM files
   - Uses pysam's `get_aligned_pairs(with_seq=True)` for efficiency
   - Applies trim (5'=11, 3'=10) and quality filters (MAPQ≥10, BASEQ≥20)
   - Fits regression model and estimates p_e

2. `scripts/run_slam_blank_model.sh` - Shell wrapper
   - Uses Stage 0a artifacts as input
   - Configurable via environment variables

### Performance

- Processing rate: ~3M reads/minute
- Total time: ~40 minutes for 4 samples (~150M reads total)

## Files Updated

- `tests/ARTIFACTS.md` - Added blank_model output location
- Created this handoff document

## GEDI (Grand-SLAM) Comparison

GEDI was run on all samples during Stage 0a. The results confirm our findings.

### GEDI T>C Rates (from mismatch tables)

| Sample | TC Count | T Coverage | TC Rate | Interpretation |
|--------|----------|------------|---------|----------------|
| **0h** | 4,166,053 | 320,707,576 | 1.30% | NOT a blank |
| **6h** | 3,009,788 | 254,781,321 | 1.18% | Labeled |
| **24h** | 3,972,935 | 321,410,404 | 1.24% | Labeled |
| **blank** | 143,385 | 297,575,258 | 0.05% | True error rate |

**Ratio 0h/blank: 27x** - GEDI confirms elevated T>C in 0h

### STAR vs GEDI Correlation (NTR)

| Sample | Pearson (≥20 reads) | Spearman (≥20 reads) | Status |
|--------|---------------------|----------------------|--------|
| 0h | 0.948 | 0.958 | PASS |
| 6h | 0.949 | 0.960 | PASS |
| 24h | 0.939 | 0.963 | PASS |

All correlations exceed the 0.9 threshold, confirming STAR-SLAM and GEDI produce
consistent NTR estimates despite the blank selection issue.

### GEDI Estimated Conversion Rates

GEDI's `single_new` parameter (estimated fraction of newly synthesized RNA):

| Sample | single_new | single_old (error) |
|--------|------------|-------------------|
| 0h | 7.28% | 0.1% |
| 6h | 8.09% | 0.1% |
| 24h | 9.91% | 0.1% |
| blank | 8.73% | 0.1% |

Note: GEDI's conversion rate is higher than raw T>C rate because it models the
fraction of newly synthesized RNA using a binomial model. The blank showing 8.73%
"new" is likely due to GEDI using the provided error rate (0.1%) rather than
learning it from the data.

## Next Steps

- **Stage 3**: Estimate p_c from labeled data using exclusion + EM/bisection
- **Update plan**: Revise slam_blank_plan.md to note 0h is not a valid blank
- **Calibrate regression**: Fit coefficients from multiple unlabeled samples

---

## Addendum: no4su as Primary Blank (2026-02-01)

Per follow-up, we re-ran Stage 1/2 using the **no4su** BAM as the *primary blank* and
treated 0h as QC-only (labeled comparison), keeping the same trims and filters.

**Command:**
```bash
OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201
python3 slam/tools/blank_model/build_blank_model.py \
  --blank-bam "$OUT_BASE/star_internal/blank_Aligned.sortedByCoord.out.bam" \
  --qc-blank-bam "$OUT_BASE/star_internal/0h_Aligned.sortedByCoord.out.bam" \
  --labeled-bams "$OUT_BASE/star_internal/6h_Aligned.sortedByCoord.out.bam" \
                "$OUT_BASE/star_internal/24h_Aligned.sortedByCoord.out.bam" \
  --reference /storage/autoindex_110_44/bulk_index/cellranger_ref/genome.fa \
  --output-dir "$OUT_BASE/blank_model_no4su" \
  --trim5p 11 --trim3p 10 \
  --min-mapq 10 --min-baseq 20 --verbose
```

**Key results (no4su as blank):**
- **no4su (blank)** T>C = **0.000835**
- **0h (QC-only)** T>C = **0.020418** (confirms 0h is labeled)
- **6h** T>C = **0.019958**
- **24h** T>C = **0.023351**

**p_e regression (blank-based):**
- p_e_blank = **0.000835**
- p_e_regression ≈ **0.000929** (6h) / **0.000888** (24h)
- Signal above background:
  - 6h: **0.019122** (≈2289× p_e)
  - 24h: **0.022515** (≈2695× p_e)

**Outputs:**  
`/mnt/pikachu/slam_blank_artifacts_20260201/blank_model_no4su/`

---

## SLAM BAM Plan Implementation (2026-02-02)

### Completed Stages

| Stage | Status | Output |
|-------|--------|--------|
| 0: SNP mask + inventory | ✅ Done | SNP BED + 283 FASTQs enumerated |
| 1: Blank model (no4su) | ✅ Done | `blank_model_no4su/` |
| 2: 1M FASTQ alignment | ✅ Done | `fastq_1M_runs/` (prior work) |
| 3: GEDI parity | ✅ Done | Pearson r=0.978 (prior work) |
| 4: STAR-SLAM with trim | Pending | `fastq_1M_runs_trimmed/` |
| 5: Documentation | ✅ Done | This document |

### 1M FASTQ Analysis Summary

From prior work in `/mnt/pikachu/slam_blank_artifacts_20260201/fastq_1M_runs/`:

- **Total samples:** 283 FASTQs (1M reads each)
- **Timepoints:** 0h, 6h, 24h
- **T>C rates by timepoint:**
  - 0h: 1.033% ± 0.366%
  - 6h: 0.962% ± 0.323%
  - 24h: 0.910% ± 0.375%
- **STAR vs GrandSLAM correlation:** Pearson r = 0.978
- **Three-method comparison:** STAR, PySAM, GEDI all show consistent T>C rates

See `fastq_1M_runs/SLAM_1M_ANALYSIS_SUMMARY.md` for full details.

### Scripts Created

| Script | Purpose |
|--------|---------|
| `scripts/run_slam_1M_alignments.sh` | STAR-SLAM on 1M FASTQs (no trim) |
| `scripts/run_slam_gedi_parity.sh` | GEDI comparison |
| `scripts/run_slam_1M_with_trim.sh` | STAR-SLAM with trimming |

### Reproducible Commands (Publication-Ready)

**Inputs**
- SNP mask: `/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz`
- FASTQ set: `/mnt/pikachu/NW-5-21/SLAM-Seq-1M/` (1M read subsets)
- STAR index: `/storage/autoindex_110_44/bulk_index`
- GEDI genome: `/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml`

**1) STAR-SLAM 1M alignments (compat mode, no trim)**
```bash
OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201 \
FASTQ_DIR=/mnt/pikachu/NW-5-21/SLAM-Seq-1M \
THREADS=24 \
./scripts/run_slam_1M_alignments.sh
```

**2) GEDI parity (compat mode, no trim; GEDI as reference)**
```bash
TRIM5_OVERRIDE=0 TRIM3_OVERRIDE=0 \
OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201 \
./scripts/run_slam_gedi_parity.sh
```

**3) STAR-SLAM trimmed run (blank-derived constants)**
```bash
OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201 \
./scripts/run_slam_1M_with_trim.sh
```

**Notes**
- GEDI runs now consume **SNP‑masked BAMs** (bedtools intersect) so STAR vs GEDI
  parity is measured on identical masked positions.
- Compat parity uses **no trims**; blank affects only trimmed mode.

### Parity Workflow (GEDI as reference)

**Important intent:** We treat **GEDI as the reference**, and only adjust **STAR** if parity fails.
The following steps **prepare GEDI artifacts** (SNP-masked BAMs and consistent trims) so parity
comparisons are fair. Blank-derived trims **must not** affect compat-mode parity.

#### 1) Compat parity (no trim, blank not involved)
This isolates STAR vs GEDI without blank influence:

```bash
TRIM5_OVERRIDE=0 TRIM3_OVERRIDE=0 ./scripts/run_slam_gedi_parity.sh
```

**Run result (2026‑02‑02):**
- 0h: Pearson 0.9777, Spearman 0.9902 (PASS)
- 24h: Pearson 0.9667, Spearman 0.9849 (PASS)
- 6h: Pearson 0.9722, Spearman 0.9861 (PASS)
- no4su: GEDI output missing MAP column → parity compare skipped
  - GEDI `ARID1A_no4su.tsv.gz` header lacks MAP; compare script skips it.

#### 2) Trimmed parity (blank-derived constants)
Uses the QC-derived trim constants from `qc/trim_6h.slam_qc.json`:

```bash
./scripts/run_slam_gedi_parity.sh
```

**Trimmed STAR outputs (auto-trim constants 5′=11, 3′=10):**
- Stored in: `/mnt/pikachu/slam_blank_artifacts_20260201/fastq_1M_runs_trimmed/`
- Blank (no4su) T>C reduced: **0.2174% → 0.0573%** (Δ −0.1601%)

**Trimmed GEDI parity (2026‑02‑02, SNP‑masked BAMs):**
- Outputs: `/mnt/pikachu/slam_blank_artifacts_20260201_trimmed_parity/`
- 0h: Pearson 0.9550, Spearman 0.9549 (PASS)
- 24h: Pearson 0.9461, Spearman 0.9572 (PASS)
- 6h: Pearson 0.9557, Spearman 0.9608 (PASS)
- no4su: GEDI output missing MAP column → parity compare skipped

**Trimmed vs untrimmed STAR (MAP correlation, readcount ≥20):**
- 0h: 0.9497
- 6h: 0.9559
- 24h: 0.9570
- no4su: 0.8110

#### 3) Trimmed vs compat STAR-only sanity check

```bash
./scripts/run_slam_1M_with_trim.sh
```

**Note:** GEDI runs now use **SNP‑masked BAMs** (via bedtools intersect) to ensure
the same masked positions as STAR’s `--slamSnpMaskIn` path.
| `slam/tools/blank_model/build_blank_model.py` | Blank transition model |

### Key Finding

The no4su sample (T>C = 0.08%) is the proper blank for this dataset. The 0h timepoint
(T>C = 2.04%) shows significant labeling and should NOT be used as a blank.

---

## SLAM BAM Plan Implementation (2026-02-02)

### Completed Stages

| Stage | Status | Output |
|-------|--------|--------|
| 0: SNP mask + inventory | ✅ Done | SNP BED + 283 FASTQs enumerated |
| 1: Blank model (no4su) | ✅ Done | `blank_model_no4su/` |
| 2: 1M FASTQ alignment | ✅ Done | `fastq_1M_runs/` (prior work) |
| 3: GEDI parity | ✅ Done | Pearson r=0.978 (prior work) |
| 4: STAR-SLAM with trim | Pending | `fastq_1M_runs_trimmed/` |
| 5: Documentation | ✅ Done | This document |

### 1M FASTQ Analysis Summary

From prior work in `/mnt/pikachu/slam_blank_artifacts_20260201/fastq_1M_runs/`:

- **Total samples:** 283 FASTQs (1M reads each)
- **Timepoints:** 0h, 6h, 24h
- **T>C rates by timepoint:**
  - 0h: 1.033% ± 0.366%
  - 6h: 0.962% ± 0.323%
  - 24h: 0.910% ± 0.375%
- **STAR vs GrandSLAM correlation:** Pearson r = 0.978

See `fastq_1M_runs/SLAM_1M_ANALYSIS_SUMMARY.md` for full details.

### Scripts Created

| Script | Purpose |
|--------|---------|
| `scripts/run_slam_1M_alignments.sh` | STAR-SLAM on 1M FASTQs (no trim) |
| `scripts/run_slam_gedi_parity.sh` | GEDI comparison |
| `scripts/run_slam_1M_with_trim.sh` | STAR-SLAM with trimming |
