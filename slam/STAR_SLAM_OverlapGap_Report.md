# STAR vs GEDI Overlap Gap Analysis

**Date:** Generated automatically  
**Purpose:** Test whether gene overlaps (especially MT genes) explain the ~0.92 correlation ceiling

---

## Inputs

- **STAR output:** `test/tmp_prod_compat/default_SlamQuant.out`
- **GEDI output:** `/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_rerun.tsv.gz`
- **GTF annotation:** `test/fixtures/slam/ref/genes.gtf`

---

## NTR Metric Comparison

GEDI provides both MAP (maximum a posteriori) and Mean (posterior mean) estimates. Comparing to STAR's NTR:

| GEDI Metric | NTR Pearson (>=20) | k/nT Pearson (>=20) |
|-------------|-------------------|---------------------|
| MAP | 0.920467 | 0.905483 |
| Mean | 0.872975 | 0.905483 |

**Selected:** GEDI MAP for NTR comparison (based on higher correlation)

---

## Gene Overlap Statistics

| Category | Count |
|----------|-------|
| Total genes in GTF | 38592 |
| Genes with exon overlaps | 2050 |
| Overlap clusters | 966 |
| MT genes (chrM) | 13 |
| MT genes with overlaps | 5 |

### Overlap Clusters (Top 10 by size)

| Cluster | Size | Genes |
|---------|------|-------|
| 1 | 22 | PCDHGA7, PCDHGA8, PCDHGB6, PCDHGA2, PCDHGC3, ... (+17 more) |
| 2 | 15 | PCDHA8, PCDHA5, PCDHAC2, PCDHAC1, PCDHA1, ... (+10 more) |
| 3 | 9 | UGT1A7, UGT1A1, UGT1A10, UGT1A4, UGT1A3, ... (+4 more) |
| 4 | 5 | GPRASP2, GPRASP3, ARMCX5-GPRASP2, ARMCX5-GPRASP2, ARMCX5 |
| 5 | 4 | CDK14, PTTG1IP2, GTPBP10, CLDN12 |
| 6 | 4 | ENSG00000291287, FAM95B1, ENSG00000240240, ENSG00000290969 |
| 7 | 4 | ENSG00000275778, TAS2R14, PRR4, PRH1 |
| 8 | 4 | SNRPN, SNHG14, ENSG00000261069, SNURF |
| 9 | 3 | ZC3H11A, ZBED6, ENSG00000288644 |
| 10 | 3 | ENSG00000285622, DPT, LINC00970 |

---

## Test A: Overlap-Gene Exclusion

### Baseline (All Genes)

**All common genes (n=11079)**

| Threshold | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|---------|-------------|--------------|--------------|---------------|
| >=20 | 4255 | 0.920467 | 0.924726 | 0.905483 | 0.931271 |
| >=50 | 2166 | 0.926568 | 0.929443 | 0.930352 | 0.936203 |
| >=100 | 1149 | 0.906414 | 0.910138 | 0.824392 | 0.917799 |

### A1: Exclude All Overlapping Genes

**Non-overlapping genes (n=10288)**

| Threshold | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|---------|-------------|--------------|--------------|---------------|
| >=20 | 3936 | 0.925057 | 0.930698 | 0.910803 | 0.936680 |
| >=50 | 2011 | 0.925550 | 0.928811 | 0.929002 | 0.935291 |
| >=100 | 1073 | 0.904206 | 0.908417 | 0.806487 | 0.916061 |

**A1 (Delta vs Baseline)**

| Threshold | N Genes | ΔNTR Pearson | ΔNTR Spearman | Δk/nT Pearson | Δk/nT Spearman |
|-----------|---------|--------------|---------------|---------------|----------------|
| >=20 | 3936 | +0.004590 | +0.005972 | +0.005321 | +0.005410 |
| >=50 | 2011 | -0.001018 | -0.000632 | -0.001350 | -0.000913 |
| >=100 | 1073 | -0.002208 | -0.001721 | -0.017905 | -0.001738 |

### A2: Exclude MT Genes Only

**Non-MT genes (n=11066)**

| Threshold | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|---------|-------------|--------------|--------------|---------------|
| >=20 | 4242 | 0.920248 | 0.924587 | 0.905804 | 0.931170 |
| >=50 | 2154 | 0.926183 | 0.929078 | 0.930266 | 0.935874 |
| >=100 | 1137 | 0.905191 | 0.908938 | 0.823161 | 0.916709 |

**A2 (Delta vs Baseline)**

| Threshold | N Genes | ΔNTR Pearson | ΔNTR Spearman | Δk/nT Pearson | Δk/nT Spearman |
|-----------|---------|--------------|---------------|---------------|----------------|
| >=20 | 4242 | -0.000219 | -0.000139 | +0.000321 | -0.000101 |
| >=50 | 2154 | -0.000385 | -0.000365 | -0.000086 | -0.000329 |
| >=100 | 1137 | -0.001223 | -0.001200 | -0.001231 | -0.001089 |

### A3: Exclude Both MT and Overlapping Genes

**Non-MT, non-overlapping genes (n=10280)**

| Threshold | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|---------|-------------|--------------|--------------|---------------|
| >=20 | 3928 | 0.924861 | 0.930610 | 0.911159 | 0.936617 |
| >=50 | 2004 | 0.925187 | 0.928530 | 0.928910 | 0.935023 |
| >=100 | 1066 | 0.903050 | 0.907500 | 0.805200 | 0.915197 |

**A3 (Delta vs Baseline)**

| Threshold | N Genes | ΔNTR Pearson | ΔNTR Spearman | Δk/nT Pearson | Δk/nT Spearman |
|-----------|---------|--------------|---------------|---------------|----------------|
| >=20 | 3928 | +0.004394 | +0.005884 | +0.005676 | +0.005346 |
| >=50 | 2004 | -0.001381 | -0.000913 | -0.001442 | -0.001180 |
| >=100 | 1066 | -0.003363 | -0.002638 | -0.019191 | -0.002602 |

---

## Test B: Locus-Collapsed Aggregation

Overlapping genes are aggregated into single "locus" units by summing read counts and taking weighted-average NTR.

**Locus-collapsed (clusters merged)**

| Threshold | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|---------|-------------|--------------|--------------|---------------|
| >=20 | 4228 | 0.924070 | 0.929459 | 0.905101 | 0.935557 |
| >=50 | 2176 | 0.925855 | 0.928506 | 0.929794 | 0.934819 |
| >=100 | 1157 | 0.903344 | 0.905612 | 0.822071 | 0.913367 |

**B (Delta vs Baseline)**

| Threshold | N Genes | ΔNTR Pearson | ΔNTR Spearman | Δk/nT Pearson | Δk/nT Spearman |
|-----------|---------|--------------|---------------|---------------|----------------|
| >=20 | 4228 | +0.003603 | +0.004733 | -0.000381 | +0.004287 |
| >=50 | 2176 | -0.000713 | -0.000938 | -0.000558 | -0.001385 |
| >=100 | 1157 | -0.003069 | -0.004526 | -0.002321 | -0.004432 |

---

## Summary and Interpretation

- **Excluding overlapping genes**: NTR Pearson = 0.925057 (Δ = +0.004590)
- **Excluding MT genes**: NTR Pearson = 0.920248 (Δ = -0.000219)
- **Excluding MT + overlapping**: NTR Pearson = 0.924861 (Δ = +0.004394)
- **Locus-collapsed aggregation**: NTR Pearson = 0.924070 (Δ = +0.003603)

### Key Findings

1. **Overlapping genes are NOT the main driver**: Excluding them changes NTR Pearson by only +0.0046
2. **MT genes have minimal impact on NTR correlation**: Excluding them changes NTR Pearson by only -0.0002
3. **Locus aggregation has limited impact**: NTR Pearson changes by only +0.0036

### Conclusion

The ~0.92 correlation ceiling is **NOT primarily driven by gene overlap handling**. 
The multi-gene assignment hypothesis is **not supported** by this analysis.

---

## Deep Dive: The True Source of Divergence

### Critical Finding: k/nT is SIMILAR, NTR is DIFFERENT

Comparing raw conversion rates (k/nT = conversions/coverage) vs NTR estimates for high-divergence genes:

| Gene | Symbol | GEDI NTR | STAR NTR | NTR Δ | GEDI k/nT | STAR k/nT | k/nT Δ |
|------|--------|----------|----------|-------|-----------|-----------|--------|
| QARS1 | — | 0.087 | 1.000 | **+0.91** | 0.0027 | 0.0189 | +0.016 |
| H4C3 | Histone | 0.216 | 1.000 | **+0.78** | 0.0053 | 0.0207 | +0.015 |
| FANCI | — | 1.000 | 0.238 | **-0.76** | 0.0194 | 0.0053 | -0.014 |
| LONP1 | — | 0.254 | 1.000 | **+0.75** | 0.0054 | 0.0209 | +0.016 |
| RPL10 | Ribosomal | 0.212 | 0.874 | **+0.66** | 0.0048 | 0.0162 | +0.011 |
| RPL17 | Ribosomal | 0.001 | 0.645 | **+0.64** | 0.0005 | 0.0124 | +0.012 |
| DCD | — | 0.726 | 0.133 | **-0.59** | 0.0144 | 0.0035 | -0.011 |

**Key Observation:** 73/73 genes with |NTR Δ| > 0.3 have |k/nT Δ| < 0.02

### Component-Level Correlation Analysis

Script: `tests/slam/test_ntr_model_gap.py`

| Component | Pearson | Interpretation |
|-----------|---------|----------------|
| ReadCount | 0.976 | High - similar gene assignment |
| Coverage | **0.972** | High - similar position counting |
| Conversions | **0.910** | Lower - some detection differences |
| k/nT (raw) | 0.905 | Driven by conversion/coverage ratio |
| NTR (original) | 0.920 | Model smoothing improves slightly |

**Observation:** Coverage corr (0.972) > Conversions corr (0.910)
- Suggests T→C detection differences contribute to the gap
- Same positions are being counted, but different T→C calls

### NTR Model Test (Approximate)

**Test:** Recompute NTR from k/nT using analytical formula π = (k/n - err) / (conv - err)

| Comparison | Pearson |
|------------|---------|
| STAR NTR vs GEDI NTR (original) | 0.920 |
| STAR NTR* vs GEDI NTR* (same formula) | 0.925 |

**Improvement:** +0.004 when using same formula

**Important Caveats:**
- This test uses gene-level aggregate (k, n), NOT STAR's histogram-based EM
- STAR's actual solver uses per-read mismatch histograms with priors
- A small improvement does NOT definitively rule out model/priors as a driver
- The analytical formula is a simplification of the full model

### Key Insight: NTR > k/nT Correlation

The NTR correlation (0.920) is **higher** than k/nT correlation (0.905).

This means NTR models are **smoothing/improving** the raw k/nT divergence, not adding noise.
The k/nT -> NTR mapping partially compensates for input differences.

### Evidence Summary

| Hypothesis | Evidence | Verdict |
|------------|----------|---------|
| Multi-gene assignment | Excluding overlaps: +0.005 | ❌ NOT the cause |
| MT gene overlaps | Excluding MT: -0.000 | ❌ NOT the cause |
| NTR model (approximate) | Same formula: +0.004 | ⚠️ Unclear (test is approximate) |
| k/nT input divergence | Conv corr 0.91 < Cov corr 0.97 | ⚠️ Contributes to gap |

### Conclusion

The gap is **NOT** due to gene assignment or overlap handling.

The gap appears to be due to a combination of:
1. **T→C conversion detection differences** (Conv corr < Cov corr)
2. **Possibly NTR model/priors** (can't rule out with approximate test)

### To Definitively Determine Root Cause

Would need:
1. **Per-position mismatch comparison** - compare T→C calls at each position
2. **Run STAR's EM on GEDI's histogram data** - true model equivalence test
3. **Mismatch position histograms** - distribution of conversions along reads
