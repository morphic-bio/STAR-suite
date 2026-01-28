# NB-EM / Poisson-EM / SCEPTRE-Parity Results (A375 + vignettes)

**Date:** 2026-01-28

This note summarizes the NB‑EM / Poisson‑EM / SCEPTRE‑parity experiments and the datasets used.
All outputs are **untracked** and listed in `tests/ARTIFACTS.md`.

## Test Sets

### A375 CRISPR-only MEX
- **MEX:** `/storage/A375/star_gex_features_cr_parity_genefull_20260124_080640_crispr_only/crispr_mex`
- **Cell Ranger reference:**
  `/storage/A375-CR-9.01/1k_CRISPR_5p_gemx_count_refmatch_2024a_fullraw/outs/per_sample_outs/1k_CRISPR_5p_gemx_count_refmatch_2024a_fullraw/count/crispr_analysis/protospacer_calls_per_cell.csv`
- **Notes:** 11 CRISPR features, 1,167 barcodes

### SCEPTRE example (low)
- **SCEPTRE output:** `tests/sceptre_example_output_low_20260128_040725/`
- **MEX:** `tests/sceptre_example_low_mex_20260128_040725/`
- **SCEPTRE assignments:** `grna_assignments_per_cell.tsv.gz`

### crispat vignette
- **MEX (0–4):** `tests/crispat_vignette_mex_0_4/`
- **MEX (0–9):** `tests/crispat_vignette_mex_0_9/`
- **crispat assignments:**
  - `tests/crispat_vignette/negative_binomial/assignments_0-4.csv`
  - `tests/crispat_vignette/negative_binomial/assignments_5-9.csv`

## Results Summary

### A375 vs Cell Ranger (CR)
| Method | Output dir | Calls summary | Concordance vs CR |
|---|---|---|---|
| **NB‑EM (NB)** | `/storage/A375/nbem_a375_20260128_075248/` | 0 calls (no features) | 7.2% unassigned match, 92.8% mismatch |
| **Poisson‑EM** (`--poisson-em`) | `/storage/A375/poisson_em_a375_20260128_075938/` | 823 single, 17 multi | **65.8% match_assigned**, 29.6% mismatch |
| **SCEPTRE‑parity** (`--sceptre-parity`) | `/storage/A375/sceptre_parity_a375_20260128_081139/` | 355 single, 0 multi | 26.3% match_assigned, 69.2% mismatch |

**Takeaway:** Poisson‑EM is substantially better than NB‑EM on A375, but still far from CR parity. GMM in CR‑compat mode remains the only method that matches CR on A375 (see `tests/crispr_feature_calling_comparison_report.md`).

### SCEPTRE example (low)
- **SCEPTRE‑parity** run: `tests/nbem_sceptre_parity_example_low_20260128_081511/`
- **Exact match vs SCEPTRE:** 98.6% (986/1000), 1.1% mismatch

**Takeaway:** SCEPTRE‑parity closely matches SCEPTRE on its own example set.

### crispat vignette
- **0–4 subset (current code):** `tests/nbem_crispat_vignette_output_0_4_current/`
  - Assigned‑cell concordance: **~95.1%**
  - Mismatches are almost entirely **None/Multi**, not wrong‑guide calls.
- **0–9 subset (current code):** `tests/nbem_crispat_vignette_output_0_9_current/`
  - Assigned‑cell concordance: **~93.6%**
  - Divergence dominated by **None** and **Multiplet** (not wrong guides).
- **Lower fallback UMI (0–9):**
  - `--backup-threshold 3`: `tests/nbem_crispat_vignette_output_0_9_bk3/` → **94.4%** assigned concordance
  - `--backup-threshold 2`: `tests/nbem_crispat_vignette_output_0_9_bk2/` → **94.2%** assigned concordance

**Takeaway:** Reducing fallback threshold slightly reduces “None” and improves assigned concordance, but gains are modest.

## Current SCEPTRE‑Parity Mode
- **Flag:** `--sceptre-parity`
- Behavior: **Poisson likelihood + per‑guide calls (no low‑MOI collapse), no fallback**
- Intended for *algorithmic parity*, not CR‑compat output parity

## Next Steps (if desired)
1) Try **Poisson‑EM + fallback threshold tuning** on A375 to approach CR calls.
2) Consider **per‑guide independent outputs** (no collapse) for SCEPTRE‑style comparison only.
3) Keep **GMM min‑UMI=10** for CR parity on A375 (documented in `tests/crispr_feature_calling_comparison_report.md`).
