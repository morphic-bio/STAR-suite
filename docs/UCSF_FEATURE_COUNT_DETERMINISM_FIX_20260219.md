# UCSF Feature-Count Determinism Fix (2026-02-19)

Date: 2026-02-19

## Summary
A real non-deterministic behavior was identified in feature matching and fixed.

- File fixed: `core/features/process_features/src/assignBarcodes.c`
- Function: `find_feature_match_parallel`
- Symptom: repeated runs with identical inputs produced small A/B guide flips
  (same total UMIs, different per-feature pair assignment in a few rows).

## Root Cause
`find_feature_match_parallel` stored per-search-slot state into arrays indexed by
`omp_get_thread_num()`, while the search itself is over fixed slots `i=0..3`.

Because OpenMP worker assignment to `i` can vary, per-slot results could be
overwritten/non-deterministically mixed, causing occasional feature-pair flips.

## Code Fix
Implemented deterministic per-slot aggregation:

1. Keep local match state per `i` inside the parallel loop.
2. Write results back into arrays indexed by `i` (not thread id).
3. Resolve exact matches and best hamming matches deterministically by `i` order.
4. Use static schedule for stable slot partitioning (`schedule(static)`).

## Validation
### 1) Repro before fix (full STAR path)
Comparing two identical reruns (same flags/inputs) on CRISPR Guide Capture raw MEX:

- `changed=4`, `a_only=1`, `b_only=1` pair-level differences
- Total UMIs unchanged

Observed flips were small A/B partner swaps (for example `RUNX2_P2_A` vs
`RUNX2_P2_B`, `CREB1_P1P2_A` vs `CREB1_P1P2_B`).

### 2) Determinism after fix (full STAR path)
Runs:

- `star_baseline_iPSC2_1_AALG2_1M_nxt_fixgex_patchA_20260219_052317`
- `star_baseline_iPSC2_1_AALG2_1M_nxt_fixgex_patchB_20260219_052729`

Pair-level comparison of CRISPR Guide Capture raw MEX:

- `all_raw`: `pairs_a=18591`, `pairs_b=18591`, `changed=0`
- `cr_filtered`: `pairs_a=16790`, `pairs_b=16790`, `changed=0`

File-level identity (`outs/raw_feature_bc_matrix`):

- `barcodes.tsv.gz` md5 identical
- `features.tsv.gz` md5 identical
- `matrix.mtx.gz` md5 identical

CRISPR call outputs are also identical (`outs/crispr_analysis/*` md5 match).

## Parity Impact
This fix removes a true internal nondeterministic source. Remaining STAR-vs-CR
feature differences are not due to stochastic feature assignment within STAR.

## Artifacts
Key run outputs used for validation:

- `/storage/ucsf-2M/star_runs/star_baseline_iPSC2_1_AALG2_1M_nxt_fixgex_patchA_20260219_052317`
- `/storage/ucsf-2M/star_runs/star_baseline_iPSC2_1_AALG2_1M_nxt_fixgex_patchB_20260219_052729`
- `/tmp/ucsf2m_gex_feature_parity_patchA_20260219_053151/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
