# Handoff: CRISPR Compat + Flex Modular (Jan 24, 2026)

## Goals
- Continue CR compatibility work on a **more challenging dataset** (beyond A375).
- Once compatibility is acceptable, **move Flex** to the new modular components.

## Current State (Branch: perturb)
- Latest pushed commit: `7757e96` (compat defaults + docs).
- Compat path uses **anchor search + per-feature mode bootstrap**:
  - `--use_feature_anchor_search`, `--require_feature_anchor_match`
  - `--feature_mode_bootstrap_reads 100000`
  - After bootstrap, uses **mode ±1 only** (no anchor scan).
- Defaults updated:
  - `min_counts` default **0** (inclusive single-read UMIs).
  - `feature_n` default **1** (Max Ns in feature).
- CRISPR compat runs + notes are summarized in:
  - `plans/STAR_pf_modular_integration_summary.md`

## Latest A375 Results (CRISPR Guide Capture, filtered barcodes)
- **Hamming ≤3, min_counts=0**
  - Output: `/storage/A375/assign_anchor_mode_h3_min0_20260123_233445`
  - STAR 3,058,669 vs CR 3,002,498 (+1.87%)

- **Hamming ≤1, min_counts=0**
  - Output: `/storage/A375/assign_anchor_mode_h1_min0_20260123_235424`
  - STAR 3,041,871 vs CR 3,002,498 (+1.31%)

- **Hamming ≤1, min_counts=0, feature_n=0**
  - Output: `/storage/A375/assign_anchor_mode_h1_min0_n0_20260124_000205`
  - STAR 3,041,865 vs CR 3,002,498 (+1.31%)

Artifacts are tracked in `tests/ARTIFACTS.md`.

## CR Molecule Validation (Barcode GTGTGATGTGTGGTCA)
- Source: `/storage/A375/outputs/1k_CRISPR_5p_gemx_count_sample_molecule_info.h5` + raw FASTQs.
- CR UMIs: 3068.
- UMI min‑Hamming distribution: **0: 2925, 1: 117, ≥2: 0**.
- Takeaway: CR‑counted UMIs do **not** require >1 mismatches for this barcode.

## Repro / Compare (A375)

### Run assignBarcodes (process_features)
```bash
OUT=/storage/A375/assign_anchor_mode_h1_min0_<stamp>
FASTQ_DIR=/storage/A375/fastqs/1k_CRISPR_5p_gemx_fastqs/crispr
WHITELIST=/storage/A375/3M-5pgex-jan-2023.txt
FEATURE_REF=/storage/A375/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv
FILTERED=/storage/A375/assign_anchor_mode_h3_20260123_230407/filtered_barcodes_nosuffix.tsv

core/features/process_features/assignBarcodes \
  --whitelist "$WHITELIST" \
  --featurelist "$FEATURE_REF" \
  --directory "$OUT" \
  --maxHammingDistance 1 \
  --feature_n 1 \
  --feature_constant_offset 0 \
  --limit_search -1 \
  --stringency 1 \
  --use_feature_anchor_search \
  --require_feature_anchor_match \
  --feature_mode_bootstrap_reads 100000 \
  --filtered_barcodes "$FILTERED" \
  --threads 1 \
  "$FASTQ_DIR"
```

### Create features.tsv/barcodes.tsv for compare
```bash
python3 tools/feature_barcodes/assignbarcodes_mex_stub.py \
  --assign-out "$OUT/crispr/filtered" \
  --feature-csv "$FEATURE_REF" \
  --default-feature-type "CRISPR Guide Capture"
```

### Compare vs CR MEX
```bash
python3 tests/compare_a375_star_mex.py \
  /storage/A375/outputs/unpacked/sample_filtered_feature_bc_matrix \
  "$OUT/crispr/filtered" \
  --feature-types "CRISPR Guide Capture"
```

## Next Steps (Compatibility)
1. **Run the same pipeline on a more challenging dataset** (identify path + CR MEX).
2. If deltas persist, inspect top‑delta barcodes with `PF_TRACE_*` (counts/dedup/reads/anchor).
3. Consider bracketing with `maxHammingDistance 0` and `barcode_n 0` if needed.

## Next Steps (Flex Modular)
- After compat is good, follow:
  - `plans/scrna_flex_plan.md`
  - `plans/agent_handoff_assignBarcodes_flex.md`
- Goal: Flex should use the **modular** components from process_features.

## Notes / Hygiene
- Build artifacts (binaries) are **not committed**; keep them untracked.
- If you need to clean: `git clean -fd` (confirm before use).
