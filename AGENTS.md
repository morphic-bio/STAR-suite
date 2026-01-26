# STAR-suite Agent Guide

This document is for coding agents working in this repo. Keep it short and
actionable; link to deeper docs rather than copying them.

## Scope and Goals

- STAR-suite reorganizes STAR into modules while keeping `core/legacy` as the
  single source of truth for the STAR core.
- Changes should preserve STAR CLI compatibility unless explicitly planned.
- Large data and test outputs stay untracked; document their locations in
  `tests/ARTIFACTS.md`.

## Repository Map

- `core/legacy/` - upstream STAR layout (single source of truth).
- `core/features/` - shared feature overlays (vbem, yremove, bamsort, etc.).
- `core/features/process_features/` - vendored `process_features` toolchain.
- `core/features/feature_barcodes/` - standalone tools (`assignBarcodes`, etc.).
- `core/features/libscrna/` - shared EmptyDrops/OrdMag/Occupancy implementations.
- `flex/` - STAR-Flex code and tools.
- `slam/` - SLAM-seq code and tools.
- `tests/` - smoke/regression scripts; keep outputs untracked.
- `docs/` - technical summaries and runbooks.

## Build and Smoke Tests

- Core build: `make core` (binary: `core/legacy/source/STAR`).
- Flex tools: `make flex` or `make flex-tools`.
- Feature tools: `make feature-barcodes-tools`.
- CB/UB regression: `tests/run_cbub_regression_test.sh`.
- CRISPR calling test: `tests/test_cr_compat_crispr_calling.sh`.
- Flex smoke: `tests/run_flex_smoke.sh` (if fixtures available).

## Data and Artifacts (Do Not Commit)

- External datasets live under `/storage/` (e.g., `/storage/A375/...`).
- Test outputs should be in `tests/*_output*/` or `/tmp/`.
- Add new artifact locations to `tests/ARTIFACTS.md`.

## Key Technical Defaults and Recent Changes

### Unsorted BAM CB/UB injection

- Unsorted BAM CB/UB tags are injected via the same path as sorted BAM.
- Buffered mode uses `g_unsortedTagBuffer` and a `SamtoolsSorter` "noSort"
  mode; output order is deterministic (spill files first, then in-memory).
- Legacy two-pass `--soloAddTagsToUnsorted` path is removed.
- Summary: `docs/UNSORTED_BAM_CBUB_IMPLEMENTATION_SUMMARY.md`.

### CB/UB Independence (Flex)

- CB and UB are handled independently in Flex.
- `status==2` means CB valid, UMI invalid; CB tags still emitted.
- Summary: `docs/CB_UB_INDEPENDENCE_IMPLEMENTATION_SUMMARY.md`.

### CRISPR Feature Calling (CR-compat)

- CR-compat mode runs GMM-based CRISPR calling automatically when Guide Capture
  features are present.
- `--crMinUmi` default is 10 (CRISPR guides); use 2-3 for lineage barcodes.
- Summary: `docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md`.

### Heatmaps (process_features)

- Cairo PNGs removed; Plotly HTML+JSON outputs only.
- Outputs: `Feature_counts_heatmap.html` + `.json`, `Feature_types_heatmap.html`
  + `.json`.
- Summary: `docs/HEATMAP_REFACTOR_SUMMARY.md`.

### Feature Offset Detection (assignBarcodes / pf_api)

- Auto-detects a global offset from the `pattern` column.
- Errors if multiple offsets are detected (heterogeneity threshold 5%).
- Use `--feature_constant_offset N` for a fixed global offset.
- Use `--force-individual-offsets` for per-feature offsets.
- Docs: `docs/feature_barcodes.md`.

## Flex Integration Notes

- Flex now uses `libscrna` for EmptyDrops/OrdMag/Occupancy (no duplicate
  implementations).
- Ensure Flex builds link `libscrna` and include `core/features/libscrna/include`.
- Plan: `plans/refactor_flex_plan.md`.

## CR-compat GEX Parity Notes

- Use `GeneFull` for GEX parity (CR includes introns since v7).
- `--soloCrGexFeature` controls which GEX MEX is merged in CR-compat mode.
- `--soloCbUbRequireTogether` controls CB/UB pairing for tag injection.

## Docs to Check First

- `docs/UNSORTED_BAM_CBUB_IMPLEMENTATION_SUMMARY.md`
- `docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md`
- `docs/HEATMAP_REFACTOR_SUMMARY.md`
- `docs/feature_barcodes.md`
- `docs/todos`
- `tests/ARTIFACTS.md`

## Branching and Merges

- Feature branches merge into `perturb`, then squash-merge into `master`.
- Keep large binaries and datasets untracked; update `.gitignore` if needed.

## Output Hygiene

- Do not commit generated binaries or test outputs.
- If a new test creates outputs, add its location to `tests/ARTIFACTS.md`.
