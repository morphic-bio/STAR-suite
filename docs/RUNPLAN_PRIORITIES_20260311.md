# STAR-suite Prioritized Run Plan (2026-03-11)

This run plan converts the active top-level TODOs into an execution order that is suitable for the paper, regression safety, and release stability.

## Priority 0: Finish the Paper-Critical Validation

### 1. Close the remaining PE bulk quant residual

- Reuse the corrected PE benchmark setup documented in `docs/PAPER_PE_BULK_BENCHMARK_RESULTS_20260310.md`.
- Compare integrated TranscriptVB against Salmon on the same corrected inputs.
- Produce a top-transcript and top-gene delta report.
- Classify the remaining difference into one of:
  - expected model/allocation difference
  - filter/input discrepancy
  - effective-length / compatibility-policy difference
- Exit criterion:
  - either explain the current `~0.06%` residual clearly enough for the paper, or reduce it further with a validated code fix.

### 2. Run the real multi-sample bulk PE benchmark

- Move beyond the current 2-sample smoke and run a real multi-sample set.
- Compare:
  - one `--batchMode 1` STAR run
  - sequential one-sample-at-a-time STAR runs
- Require per-sample parity for:
  - trimmed FASTQs
  - `Y/noY` FASTQs
  - `Y/noY` BAMs
  - transcriptome BAMs
  - `TranscriptVB` outputs
- Exit criterion:
  - byte-identical or otherwise explained deterministic differences with documented policy.

### 3. Revisit SLAM VB over-dispersion

- Use technical replicates and count-stratified slices.
- Sweep `slamVbOverdispPhi`.
- Decide whether the current default remains justified.
- Exit criterion:
  - a short recommendation note with the chosen default and why.

## Priority 1: Close the Most Important Correctness Gaps

### 4. Finish perturb mixed-chemistry support

- Implement per-library `NXT/TRU` auto-detection in `process_features` bootstrap.
- Add automated tests for:
  - mixed `NXT + TRU` feature-library merge
  - 2-column mapping edge cases
  - chemistry precedence matrix
- Exit criterion:
  - mixed-chemistry runs do not depend on a single global chemistry assumption.

### 5. Resolve remaining cell-calling and intronic-parity questions

- Test ambient-window sizing sensitivity for smaller datasets.
- Investigate GeneFull intronic-count inflation vs Cell Ranger.
- Exit criterion:
  - either validated defaults, or explicit caveats/documented behavior.

### 6. Finish poly-G auto mode

- Replace the current heuristic aliasing of `auto` to CellRanger4 mode with data-driven detection.
- Exit criterion:
  - `auto` is dataset-driven and tested.

### 7. Clean up CB/UB and FLEX raw-tag behavior

- Relax CB/UB coupling for non-Flex only.
- Confirm Flex policy remains strict.
- Re-enable and validate FLEX `ZG/ZX` tags once the tag path is stable.
- Exit criterion:
  - tag emission policy is documented and regression-tested.

## Priority 2: Throughput and Architecture

### 8. Implement perturb dynamic permit scheduling

- Phase 1: STAR-only dynamic active-thread limit behind a feature flag.
- Phase 2: queue-aware producer/consumer permit rebalancing.
- Validate on A375 first, then UCSF full data.
- Exit criterion:
  - measurable throughput gain with parity preserved.

### 9. Finish the QC graphics consolidation

- Define the shared QC graphics core.
- Keep Plotly HTML+JSON as the common format.
- Audit downstream scripts for removed PNG expectations.
- Exit criterion:
  - one reusable QC graphics path instead of duplicated per-module outputs.

### 10. Harden feature-offset and decompression behavior

- Add the read-level feature-offset preflight if the false-positive rate is acceptable.
- Confirm native-gz and `--readFilesCommand` paths remain behaviorally identical on representative inputs.
- Exit criterion:
  - no hidden behavior split between decompression modes.

## Priority 3: Coverage Expansion

### 11. Finish missing fixture-dependent smokes

- Run the outstanding memory-leak repair smokes once fixtures are available.

### 12. Expand CRISPR feature-calling threshold validation

- Validate `--crMinUmi` on lineage, FLEX, and antibody datasets.
- Convert the resulting recommendations into tests and user guidance.

## Suggested Execution Order

1. PE bulk residual analysis
2. Real multi-sample bulk PE benchmark
3. SLAM over-dispersion retest
4. Mixed-chemistry perturb support and tests
5. Cell-calling and GeneFull parity investigations
6. Poly-G auto detection
7. CB/UB and FLEX raw-tag cleanup
8. Dynamic permit scheduling
9. QC graphics consolidation
10. Fixture and assay coverage expansion

## Explicitly Out of This Plan

- Release-compatibility tarball work
- Installer-bundle environment selection
- Published `v0.50.3` release smoke

Those are complete enough for current prerelease use and should not block the paper/debugging work.
