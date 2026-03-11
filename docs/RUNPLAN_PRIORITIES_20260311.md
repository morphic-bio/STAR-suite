# STAR-suite Prioritized Run Plan (2026-03-11)

This run plan converts the active top-level TODOs into an execution order that is suitable for the paper, regression safety, and release stability.

## Priority 0: Finish the Paper-Critical Validation

- No active paper-critical blockers at present.

## Priority 1: Close the Most Important Correctness Gaps

### 1. Establish FLEX parity versus Cell Ranger 9.0

- Run a representative FLEX dataset against Cell Ranger 9.0.
- Confirm output parity at the level needed for the paper and for product confidence.
- Include runtime comparison because FLEX is one of the few workflows still slightly slower than CR.
- Evaluate whether targeted acceleration in feature processing or related stages is worth the extra implementation complexity.
- Exit criterion:
  - clear parity status against CR9 and a yes/no recommendation on performance optimization scope.

### 2. Harden perturb mixed-chemistry support

- Treat the core behavior as implemented and E2E-validated rather than missing.
- Focus this work on updating regression coverage:
  - mixed `NXT + TRU` feature-library merge
  - 2-column mapping edge cases
  - chemistry precedence matrix
  - remaining filtered-cell / multi-library validation gaps
- Exit criterion:
  - mixed-chemistry support is covered by current unit/smoke tests strongly enough
    that future regressions are unlikely to slip through unnoticed.

### 3. Resolve remaining cell-calling and intronic-parity questions

- Test ambient-window sizing sensitivity for smaller datasets.
- Investigate GeneFull intronic-count inflation vs Cell Ranger.
- Exit criterion:
  - either validated defaults, or explicit caveats/documented behavior.

### 4. Finish poly-G auto mode

- Replace the current heuristic aliasing of `auto` to CellRanger4 mode with data-driven detection.
- Exit criterion:
  - `auto` is dataset-driven and tested.

### 5. Clean up CB/UB and FLEX raw-tag behavior

- Relax CB/UB coupling for non-Flex only.
- Confirm Flex policy remains strict.
- Re-enable and validate FLEX `ZG/ZX` tags once the tag path is stable.
- Exit criterion:
  - tag emission policy is documented and regression-tested.

## Priority 2: Throughput and Architecture

### 6. Implement perturb dynamic permit scheduling

- Phase 1: STAR-only dynamic active-thread limit behind a feature flag.
- Phase 2: queue-aware producer/consumer permit rebalancing.
- Validate on A375 first, then UCSF full data.
- Exit criterion:
  - measurable throughput gain with parity preserved.

### 7. Finish the QC graphics consolidation

- Define the shared QC graphics core.
- Keep Plotly HTML+JSON as the common format.
- Audit downstream scripts for removed PNG expectations.
- Exit criterion:
  - one reusable QC graphics path instead of duplicated per-module outputs.

### 8. Harden feature-offset and decompression behavior

- Add the read-level feature-offset preflight if the false-positive rate is acceptable.
- Confirm native-gz and `--readFilesCommand` paths remain behaviorally identical on representative inputs.
- Exit criterion:
  - no hidden behavior split between decompression modes.

## Priority 3: Coverage Expansion

### 9. Finish missing fixture-dependent smokes

- Run the outstanding memory-leak repair smokes once fixtures are available.

### 10. Expand CRISPR feature-calling threshold validation

- Validate `--crMinUmi` on lineage, FLEX, and antibody datasets.
- Convert the resulting recommendations into tests and user guidance.

## Suggested Execution Order

1. FLEX parity versus Cell Ranger 9.0
2. Mixed-chemistry perturb support and tests
3. Cell-calling and GeneFull parity investigations
4. Poly-G auto detection
5. CB/UB and FLEX raw-tag cleanup
6. Dynamic permit scheduling
7. QC graphics consolidation
8. Fixture and assay coverage expansion

## Explicitly Out of This Plan

- Release-compatibility tarball work
- Installer-bundle environment selection
- Published `v0.50.3` release smoke
- PE bulk `TranscriptVB` vs Salmon residual deep dive (`~0.06%`)
- Expanded real-dataset PE multisample benchmark beyond the current regression smoke
- SLAM over-dispersion parameter sweep until another SLAM-seq dataset is available

Those are complete enough for current prerelease use and should not block the paper/debugging work. The PE bulk residual is currently treated as an acceptable small policy/model difference unless a materially larger dataset-specific gap appears. The PE multisample batch-vs-sequential path remains covered by the checked-in smoke test and can be revisited later if a broader formal benchmark is needed. The SLAM over-dispersion work is also deferred because the initial conservative implementation made the term effectively negligible and there is no urgency to complicate the analysis without an additional independent dataset to justify the extra effort.
