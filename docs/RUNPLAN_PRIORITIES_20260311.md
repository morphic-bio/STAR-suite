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
  - any remaining multi-library validation gaps beyond the filtered-cell smoke
- Exit criterion:
  - mixed-chemistry support is covered by current unit/smoke tests strongly enough
    that future regressions are unlikely to slip through unnoticed.

### 3. Finish poly-G auto mode

- Replace the current heuristic aliasing of `auto` to CellRanger4 mode with data-driven detection.
- Exit criterion:
  - `auto` is dataset-driven and tested.

## Priority 3: Coverage Expansion

### 4. Finish missing fixture-dependent smokes

- Run the outstanding memory-leak repair smokes once fixtures are available.

### 5. Expand CRISPR feature-calling threshold validation

- Validate `--crMinUmi` on lineage, FLEX, and antibody datasets.
- Convert the resulting recommendations into tests and user guidance.

## Suggested Execution Order

1. FLEX parity versus Cell Ranger 9.0
2. Mixed-chemistry perturb support and tests
3. Poly-G auto detection
4. Fixture and assay coverage expansion

## Explicitly Out of This Plan

- Release-compatibility tarball work
- Installer-bundle environment selection
- Published `v0.50.3` release smoke
- PE bulk `TranscriptVB` vs Salmon residual deep dive (`~0.06%`)
- Expanded real-dataset PE multisample benchmark beyond the current regression smoke
- SLAM over-dispersion parameter sweep until another SLAM-seq dataset is available
- Further CR-compatibility cleanup for ambient-window sizing and GeneFull intronic edge cases
- Additional feature-offset preflight work until there is a broader formal preflight framework
- `--readFilesCommand` compatibility work beyond the current legacy support path
- QC graphics consolidation until we are actively working on new graphical output again
- CB/UB semantics cleanup unless stats/reporting problems make it worth revisiting
- Legacy FLEX `ZG/ZX` cleanup until there is time to remove or replace the old workaround path
- Perturb dynamic permit scheduling unless a new workload exposes a fresh scheduling bottleneck

Those are complete enough for current prerelease use and should not block the paper/debugging work. The PE bulk residual is currently treated as an acceptable small policy/model difference unless a materially larger dataset-specific gap appears. The PE multisample batch-vs-sequential path remains covered by the checked-in smoke test and can be revisited later if a broader formal benchmark is needed. The SLAM over-dispersion work is also deferred because the initial conservative implementation made the term effectively negligible and there is no urgency to complicate the analysis without an additional independent dataset to justify the extra effort. The remaining CR-compatibility cleanup is likewise deferred because compatibility is already much stronger than before and the expected return is low unless a concrete accuracy problem points back to those areas. Additional feature-offset preflight work is deferred because there is already a working offset-detection phase, and any deeper consistency check belongs in a more formal preflight framework. `--readFilesCommand` remains for flexibility and legacy support, but active development should target the direct internal gzip path rather than chasing full compatibility for every command-stream combination. QC graphics consolidation is deferred because only the `process_features` heatmap consolidation/optimization remains, and that work is better done the next time we add or revise graphical output. CB/UB and legacy `ZG/ZX` cleanup are also deferred unless they create real stats or maintenance pain. Perturb dynamic permit scheduling is likewise considered complete unless new workloads show another clear bottleneck.
