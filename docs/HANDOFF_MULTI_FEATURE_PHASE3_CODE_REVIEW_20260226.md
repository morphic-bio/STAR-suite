# Handoff: Multi-Feature Phase 3 Code Review

Date: 2026-02-26  
Branch: `multi-feature`  
Reviewer: Codex  
Scope: Review current Phase 3 implementation status and hand off actionable next steps.

## Executive Summary

Phase 3 is **not implemented yet** in the current diff. The active code changes are still Phase 0-2 focused (config parsing, per-library refs/IDs, and assign wiring).  
Merge/provenance hardening objectives from `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md` remain open.

## Files Reviewed

- `core/legacy/source/PfMultiConfig.h`
- `core/legacy/source/PfMultiConfig.cpp`
- `core/legacy/source/PfMultiProcess.cpp`
- `core/legacy/source/PfMultiFeatureSpecs.h`
- `core/legacy/source/PfMultiMexStub.cpp`
- `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md`
- `docs/HANDOFF_MULTI_FEATURE_PHASE3_READY_20260226.md`
- `tests/multi_feature/test_multi_feature_config.sh`
- `tests/multi_feature/test_data_driven_specs.sh`

## Findings (Ordered by Severity)

### High: Phase 3 provenance outputs are missing

- No per-library manifest writer exists in PF merge path (`core/legacy/source/PfMultiProcess.cpp`).
- `FeatureRun` does not carry library provenance fields needed for machine-readable manifests (`core/legacy/source/PfMultiProcess.cpp:1308`).
- No top-level Phase 3 summary table is emitted to `Log.out` after merge.

Impact: Runbook Phase 3 acceptance criteria are not met, and downstream audits cannot reconstruct exact per-library inputs/settings.

### High: Feature-type labeling risk still exists at merge stub boundary

- Merge path still calls `PfMultiMexStub::processAssignOutput(run.assignOut, run.featureRefPath, run.featureType, ...)` (`core/legacy/source/PfMultiProcess.cpp:1747`).
- In stub generation, `features.tsv` keeps CSV `feature_type` when present, only falling back to `defaultType` when empty (`core/legacy/source/PfMultiMexStub.cpp:247`).

Impact: Libraries declared as CRISPR can be mislabeled as `Custom` if their input CSV carries `feature_type=Custom`, breaking CRISPR feature-type-dependent downstream logic.

### High: Partial-merge behavior can silently degrade outputs

- Failed feature MEX reads are logged and skipped (`core/legacy/source/PfMultiProcess.cpp:1894`), with merge continuing if at least one feature MEX remains.
- If all feature MEX reads fail, code logs warning and returns success (`core/legacy/source/PfMultiProcess.cpp:1900`).

Impact: Pipeline can complete with incomplete feature coverage without hard failure, violating strict fail-fast expectations for this phase.

### Medium: Assign failure error is not library-identifying

- On assign failure, error message includes only library type (`core/legacy/source/PfMultiProcess.cpp:1700`), not `star_library_id`, sample, FASTQ, or reference path.

Impact: Troubleshooting multi-library failures is slower and ambiguous when multiple libraries share a type.

## Validation Performed

- Build:
  - `make -C core/legacy/source -j8 STAR` (pass)
- Fast multi-feature tests:
  - `bash tests/multi_feature/test_multi_feature_config.sh` (21/21 pass)
  - `bash tests/multi_feature/test_data_driven_specs.sh` (16/16 pass)

Note: These tests validate Phase 0-2 behavior. They do not close Phase 3 provenance/merge-hardening criteria.

## Phase 3 Implementation Plan for Next Agent

1. Add `FeatureRun` provenance fields:
   - `libraryId`, `sampleName`, `resolvedFastq`, `resolvedWhitelist`, `resolvedFeatureRef`, `assignStatus`, `assignReturnCode`.
2. Emit per-library manifest files under each assign root:
   - Suggested path: `{assignOut}/pf_library_manifest.tsv` (or JSON if preferred).
3. Emit top-level summary table in `Log.out`:
   - One row per configured library, including `status=SUCCESS|FAILED|SKIPPED`.
4. Harden merge policy:
   - Default fail-fast if any configured feature library fails assign or MEX read.
   - If controlled-skip mode is added later, gate it explicitly and log it loudly.
5. Fix feature-type label authority:
   - Enforce declared library type (config-driven) at merge stub generation, or add explicit remap before merge so `CRISPR Guide Capture` remains CRISPR in final MEX.
6. Add Phase 3 tests:
   - Provenance manifest existence/content checks.
   - Failure-in-one-library behavior (assert fail-fast default).
   - CRISPR+Custom mixed refs where CSV `feature_type` conflicts with declared type.

## Recommended Acceptance Checks (After Fixes)

- Multi-feature E2E smoke includes both perturb and lineage features with correct types.
- `Log.out` contains a deterministic per-library summary table.
- Per-library manifest captures all resolved inputs and chemistry/match mode.
- A375 single-library CRISPR regression still passes unchanged.

