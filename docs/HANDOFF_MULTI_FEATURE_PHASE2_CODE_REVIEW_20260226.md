# Handoff: Multi-Feature Libraries Phase 2 Code Review

Date: 2026-02-26  
Branch: `multi-feature`  
Reviewer: Codex  
Scope: Review **Phase 2** (per-library assign wiring + chemistry anchoring), not Phase 3.

## Executive Summary

Phase 2 implementation is largely correct and test-backed in this workspace.  
Per-library `star_feature_ref`, `star_library_id`, and per-library chemistry composition are wired through prep and assign paths.

No blocking correctness bug was found for the intended Phase 2 contract. Three edge-case gaps were identified and have all been resolved (see Findings section below).

## What Was Reviewed

- `core/legacy/source/PfMultiConfig.cpp`
- `core/legacy/source/PfMultiConfig.h`
- `core/legacy/source/PfMultiProcess.cpp`
- `core/legacy/source/PfMultiFeatureSpecs.h`
- `tests/multi_feature/test_multi_feature_config.sh`
- `tests/multi_feature/test_data_driven_specs.sh`
- `tests/multi_feature/test_no_global_ref_guard.sh`
- `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md`

## Phase 2 Contract Check

### 1) Per-library feature reference wiring

Pass. `star_feature_ref` takes precedence and bypasses global filtering in prep:
- `core/legacy/source/PfMultiProcess.cpp:989`
- `core/legacy/source/PfMultiProcess.cpp:991`
- `core/legacy/source/PfMultiProcess.cpp:999`

Global ref guard is correctly conditional (only required when needed):
- `core/legacy/source/PfMultiProcess.cpp:861`
- `core/legacy/source/PfMultiProcess.cpp:867`
- `core/legacy/source/PfMultiProcess.cpp:873`

### 2) Per-library assign invocation

Pass. Each prepared library carries its own resolved inputs and output root, then runs assign independently:
- `core/legacy/source/PfMultiProcess.cpp:980`
- `core/legacy/source/PfMultiProcess.cpp:982`
- `core/legacy/source/PfMultiProcess.cpp:1321`
- `core/legacy/source/PfMultiProcess.cpp:1611`

### 3) Per-library chemistry composition and autodetect anchoring

Pass. Explicit chemistries bypass autodetect, `auto` libraries use detected match mode with anchored namespace logic:
- `core/legacy/source/PfMultiProcess.cpp:1003`
- `core/legacy/source/PfMultiProcess.cpp:1015`
- `core/legacy/source/PfMultiProcess.cpp:1372`
- `core/legacy/source/PfMultiProcess.cpp:1711`
- `core/legacy/source/PfMultiProcess.cpp:1717`

### 4) Config parsing + validation for new columns

Pass. Parser supports aliases, auto-generates IDs, validates duplicate/sanitized-collision IDs, resolves relative ref paths:
- `core/legacy/source/PfMultiConfig.cpp:170`
- `core/legacy/source/PfMultiConfig.cpp:172`
- `core/legacy/source/PfMultiConfig.cpp:359`
- `core/legacy/source/PfMultiConfig.cpp:376`
- `core/legacy/source/PfMultiConfig.cpp:394`
- `core/legacy/source/PfMultiConfig.cpp:398`

## Findings (Ordered by Severity)

### Medium: Empty `feature_types` can be silently skipped — RESOLVED

- Data-driven spec builder drops libraries whose normalized type is empty (`norm.empty() -> continue`).
- This means a malformed library row can be silently ignored instead of failing fast.

Evidence:
- `core/legacy/source/PfMultiFeatureSpecs.h` (`if (norm.empty()) continue;`)

Recommendation:
- In prep, fail fast when a non-GEX library has empty `feature_types`, and include row/sample/library_id in the error.

**Fix**: Changed `continue` to `throw runtime_error(...)` with sample and starLibraryId in the message.

### Medium: Assign failure message lacks library identity context — RESOLVED

- Current throw only reports `libraryType`, which is insufficient when multiple libs share type.

Evidence:
- `core/legacy/source/PfMultiProcess.cpp:1703`

Recommendation:
- Include `libraryId`, `sampleName`, `refPath`, and `resolvedFastq` in failure message.

**Fix**: Error message now includes `type`, `libraryId`, `sample`, `featureRef`, and `fastq`.

### Low: Expected per-library reference path emits warning-like language — RESOLVED

- For `star_feature_ref`, log says “feature reference not filtered … using full reference,” which is expected behavior in this mode and can look like a problem.

Evidence:
- `core/legacy/source/PfMultiProcess.cpp:1028`

Recommendation:
- Downgrade wording to `NOTICE` when `star_feature_ref` is explicitly provided.

**Fix**: When `star_feature_ref` is set, emits `NOTICE: per-library star_feature_ref provided ... type-based filtering skipped as expected`. The `WARNING` is preserved only when no per-library ref was set and filtering genuinely failed.

## Validation Run (This Review)

- Build: `make -C core/legacy/source -j8 STAR` (pass)
- `bash tests/multi_feature/test_multi_feature_config.sh` (21/21 pass)
- `bash tests/multi_feature/test_data_driven_specs.sh` (16/16 pass)
- `bash tests/multi_feature/test_no_global_ref_guard.sh` (4/4 pass)

## Resolution Summary

All three findings have been addressed:

| Finding | Severity | Status | File |
|---------|----------|--------|------|
| Empty `feature_types` silently skipped | Medium | **Resolved** | `PfMultiFeatureSpecs.h` |
| Assign failure message lacks context | Medium | **Resolved** | `PfMultiProcess.cpp` |
| Per-library ref emits warning-like log | Low | **Resolved** | `PfMultiProcess.cpp` |
| GEX detection is substring-based, misclassifies types containing `gex` | Medium | **Resolved** | `PfMultiFeatureSpecs.h`, `PfMultiConfig.cpp`, `PfMultiProcess.cpp` |

Build: clean (exit 0). All 64 fast unit tests pass (21 + 16 + 4 + 8 + 15).

## Next-Agent Action List

1. ~~Add strict validation for empty `feature_types` on feature libraries.~~ Done.
2. ~~Improve assign failure diagnostics with `libraryId` and concrete paths.~~ Done.
3. ~~Tidy log wording for expected per-library-ref path behavior.~~ Done.
4. ~~Proceed to Phase 3 merge/provenance hardening.~~ Done. See `HANDOFF_MULTI_FEATURE_PHASE3_MERGE_PROVENANCE_20260226.md`.

