# Handoff: NXT/TRU Namespace-Agnostic Merge

**Date:** 2026-02-24  
**Branch:** `autodetect-NXT-TRU` (based on `core-alignment-threads-integration`)  
**Scope:** Make the pf-multi merge pipeline agnostic to NXT/TRU barcode chemistry  
**Status:** Implemented, tested, reviewed and revised

---

## Problem Statement

The pf-multi pipeline assumed a **single global NXT/TRU chemistry** for all libraries.
For 3' HT experiments where GEX libraries use TRU barcodes and feature libraries use
NXT barcodes (or vice versa), this caused:

1. **Merge barcode mismatch**: `mergeMex` joins barcodes by exact string match. NXT
   feature barcodes never match TRU GEX barcodes → feature counts silently dropped.
2. **Uniform translation corruption**: `writeCombinedMex` applied the same NXT→TRU
   translation to ALL barcodes, corrupting GEX barcodes that were already in TRU.
3. **No per-library chemistry tracking**: `PfPreparedFeatureLibrary`, `FeatureRun`,
   and `AssignOptions` had no chemistry field, making per-library handling impossible.
4. **Unsafe auto-detection**: `detectChemistryFromWhitelistPath` returned NXT for any
   2-column whitelist, even when zero rows matched the complement rule. With boundary
   normalization, this could silently flip barcodes incorrectly.

A full audit of all affected locations is in
`docs/AUDIT_NXT_TRU_PER_LIBRARY_CHEMISTRY_20260224.md`.

---

## Design Decision

**Normalize all barcode namespaces to TRU at the merge boundary.**

- Each library's `assignBarcodes` works in read-space (unchanged).
- Before `mergeMex`, a new helper `normalizeMexBarcodesToTru()` translates NXT→TRU
  for any NXT-sourced MexData (complementing positions 7-8, 0-indexed).
- After normalization, all inputs are in TRU → merge joins work correctly.
- `writeCombinedMex` receives `inputChemistry="TRU"` and `outputChemistry` from the
  user's `--crOutputChemistry` setting (default TRU). When output is TRU this is a
  no-op; when output is NXT, the final write translates TRU→NXT.

This design was chosen because:
- It puts all namespace logic at **one point** (the merge boundary).
- It handles same-chemistry runs (all NXT or all TRU) correctly.
- TRU is the default output namespace, simplifying the common path.
- It is forward-compatible with future per-library auto-detection.

---

## Current Scope and Known Limitations

This change addresses **same-whitelist** scenarios where all libraries share a single
whitelist (the common case). Specifically:

- **Per-library `effectiveChem` is currently set from the global context** (line 919).
  True mixed-chemistry (different whitelists per library) is not yet supported — it
  requires per-library whitelist fields in `PfMultiConfig` and per-library whitelist
  resolution in `PfMultiProcess`.
- **`AssignOptions.translateNxt` is wired through but never activated** from the
  pf-multi caller path. It exists as infrastructure for the future auto-detection
  phase, where `assignBarcodes` will set `translate_NXT` per-library based on
  bootstrap-phase detection. Currently the field defaults to `false` and the
  normalization happens entirely at the merge boundary.
- **`PfMultiMexStub::processAssignOutput`** still uses a single global whitelist path.

---

## Files Changed

### `core/legacy/source/PfMultiAssign.h`

Added `bool translateNxt = false` to `AssignOptions`. This is infrastructure for
future per-library auto-detection. Not activated from the pf-multi caller path.

### `core/legacy/source/PfMultiAssign.cpp`

1. **`applyAssignOptions`**: Calls `pf_config_set_translate_nxt(cfg, 1)` when
   `options.translateNxt` is true.
2. **`writeApiRunSummary`**: Logs `translateNxt=` value for debugging.

### `core/legacy/source/PfMultiProcess.cpp`

Seven changes:

1. **`PfPreparedFeatureLibrary`**: Added `string effectiveChem` field. Set from the
   global `context.effectiveChem` during `buildPfMultiPreparedContext`.

2. **`FeatureRun`**: Added `string effectiveChem` field. Populated from
   `preparedLib.effectiveChem` when building the feature run list.

3. **`normalizeMexBarcodesToTru()`**: New static helper that complements positions
   7-8 for NXT→TRU. No-op when `fromChem != "NXT"`.

4. **`detectChemistryFromWhitelistPath()`**: Fixed a bug where any 2-column whitelist
   was classified as NXT. Now requires ≥80% of sampled rows to match the center-2bp
   complement rule before returning NXT. Below this threshold, falls through to
   filename/default detection.

5. **`FeatureMexEntry` struct**: Feature MEX data is now paired with its chemistry
   and feature type in a struct, so a failed read mid-list does not cause index
   misalignment between `featureDataVec` and `featureRuns`.

6. **Merge boundary normalization**: Before `mergeMex`, normalizes all inputs
   (GEX raw, GEX filtered, and each feature entry) to TRU.

7. **`writeCombinedMex` calls**: Pass `inputChemistry="TRU"` and
   `outputChemistry=outputChem` (from `--crOutputChemistry`, default TRU).
   `--crOutputChemistry NXT` correctly triggers TRU→NXT translation at write time.

8. **`normalizeBarcodeVecToTru()`**: New static helper (parallel to
   `normalizeMexBarcodesToTru`) for plain `vector<string>` barcode lists.
   `getFilteredBarcodesFromSolo` is now called with `useOutputNamespace=false`
   to always fetch read-space barcodes (unambiguous namespace regardless of
   1-column vs 2-column whitelist). These are normalized to TRU immediately
   after retrieval. MEX-derived fallback barcodes (`computeObservedGexBarcodes`
   from already-normalized data) are already in TRU and skip normalization,
   avoiding double-translation.

---

## What Did NOT Change

- **`PfMultiMerge.cpp`**: No changes. `needsNamespaceTranslation` and
  `translateNxtMiddleTwoBasesInplace` remain and are correctly activated when
  `outputChem` differs from `"TRU"`.
- **`assignBarcodes.c` / `pf_api.c`**: No changes. `translate_NXT` global is still
  set to 0 by default. Per-library translation happens at the STAR merge boundary.
- **`call_features.c` / `gmm.c`**: Barcode-agnostic — no changes needed.
- **`PfMultiConfig`**: No per-library chemistry fields yet — deferred to the
  auto-detection phase.

---

## Data Flow (Before vs After)

### Before (global chemistry)

```
assignBarcodes(NXT whitelist) → Feature MEX [NXT barcodes]
STAR Solo(NXT whitelist)      → GEX MEX    [NXT barcodes]
                                    ↓
                              mergeMex (NXT=NXT, match OK)
                                    ↓
                        writeCombinedMex(NXT→outputChem)
                                    ↓
                           output [outputChem barcodes]
```

Problem: fails when GEX=TRU and Feature=NXT (3' HT).
Problem: any 2-column whitelist could be mis-detected as NXT.

### After (normalize at boundary)

```
assignBarcodes(any whitelist) → Feature MEX [read-space barcodes]
STAR Solo(any whitelist)      → GEX MEX    [read-space barcodes]
                                    ↓
                    normalizeMexBarcodesToTru(each source, per-lib chem)
                                    ↓
                              mergeMex (TRU=TRU, match OK)
                                    ↓
                        writeCombinedMex(TRU→outputChem)
                                    ↓
                           output [outputChem barcodes]
```

Handles same-whitelist runs correctly (all-TRU or all-NXT).
Provides the per-library chemistry plumbing for future mixed-whitelist support.

---

## Validation

### A375 (TRU path — full end-to-end)

- Ran `tests/test_cr_compat_crispr_calling.sh` twice (pre- and post-review fixes)
- **Result**: PASSED both times
- 1,186 cells, 1,051 single-feature calls, 35 multi-feature
- `normalizeMexBarcodesToTru` is a no-op for TRU inputs → zero behavioral change

### UCSF iPSC2 (NXT path — barcode translation parity)

- Compared against prior bootstrap run output (7,286 filtered barcodes)
- Applied `normalizeMexBarcodesToTru` logic (complement positions 7-8) to native
  NXT barcodes from prior run
- **Result**: 7,286/7,286 (100%) match with the TRU output barcodes
- All output barcodes are in the translated GEX whitelist (100%)

### Build

- Clean build with `make -C core/legacy/source -j8 STAR`
- No new warnings in modified files

---

## Behavioral Differences

| Aspect | Before | After |
|--------|--------|-------|
| Merge namespace | Read-space (NXT or TRU) | Always TRU |
| `writeCombinedMex` inputChem | `effectiveChem` (NXT or TRU) | Always TRU |
| `writeCombinedMex` outputChem | `outputChem` (user setting) | `outputChem` (user setting, unchanged) |
| `--crOutputChemistry NXT` | Honored | Honored (TRU→NXT at write) |
| `barcodes.native.tsv.gz` | Emitted when input≠output chem | Emitted when outputChem≠TRU |
| 2-col whitelist auto-detection | Any 2-col → NXT | Requires ≥80% complement-rule match |
| Feature MEX read failure | Index drift risk | Chemistry paired with data, no drift |

---

## Review Fixes Applied

These issues were identified during code review and fixed in this revision:

### High: `--crOutputChemistry` was ignored (hardcoded TRU/TRU)

`writeCombinedMex` was called with `"TRU", "TRU"`, ignoring the user's
`--crOutputChemistry` setting. Fixed to pass `"TRU", outputChem`, where `outputChem`
comes from `prepared.outputChem` (which respects `--crOutputChemistry`).

### High: Auto-detection misclassified non-NXT 2-column whitelists

`detectChemistryFromWhitelistPath` returned "NXT" for any 2-column whitelist, even
if zero rows matched the center-2bp complement rule. With boundary normalization this
could silently corrupt barcodes. Fixed to require ≥80% complement-rule agreement
before classifying as NXT; below that threshold, falls through to filename/default.

### Medium: Index alignment drift on partial feature MEX read failure

Feature MEX data was pushed to a `vector<MexData>` and chemistry was looked up by
index from the parallel `featureRuns` vector. If a MEX read failed mid-list,
`featureDataVec` would have fewer entries than `featureRuns`, causing the wrong
chemistry to be applied. Fixed by introducing `FeatureMexEntry` struct that pairs
each MexData with its chemistry and feature type at read time.

### Medium: Filtered-barcode namespace mismatch for 1-column NXT whitelists

`getFilteredBarcodesFromSolo(..., useOutputNamespace=true)` pulls from `cbWLstrOut`.
For 1-column whitelists, `cbWLstrOut == cbWLstr` (read-space), so these remain in NXT.
After merge-boundary normalization, matrix barcodes are TRU, but the filter set passed
to `writeCombinedMex` was NXT → no barcodes match → empty filtered output. Initial fix
unconditionally normalized all filtered barcodes, which caused double-translation for
2-column whitelists and MEX-derived fallbacks (already in TRU). Final fix: changed
`getFilteredBarcodesFromSolo` to `useOutputNamespace=false` (read-space), normalize
immediately after, and leave MEX-derived fallbacks (already TRU) untouched.

### Low: `translateNxt` not activated from pf-multi path

Noted in the handoff as infrastructure for future use. Not a bug — the normalization
currently happens at the merge boundary, not in `assignBarcodes`.

---

## Remaining Work (Future PRs)

1. **Per-library whitelist and chemistry in `PfMultiConfig`**:
   - Add `whitelist` and `chemistry` fields to `LibraryEntry`.
   - Support different whitelists for GEX vs feature libraries (3' HT pattern).
   - This is a prerequisite for true mixed-chemistry support.

2. **Auto-detect NXT/TRU per library during bootstrap phase**:
   - During `feature_mode_bootstrap_reads` window, trial-match barcodes both as-is
     and after NXT translation; tally hit rates.
   - In `feature_mode_finalize()`, decide NXT vs TRU from match ratio.
   - Export via new `pf_get_detected_chemistry()` API.

3. **Synthetic test harness for mixed-chemistry**:
   - Generate synthetic FASTQs with NXT and TRU barcodes in separate libraries.
   - Verify merged MEX has correct counts in TRU namespace.

4. **Edge-case tests to add**:
   - `--crOutputChemistry NXT` should produce NXT barcodes in output.
   - 2-column whitelist with low/zero NXT-rule agreement should not auto-force NXT.
   - Feature MEX partial-read failure should preserve correct chemistry mapping.

---

## Review Checklist

- [x] `normalizeMexBarcodesToTru` correctly complements positions 7-8 (0-indexed)
- [x] Per-library `effectiveChem` is set from global context (correct for same-whitelist runs)
- [x] GEX raw, GEX filtered, and feature MexData are all normalized before merge
- [x] `writeCombinedMex` respects `--crOutputChemistry` (passes `outputChem`, not hardcoded TRU)
- [x] Auto-detection requires ≥80% complement-rule agreement for NXT classification
- [x] Feature MEX reads paired with chemistry (no index drift on partial failure)
- [x] Solo filtered barcodes normalized to TRU before `writeCombinedMex` (1-column NXT safe)
- [x] No new warnings in build
- [x] A375 TRU path passes (full end-to-end, no behavioral change)
- [x] UCSF iPSC2 NXT path produces identical TRU barcodes (7,286/7,286)
- [x] `docs/todos` updated with DONE items and remaining auto-detection TODOs
- [x] Audit document updated with resolution section
- [x] Handoff accurately describes current scope and known limitations
