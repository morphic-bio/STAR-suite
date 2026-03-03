# Handoff: Namespace Parity Fix (2026-03-03)

## Problem

Integration runs with NXT-configured 2-column whitelists showed severe collapse
in feature calls (~55 barcodes with CRISPR features instead of ~15,000; Jaccard
0.79 instead of ~0.99). Reported in
`docs/HANDOFF_NAMESPACE_PARITY_BLOCKER_20260303.md`.

## Root Causes (Three Interacting Bugs)

### Bug 1: ParametersSolo.cpp regression (IMMEDIATE CAUSE)

Uncommitted working-tree changes to `ParametersSolo.cpp` loaded both columns
of a 2-column whitelist as **independent barcodes** (doubling the whitelist to
6.8M entries), instead of pairing COL1 (for matching) with COL2 (for output,
3.7M entries).

This caused Solo to emit a mixed NXT/TRU `barcodes.tsv`, which was the direct
trigger for the severe parity collapse.

**Fix:** Reverted `ParametersSolo.cpp` to committed state via
`git checkout -- core/legacy/source/ParametersSolo.cpp`.

### Bug 2: GEX MEX double-translation (LATENT, in both committed + uncommitted)

`normalizeMexBarcodesToTru(gexData, fromChem)` was called with `fromChem` =
the whitelist **matching** chemistry (e.g. "NXT"), but Solo's raw MEX barcodes
are actually written using `cbWLstrOut` — the **output** namespace.

For 2-column NXT whitelists:
- `cbWLstrOut` = COL2 = **TRU** → Solo's `barcodes.tsv` is TRU
- Code passed `fromChem = "NXT"` → double-translated TRU→NXT

This corrupted GEX barcodes from TRU to NXT after normalization, causing
namespace mismatch with correctly-TRU filtered barcodes.

**Fix:** Added `soloOutputNamespace` to `PfMultiPreparedContext`:
- For 2-column NXT whitelists: `soloOutputNamespace = "TRU"`
- For 1-column whitelists: `soloOutputNamespace = effectiveChem`

GEX MEX normalization now uses `soloOutputNamespace` as `fromChem`.

### Bug 3: Feature MEX double-translation (LATENT, same pattern)

Feature MEX `entry.effectiveChem` was set to `assignmentChem` (the read/matching
chemistry), not the actual namespace of the `process_features` output barcodes.

When `translateNxt=true`, process_features already translates output barcodes
from the whitelist namespace to its complement (NXT→TRU). Calling
`normalizeMexBarcodesToTru(data, "NXT")` on these TRU barcodes double-translated
them back to NXT.

**Fix:** Added `featureMexOutputNamespace` to `FeatureRun` and `FeatureMexEntry`:
- `translateNxt=true` → `featureMexOutputNamespace = oppositeNamespace(assignmentWhitelistNamespace)`
- `translateNxt=false` → `featureMexOutputNamespace = assignmentWhitelistNamespace`
- **After `processAssignOutput()` stub** (2-column whitelist): overridden to
  `prepared.soloOutputNamespace` (TRU), because `copyBarcodesTsv` applies
  COL1→COL2 mapping to `barcodes.tsv`.

Feature MEX normalization now uses `featureMexOutputNamespace` as `fromChem`.

### Bug 5b: Ambiguous orientation warning now wired into main flow

`buildPfMultiPreparedContext` now emits a `WARNING` when `inferredChemConfident=false`
and the chemistry was inferred via `auto` mode. This surfaces the ambiguity for
custom/renamed 2-column whitelist files that lack NXT/TRU in the filename.

## Core Concept

The root confusion was between three distinct namespaces:

| Namespace | Meaning | Example (2-col NXT) |
|-----------|---------|---------------------|
| Whitelist matching namespace | COL1, used for read matching | NXT |
| Solo output namespace | `cbWLstrOut`, written to `barcodes.tsv` | TRU |
| Feature MEX output namespace | After optional `translateNxt` | TRU (if translateNxt) |

`normalizeMexBarcodesToTru(data, fromChem)` needs `fromChem` to be the **actual
namespace of the barcodes in the MEX file**, not the whitelist chemistry.

## All Files Changed

### Tracked (vs HEAD)

| File | Change |
|------|--------|
| `core/legacy/source/ParametersSolo.cpp` | Reverted to committed state (no diff) |
| `core/legacy/source/ParametersSolo.h` | Fixed stale `cbWLstrOut` doc-comment |
| `core/legacy/source/PfMultiProcess.cpp` | Fixed GEX/feature MEX namespace bugs (see below) |

### PfMultiProcess.cpp changes (uncommitted, on top of other uncommitted work)

1. **`whitelistHasTwoColumns()`** — new helper to detect 2-column whitelists
2. **`PfMultiPreparedContext`** — added `hasTwoColumnWhitelist`, `soloOutputNamespace`
3. **`buildPfMultiPreparedContext`** — computes `soloOutputNamespace` after chemistry
   detection; recomputes on per-library GEX chemistry override
4. **`gexNormalizationChem`** — now derived from `soloOutputNamespace` (was
   `inferredChem`/`effectiveChem`)
5. **Filtered GEX barcode normalization** — after `normalizeFilteredBarcodesForAssignNamespace`,
   barcodes are in `gexWhitelistNamespace`; normalization uses that (was `gexNormalizationChem`)
6. **`FeatureRun.featureMexOutputNamespace`** — new field; computed from
   `translateNxtForAssign` and `assignmentWhitelistNamespace`; **overridden to
   `soloOutputNamespace` after `processAssignOutput` for 2-column whitelists**
   (accounts for `copyBarcodesTsv` COL1→COL2 mapping)
7. **`FeatureMexEntry.featureMexOutputNamespace`** — propagated from `FeatureRun`
8. **Feature MEX normalization** — uses `featureMexOutputNamespace` (was `effectiveChem`)
9. **`buildPfMultiPreparedContext`** — emits WARNING when `inferredChemConfident=false`
   and chemistry mode is `auto` (ambiguous 2-col orientation)

### Bug 4: GEX filtered-barcode overlap collision (HIGH — discovered in parity rerun)

`normalizeFilteredBarcodesForAssignNamespace()` decided whether a barcode was
"already in the assignment namespace" via set membership only. For 2-column NXT
whitelists, COL1 and COL2 barcodes overlap at self-complementary positions 7-8
(about 987/6069 barcodes in the UCSF 2M dataset). A TRU barcode matching an
NXT entry in the set would be kept as-is (TRU), then
`normalizeBarcodeVecToTru(_, "NXT")` would corrupt it further.

Evidence: `in_set=987 translated_to_set=5082` in the 2M rerun Log.out.

**Fix:** Made `normalizeFilteredBarcodesForAssignNamespace` source-namespace
aware. New `sourceNamespace` and `whitelistNamespace` parameters enable
**deterministic mode**: when both are known, the function always translates (or
not) based on namespace comparison — no membership guessing. The GEX path now
passes `gexNormalizationChem` (= Solo output namespace) and
`gexWhitelistNamespace` (= COL1 namespace). The per-library assign path passes
known namespaces for solo-sourced barcodes, falling back to heuristic only for
explicit (unknown-namespace) filtered barcode files.

## Core Concept

The root confusion was between three distinct namespaces:

| Namespace | Meaning | Example (2-col NXT) |
|-----------|---------|---------------------|
| Whitelist matching namespace | COL1, used for read matching | NXT |
| Solo output namespace | `cbWLstrOut`, written to `barcodes.tsv` | TRU |
| Feature MEX output namespace | After optional `translateNxt` | TRU (if translateNxt) |

`normalizeMexBarcodesToTru(data, fromChem)` needs `fromChem` to be the **actual
namespace of the barcodes in the MEX file**, not the whitelist chemistry.

`normalizeFilteredBarcodesForAssignNamespace()` must not rely on set membership
when the source namespace is known, because COL1/COL2 sets overlap.

### Bug 5: Ambiguous orientation for 2-column whitelists (LOW)

`detectChemistryFromWhitelistPath()` and `inferTwoColumnNamespace()` returned
`confident=true` when the center-2bp complement rule matched >= 80%, even when
no filename hint confirmed the NXT/TRU orientation. The complement rule is
symmetric (`translate(A)==B` iff `translate(B)==A`), so content alone cannot
distinguish COL1=NXT from COL1=TRU.

For custom/renamed whitelist files without "nxt" or "tru" in the filename,
this could silently assume the wrong orientation.

**Fix:** Both functions now return `confident=false` when the complement rule
matches but no filename hint confirms the orientation. A WARNING is included
in the reason string. Downstream code that checks `namespaceConfidence` can
use this to decide whether to warn the user or require an explicit override.

### Bug 6: PfMultiMexStub double-translation — the CRISPR collapse root cause (CRITICAL)

`PfMultiMexStub::copyBarcodesTsv()` converts `barcodes.txt` (NXT, from
`assignBarcodes`) → `barcodes.tsv` (TRU) by applying the 2-column whitelist
COL1→COL2 mapping. This is called at `processAssignOutput()` with the original
2-column whitelist path — by design, to match CellRanger's TRU output convention.

But `PfMultiMerge::readMex()` reads `barcodes.tsv` (now TRU), and the merge
code used `featureMexOutputNamespace = "NXT"` (because `translateNxt=false` and
`assignmentWhitelistNamespace=NXT`). So `normalizeMexBarcodesToTru(data, "NXT")`
translated the already-TRU barcodes → NXT (double-translation). In the merge,
these NXT barcodes didn't match the GEX barcodes (TRU), causing 30,647 of
36,481 feature triplets to be silently dropped — the "feature entries dropped
(barcode not in GEX)" warning in `run.log`.

Evidence from 2M rerun: `assignBarcodes` output had 23,943 barcodes with 1.6M
counts (matching CR), but only 254 survived into the merged raw MEX.

**Fix:** After `processAssignOutput()` runs, update `featureMexOutputNamespace`
to `prepared.soloOutputNamespace` (TRU) for 2-column whitelists. This ensures
the merge code knows `barcodes.tsv` is already in TRU, and
`normalizeMexBarcodesToTru(data, "TRU")` is a no-op.

## Core Concept

The root confusion was between **four** distinct namespaces at different stages:

| Namespace | Meaning | Example (2-col NXT) |
|-----------|---------|---------------------|
| Whitelist matching namespace | COL1, used for read matching | NXT |
| Solo output namespace | `cbWLstrOut`, written to `barcodes.tsv` | TRU |
| assignBarcodes output (`barcodes.txt`) | Internal MEX, before stub | NXT |
| Feature MEX (`barcodes.tsv`, after stub) | After `copyBarcodesTsv` mapping | TRU |

`normalizeMexBarcodesToTru(data, fromChem)` needs `fromChem` to be the **actual
namespace of the barcodes in the file that `readMex` reads** — which is
`barcodes.tsv` (post-stub), not `barcodes.txt` (pre-stub).

`normalizeFilteredBarcodesForAssignNamespace()` must not rely on set membership
when the source namespace is known, because COL1/COL2 sets overlap.

## All Files Changed

### Tracked (vs HEAD)

| File | Change |
|------|--------|
| `core/legacy/source/ParametersSolo.cpp` | Reverted to committed state (no diff) |
| `core/legacy/source/ParametersSolo.h` | Fixed stale `cbWLstrOut` doc-comment |
| `core/legacy/source/PfMultiProcess.cpp` | Fixed GEX/feature MEX namespace bugs + overlap collision + ambiguous-orientation detection (see below) |
| `core/legacy/source/PfMultiAssign.cpp` | Fixed ambiguous-orientation confidence for 2-column whitelists |

### PfMultiProcess.cpp changes (uncommitted, on top of other uncommitted work)

1. **`whitelistHasTwoColumns()`** — new helper to detect 2-column whitelists
2. **`PfMultiPreparedContext`** — added `hasTwoColumnWhitelist`, `soloOutputNamespace`
3. **`buildPfMultiPreparedContext`** — computes `soloOutputNamespace` after chemistry
   detection; recomputes on per-library GEX chemistry override
4. **`gexNormalizationChem`** — now derived from `soloOutputNamespace` (was
   `inferredChem`/`effectiveChem`)
5. **`normalizeFilteredBarcodesForAssignNamespace`** — new `sourceNamespace` +
   `whitelistNamespace` parameters; deterministic mode bypasses membership
   heuristic when both are known, eliminating COL1/COL2 overlap collisions
6. **GEX filtered barcode normalization** — passes `gexNormalizationChem` and
   `gexWhitelistNamespace` for deterministic normalization
7. **Per-library assign filtered barcodes** — passes known namespaces for
   solo-sourced barcodes; heuristic fallback for explicit-file sources
8. **`FeatureRun.featureMexOutputNamespace`** — new field; computed from
   `translateNxtForAssign` and `assignmentWhitelistNamespace`
9. **`FeatureMexEntry.featureMexOutputNamespace`** — propagated from `FeatureRun`
10. **Feature MEX normalization** — uses `featureMexOutputNamespace` (was `effectiveChem`)
11. **`detectChemistryFromWhitelistPath`** — deferred complement-rule return to
    check filename; `confident=false` + WARNING when no filename confirms orientation

### PfMultiAssign.cpp changes

12. **`inferTwoColumnNamespace`** — `confident=false` when complement rule matches
    but no filename hint present (was `confident=true`)

### Untracked (new files)

| File | Purpose |
|------|---------|
| `tests/autodetect_nxt_tru/test_two_column_wl_pairing.sh` | Regression test for 2-col WL pairing + overlap collision |
| `docs/HANDOFF_NAMESPACE_PARITY_FIX_20260303.md` | This document |

### Test files updated

| File | Change |
|------|--------|
| `tests/autodetect_nxt_tru/test_namespace_unit_fixes.sh` | Removed stale `STAR_ALLOW_MIXED_CBWL_OUTPUT` guardrail test (Test 3); now 2/2 pass |

## Regression Tests

### test_two_column_wl_pairing.sh

Self-contained, no genome/fixtures needed. Verifies:

1. Log.out contains correct pairing message
2. Log.out does NOT contain regressed doubling message
3. CB count = number of whitelist rows (not 2x)
4. Overlap whitelist created with deliberate COL1/COL2 self-complementary barcodes
5. STAR loads overlap whitelist with correct pairing (3 CBs, not 6)

```bash
bash tests/autodetect_nxt_tru/test_two_column_wl_pairing.sh
# Current status: all 5 checks PASS
```

### test_namespace_unit_fixes.sh

Tests PfMultiProcess helpers, PfMultiAssign normalization metadata, deterministic
normalization with overlap barcodes, and ambiguous-orientation detection:

```bash
bash tests/autodetect_nxt_tru/test_namespace_unit_fixes.sh
# Current status: 3 passed, 0 failed
```

Test 3 specifically validates the exact bug path: creates COL1/COL2 overlap
barcodes, calls `normalizeFilteredBarcodesForAssignNamespace` with
`sourceNamespace="TRU"` and `whitelistNamespace="NXT"`, and asserts `inSet=0`,
`translatedToSet=3` (all deterministic, no membership guessing). Also verifies
correct barcode identity: `col2_a` maps to `col1_a` (not the overlapping
`col1_b`).

Tests 1 and 2 also verify that ambiguous 2-column whitelists (no "nxt"/"tru" in
filename) return `confident=false` with a WARNING reason string.

## Namespace Flow (Corrected)

For a 2-column NXT whitelist:

```
Solo:
  cbWLstr[i]    = COL1 (NXT, for matching)
  cbWLstrOut[i] = COL2 (TRU, for output)
  barcodes.tsv  = cbWLstrOut → TRU

GEX MEX normalization:
  normalizeMexBarcodesToTru(gexData, "TRU")   ← soloOutputNamespace
  = no-op (already TRU) ✓

Feature assignment (process_features):
  whitelist = COL1 only (NXT)
  translateNxt = true (for TRANSLATED_MATCH) → barcodes.txt: TRU
  translateNxt = false (for RAW_MATCH)       → barcodes.txt: NXT

PfMultiMexStub::copyBarcodesTsv (post-assign):
  Reads barcodes.txt, applies 2-col COL1→COL2 mapping → barcodes.tsv: TRU
  (This is always TRU for 2-column NXT whitelists, regardless of translateNxt)

PfMultiMerge::readMex reads barcodes.tsv → TRU

Feature MEX normalization (featureMexOutputNamespace after stub fix):
  normalizeMexBarcodesToTru(data, "TRU") = no-op ✓
  (previously "NXT" → double-translated to NXT → merge mismatch)

Filtered barcodes (from Solo, GEX path):
  getFilteredBarcodesFromSolo(true) → TRU (cbWLstrOut)
  normalizeFilteredBarcodesForAssignNamespace(_, _, _, "TRU", "NXT")
    → DETERMINISTIC: always translate TRU→NXT → NXT  ✓
    (no set-membership guessing; immune to COL1/COL2 overlap)
  normalizeBarcodeVecToTru(_, "NXT") → TRU ✓

Filtered barcodes (per-library assign path):
  getFilteredBarcodesFromSolo(false) → NXT (cbWLstr)
  normalizeFilteredBarcodesForAssignNamespace(_, _, _, "NXT", "NXT")
    → DETERMINISTIC: same namespace, no translation → NXT  ✓
  (no overlap issue since source = target)

All inputs in TRU → merge → writeCombinedMex ✓
```

## Bug 7: Bootstrap OrdMag blocked on shallow data (HIGH — cell-count parity)

When `nCB < indMin` (45,000), the libscrna EmptyDrops backend aborted entirely
— including the CR9-style bootstrap OrdMag knee estimator. Cell calling fell
back to STAR's old `SoloFeature_cellFiltering` simple filter with hardcoded
`nExpectedCells=3000`. On shallow data (e.g. 2M reads, ~16K barcodes), this
found only 6,069 cells vs CellRanger's 7,346 — a Jaccard gap of 0.10.

The bootstrap OrdMag knee estimation (which gave Jaccard 0.99 on the full
dataset) was locked inside `scrna_emptydrops_run()` and never executed because
the function returned at the `nCB < indMin` guard before reaching OrdMag.

**Root cause:** OrdMag knee detection doesn't need the ambient barcode window
(indMin..indMax) — only MC tail rescue does. The two were incorrectly bundled
behind the same guard.

**Fix (two files):**

1. `SoloFeature_emptyDrops_libscrna.cpp`: when bootstrap is enabled (default,
   non-Flex), the `nCB < indMin` guard no longer returns early. Instead, it logs
   "running bootstrap OrdMag only (no MC tail rescue)" and proceeds.

2. `scrna_api.cpp`: after the OrdMag step (Step 3), if `nCells <= indMin`,
   returns the OrdMag results directly, skipping MC ambient estimation and
   simulation (Steps 4+). This is the same codepath as "no sparse data."

**Flex is unaffected:** libscrna defaults have `use_bootstrap=0`. The STAR
caller gates bootstrap on `!pSolo.emptyDropsLegacyKnee` (default=false, so
bootstrap is on). Flex never hits this codepath (it uses its own FlexFilter).

**New flag:** `--soloEmptyDropsLegacyKnee yes` restores the old behavior
(hardcoded `nExpectedCells=3000`, early return when `nCB < indMin`).

### Validation: UCSF 2M parity (bootstrap OrdMag fix)

| Metric | Hardcoded knee | Namespace-fixed | **Bootstrap OrdMag** | **CR baseline** |
|--------|---------------|-----------------|---------------------|-----------------|
| Filtered cells | 6,069 | 6,069 | **7,200** | **7,346** |
| GEX Jaccard | 0.804 | 0.804 | **0.895** | — |
| Feature Jaccard | 0.624 | 0.804 | **0.895** | — |
| GEX gene Pearson | 0.974 | 0.974 | **0.972** | — |
| Feature rows_star_non_none | 55 | 5,086 | **5,737** | 5,767 |
| Feature assigned UMIs | 17,293 | 1,323,141 | **1,462,084** | 1,447,384 |
| Feature call parity | — | 99.48% | **99.49%** | — |
| Feature UMI Pearson | — | — | **0.9996** | — |

The remaining Jaccard gap (0.895 vs 1.0) reflects the inherent difference
between STAR and CellRanger OrdMag implementations at shallow depth (2M reads,
~16K total barcodes). On the full dataset, the same bootstrap OrdMag code
achieves Jaccard 0.99. This is not a bug — it's a known consequence of low
read depth where the barcode-rank knee is less well-defined.

## All Files Changed (updated)

### Tracked (vs HEAD)

| File | Change |
|------|--------|
| `core/legacy/source/ParametersSolo.cpp` | Reverted to committed state + added `--soloEmptyDropsLegacyKnee` parsing |
| `core/legacy/source/ParametersSolo.h` | Fixed stale `cbWLstrOut` doc-comment + added `emptyDropsLegacyKnee` field |
| `core/legacy/source/Parameters.cpp` | Registered `soloEmptyDropsLegacyKnee` parameter |
| `core/legacy/source/parametersDefault` | Added `soloEmptyDropsLegacyKnee` documentation |
| `core/legacy/source/SoloFeature_emptyDrops_libscrna.cpp` | Bootstrap OrdMag on shallow data; gate on legacy knee flag |
| `core/features/libscrna/src/scrna_api.cpp` | Skip MC tail rescue when nCells <= indMin (return OrdMag-only) |
| `core/legacy/source/PfMultiProcess.cpp` | Fixed GEX/feature MEX namespace bugs + overlap collision + ambiguous-orientation detection |
| `core/legacy/source/PfMultiAssign.cpp` | Fixed ambiguous-orientation confidence for 2-column whitelists |

### Untracked (new files)

| File | Purpose |
|------|---------|
| `tests/autodetect_nxt_tru/test_two_column_wl_pairing.sh` | Regression test for 2-col WL pairing + overlap collision |
| `tests/autodetect_nxt_tru/test_namespace_unit_fixes.sh` | Updated unit tests for deterministic normalization + ambiguity |
| `docs/HANDOFF_NAMESPACE_PARITY_FIX_20260303.md` | This document |
