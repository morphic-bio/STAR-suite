# Union Whitelist Support & Strict Namespace Normalization — Handoff Document

**Date**: 2026-03-03  
**Branch**: CR-Larry-perturb  
**Status**: Implementation complete, tests passing

## Problem

Legacy workflows using the raw `3M-february-2018.txt` barcode whitelist (or
early STARsolo runs before NXT/TRU namespace separation) produce filtered
barcode sets containing a mix of NXT and TRU namespace barcodes. When
`assignBarcodes` receives these mixed-namespace filtered barcode files, it can
miss cells because the filtered barcode hash only contains exact-form entries.

## Solution: Strict Namespace Normalization

This implements the full strict-namespace normalization runbook
(`docs/RUNBOOK_STRICT_BARCODE_NAMESPACE_NORMALIZATION_20260303.md`):

1. **Exact-only matching** — The runtime translated fallback in
   `filtered_barcode_hash_contains()` has been removed. Both the
   `process_features` and `feature_barcodes` implementations now use exact-only
   hash lookup. This eliminates the heuristic translation path and makes
   namespace handling deterministic.

2. **Ingress normalization** — New `--source_namespace` and `--target_namespace`
   flags allow callers to declare the namespace of their filtered barcode file
   and desired output. When both are known single-namespace types (NXT vs TRU)
   and differ, `pf_load_filtered_barcodes()` translates all barcodes at ingress.

3. **Union whitelist escape hatch** — `--allow_union_whitelist` expands the
   filtered barcode hash at ingress time by inserting each barcode's NXT/TRU
   translation. This ensures both namespace forms are present for exact-only
   lookup. Legacy workflows that used mixed-namespace files will work with this
   flag.

4. **Namespace API contract** — `pf_namespace_t` enum (`PF_NS_NXT`, `PF_NS_TRU`,
   `PF_NS_UNION`, `PF_NS_UNKNOWN`), conversion helpers, and
   `pf_translate_barcode_inplace()` are now public API.

5. **Telemetry** — `api_run.txt` now records `namespace_policy.strict_exact_only=1`
   and `namespace_policy.allow_union_whitelist=` for provenance.

6. **Hard error on incomplete namespace** — When filtered barcodes are loaded
   without both `source_namespace` and `target_namespace` resolved (NXT/TRU) and
   without union mode, a hard error is returned (`PF_ERR_NAMESPACE` / `EXIT_FAILURE`).
   This covers partial specification (only one of source/target set) as well as
   completely missing namespace metadata.

### Key Properties

- **Exact-only matching**: No runtime fallback. Hash must contain the exact
  barcode form at lookup time.
- **Ingress-only changes**: Normalization and expansion happen at load time.
  Feature matching, Hamming distance scoring, and assignment logic are unchanged.
- **Allocation-safe**: `expand_hash_union_namespace()` and
  `pf_normalize_hash_namespace()` check every allocation and return `-1` on OOM.
- **Involution**: NXT/TRU translation is its own inverse. Applying it twice
  yields the original barcode.

## Files Changed

### process_features (standalone tool)

| File | Change |
|------|--------|
| `include/pf_api.h` | Added `pf_namespace_t` enum, `pf_translate_barcode_inplace()`, `pf_namespace_to_string/from_string()`, `pf_config_set_source/target_namespace()`, `PF_ERR_NAMESPACE`. |
| `src/pf_api.c` | Namespace API implementation. Added `source_namespace`/`target_namespace` to config. Ingress normalization in `pf_load_filtered_barcodes()`. Returns `PF_ERR_NAMESPACE` hard error when namespace metadata is missing and union mode is off. |
| `src/barcode_match.c` | Removed runtime fallback from `filtered_barcode_hash_contains()`. Added `pf_normalize_hash_namespace()` (with OOM-safe key cleanup on `kh_put` failure). |
| `include/prototypes.h` | Declared `pf_normalize_hash_namespace()`. |
| `src/main.c` | Added `--source_namespace` (required with `--filtered_barcodes`) and `--target_namespace` CLI flags with ingress normalization. Hard error when filtered barcodes loaded with incomplete namespace (either or both unknown) and union mode is off. |
| `tests/test_union_whitelist.c` | Added involution and namespace enum helper tests. Updated `test_exact_only_without_expansion` to assert no-fallback. |
| `tests/test_filtered_barcode_namespace.c` | Updated: `test_translated_miss_exact_only` (asserts miss), `test_translated_found_after_expansion` (asserts hit after expand). |
| `README.md` | Updated namespace flags table, removed fallback note, added exact-only warning. |

### feature_barcodes

| File | Change |
|------|--------|
| `src/assignBarcodes.c` | Removed runtime fallback from `filtered_barcode_hash_contains()`. Added `expand_filtered_hash_union_namespace()` (non-static, GLib equivalent). Per-sample path also expands conditionally (fallback when hash is loaded per-sample rather than preloaded). |
| `include/common.h` | Added `allow_union_whitelist` field to `sample_args` struct. |
| `include/prototypes.h` | Declared `expand_filtered_hash_union_namespace()`. |
| `src/main.c` | Added `--allow_union_whitelist` CLI flag (case 27). Union expansion runs on the preloaded hash in `main()` (the normal path), not only in the per-sample fallback. Help text documents the new flag under "Namespace & Compatibility". |

### STAR integration (core/legacy)

| File | Change |
|------|--------|
| `PfMultiAssign.h` | Added `sourceNamespace`/`targetNamespace` to `AssignOptions`. |
| `PfMultiAssign.cpp` | Wired `pf_config_set_source/target_namespace()` in `applyAssignOptions()`. Added `PF_ERR_NAMESPACE` to `pfErrorCodeString()`. |
| `PfMultiProcess.cpp` | Removed heuristic membership fallback from `normalizeFilteredBarcodesForAssignNamespace()` — now requires resolved namespaces or throws. Hard error when assignment whitelist namespace is unresolved. Added `allowUnionWhitelist` to `PfLibraryNamespaceContext`. Set `sourceNamespace`/`targetNamespace` from assignment whitelist namespace. Added `namespace_policy` telemetry to `api_run.txt`. Union mode with explicit filtered barcodes now bypasses normalization (pass-through) to avoid dropping opposite-namespace barcodes before downstream expansion. |

## Behavioral Impact

| Scenario | Before | After |
|----------|--------|-------|
| Single-namespace filtered barcodes, correct namespace | Works | Works (with namespace flags or union mode) |
| Filtered barcodes with incomplete/missing namespace, no union | Runtime fallback masks mismatch | **Hard error** (`PF_ERR_NAMESPACE` / `EXIT_FAILURE`) — covers both-unknown and partial-spec |
| Mixed-namespace + `--allow_union_whitelist` | N/A | Both forms in hash, all cells found |
| Wrong namespace + `--source/target_namespace` | Runtime fallback masks mismatch | Ingress normalization fixes it |
| STAR explicit filtered barcodes + union mode | Heuristic membership fallback | Pass-through to downstream union expansion (no prefilter normalization drop) |
| STAR explicit filtered barcodes, unknown namespace, no union | Heuristic membership fallback | **Hard error** (requires resolved assignment whitelist namespace) |
| STAR integration path (Solo-sourced) | normalizeFilteredBarcodesForAssignNamespace handles it | Same, deterministic only (heuristic removed) + telemetry |
| feature_barcodes path | Runtime fallback in hash lookup | Expansion conditional on `--allow_union_whitelist` (explicit opt-in) |

## Test Results

All tests pass:
- `test_filtered_barcode_namespace`: 6/6 (exact-only, expansion, miss, null, short)
- `test_union_whitelist`: 8 test groups, all pass (expansion, involution, enum helpers)
- `test_prehash_ambiguity`: PASS
- `test_barcode_filter`: 9/9
- `test_offset_detection`: 0 failures
- `test_pf_api_context_isolation`: PASS
- Regression test: 6/6 files match baseline

## Risk Assessment

**Medium risk** — the removal of both the runtime fallback and the heuristic
membership fallback are behavioral changes that will surface as hard errors for
workflows that previously relied on silent correction.

Mitigations:
1. STAR path already normalizes filtered barcodes via `normalizeFilteredBarcodesForAssignNamespace()` (now deterministic-only) before they reach `assignBarcodes`. No silent data loss for STAR users — unresolved namespace is a hard error instead.
2. Hard error (`PF_ERR_NAMESPACE` / `EXIT_FAILURE`) when filtered barcodes are loaded without namespace metadata and without union mode — fails fast instead of silently dropping.
3. `--allow_union_whitelist` escape hatch for legacy mixed-namespace files.
4. `--source_namespace` / `--target_namespace` for explicit ingress normalization.
5. `feature_barcodes` path now has explicit `--allow_union_whitelist` flag (consistent with `process_features` policy).
6. All existing regression tests pass (namespace-specific tests updated to match exact-only policy).

## STAR CLI Usage

```
--crAssignAllowUnionWhitelist 1
```

## Standalone CLI Usage

```bash
# For mixed-namespace legacy files:
assignBarcodes --allow_union_whitelist \
    --filtered_barcodes /path/to/mixed_barcodes.txt ...

# For known-namespace files that need conversion:
assignBarcodes --source_namespace NXT --target_namespace TRU \
    --filtered_barcodes /path/to/nxt_barcodes.txt ...
```
