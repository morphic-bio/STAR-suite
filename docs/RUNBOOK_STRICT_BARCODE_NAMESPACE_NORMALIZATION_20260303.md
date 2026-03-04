# Runbook: Strict Barcode Namespace Normalization

Date: 2026-03-03 (updated 2026-03-03)  
Branch: `CR-Larry-perturb`  
Status: Phases 0-5 complete; runtime fallback removed, ingress normalization + namespace API + telemetry implemented

## Objective

Prevent further NXT/TRU namespace regressions by enforcing one policy everywhere:

1. Every barcode is normalized at ingress to a declared target namespace.
2. Unknown namespace triggers auto-detection.
3. Ambiguous/low-confidence detection is a hard error.
4. After normalization, all matching/filtering is exact only (no fallback dual-namespace matching).

This runbook is implementation-first and test-heavy.

## Non-Negotiable Policy

1. `UNKNOWN` namespace must not enter core matching/filter paths.
2. Any component reading barcodes must either:
   - receive explicit source namespace, or
   - auto-detect source namespace and require confidence.
3. If source namespace is not confidently determined, fail with actionable error.
4. Normalize exactly once at each boundary; avoid double-translation.
5. Post-normalization operations are exact set/hash lookups only.

### Legacy Exception: Union Whitelist Mode

Some legacy workflows (raw `3M-february-2018.txt` union whitelist, early
STARsolo runs before namespace separation) produce filtered barcode sets
containing mixed NXT+TRU barcodes.  These cannot be assigned a single source
namespace and would fail under rule 3 above.

The `--allow_union_whitelist` flag (standalone CLI) /
`--crAssignAllowUnionWhitelist 1` (STAR parameter) provides an explicit
opt-in escape hatch:

- At filtered barcode ingress, both namespace forms of each barcode are
  inserted into the hash (`expand_hash_union_namespace()`).
- This makes the hash self-contained so Phase 3 fallback removal does not
  regress union-whitelist workflows.
- A `WARNING` is emitted to stderr.  On allocation failure the tool aborts
  rather than silently proceeding with partial expansion.
- This mode is **not** the default and requires explicit opt-in.

Implementation: `HANDOFF_UNION_WHITELIST_SUPPORT_20260303.md`.

## Canonical Internal Conventions

Use explicit conventions to avoid implicit behavior:

1. `assignBarcodes` stage:
   - canonical namespace = assignment whitelist namespace.
   - all filtered barcodes loaded for assign are normalized to assignment
     whitelist namespace before `pf_load_filtered_barcodes`.
   - **Exception — union mode**: when `--allow_union_whitelist` is active,
     filtered barcodes are *not* normalized to a single namespace.  Instead,
     both NXT and TRU forms of each barcode are inserted at ingress via
     `expand_hash_union_namespace()`.  This is the only case where the
     filtered barcode hash intentionally contains two namespace forms of
     the same logical barcode.  Downstream code must not double-normalize
     a union-expanded set (the expansion is idempotent as a guard).
2. `pf-multi` integration boundary:
   - canonical namespace = `TRU`.
   - all MEX barcode vectors converted to `TRU` immediately after read,
     before merge.
   - Union mode does **not** apply to MEX vectors — they are always
     single-namespace-normalized to `TRU` regardless of flag state.
3. final output:
   - convert from canonical `TRU` to `--crOutputChemistry` only in writer.

## Code Touchpoints (Must Be Audited)

### Orchestration / namespace composition

1. `core/legacy/source/PfMultiProcess.cpp`
   - `normalizeFilteredBarcodesForAssignNamespace(...)`
   - `getFilteredBarcodesFromSolo(...)`
   - chemistry composition (`effectiveReadNamespace`, `assignmentWhitelistNamespace`, `translateNxt`)
   - TRU merge-boundary normalization (`normalizeMexBarcodesToTru`, `normalizeBarcodeVecToTru`)
2. `core/legacy/source/PfMultiAssign.cpp`
   - whitelist namespace inference and normalization (`normalizeWhitelistForAssign`)

### Process-features ingress + filtering

1. `core/features/process_features/src/io.c`
   - `read_barcodes_into_hash(...)`
2. `core/features/process_features/src/pf_api.c`
   - `pf_load_filtered_barcodes(...)`
3. `core/features/process_features/src/barcode_match.c`
   - `filtered_barcode_hash_contains(...)` (must become exact-only after ingress normalization)
4. `core/features/process_features/src/assignBarcodes.c`
   - any direct barcode filter membership checks
5. `core/features/process_features/src/mex_writer.c`
   - filtered membership path

### Legacy feature_barcodes parity path

1. `core/features/feature_barcodes/src/assignBarcodes.c`
   - local `filtered_barcode_hash_contains(...)`
   - whitelist/filtered ingestion and filter checks
2. `core/features/feature_barcodes/src/io.c`
   - `read_barcodes_into_hash(...)`

## Implementation Plan

## Phase 0: Freeze and Safety Rails

1. Add this runbook to repo and reference it from active handoff docs.
2. Require clean rebuild before any debugging/benchmark run:
   - `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`
3. Add a short design note in code comments at each normalization boundary:
   - source namespace
   - target namespace
   - confidence requirement

Exit criteria:
1. Team agrees canonical conventions above.
2. All affected owners confirm no heuristic fallback allowed.

## Phase 1: Namespace API Contract

Introduce a shared internal contract (in C++ + C wrappers as needed):

1. Namespace enum:
   - `NXT`, `TRU`, `UNION`, `UNKNOWN`
   - `UNION` represents an explicitly opted-in mixed-namespace set
     (only valid when `allow_union_whitelist` is active).
2. Detection result:
   - `{namespace, confident, reason}`
3. Normalization helper:
   - `normalize_barcode_to_namespace(barcode, from, to)` (idempotent).
4. Involution tests:
   - `translate(translate(x)) == x`.

Exit criteria:
1. All normalization/detection call sites use the same contract shape.
2. No stringly-typed namespace branching in new code paths.

## Phase 2: Enforce Ingress Normalization

For each ingress point (file/memory API), enforce:

1. validate barcode syntax
2. resolve source namespace (explicit or auto-detect)
3. if unresolved/ambiguous → hard error
4. normalize to target namespace
5. store normalized token only

**Filtered-barcode-only exception**: for ingress paths 1–4 below (filtered
barcode sets), if `allow_union_whitelist` is active *and* namespace is
unresolved/ambiguous, classify as `UNION` and expand both forms instead of
erroring.  This exception does **not** apply to ingress path 5 (MEX barcode
vectors) or any whitelist ingress — those always require single-namespace
resolution.

Required ingress paths:

1. filtered barcode files for assign (`pf-multi` path) — union exception eligible
2. filtered barcode files in process_features CLI/API — union exception eligible
3. filtered barcode files in feature_barcodes CLI — union exception eligible
4. Solo in-memory filtered barcodes pulled by `PfMultiProcess` — union exception eligible
5. MEX barcode vectors immediately after `readMex` — **always single-namespace; union mode does not apply**

Exit criteria:
1. `filtered_barcode_hash_contains` can be exact-only in both toolchains.
2. No caller relies on “try translated variant if exact miss”.

## Phase 3: Remove Heuristic Matching

**Prerequisite**: `--allow_union_whitelist` escape hatch is implemented and
tested (done — see `HANDOFF_UNION_WHITELIST_SUPPORT_20260303.md`).  This
ensures legacy union-whitelist workflows survive fallback removal by having
both namespace forms in the hash at ingress.

1. Delete/disable namespace-agnostic fallback behavior in:
   - `process_features/src/barcode_match.c::filtered_barcode_hash_contains`
     (lines 237–252: the translated-query fallback block)
   - `feature_barcodes/src/assignBarcodes.c::filtered_barcode_hash_contains`
2. Replace with exact membership on already-normalized sets.
3. Add fail-fast assertions/log checks if non-canonical namespace reaches lookup.

Exit criteria:
1. Filter behavior is deterministic and exact.
2. No hidden dual-namespace acceptance.
3. `--allow_union_whitelist` workflows pass all regression tests.

## Phase 4: Ambiguity Hard-Error Path

1. For unknown source namespace:
   - run detection
   - if confidence below threshold -> fail
2. Error must include:
   - file/path or source type
   - detection evidence summary
   - exact remediation (set explicit `star_chemistry` / `--crChemistry` / explicit namespace option)

Exit criteria:
1. No silent fallback to guessed namespace.
2. Ambiguous inputs fail early before assignment/merge.

## Phase 5: Telemetry + Provenance

Per library/run output (`assignBarcodes.api_run.txt` + manifest) must record:

1. source namespace
2. target namespace
3. detection mode and confidence
4. number translated / in-set / unmatched
5. strict mode enabled marker

Exit criteria:
1. A reviewer can reconstruct namespace decisions from artifacts alone.

## Phase 6: Regression Suite (Extensive)

### A. Unit tests (fast, deterministic)

1. Translation involution and edge lengths.
2. Namespace detection:
   - explicit NXT/TRU
   - auto-detect confident
   - ambiguous -> error
3. Normalization idempotence:
   - normalize twice yields same output.
4. Exact filter policy:
   - exact hit passes
   - translated-only hit fails unless normalized at ingress beforehand.

### B. Component tests (process_features / feature_barcodes)

1. `read_barcodes_into_hash` stores canonical namespace only.
2. `pf_load_filtered_barcodes` rejects ambiguous source namespace.
3. `filtered_barcode_hash_contains` exact-only semantics.
4. MEX writer/filter behavior unchanged after explicit pre-normalization.

### C. Integration tests (pf-multi)

1. 1-column TRU whitelist, TRU reads.
2. 1-column NXT whitelist, NXT reads.
3. 2-column whitelist, mixed libraries:
   - GEX TRU, feature NXT
   - GEX TRU, feature TRU
4. Explicit filtered barcodes in opposite namespace are normalized correctly.
5. Unknown/ambiguous namespace input fails before assign.
6. STAR flag wiring (`--crAssignAllowUnionWhitelist`):
   - With mixed-namespace filtered barcodes and flag set to `0`:
     verify behavior matches default (no expansion, post-Phase-3 would
     fail on exact-only mismatch).
   - With same input and flag set to `1`: verify expansion occurs,
     both namespace forms are in the hash, and assignment succeeds.
   - This locks the Parameters → PfMultiProcess → PfMultiAssign →
     `applyAssignOptions` → `pf_config_set_allow_union_whitelist` →
     `pf_load_filtered_barcodes` wiring end-to-end.

### D. Real-data regressions

1. UCSF 2M known parity fixture:
   - verify no barcode collapse
   - filtered/raw overlap and Pearson/Spearman remain within prior acceptance bands.
2. MSK 30KO 2-feature benchmark fixtures:
   - ensure no namespace-induced drops in LARRY or PolyIII.
3. A375 parity smoke (existing baseline path).

### E. Negative tests (must fail)

1. Ambiguous namespace source without explicit override.
2. Mixed-namespace filtered file without resolvable source metadata
   (unless `--allow_union_whitelist` is set — then it must succeed).
3. Contradictory explicit namespace vs observed evidence.

### F. Union whitelist escape hatch (implemented)

1. `test_union_whitelist.c` — 6 test cases, 22 assertions:
   - Single-barcode expansion adds correct NXT/TRU translation.
   - Both forms present → no-op (0 additions).
   - Short barcodes (< 9 bases) skipped.
   - Multiple mixed barcodes expand correctly.
   - Runtime fallback works without expansion (pre-Phase-3 compat).
   - Idempotent double expansion.
2. Allocation failure returns `-1`; API path returns `PF_ERR_ALLOC`;
   CLI path exits with `EXIT_FAILURE`.
3. Fault-injection method for OOM validation:
   - Preferred: add a `test_union_whitelist_oom.c` harness that uses
     `LD_PRELOAD` with a `malloc`-fail shim (e.g., `libfailmalloc` or a
     minimal wrapper that fails after N allocations) to force `strdup` or
     `kh_put` failure inside `expand_hash_union_namespace`.
   - Alternative (no LD_PRELOAD): add a compile-time `#ifdef PF_TEST_OOM`
     guard in `expand_hash_union_namespace` that injects a failure after
     the first successful `strdup`, then assert the function returns `-1`
     and the hash is not corrupted (existing entries still valid).
   - Either approach must verify: return value is `-1`, no null keys in
     hash, no leaked pending keys, and caller propagates the error.

## Acceptance Gates

All must pass:

1. No heuristic dual-namespace matching in filter lookups.
2. Ambiguous namespace always errors before compute-heavy stages.
3. MEX merge boundary remains canonical TRU internally.
4. Existing parity metrics do not regress versus known baselines.
5. New tests are wired into `make` targets (not ad-hoc only).

## Execution Checklist

1. Implement Phases 1-3 in small commits (one concern per commit).
2. Run fast unit/component suite on each commit.
3. Run integration suite after Phase 4.
4. Run real-data regression pack after Phase 6.
5. Publish one summary handoff with:
   - commands
   - artifact paths
   - pass/fail table per gate
   - any remaining risk.

## Build + Test Commands

Core rebuild:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
make feature-barcodes-tools
```

Process-features focused tests:

```bash
make -C core/features/process_features -j8 test
```

Targeted namespace tests (expected to expand in this effort):

```bash
core/features/process_features/tests/test_filtered_barcode_namespace
core/features/process_features/tests/test_pf_api_context_isolation
core/features/process_features/tests/test_union_whitelist
```

## Risk Controls

1. Keep legacy flag aliases for user-facing compatibility, but do not allow legacy semantics that violate strict normalization policy.
2. For any temporary compatibility mode, require explicit opt-in and loud warning in logs.
   - `--allow_union_whitelist` / `--crAssignAllowUnionWhitelist` is the only
     current exception; it requires explicit opt-in, emits a WARNING, and
     aborts on allocation failure rather than silently degrading.
3. Do not merge partial namespace work without full regression matrix results.
4. Phase 3 (runtime fallback removal) must not land without confirming that
   all known union-whitelist workflows pass with `--allow_union_whitelist`.

## Deliverables

1. Code changes implementing strict normalization-at-ingress policy.
2. Union whitelist escape hatch for legacy workflows (done —
   `HANDOFF_UNION_WHITELIST_SUPPORT_20260303.md`).
3. Expanded regression tests covering unit/component/integration/E2E.
4. Final validation report with artifact links and gate outcomes.
