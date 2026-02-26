# Handoff: Multi-Feature Libraries — Phases 0–2

Date: 2026-02-26  
Branch: `multi-feature`  
Base: `master` at `8c402fc` (includes NXT/TRU auto-detection and `star_chemistry` column)  
Status: Phases 0–2 complete; two review rounds applied; E2E validated  
Runbook: `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md`

## What Was Done

### Problem

STAR's `pf-multi` pipeline assumed a single global feature reference for all
feature libraries, and only recognized a hardcoded list of 10x feature types
(CRISPR Guide Capture, Antibody Capture, Multiplexing Capture). This prevented
running experiments with multiple feature library types in one run — for
example, MSK 30polyKO experiments that combine:

- **mRNA** (Gene Expression, TRU)
- **PolyIII** (gRNA / CRISPR Guide Capture, NXT)
- **LARRY** (lineage barcodes / Custom, TRU)

The LARRY library used a `feature_types=Custom` value that the hardcoded
list silently ignored, and there was no mechanism to specify independent
feature reference CSVs per library.

### Changes

**Phase 0: Data-driven feature specs and fixtures**

1. Replaced `defaultPfFeatureSpecs()` (hardcoded 4-entry list) with
   `buildFeatureSpecsFromConfig()` that scans all non-GEX `feature_types`
   values from the config. Any `feature_types` not matching "Gene Expression"
   is now automatically routed to `assignBarcodes`.
2. Known 10x types map to their canonical `featureRefType` via
   `knownFeatureRefTypeMap()`; unknown types (Custom, Lineage, etc.) use
   their `feature_types` value verbatim and emit a `NOTICE` in the prep log.
3. Created fixture tooling: `tests/multi_feature/create_fixture.sh` —
   downsamples MSK 30polyKO data (mRNA + PolyIII + LARRY, 100K reads each).
4. Verified no regression: 8/8 autodetect + 15/15 star_chemistry tests pass.

**Phase 1: Config parsing and feature ref precedence**

1. Extended `PfMultiConfig::LibraryEntry` with two new fields:
   - `starFeatureRef` — per-library feature reference CSV path
   - `starLibraryId` — stable output/provenance identifier
2. Parser handles columns `star_feature_ref` / `starfeatureref` and
   `star_library_id` / `starlibraryid` (case-insensitive, alias support).
3. Feature ref precedence: `star_feature_ref` (per-library) > `--crFeatureRef`
   (CLI flag) > `[feature] ref` (config global).
4. When `star_feature_ref` is set, `filterFeatureRefCsv` is skipped for that
   library (the per-library ref is assumed to already be specific).
5. Auto-generates `star_library_id` when absent:
   `{sample}_{feature_types}_{index}` (spaces replaced with `_`).
6. Post-parse validation:
   - Duplicate `star_library_id` → hard `runtime_error`
   - Nonexistent `star_feature_ref` path → hard `runtime_error`
7. `PfPreparedFeatureLibrary` now carries `libraryId` and logs it in prep output.

**Phase 2: Per-library assign wiring and E2E validation**

1. The per-library assign wiring was already in place from Phase 1: each
   `PfPreparedFeatureLibrary` carries its own `resolvedFastq`, `assignOut`
   (using `libraryId`), `featureRefPath`, `effectiveChem`, `explicitChem`,
   and `resolvedChemRequest`. The assign loop (line 1317) calls
   `runAssignBarcodes` with per-library arguments.
2. Chemistry composition is per-library: `useAutodetect` and `explicitChem`
   flags are computed per `preparedLib`. Auto-detect results are anchored
   per library independently.
3. **E2E smoke test** (`run_msk_multifeature_smoke.sh`):
   - Config: mRNA (GEX/TRU) + PolyIII (CRISPR/NXT) + LARRY (Custom/TRU)
   - No global `[feature] ref` section (all feature libs provide `star_feature_ref`)
   - Explicit `star_library_id` for output provenance
   - **9/9 validation checks pass**: output dirs exist, chemistry logged per lib,
     per-library feature ref used, CRISPR recognized as known type, Custom flagged
     as unknown, no global ref error
4. **Chemistry validation** from `Log.out`:
   - GEX: `star_chemistry=tru → effectiveChem=TRU` (re-anchors namespace)
   - CRISPR: `star_chemistry=nxt → effectiveChem=NXT` (auto-detect skipped)
   - LARRY: `star_chemistry=tru → effectiveChem=TRU` (auto-detect skipped)
   - No cross-library chemistry leakage
5. **Regression test**: A375 single-library CRISPR test passes (1,186 cells,
   1,051 with single feature call, full GMM pipeline working). The auto-generated
   `libraryId` (`A375_CRISPR_Guide_Capture_1`) is backward-compatible.
6. Fixed fixture naming: `create_fixture.sh` now produces `_R1_001.fastq.gz`
   (not `_R1.fastq.gz`) to match the `_R1_` pattern expected by `assignBarcodes`.

## Files Changed

| File | Change |
|------|--------|
| `core/legacy/source/PfMultiConfig.h` | Added `starFeatureRef`, `starLibraryId` to `LibraryEntry` |
| `core/legacy/source/PfMultiConfig.cpp` | Parse `star_feature_ref`/`star_library_id` in both parsing blocks; post-parse auto-ID generation + duplicate/path validation; relative-path resolution against config dir; exact matching in `getFeatureLibraries`; added `<set>` include |
| `core/legacy/source/PfMultiProcess.cpp` | Extracted `buildFeatureSpecsFromConfig()` into `PfMultiFeatureSpecs.h`; relaxed global featureRef requirement (only if a library needs it); `assignOut` uses `libraryId` instead of `sampleName`; added feature ref precedence logic (per-lib > CLI > config global); enriched prep log |
| `core/legacy/source/PfMultiFeatureSpecs.h` | Extracted `buildFeatureSpecsFromConfig()` + `knownFeatureRefTypeMap()` into a testable header |
| `docs/feature_barcodes.md` | New sections: "Per-Library Feature Reference and Library ID", "Data-Driven Feature Routing"; corrected precedence order |
| `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md` | Post-review updates: dropped `star_whitelist` (deferred), updated phases, resolved open questions, corrected precedence |
| `tests/multi_feature/create_fixture.sh` | MSK 30polyKO fixture creation script; proper SIGPIPE handling + non-empty output validation |
| `tests/multi_feature/multi_config_msk.csv.template` | Template config with `star_feature_ref` and `star_library_id` for the MSK fixture |
| `tests/multi_feature/test_multi_feature_config.sh` | 20-case unit test for config parsing (incl. relative path resolution, no-global-ref) |
| `tests/multi_feature/test_data_driven_specs.sh` | 16-case test exercising `buildFeatureSpecsFromConfig` + `getFeatureLibraries` (exact matching, CellPlex, dedup) |
| `tests/multi_feature/run_msk_multifeature_smoke.sh` | E2E smoke test: mRNA + PolyIII + LARRY with per-library feature ref, chemistry, and library ID (9 checks) |

## Architecture

```
pfMultiConfig CSV
 ├─ [libraries] header: ..., star_feature_ref, star_library_id
 │
 └─ parseConfig()
      ├─ parse star_feature_ref / star_library_id per row
      ├─ auto-generate star_library_id when absent
      ├─ validate: no duplicate raw IDs, no sanitized-ID collisions,
      │            star_feature_ref paths exist
      │
      └─ PfMultiProcess::buildPfMultiPreparedContext()
           │
           ├─ buildFeatureSpecsFromConfig(config)
           │    ├─ scan all non-GEX feature_types
           │    ├─ dedup by normalized type
           │    └─ map known→canonical, unknown→verbatim + NOTICE
           │
           ├─ for each spec × matching library:
           │    ├─ if star_feature_ref set:
           │    │    use directly, skip filterFeatureRefCsv
           │    ├─ else:
           │    │    filter global ref by featureRefType
           │    ├─ resolve chemistry (star_chemistry > --crChemistry > auto)
           │    ├─ assignOut = cr_assign/{type}/{libraryId}
           │    └─ log: library_id, feature_types, featureRef, effectiveChem
           │
           └─ runPfMulti()  [Phase 2]
                │
                └─ for each preparedLib:
                     ├─ resolve whitelist (global)
                     ├─ set per-library chemistry (useAutodetect / explicitChem)
                     ├─ runAssignBarcodes(whitelist, refPath, fastq, assignOut)
                     ├─ post-assign: compose effectiveChem from autodetect result
                     │    (RAW_MATCH → same namespace, TRANSLATED → opposite)
                     └─ record FeatureRun for merge stage
```

## Test Results

| Suite | Result |
|-------|--------|
| `tests/autodetect_nxt_tru/test_autodetect.sh` | 8/8 pass (regression) |
| `tests/autodetect_nxt_tru/test_star_chemistry_column.sh` | 15/15 pass (regression) |
| `tests/multi_feature/test_multi_feature_config.sh` | 21/21 pass (expanded after reviews) |
| `tests/multi_feature/test_data_driven_specs.sh` | 16/16 pass (rewritten after review 1) |
| `tests/multi_feature/test_no_global_ref_guard.sh` | 4/4 pass (new, review 2) |
| `tests/multi_feature/run_msk_multifeature_smoke.sh` | 9/9 pass (E2E, Phase 2) |
| `tests/test_cr_compat_crispr_calling.sh` | pass (A375 regression, Phase 2) |
| **Total** | **73/73 pass + 1 regression pass** |

### New Test Coverage (Phase 2)

**E2E smoke test** (`run_msk_multifeature_smoke.sh`):
- Output directories exist for each library (grna_de, larry_de) (2 assertions)
- Per-library chemistry logged in Log.out (2 assertions)
- Per-library feature ref logged in Log.out (2 assertions)
- CRISPR recognized as known type / Custom flagged as unknown (2 assertions)
- No global feature ref error when all feature libs provide star_feature_ref (1 assertion)

**Regression test** (`test_cr_compat_crispr_calling.sh`):
- A375 single-library CRISPR test: 1,186 cells, full GMM pipeline, protospacer calls

### New Test Coverage (Phase 1)

**Config parser tests** (`test_multi_feature_config.sh`):
- Per-library `star_feature_ref` parsing (3 assertions)
- Explicit `star_library_id` preservation (2 assertions)
- Auto-generated `star_library_id` format (2 assertions)
- Duplicate `star_library_id` rejection (1 assertion)
- Nonexistent `star_feature_ref` rejection (1 assertion)
- Backward compatibility without new columns (2 assertions)
- Column aliases `starfeatureref` / `starlibraryid` (2 assertions)
- Trailing field padding with new columns (2 assertions)
- Full multi-feature config with all columns (3 assertions)

**Data-driven routing tests** (`test_data_driven_specs.sh`):
- Known types (CRISPR + Antibody) route correctly (4 assertions)
- Custom type routed as non-GEX feature library (2 assertions)
- Multiple Custom libraries with different refs (2 assertions)
- Non-standard "Lineage" type routed as non-GEX (1 assertion)

## Design Decisions

1. **`star_whitelist` deferred**: All libraries in one GEM well share cell
   barcodes. NXT/TRU barcode variants are already handled by `star_chemistry`
   plus merge normalization. Per-library CB whitelists would require
   barcode-set union semantics at merge, which is out of scope for this slice.

2. **Data-driven over hardcoded**: Rather than extending the hardcoded list
   each time a new feature type appears, the config itself declares what
   feature types exist. Known types still map to canonical `featureRefType`
   for backward compatibility. The `NOTICE` log for unknown types provides
   visibility without blocking.

3. **Filter skip for per-library refs**: When a library provides its own
   `star_feature_ref`, type-based filtering is skipped. The assumption is
   that a per-library ref is already specific to the library's feature type.
   This avoids misfiltering when the ref's `feature_type` column doesn't
   match the 10x canonical name.

4. **Auto-generated library IDs**: Format `{sample}_{feature_types}_{index}`
   is deterministic and collision-free within a config. Spaces in
   `feature_types` are replaced with `_`. Explicit `star_library_id` overrides
   this for stable provenance across config edits.

## Review Fixes

| Severity | Finding | Fix |
|----------|---------|-----|
| **High** | Global `context.featureRef` hard-fails at line 910 before per-library refs are applied, blocking configs where all feature libs provide `star_feature_ref` | Deferred check: only require global ref if at least one non-GEX library lacks `star_feature_ref` |
| **High** | `getFeatureLibraries` uses `find()` substring match (Custom matches Custom2) and only strips spaces from target, not punctuation (CellPlex (CMO) can't match itself) | Replaced with exact `==` comparison on fully normalized strings (strip all non-alphanumeric, lowercase), same normalization for both target and library |
| **High** | `assignOut` uses only `sampleName`, so parallel libraries with same sample overwrite each other | Changed `assignOut` to use `libraryId` (unique by construction) |
| **Medium** | `test_data_driven_specs.sh` only linked `PfMultiConfig.o`, never exercised `buildFeatureSpecsFromConfig` | Extracted `buildFeatureSpecsFromConfig` + `knownFeatureRefTypeMap` into testable `PfMultiFeatureSpecs.h`; rewrote test to link and exercise both functions; added tests for exact matching, CellPlex punctuation, dedup |
| **Medium** | Docs claim precedence `star_feature_ref > [feature] ref > --crFeatureRef` but code does `--crFeatureRef` first | Corrected docs/handoff/runbook to `star_feature_ref > --crFeatureRef > [feature] ref` |
| **Medium** | `star_feature_ref` validation uses `stat()` at parse time with no CWD-relative resolution | Relative paths now resolved against the config file's directory before `stat()` |
| **Low** | Fixture script `|| true` masks extraction failures | Replaced with `set +o pipefail` around the pipe (SIGPIPE is expected) + explicit non-empty output check |

**Additional fix found during review 1**: `knownFeatureRefTypeMap` had keys with spaces
(`"crisprguide capture"`) that could never match the space-stripped normalized output
of `normalizedFeatureType()`. Fixed keys to match normalization format (`"crisprguidecapture"`).

### Review 2

| Severity | Finding | Fix |
|----------|---------|-----|
| **High** | Sanitized library ID collision: distinct raw IDs like `a:b` and `a/b` both sanitize to `a_b` via `sanitizeDirName`, causing output directory overwrites | Added second uniqueness check on sanitized IDs at parse time; replicates `sanitizeDirName` logic (non-alnum/dash/underscore → `_`) and throws with a clear message naming both the offending ID and the collision result |
| **Medium** | "No global ref needed" path only parser-tested, not pipeline-tested | Added `test_no_global_ref_guard.sh` (4 cases) that reproduces the exact guard logic from `buildPfMultiPreparedContext` and exercises all branches: all-covered, one-missing, global-set, GEX-only |
| **Low** | Handoff verify block listed stale test counts (18/18, 9/9) vs actual (20/20, 16/16) | Updated to current counts (21/21, 16/16, 4/4) |

## What Remains (Phases 3–5)

| Phase | Description | Status |
|-------|-------------|--------|
| **Phase 2** | Per-library assign wiring + E2E validation | **Complete** |
| **Phase 3** | Merge and provenance hardening: merge all library MEX outputs; per-library provenance manifests | Not started |
| **Phase 4** | PF intra-domain scheduler: library-aware permit allocation | Not started |
| **Phase 5** | Multi-library validation suite + E2E parity tests + performance benchmarks | Not started |

## How to Verify

```bash
# Build
make -C core/legacy/source -j8 STAR

# Run unit tests (fast, no fixtures needed)
bash tests/autodetect_nxt_tru/test_autodetect.sh          # 8/8
bash tests/autodetect_nxt_tru/test_star_chemistry_column.sh # 15/15
bash tests/multi_feature/test_multi_feature_config.sh       # 21/21
bash tests/multi_feature/test_data_driven_specs.sh          # 16/16
bash tests/multi_feature/test_no_global_ref_guard.sh        # 4/4

# Create fixture (requires MSK data at /mnt/pikachu/MSK-perturb/)
bash tests/multi_feature/create_fixture.sh /tmp/msk_multi_fixture

# E2E smoke test: multi-feature (requires fixture + genome index; ~25 min)
bash tests/multi_feature/run_msk_multifeature_smoke.sh      # 9/9

# Regression: A375 single-library CRISPR (requires /storage/A375 fixtures; ~6 min)
bash tests/test_cr_compat_crispr_calling.sh
```

## Review Focus Areas (Post-Fix)

1. **GEX exclusion edge cases** (`PfMultiFeatureSpecs.h`): The GEX check uses
   `find("geneexpression")`. A type like `"Gene Expression Plus"` would
   normalize to `"geneexpressionplus"` which still contains `"geneexpression"`,
   so it would be excluded from feature specs. Is this the correct behavior?

2. **Filter-skip safety**: When `star_feature_ref` is set, the per-library ref
   is used directly without type-based filtering. Should a warning be emitted
   if the ref CSV contains features of multiple types?

3. **Auto-generated `star_library_id` stability**: The format uses row index
   (`_0`, `_1`), making it unstable if rows are reordered. Is this acceptable
   for auto-generated IDs, or should the provenance docs recommend explicit
   `star_library_id` for reproducibility?

4. **`assignOut` path change**: Switching from `sampleName` to `libraryId`
   changes the output directory structure. Existing single-library configs
   get auto-generated IDs like `S1_CRISPR_Guide_Capture_1` instead of `S1`.
   Is this backward-compatible for downstream tools that read `cr_assign/` paths?
