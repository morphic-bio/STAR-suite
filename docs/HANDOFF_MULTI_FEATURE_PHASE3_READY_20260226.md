# Handoff: Multi-Feature Libraries — Ready for Phase 3

Date: 2026-02-26
Author: AI agent (Phases 0–2)
Branch: `core-alignment-threads-integration` (working tree, uncommitted)
Base: `master` at `8c402fc`
Runbook: `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md`
Prior handoff: `docs/HANDOFF_MULTI_FEATURE_LIBRARIES_PHASE012_20260226.md`
Prior conversation: [Multi-feature phases 0-2](04048633-e62b-447e-8b6e-215ee73b03b4)

## Current State

Phases 0–2 are **complete and validated**. The next agent should begin
**Phase 3: Merge and Provenance Hardening**.

| Phase | Status | Tests |
|-------|--------|-------|
| Phase 0: Data-driven feature specs + fixtures | Complete | 8+15 regression pass |
| Phase 1: Config parsing + feature ref precedence | Complete, 2 review rounds | 21+16+4 = 41 pass |
| Phase 2: Per-library assign wiring + E2E | Complete | 9/9 E2E + A375 regression |
| **Phase 3: Merge and provenance hardening** | **Not started** | — |
| Phase 4: PF intra-domain scheduler | Not started | — |
| Phase 5: Validation suite + parity | Not started | — |

## What Was Built (Phases 0–2)

### Core capability

STAR's `pf-multi` pipeline now supports multiple feature library types in a
single run. Each library can independently specify:

- **`star_feature_ref`** — per-library feature reference CSV (skips global filtering)
- **`star_library_id`** — stable output/provenance identifier (auto-generated if absent)
- **`star_chemistry`** — per-library chemistry override (TRU, NXT, auto)

Feature routing is **data-driven**: any `feature_types` value in the config
that doesn't normalize to "geneexpression" or "gex" is automatically routed
to `assignBarcodes`. Known 10x types map to canonical `featureRefType`;
unknown types (Custom, Lineage) use their value verbatim with a NOTICE log.

### Key files modified

| File | Role |
|------|------|
| `core/legacy/source/PfMultiConfig.h` | `LibraryEntry` struct with `starFeatureRef`, `starLibraryId` |
| `core/legacy/source/PfMultiConfig.cpp` | Parser for new columns + validation (duplicate ID, sanitized collision, path existence, relative path resolution) |
| `core/legacy/source/PfMultiProcess.cpp` | Feature ref precedence, conditional global ref requirement, `assignOut` uses `libraryId`, per-library chemistry composition |
| `core/legacy/source/PfMultiFeatureSpecs.h` | Extracted `buildFeatureSpecsFromConfig()` + `knownFeatureRefTypeMap()` |
| `docs/feature_barcodes.md` | User-facing docs for new columns and data-driven routing |
| `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md` | Master runbook for all phases |

### Key test files

| File | What it tests | Count |
|------|--------------|-------|
| `tests/multi_feature/test_multi_feature_config.sh` | Config parsing: new columns, aliases, validation, auto-ID | 21 |
| `tests/multi_feature/test_data_driven_specs.sh` | Feature routing: known types, custom, exact matching, dedup | 16 |
| `tests/multi_feature/test_no_global_ref_guard.sh` | Conditional global ref requirement (pipeline-level) | 4 |
| `tests/multi_feature/run_msk_multifeature_smoke.sh` | E2E: mRNA+PolyIII+LARRY with per-lib ref/chem/ID | 9 |
| `tests/multi_feature/create_fixture.sh` | Fixture creation (100K reads from MSK 30polyKO) | — |
| `tests/test_cr_compat_crispr_calling.sh` | Regression: A375 single-library CRISPR calling | — |

## How to Build and Run Tests

```bash
# Build STAR
make -C core/legacy/source -j8 STAR

# Fast unit tests (no fixtures needed, <30s total)
bash tests/autodetect_nxt_tru/test_autodetect.sh            # 8/8
bash tests/autodetect_nxt_tru/test_star_chemistry_column.sh  # 15/15
bash tests/multi_feature/test_multi_feature_config.sh        # 21/21
bash tests/multi_feature/test_data_driven_specs.sh           # 16/16
bash tests/multi_feature/test_no_global_ref_guard.sh         # 4/4

# Create fixture (requires MSK data at /mnt/pikachu/MSK-perturb/)
bash tests/multi_feature/create_fixture.sh /tmp/msk_multi_fixture

# E2E smoke test (~25 min, LARRY library with 245K features is slow)
bash tests/multi_feature/run_msk_multifeature_smoke.sh       # 9/9

# Single-library regression (~6 min, requires /storage/A375)
bash tests/test_cr_compat_crispr_calling.sh
```

## Phase 2 E2E Results (Reference)

The MSK multi-feature smoke test ran with this config:

| Library | feature_types | star_chemistry | star_feature_ref | star_library_id |
|---------|--------------|----------------|------------------|-----------------|
| mRNA | Gene Expression | TRU | (none) | gex_de |
| PolyIII | CRISPR Guide Capture | NXT | ref_feature_geneBC.csv | grna_de |
| LARRY | Custom | TRU | ref_feature_larryBC.csv | larry_de |

Chemistry log from `Log.out`:
```
GEX star_chemistry=tru → effectiveChem=TRU, inferredChem=TRU (re-anchored)
CRISPR star_chemistry=nxt → effectiveChem=NXT (auto-detect skipped)
Custom star_chemistry=tru → effectiveChem=TRU (auto-detect skipped)
```

Output directory structure:
```
/tmp/msk_smoke_p2/
├── cr_assign/
│   ├── CRISPR_Guide_Capture/grna_de/PolyIII/   ← CRISPR assign output
│   └── Custom/larry_de/LARRY/                    ← LARRY assign output
├── outs/
│   ├── raw_feature_bc_matrix/                    ← Merged MEX
│   └── filtered_feature_bc_matrix/               ← Merged MEX (filtered)
├── Solo.out/                                     ← GEX Solo output
└── Log.out
```

Merged MEX feature breakdown:
- Gene Expression: 38,606 features
- Custom: 246,009 features (LARRY barcodes)
- **CRISPR Guide Capture: 0 features** (see Known Issue #1 below)

## Known Issues for Phase 3

### Issue 1 (High): CRISPR features missing from merged MEX

The CRISPR Guide Capture library's assign output (`grna_de/PolyIII/features.tsv`)
labels all features as `Custom` because that's the `feature_type` column value
in the upstream feature reference CSV (`ref_feature_geneBC.csv`). The merge
*does* include these features, but they appear as "Custom" type rather than
"CRISPR Guide Capture". As a result:

- `runCrisprFeatureCalling()` finds zero CRISPR features and skips calling
- The features are merged into the output but mislabeled

**Root cause**: `PfMultiMexStub::processAssignOutput()` preserves the
`feature_type` from the feature reference CSV verbatim. When the ref uses
"Custom" instead of "CRISPR Guide Capture", the label is wrong.

**Suggested fix**: In the merge stage, either:
1. Override the feature type in the MEX stub to match the library's declared
   `feature_types` from the config (preferred — the config is the source of truth)
2. Or: add a per-library provenance annotation that maps assign-output feature
   types to declared library types

This is the most important issue for Phase 3 to resolve.

### Issue 2 (Medium): No per-library provenance manifests

The runbook (Phase 3, item 2) calls for per-library provenance manifests under
each assign output directory. Currently, the only provenance is the prep log
in `Log.out`. Each library should get a machine-readable manifest with:
- resolved FASTQ path
- resolved feature ref path
- resolved whitelist path
- detected match mode
- effective chemistry
- library ID

### Issue 3 (Low): Merged MEX has no CRISPR calling for multi-library

The E2E smoke test skipped CRISPR calling because the features were
mislabeled (Issue #1). Once feature types are correctly mapped, CRISPR
calling should work for multi-library configs. However, verify that the
`runCrisprFeatureCalling` function handles merged MEX with multiple non-GEX
feature types gracefully (it should only use CRISPR features).

### Issue 4 (Low): Double-directory nesting in assign output

The assign output path is `cr_assign/{type}/{libraryId}/{sampleName}/` —
the `sampleName` subdirectory is created by `assignBarcodes` when given a
directory as input. For the MSK fixture, this produces paths like
`cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII/`. The `findAssignOutputDir`
function handles this by searching one level of subdirectories for `matrix.mtx`,
but it's fragile. Consider whether the sample subdirectory is needed when
`libraryId` already provides uniqueness.

## Phase 3 Scope (from Runbook)

> 1. Merge all successful library MEX outputs after normalization to output namespace.
> 2. Add per-library provenance manifests and top-level summary table in `Log.out`.
> 3. Guard against silent partial merge:
>    - if a configured library fails assign, record explicit failure and apply
>      existing run policy (fail-fast or controlled skip).

### Concrete tasks

1. **Fix feature type labeling in merge** (Issue #1): When building MEX stubs
   from assign output, use the library's `featureRefType` (from config routing)
   instead of the raw `feature_type` from the feature ref CSV. The mapping is
   available in `FeatureRun.featureType`.

2. **Add per-library provenance manifest**: Write a JSON or TSV file under
   each library's assign output directory with resolved inputs, detected
   match mode, and effective chemistry.

3. **Add top-level summary table**: Log a structured summary in `Log.out`
   after all libraries complete, showing library_id, feature_type, reads
   processed, barcodes found, effective chemistry, and status.

4. **Guard against partial merge**: If a library fails `runAssignBarcodes`,
   currently it throws and aborts. Phase 3 should (at minimum) log which
   library failed and optionally support controlled-skip mode (deferred per
   runbook).

5. **Verify CRISPR calling works with corrected labels**: After fixing Issue #1,
   re-run the E2E smoke and confirm `crispr_analysis/` is generated when
   CRISPR Guide Capture features are present.

### Exit criteria (from Runbook)

1. Output includes all configured feature types in merged matrices with correct labels.
2. Provenance shows exact resolved inputs per library.
3. A375 single-library regression still passes.

## Architecture Reference

```
pfMultiConfig CSV
 ├─ [libraries]: ..., star_feature_ref, star_library_id
 │
 └─ parseConfig()
      ├─ parse + validate (duplicate IDs, sanitized collisions, paths)
      │
      └─ buildPfMultiPreparedContext()
           ├─ buildFeatureSpecsFromConfig(config)     [PfMultiFeatureSpecs.h]
           │    ├─ scan non-GEX feature_types, dedup
           │    └─ map known→canonical, unknown→verbatim
           │
           ├─ for each spec × matching library:
           │    ├─ star_feature_ref set? → use directly
           │    ├─ else: filter global ref by featureRefType
           │    ├─ resolve chemistry per library
           │    └─ assignOut = cr_assign/{type}/{libraryId}
           │
           └─ runPfMulti()
                ├─ for each preparedLib:
                │    ├─ runAssignBarcodes(whitelist, refPath, fastq, assignOut)
                │    └─ compose effectiveChem from autodetect result
                │
                ├─ for each featureRun:                  ← PHASE 3 FIX HERE
                │    ├─ findAssignOutputDir()
                │    └─ processAssignOutput()             (feature type label issue)
                │
                ├─ normalize all MEX to TRU namespace
                ├─ mergeMex(gexData, featureDataVec)
                ├─ writeCombinedMex(raw + filtered)
                └─ runCrisprFeatureCalling() if CRISPR present
```

## Data Locations

| Item | Path |
|------|------|
| STAR binary | `core/legacy/source/STAR` |
| Genome index | `/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star` |
| MSK fixture | `/tmp/msk_multi_fixture` (recreate with `create_fixture.sh`) |
| MSK source FASTQs | `/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM` |
| gRNA feature ref | `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv` |
| LARRY feature ref | `/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv` |
| TRU whitelist | `/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt` |
| A375 fixtures | `/storage/A375/` |
| Phase 2 E2E output | `/tmp/msk_smoke_p2/` |

## Precedence Rules (Implemented)

- **Feature ref**: `star_feature_ref` (per-library) > `--crFeatureRef` (CLI) > `[feature] ref` (config global)
- **Chemistry**: `star_chemistry` (per-library) > `--crChemistry` (CLI) > `auto`
- **Global ref requirement**: Only required if at least one non-GEX library lacks `star_feature_ref`

## Open Review Questions (Inherited from Phase 1)

1. **GEX exclusion edge cases**: `find("geneexpression")` in normalized type
   would also exclude "Gene Expression Plus". Acceptable?
2. **Filter-skip safety**: Per-library ref used without type validation.
   Should we warn if it contains multiple feature types?
3. **Auto-generated library ID stability**: Row-index-based (`_0`, `_1`) is
   unstable across reorders. Recommend explicit `star_library_id` in docs?
4. **`assignOut` path change**: Auto-generated IDs like `A375_CRISPR_Guide_Capture_1`
   replace old `sampleName`-based paths. Backward compatible for downstream tools?

## Commit Status

All changes are in the working tree (uncommitted) on branch
`core-alignment-threads-integration`. The branch is 15 commits ahead of
origin. No changes have been committed for the multi-feature work yet;
consider squashing into a single commit before merge.
