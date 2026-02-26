# Handoff: Multi-Feature Libraries — Phase 3: Merge and Provenance Hardening

Date: 2026-02-26  
Branch: `multi-feature`  
Prior: `docs/HANDOFF_MULTI_FEATURE_LIBRARIES_PHASE012_20260226.md`  
Runbook: `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md`

## What Was Done

### Problem

Phase 2 left three issues in the merge pipeline:

1. **Feature type mislabeling (High)**: `PfMultiMexStub::processAssignOutput` preserved the
   `feature_type` column from the feature reference CSV verbatim. When the CSV used "Custom"
   for gRNA guides, the merged MEX labeled those features as "Custom" instead of
   "CRISPR Guide Capture", causing `runCrisprFeatureCalling()` to find zero CRISPR features.

2. **No per-library provenance (Medium)**: The only provenance was the prep log in `Log.out`.
   No machine-readable manifest existed per library with resolved inputs, detected match mode,
   and effective chemistry.

3. **Insufficient merge diagnostics (Low)**: MEX-read failures and merge summaries lacked
   library identity context, making multi-library debugging difficult.

### Changes

**1. Feature type override in MEX stub**

- Added `featureTypeOverride` parameter to `PfMultiMexStub::writeFeaturesTsv` and
  `processAssignOutput`. When non-empty, all features in `features.tsv` are labeled with
  the override value regardless of what the CSV says.
- The call site in `PfMultiProcess.cpp` passes `run.featureType` (the config-declared
  `featureRefType` from routing) as the override, making the config the source of truth.
- This ensures CRISPR Guide Capture features are correctly labeled even when the upstream
  feature reference CSV uses a generic type like "Custom".

**2. Per-library provenance manifests**

- After all libraries pass assign + stub + MEX-read, a `pf_library_provenance.tsv` file
  is written to each assign output directory with key-value pairs:
  - `library_id`, `sample`, `feature_type`, `fastq_dir`, `feature_ref`, `whitelist`
  - `chemistry_request`, `chemistry_explicit`, `effective_chemistry`
  - `detected_match_mode`, `return_code`, `status`, `assign_output_dir`
- Manifests are written only after all validation stages pass (no partial provenance
  on failure). The `status` field is always `OK` in written manifests; a FAIL condition
  triggers a `runtime_error` before manifest writing.

**3. Top-level multi-library summary table**

- After all library assigns complete, a structured table is logged to `Log.out` showing
  `library_id`, `feature_type`, `sample`, `effective_chem`, `match_mode`, `status`,
  and `assign_out` for each library. All entries show `OK` because fail-fast throws
  before the table is emitted if any library fails.

**4. Merged feature type breakdown**

- After `mergeMex()`, a per-type feature count breakdown is logged showing how many
  features of each type (Gene Expression, CRISPR Guide Capture, Custom, etc.) ended up
  in the merged matrix, along with total barcodes and triplets.

**5. Enriched FeatureRun struct**

- `FeatureRun` now carries `libraryId`, `sampleName`, `resolvedFastq`,
  `resolvedChemRequest`, `explicitChem`, and `returnCode` for provenance and diagnostics.

**6. Fail-fast on stub and MEX-read failures**

- MEX stub generation failure throws `runtime_error` with library type, ID, and output
  path. Previously was warning-only and continued.
- MEX-read failure throws `runtime_error` with library type, ID, output path, and the
  underlying exception. Previously was warning-only and continued with partial data.
- Fail-fast throws before manifest writing and before the summary table, so neither
  artifact is emitted on failure.

**7. Stub always force-writes (rerun freshness)**

- Production now calls `processAssignOutput` with `force=true`. This ensures that if a
  rerun uses a different config (e.g., changed `featureRefType`), `features.tsv` and
  `barcodes.tsv` are rewritten with the current override labels instead of silently
  keeping stale outputs from a previous run.

**8. Stale manifest cleanup on rerun**

- Before stub/MEX validation, existing `pf_library_provenance.tsv` and `.tmp` files
  are deleted from each assign output directory. `ENOENT` is ignored; other errors
  (permission, I/O) throw. This ensures a failed rerun never leaves behind a stale
  `status=OK` manifest from a previous successful run.

**9. Atomic provenance manifests**

- Manifests are written to `.tmp` files first, then atomically renamed. If any manifest
  write or rename fails, all `.tmp` files are cleaned up before re-throwing. This prevents
  partial or truncated `status=OK` manifests from surviving I/O failures.

## Files Changed

| File | Change |
|------|--------|
| `core/legacy/source/PfMultiMexStub.h` | Added `featureTypeOverride` parameter to `writeFeaturesTsv` and `processAssignOutput` |
| `core/legacy/source/PfMultiMexStub.cpp` | `writeFeaturesTsv` applies override when non-empty; `processAssignOutput` passes it through; idempotent `!wroteAny` path returns 0 when outputs exist and `force=false` |
| `core/legacy/source/PfMultiProcess.cpp` | Enriched `FeatureRun` struct; clears stale manifests before validation; calls stub with `force=true` (rerun freshness); writes `pf_library_provenance.tsv` atomically (tmp+rename, cleanup on failure); logs summary table and feature breakdown; fail-fast `runtime_error` on stub/MEX-read failure |
| `tests/multi_feature/test_phase3_feature_type_override.sh` | 5-case test: no override, CRISPR override, Antibody override, empty CSV + override, empty CSV + no override |
| `tests/multi_feature/test_phase3_provenance_and_failfast.sh` | 8-case test: manifest key completeness (force=true), status=OK, feature type override, rerun-with-changed-config overwrites stale labels, stub failure is hard error, MEX-read failure is hard error, manifest not written on failure, failed rerun clears stale manifest |

## Test Results

| Suite | Result |
|-------|--------|
| `test_multi_feature_config.sh` | 21/21 pass |
| `test_data_driven_specs.sh` | 16/16 pass |
| `test_no_global_ref_guard.sh` | 4/4 pass |
| `test_autodetect.sh` | 8/8 pass |
| `test_star_chemistry_column.sh` | 15/15 pass |
| `test_phase3_feature_type_override.sh` | 5/5 pass |
| `test_phase3_provenance_and_failfast.sh` | 8/8 pass |
| **Total** | **77/77 pass** |

## Exit Criteria Check (from Runbook Phase 3)

| Criterion | Status |
|-----------|--------|
| Output includes all configured feature types in merged matrices with correct labels | **Pass** — `featureTypeOverride` ensures config-declared types are used |
| Provenance shows exact resolved inputs per library | **Pass** — `pf_library_provenance.tsv` per library + summary table in `Log.out` |
| Guard against silent partial merge | **Pass** — Stub failures and MEX-read failures are now `runtime_error` (fail-fast); provenance manifests deferred until all validations pass; no partial/misleading status |

## What Remains (Phases 4–5)

| Phase | Description | Status |
|-------|-------------|--------|
| **Phase 3** | Merge and provenance hardening | **Complete** |
| **Phase 4** | PF intra-domain scheduler: library-aware permit allocation | **Complete** (see `HANDOFF_MULTI_FEATURE_PHASE4_LIBRARY_SCHEDULER_20260226.md`) |
| **Phase 5** | Multi-library validation suite + E2E parity tests + benchmarks | Not started |

## How to Verify

```bash
# Build
make -C core/legacy/source -j8 STAR

# Fast unit tests (no fixtures needed)
bash tests/autodetect_nxt_tru/test_autodetect.sh                    # 8/8
bash tests/autodetect_nxt_tru/test_star_chemistry_column.sh         # 15/15
bash tests/multi_feature/test_multi_feature_config.sh               # 21/21
bash tests/multi_feature/test_data_driven_specs.sh                  # 16/16
bash tests/multi_feature/test_no_global_ref_guard.sh                # 4/4
bash tests/multi_feature/test_phase3_feature_type_override.sh       # 5/5
bash tests/multi_feature/test_phase3_provenance_and_failfast.sh     # 8/8

# E2E (requires fixtures)
bash tests/multi_feature/run_msk_multifeature_smoke.sh              # 9/9
bash tests/test_cr_compat_crispr_calling.sh                         # regression
```

## Output Structure (After Phase 3)

```
<outPrefix>/
├── cr_assign/
│   ├── CRISPR_Guide_Capture/
│   │   └── grna_de/
│   │       └── PolyIII/
│   │           ├── features.tsv          ← now labeled "CRISPR Guide Capture"
│   │           ├── barcodes.tsv
│   │           ├── matrix.mtx
│   │           └── pf_library_provenance.tsv   ← NEW
│   └── Custom/
│       └── larry_de/
│           └── LARRY/
│               ├── features.tsv          ← labeled "Custom"
│               ├── barcodes.tsv
│               ├── matrix.mtx
│               └── pf_library_provenance.tsv   ← NEW
├── outs/
│   ├── raw_feature_bc_matrix/            ← merged with correct types
│   ├── filtered_feature_bc_matrix/
│   └── crispr_analysis/                  ← now generated (features found)
├── Solo.out/
└── Log.out                               ← includes summary table + breakdown
```

## Provenance Manifest Format

`pf_library_provenance.tsv` (tab-separated key-value):

```
key	value
library_id	grna_de
sample	S1
feature_type	CRISPR Guide Capture
fastq_dir	/path/to/PolyIII
feature_ref	/path/to/ref_feature_geneBC.csv
whitelist	/path/to/3M-february-2018_TRU.txt
chemistry_request	nxt
chemistry_explicit	yes
effective_chemistry	NXT
detected_match_mode	(not applicable when explicit)
return_code	0
status	OK
assign_output_dir	/tmp/out/cr_assign/CRISPR_Guide_Capture/grna_de/PolyIII
```

Note: Before each run, existing manifests and `.tmp` files are deleted. New manifests are
written atomically (tmp file + rename) only after all stages (assign, stub generation,
MEX read) succeed. If any stage fails, a `runtime_error` is thrown and no manifests exist.
If a manifest write itself fails, all `.tmp` files are cleaned up. The `status` field is
therefore always `OK` in committed manifests, and stale manifests never survive failed reruns.
