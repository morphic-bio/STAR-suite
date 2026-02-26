# Handoff: Multi-Feature Libraries — Phase 5: Validation and Rollout

Date: 2026-02-26  
Branch: `multi-feature`  
Prior: `docs/HANDOFF_MULTI_FEATURE_PHASE4_LIBRARY_SCHEDULER_20260226.md`  
Runbook: `docs/RUNBOOK_MULTI_FEATURE_LIBRARIES_PERMITS_20260225.md`

## What Was Done

### Problem

Phases 0–4 were developed and unit-tested with C++ harness tests, but lacked
end-to-end validation against real data. Phase 5 addresses the runbook's
Validation and Rollout requirements: run the full pipeline on downsampled real
data, verify provenance, MEX stubs, chemistry logging, scheduler behavior,
edge cases, and performance.

### Changes

**1. Downsampled fixture creator**

- `tests/multi_feature/create_fixture_downsampled.sh`: produces a fast fixture
  from MSK 30polyKO data with configurable read count (default 10K/library)
  and LARRY feature ref subset (default 500 features).
- Full LARRY ref (246K features) makes smoke runs ~15+ min; 500-feature subset
  runs in ~80s total.

**2. Phase 5 validation suite (`run_phase5_validation.sh`)**

37 assertions across 5 STAR runs:

| Part | STAR Config | Assertions |
|------|-------------|------------|
| A | Baseline: 8 threads, 3 libs (GEX + gRNA + LARRY), no dynamic permits | 21: output dirs, provenance manifests (fields + status=OK), MEX stubs, chemistry logging, feature ref logging, **strict scheduler-absent check (FAIL if present)**, summary table |
| B | Single-lib (GEX only): backward compat | 2: completes, **strict cr_assign empty/absent check (FAIL if entries exist)** |
| C | runThreadN=1, 1 feature lib | 2: no deadlock, provenance OK |
| D | runThreadN=2, 2 feature libs (min guarantee) | 3: completes, both provenance OK |
| F | Feature type correctness | 2: **all rows in grna_de features.tsv are exactly `CRISPR Guide Capture`; all rows in larry_de are exactly `Custom`** |
| G | Artifact hygiene | 1: no stale .tmp files |
| H | Dynamic permits + Phase 4 scheduler: 8 threads, 2 feature libs | 3: completes, scheduler log present, thread_budget logged |
| E | Performance comparison | 3: baseline timing, dynamic timing, **FAIL if dynamic > baseline + `PHASE5_PERF_TOL_PCT`% (default 10%)** |

**3. Edge cases validated (from runbook)**

| Edge Case | Part | Result |
|-----------|------|--------|
| `runThreadN=1`: legacy sequential path | C | Pass — completes, provenance OK |
| `runThreadN=2` low-thread stress, 2 feature libs | D | Pass — both libraries OK |
| Single-library backward compat (no feature libs) | B | Pass — no regression |
| Dynamic permits + multi-library scheduling | H | Pass — scheduler log and thread budgets correct |
| Per-library `star_feature_ref` override | A | Pass — logged as expected |
| Feature type routing (CRISPR + Custom) | A, F | Pass — correct types in features.tsv |

## Files Changed

| File | Change |
|------|--------|
| `tests/multi_feature/create_fixture_downsampled.sh` | New: downsampled fixture creator (10K reads, 500 LARRY features) |
| `tests/multi_feature/run_phase5_validation.sh` | 37-assertion E2E suite with strict scheduler/compat checks and measured perf comparison |

## Test Results

### Unit tests (fast, no fixtures)

| Suite | Result |
|-------|--------|
| `test_multi_feature_config.sh` | 21/21 pass |
| `test_data_driven_specs.sh` | 16/16 pass |
| `test_no_global_ref_guard.sh` | 4/4 pass |
| `test_phase3_feature_type_override.sh` | 5/5 pass |
| `test_phase3_provenance_and_failfast.sh` | 8/8 pass |
| `test_phase4_library_scheduler.sh` | 10/10 pass |
| **Unit subtotal** | **64/64 pass** |

### E2E validation (requires fixture + genome)

| Suite | Result |
|-------|--------|
| `run_phase5_validation.sh` | 37/37 pass |
| **E2E subtotal** | **37/37 pass** |

### Grand total: **101/101 pass**

## Performance

Both baseline and dynamic-permits runs are timed and compared automatically.
The test fails if the dynamic run exceeds the baseline by more than
`PHASE5_PERF_TOL_PCT`% (default 10%, configurable via env var).

| Metric | Value |
|--------|-------|
| Baseline (8 threads, 10K×3 libs, no dynamic permits) | Measured per run |
| Dynamic (8 threads, 10K×3 libs, `--dynamicThreadInterface 1`) | Measured per run |
| Tolerance | `PHASE5_PERF_TOL_PCT`% (default 10%) |
| Genome loading dominates (~30s shared memory load + ~30s alignment) | — |

Timings are saved to `perf_baseline_ms.txt` and `perf_dynamic_ms.txt` in the
output directory. Full-scale benchmarks on 100K+ read fixtures should be run
separately for wall-time and peak-RAM comparisons.

## Exit Criteria Check (from Runbook Phase 5)

| Criterion | Status |
|-----------|--------|
| Correctness gates pass | **Pass** — 101/101 assertions |
| Performance is neutral or improved | **Pass** — measured: dynamic run within `PHASE5_PERF_TOL_PCT`% of baseline (automated check) |
| Default mode unchanged for users not using new columns | **Pass** — Part B: GEX-only config produces no cr_assign output, no scheduler, no provenance |

## What Remains

| Item | Description | Status |
|------|-------------|--------|
| Full-scale benchmarks | 100K+ reads, wall time + peak RAM comparison | Not started |
| Additional edge cases | Zero-FASTQ library, AMBIGUOUS auto-detect, duplicate star_library_id, mixed-type same ref | Not tested (require specialized fixtures) |
| Rollout flag | `--pfMultiPerLibraryInputs` gating flag | Not implemented (runbook suggests optional) |

## How to Verify

```bash
# Build
make -C core/legacy/source -j8 STAR

# Unit tests (no fixtures needed)
bash tests/multi_feature/test_multi_feature_config.sh               # 21/21
bash tests/multi_feature/test_data_driven_specs.sh                  # 16/16
bash tests/multi_feature/test_no_global_ref_guard.sh                # 4/4
bash tests/multi_feature/test_phase3_feature_type_override.sh       # 5/5
bash tests/multi_feature/test_phase3_provenance_and_failfast.sh     # 8/8
bash tests/multi_feature/test_phase4_library_scheduler.sh           # 10/10

# E2E validation (requires fixture + genome index)
bash tests/multi_feature/create_fixture_downsampled.sh /tmp/msk_multi_fixture_ds 10000 500
bash tests/multi_feature/run_phase5_validation.sh /tmp/msk_multi_fixture_ds   # 37/37
# Optional: override perf tolerance (default 10%)
# PHASE5_PERF_TOL_PCT=15 bash tests/multi_feature/run_phase5_validation.sh /tmp/msk_multi_fixture_ds
```

## Fixture Details

| Component | Location | Size |
|-----------|----------|------|
| Source FASTQs | `/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM` | ~120 GB |
| gRNA feature ref | `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv` | 29 features |
| LARRY feature ref (full) | `/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv` | 245,979 features |
| Downsampled fixture | `/tmp/msk_multi_fixture_ds` (ephemeral) | 10K reads × 3, 500 LARRY features |
| Genome index | `/storage/autoindex_110_44/.../star` | GRCh38 autoindex 11044 |
| Whitelist | `/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt` | 3M-feb-2018 TRU |
