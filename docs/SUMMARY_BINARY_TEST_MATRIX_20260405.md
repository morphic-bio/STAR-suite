# Binary Test Matrix — Implementation Summary

**Date:** 2026-04-05
**Branch:** `feature-binary-tests`
**Status:** Tier A green (15 passed, 0 failed, 0 skipped), not yet committed

## Overview

Expanded the binary validation suite from the original 5-test workflow matrix
to a 15-test Tier A + 8-test Tier B regression suite. The new tests cover
genomeGenerate, BAM output modes, Y-removal, TRU/NXT chemistry detection,
adapter trimming, trim QC, cutadapt parity, and TranscriptVB/Salmon.

Based on [RUNBOOK_BINARY_TEST_MATRIX_20260405.md](RUNBOOK_BINARY_TEST_MATRIX_20260405.md).

## New Files

| File | Purpose |
|------|---------|
| `tests/run_binary_test_matrix.sh` | Top-level Tier A/B runner with PASS/FAIL/SKIP tracking |
| `tests/run_gex_binary_smoke.sh` | scRNA GEX smoke — synthetic CB+UMI on public chr22+chrY |
| `tests/run_flex_tiny_public_binary_smoke.sh` | Flex smoke — public tiny 10x fixture |
| `tests/run_genome_generate_validation.sh` | **NEW** — genomeGenerate index structure validation |
| `tests/run_adapter_clip_synthetic_test.sh` | **NEW** — `--clip3pAdapterSeq` on synthetic adapter-contaminated reads |

## Modified Files

| File | Change |
|------|--------|
| `tests/run_binary_test_matrix.sh` | Expanded from 5 Tier A + 5 Tier B to 15 Tier A + 8 Tier B entries |
| `tests/run_y_removal_comprehensive_test.sh` | Fixed tool path (now uses `core/features/yremove_fastq/...`) |
| `tests/test_cutadapt32_parity.sh` | Fixed trimvalidate path (now uses `core/features/vbem/...`) |
| `tests/run_gex_binary_smoke.sh` | Added `--limitBAMsortRAM` for sorted BAM mode |
| `scripts/release/docker/container_check_release_smokes.sh` | Added `--profile` flag (core, binary-workflows-tier-a, all) |
| `scripts/release/docker/Dockerfile.runtime-check` | Added `curl` and `salmon` to apt-get install |
| `core/legacy/source/STAR.cpp` | **Bug fix:** added `!P.outBAMcoord` to RAchunk early-free guard (use-after-free) |
| `docs/codespaces/RUNBOOK_bulk.md` | Fixed Genome file size from ~60 MB to ~106 MB |

## Tier A Results

Full matrix run on host `pikachu` (2026-04-05T08:32:37Z):

| # | Test | Feature Covered | Time | Result |
|---|------|-----------------|------|--------|
| 1 | `gex-binary-smoke-no-bam` | scRNA GEX workflow, Solo counting, no-BAM | ~29s | PASS |
| 2 | `perturb-a375-5prime` | Perturb-seq, multi-library, public download | ~89s | PASS |
| 3 | `flex-tiny-public-binary` | Flex pipeline, probe catalog, sample detection | ~29s | PASS |
| 4 | `pe-bulk-public` | PE bulk, adapter trim, Y-removal, Salmon/VB | ~41s | PASS |
| 5 | `slam-snp-mask-build` | SLAM SNP mask construction | <1s | PASS |
| 6 | `gex-binary-smoke-bam-sorted` | BAM SortedByCoordinate + Solo, SO:coordinate header | ~29s | PASS |
| 7 | `gex-binary-smoke-bam-unsorted` | BAM Unsorted + Solo, CB/UB tags | ~29s | PASS |
| 8 | `genome-generate-validation` | genomeGenerate index structure, chr/gene validation | <1s | PASS |
| 9 | `y-removal-unit` | Y-contig detection, path derivation | <1s | PASS |
| 10 | `y-removal-comprehensive` | remove_y_reads split, ordering, multithreaded | <1s | PASS |
| 11 | `tru-nxt-autodetect` | TRU/NXT chemistry auto-detection | <1s | PASS |
| 12 | `tru-nxt-config-parsing` | star_chemistry CSV column parsing (C++ harness) | <1s | PASS |
| 13 | `trim-qc-fastq` | trim_qc_fastq JSON output, metrics validation | <1s | PASS |
| 14 | `adapter-clip-synthetic` | --clip3pAdapterSeq mapping rate recovery | <1s | PASS |
| 15 | `cutadapt-parity` | trimvalidate Cutadapt 3.2 compatibility | <1s | PASS |

Total wall time: ~4 min 8s. Logs: `/tmp/binary_test_matrix_20260405_083237_677537`

## Tier B Results

> Historical result from 2026-04-05. As of 2026-07-24, Tier B uses
> `slam-pinned-fixture-parity` / `tests/run_slam_parity_smoke.sh`. The former
> `slam-e2e-fixture` implementation is a deprecated delegating shim.

| # | Test | Result | Notes |
|---|------|--------|-------|
| 1 | `sorted-bam-cbub-nonflex` | FAIL | Hardcoded `/storage/` paths |
| 2 | `unsorted-bam-cbub` | FAIL | Hardcoded `/storage/` paths |
| 3 | `solo-yremove` | FAIL | Needs local artifacts |
| 4 | `flex-cr-config` | PASS | |
| 5 | `slam-e2e-fixture` | FAIL | Needs SLAM fixture |
| 6 | `transcriptvb-quick` | SKIP | VB genome index not built |
| 7 | `transcriptvb-edge-cases` | SKIP | VB genome index not built |
| 8 | `salmon-parity` | SKIP | VB genome + salmon prereqs |

Tier B failures are expected — these tests require host-specific artifacts that
are not yet portable.

## Feature Coverage Matrix

| Feature | Tier A Tests | Tier B Tests |
|---------|-------------|-------------|
| scRNA GEX workflow | gex-binary-smoke (no-BAM, unsorted) | — |
| Perturb-seq | perturb-a375-5prime | — |
| Flex pipeline | flex-tiny-public-binary | flex-cr-config |
| PE bulk workflow | pe-bulk-public | — |
| SLAM | slam-snp-mask-build, slam-fixture-contract | slam-pinned-fixture-parity |
| genomeGenerate | genome-generate-validation | — |
| BAM SortedByCoordinate | gex-binary-smoke-bam-sorted | sorted-bam-cbub-nonflex |
| BAM Unsorted | gex-binary-smoke-bam-unsorted | unsorted-bam-cbub |
| Y-removal | y-removal-unit, y-removal-comprehensive | solo-yremove |
| TRU/NXT detection | tru-nxt-autodetect, tru-nxt-config-parsing | — |
| Adapter trimming | adapter-clip-synthetic, cutadapt-parity | — |
| Trim QC | trim-qc-fastq | — |
| TranscriptVB/Salmon | *(via pe-bulk-public)* | transcriptvb-quick, edge-cases, salmon-parity |

## New Test Details

### genomeGenerate Validation (`run_genome_generate_validation.sh`)

- Generates synthetic 2-chromosome reference (chr1=500bp, chrY=300bp) with GTF
- Runs `STAR --runMode genomeGenerate`
- Validates 22 checks: all index files (Genome, SA, SAindex, chrName.txt,
  chrLength.txt, chrStart.txt, chrNameLength.txt, genomeParameters.txt,
  sjdbList.fromGTF.out.tab, sjdbInfo.txt, geneInfo.tab), chromosome names
  and lengths, gene annotations, genomeSAindexNbases parameter, versionGenome
- Quick alignment sanity: loads index and maps 5 synthetic reads

### Adapter Clip Synthetic (`run_adapter_clip_synthetic_test.sh`)

- Generates 50kb reference, builds index
- Creates 200 clean reads and 200 adapter-contaminated reads (30-59bp genomic +
  AGATCGGAAGAG adapter fill)
- Runs STAR three ways: no trim (baseline), with `--clip3pAdapterSeq`, clean reads
- Validates: trimming improves mapping rate (17.5% -> 100%), trimmed rate matches
  clean rate, all runs produce Log.final.out

## Bug Fix: RAchunk use-after-free in BAM SortedByCoordinate

The BAM sorted test initially crashed with a non-deterministic segfault during
`BAMbinSortByCoordinate`. Investigation traced it to a **use-after-free** in
[STAR.cpp:1767](../core/legacy/source/STAR.cpp#L1767).

**Root cause:** The early-free guard for `RAchunk` checked `!P.outSAMbool` but
not `!P.outBAMcoord`. When `--outSAMtype BAM SortedByCoordinate` is used,
`outSAMbool` is false (no SAM output) but `outBAMcoord` is true. RAchunk was
freed after Solo counting, then the BAM sorter at
[bamSortByCoordinate.cpp:217](../core/legacy/source/bamSortByCoordinate.cpp#L217)
tried to read `RAchunk[it]->chunkOutBAMcoord->binTotalN[ibin]` from freed memory.

**Fix:** Added `!P.outBAMcoord` to the guard condition:
```cpp
if (RAchunk != nullptr && !P.outSAMbool && !P.outBAMcoord && !P.quant.geCount.yes
    && !P.quant.transcriptVB.yes && !P.quant.slam.yes && !P.trimQcEnabled
    && !P.emitYReadNamesyes && !P.emitYNoYFastqyes && !batchModeActive) {
```

**Verification:** 5 consecutive runs with stdout redirect (the worst-case trigger)
all pass. Full matrix green 15/15.

## Known Issues

1. **Tier B hardcoded paths**: `run_sorted_bam_cbub_nonflex_test.sh`,
   `run_unsorted_bam_cbub_test.sh`, and `run_solo_yremove_smoke.sh` have
   hardcoded `/storage/` paths and need env-var overrides or rewrite to be
   portable.

2. **TranscriptVB genome index**: The VB tests need a pre-built genome index at
   `/tmp/star_vb_test/star_new_index` and GSE110004 reads. These are available
   on this host but not built automatically.

## Design Decisions

1. **Self-contained Tier A**: All Tier A tests generate their own data or use
   public downloads. No host-specific fixtures required.

2. **Skip vs fail**: Tests check prerequisites (binaries, fixtures, genome dirs)
   and skip gracefully rather than failing the matrix. This lets the matrix
   provide useful results even in partial environments.

3. **Tool path fixes**: `run_y_removal_comprehensive_test.sh` and
   `test_cutadapt32_parity.sh` had stale tool paths from before the repo
   reorganization. Fixed to use `core/features/` paths with env-var overrides.

4. **Adapter clip validation**: Tests mapping rate improvement rather than exact
   trim positions. Baseline (no trim) maps 17.5% of adapter-contaminated reads;
   with `--clip3pAdapterSeq` all 200 map (100%), matching the clean control.

## What Remains

### Not yet committed
All new/modified files listed above are unstaged on `feature-binary-tests`.

### Future work
- **BAM sorted**: Need a larger public fixture or STAR fix for SortedByCoordinate + Solo
- **TranscriptVB portability**: Make quick_test.sh self-contained with synthetic PE reference
- **Tier B portability**: Add env-var overrides to hardcoded `/storage/` paths
- **Container validation**: Tier A has not been run inside Docker yet
- **Ubuntu 22.04**: Not yet added to CI
