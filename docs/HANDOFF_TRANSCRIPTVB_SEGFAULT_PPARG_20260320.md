# Handoff: TranscriptVB Segfault on PPARG 35M PE Dataset (2026-03-20)

## Status: OPEN — segfault reproduces on both `master` and `benchmark-flex`

## Summary

`--quantMode TranscriptomeSAM TranscriptVB` segfaults (signal 11) on the
PPARG 35M-read-pair PE dataset. The crash occurs immediately after
"started transcript quantification" — mapping and BAM sorting complete
successfully. Removing `TranscriptVB` from `--quantMode` (leaving just
`TranscriptomeSAM`) succeeds without issue.

This is a **pre-existing bug on `master`**, not introduced by the
`benchmark-flex` branch.

## Reproduction

### Minimal command (crashes)

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR

/mnt/pikachu/STAR-suite/core/legacy/source/STAR \
  --runMode alignReads \
  --runThreadN 32 \
  --genomeDir /storage/autoindex_110_44/bulk_index \
  --readFilesIn \
    /storage/JAX_PE/PPARG_R_WT_3_GT24-02468_ATGAGGCC-CAATTAAC_S21_L001_R1_001.fastq.gz \
    /storage/JAX_PE/PPARG_R_WT_3_GT24-02468_ATGAGGCC-CAATTAAC_S21_L001_R2_001.fastq.gz \
  --readFilesCommand zcat \
  --trimCutadapt Yes \
  --trimCutadaptQuality 20 \
  --trimCutadaptMinLength 20 \
  --trimCutadaptAdapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortMethod samtools \
  --keepBAM yes \
  --quantMode TranscriptomeSAM TranscriptVB \
  --quantVBgcBias 1 \
  --outFileNamePrefix /tmp/vb_crash_repro/
```

### Same command without TranscriptVB (succeeds)

Replace `--quantMode TranscriptomeSAM TranscriptVB` with
`--quantMode TranscriptomeSAM` and remove `--quantVBgcBias 1`.
Completes in 7m17s with exit code 0.

## Crash Details

### On `benchmark-flex` (commit ccf1735)

| Metric | Value |
|--------|-------|
| Output dir | `/tmp/pe_bulk_benchmark_pparg_noy_20260320_085841/storage/integrated/` |
| Signal | 11 (SIGSEGV) |
| Wall time at crash | 7m42s |
| Max RSS | 66.9 GB (65.3 GB VmRSS) |
| Last log line | `Mar 20 09:06:22 ..... started transcript quantification` |
| Mapping completed | Yes (09:04:13) |
| BAM sort completed | Yes (09:06:22) |
| Transcriptome BAM | 6.9 GB (written successfully) |
| Coord-sorted BAM | 4.5 GB (written successfully) |
| quant.sf output | Not created (crash before completion) |

### On `master` (9 commits ahead of origin/master)

| Metric | Value |
|--------|-------|
| Output dir | `/tmp/pe_bulk_vb_crash_test_master_20260320_091913/` |
| Signal | 11 (SIGSEGV) |
| Wall time at crash | 7m54s |
| Max RSS | 66.9 GB |
| Last log line | `Mar 20 09:27:03 ..... started transcript quantification` |
| Transcriptome BAM | 6.9 GB (written successfully) |
| Coord-sorted BAM | 4.5 GB (written successfully) |
| quant.sf output | Not created |

### Successful run (no TranscriptVB, on `benchmark-flex`)

| Metric | Value |
|--------|-------|
| Output dir | `/tmp/pe_bulk_benchmark_pparg_manual_20260320_090812/integrated/` |
| Exit code | 0 |
| Wall time | 7m17s |
| Max RSS | 65.7 GB |
| `--quantMode` | `TranscriptomeSAM` (no `TranscriptVB`) |

## What Works vs What Crashes

| Dataset | Read pairs | `TranscriptVB` | Result |
|---------|-----------|----------------|--------|
| 21033 L007 (PE benchmark) | 6.5M | Yes | OK (5.7–29.2s) |
| PPARG JAX_PE | 35.1M | No | OK (7m17s) |
| PPARG JAX_PE | 35.1M | Yes | **SIGSEGV** |

The 6.5M sample works with TranscriptVB. The 35.1M sample crashes.
This suggests a size-dependent bug: either a buffer overflow, integer
overflow, or memory allocation failure in the VB engine that manifests
at scale.

## Key Observations

1. **The crash is instantaneous.** "started transcript quantification"
   appears in the log, then immediate SIGSEGV. The VB engine crashes
   during initialization, not during iterative convergence.

2. **Memory is not exhausted.** Peak RSS is ~67 GB on a 128 GB machine.
   The crash is not an OOM.

3. **The transcriptome BAM is valid.** Both the transcriptome BAM
   (6.9 GB) and coordinate-sorted BAM (4.5 GB) are fully written
   before the VB step starts.

4. **Branch-independent.** Identical crash on `master` and
   `benchmark-flex`. The `benchmark-flex` diff touches only Flex hash
   screen code and `ReadAlign_outputAlignments.cpp` (synthetic probe
   guard); no VB/quant code was modified.

5. **`--quantVBgcBias 1` may be a factor.** Both crashes used
   `--quantVBgcBias 1`. This was not tested with `--quantVBgcBias 0`
   on the PPARG dataset. The GC bias model allocates additional
   per-transcript arrays that could be the crash site.

## Likely Code Path

The TranscriptVB quantification is triggered in `STAR.cpp` after BAM
sort. The entry point is in `TranscriptQuantEC.cpp` which calls into
`core/features/vbem/source/libem/`. The VB engine:

1. Loads equivalence classes from the transcriptome BAM
2. Builds the EC matrix
3. Allocates per-transcript arrays (effective lengths, GC bias, etc.)
4. Runs VB/EM iterations

The crash at "started transcript quantification" means it's dying in
steps 1–3 (loading/allocation), not during iterative convergence.

### Files to investigate

| File | Role |
|------|------|
| `core/legacy/source/TranscriptQuantEC.cpp` | Entry point for VB quant after BAM sort |
| `core/features/vbem/source/libem/ec_loader.cpp` | EC loading from transcriptome BAM |
| `core/features/vbem/source/libem/ec_builder.cpp` | EC matrix construction |
| `core/features/vbem/source/libem/vb_engine.cpp` | VB iterations |
| `core/features/vbem/source/libem/gc_bias.cpp` | GC bias model (if `--quantVBgcBias 1`) |
| `core/features/vbem/source/libem/effective_length.cpp` | Effective length computation |

## Suggested Debugging Steps

## Root Cause and Fix

This was traced under `gdb` on a clean debug build and is **not** an EC/VB
math crash. The segfault was a null dereference in
`core/legacy/source/STAR.cpp` during the TranscriptVB merge step:

- crash site: `STAR.cpp:1848`
- `RAchunk == nullptr` at crash time

The immediate cause was a recent optimization introduced in commit
`54f1108` (`flex: Solo phase optimization — 5.1x speedup`) on
2026-03-19. That commit added an early `RAchunk` free block in
`STAR.cpp` before transcript quantification. The predicate only checked
`outSAM`, gene-count quant, Y-read emission, and batch mode, but it
missed later consumers that still read chunk-local state:

- `TranscriptVB`
- `SLAM`
- trim QC merge

As a result, PE runs with `--quantMode TranscriptomeSAM TranscriptVB`
could free `RAchunk` immediately after mapping and then crash as soon as
TranscriptVB started.

The fix is now in `core/legacy/source/STAR.cpp`: the early-free guard
keeps `RAchunk` alive whenever any downstream phase still needs it:

- `P.quant.transcriptVB.yes`
- `P.quant.slam.yes`
- `P.trimQcEnabled`

Validation:

- clean rebuild completed successfully
- patched smoke repro on the PPARG pair with
  `--quantMode TranscriptomeSAM TranscriptVB --quantVBgcBias 1`
  and `--readMapNumber 200000` completed transcript quantification
  without crashing
- validation artifact:
  `/tmp/pe_bulk_vb_patch_smoke_20260320_1009/`

1. **Compile with `-g -O0 -fsanitize=address`** and reproduce:
   ```bash
   # In core/legacy/source/Makefile, add to CXXFLAGS:
   #   -g -O0 -fsanitize=address -fno-omit-frame-pointer
   # Rebuild and run the crash command. ASAN will report the exact
   # allocation and access site.
   ```

2. **Test without GC bias**: run with `--quantVBgcBias 0` on the PPARG
   dataset to isolate whether the GC bias model is the crash site.

3. **Test with fewer threads**: try `--runThreadN 1` to rule out any
   threading issue in the VB init.

4. **Get a stack trace**: if ASAN is too slow on 35M reads, use
   `gdb --batch -ex run -ex bt` or `coredumpctl` to get a backtrace.

5. **Check for integer overflow**: the transcriptome BAM has ~35M
   alignments. Look for `int` or `uint32_t` counters in the EC loader
   that might overflow with large alignment counts. Particularly check
   array sizing based on read counts.

## Repo State Left by This Session

- **Current branch**: `master` (switched from `benchmark-flex` to test
  the crash; confirmed it reproduces on both)
- **Stash**: `stash@{0}` contains the `benchmark-flex` working tree
  (`git stash pop` to restore)
- **STAR binary**: built from `master` (clean build at 09:18:36)
- **To return to benchmark-flex**: `git checkout benchmark-flex && git stash pop`
- **To return to clean state**: clean rebuild required after any checkout

## Artifacts

| Item | Path |
|------|------|
| Crash run (benchmark-flex) | `/tmp/pe_bulk_benchmark_pparg_noy_20260320_085841/storage/integrated/` |
| Crash run (master) | `/tmp/pe_bulk_vb_crash_test_master_20260320_091913/` |
| Successful run (no VB) | `/tmp/pe_bulk_benchmark_pparg_manual_20260320_090812/integrated/` |
| PPARG FASTQs | `/storage/JAX_PE/PPARG_R_WT_3_GT24-02468_ATGAGGCC-CAATTAAC_S21_L001_{R1,R2}_001.fastq.gz` |
| Small PE FASTQs (works) | `/storage/PE/21033-09-01-13-01_S1_L007_{R1,R2}_001.fastq.gz` |
| STAR index | `/storage/autoindex_110_44/bulk_index` |
| Benchmark script | `scripts/paper/run_pe_bulk_feature_benchmark.sh` |
| Tool env (reusable) | `/tmp/pe_bulk_feature_benchmark_full_iurun_20260310_033209/tools/tool_env.sh` |

## Context: Why This Was Being Investigated

The PE bulk benchmark in the README uses a 6.5M read pair sample
(`21033-09-01-13-01_S1_L007`). To get more representative benchmark
numbers for the paper, we attempted to run the same benchmark on the
larger PPARG sample (35.1M read pairs, 5.4x larger). The benchmark
script (`scripts/paper/run_pe_bulk_feature_benchmark.sh`) uses
`--quantMode TranscriptomeSAM TranscriptVB --quantVBgcBias 1` for the
integrated arm, which hits this segfault.

**Workaround for the benchmark**: Run with `--quantMode TranscriptomeSAM`
only (no TranscriptVB) and use Salmon on both arms. This still captures
the speed advantage from eliminating decompress/trim/Y-removal steps.
The integrated STAR step completed in 7m17s with this workaround.
