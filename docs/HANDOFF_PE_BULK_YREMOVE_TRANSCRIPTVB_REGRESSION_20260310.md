# Handoff: PE Bulk Y-Removal and TranscriptVB Auto Regression

**Date**: 2026-03-10  
**Checkout**: `/mnt/pikachu/STAR-suite`  
**Branch**: `benchmark-pe`  
**Commit**: `0de2e6e`  

## Hard Restriction

Some PE files are covered by a data-use agreement. Treat all per-read
Y-classified outputs as **restricted local artifacts**.

Do not distribute, upload, or commit any of the following without explicit
user approval:

- `_Y*.fastq*`
- `_noY*.fastq*`
- `_Y.bam`
- `_noY.bam`
- `ynames.txt` or any per-read Y-name list
- any derivative file that exposes per-read Y assignment

Only aggregate metrics, logs, and summarized tables should be shared by
default. Keep restricted outputs on local storage only.

## Current Status

Two regressions were identified on the real PE downsampled fixture:

1. **Integrated Y/noY FASTQ emission bug**: fixed.
2. **TranscriptVB PE `--quantVBLibType A` autodetect bug**: fixed on the real
   downsampled case.

Benchmarking was intentionally paused. The next step is not broad benchmarking;
it is to continue debugging any remaining integrated-vs-stepwise divergence on
the downsampled PE set first.

## What Was Fixed

### 1. Integrated Y/noY FASTQ emission

Root cause:

- `RAchunk[0]` was reused across the TranscriptVB autodetect pre-pass and the
  main mapping pass.
- Y/noY FASTQ chunk streams were closed after the autodetect pass and not
  reopened before the main pass.
- Result: reads could appear correctly in BAM, but be silently missing from
  emitted Y/noY FASTQs.

Files changed:

- [ReadAlignChunk.h](/mnt/pikachu/STAR-suite/core/legacy/source/ReadAlignChunk.h)
- [ReadAlignChunk.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/ReadAlignChunk.cpp)
- [ReadAlignChunk_processChunks.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/ReadAlignChunk_processChunks.cpp)

Validation summary:

- previously missing Y read:
  `LH00150:167:22CGNHLT3:7:1106:20593:7745`
  is now present in the emitted `_Y` FASTQ and absent from `_noY`
- `bam_minus_fastq = 0` after the fix on the downsampled run
- remaining FASTQ-only reads are expected unmapped / too-many-loci classes

## 2. TranscriptVB PE auto-detect (`--quantVBLibType A`)

Problem before fix:

- the downsampled PE bulk run detected the wrong library format and then
  dropped a large number of alignments as incompatible
- old broken run:
  - [Log.out](/tmp/pe_bulk_benchmark_real_v4/downsampled/integrated/Log.out:392)
    votes: `UNKNOWN(3)=406 ISF(5)=89 UNKNOWN(11)=1 ISR(13)=307`
  - [Log.out](/tmp/pe_bulk_benchmark_real_v4/downsampled/integrated/Log.out:396)
    detected `UNKNOWN (formatID=3)`
  - [Log.out](/tmp/pe_bulk_benchmark_real_v4/downsampled/integrated/Log.out:433)
    `dropped_incompat: 265629`

Fix approach:

- rewrote observed-format derivation to use true mate geometry and true 5'
  positions
- added explicit names for outward / same-strand formats in logging
- collapsed mixed inward paired-end evidence to `IU` instead of forcing one
  exact stranded format
- updated EC compatibility logic to use 5' positions when present

Files changed:

- [LibFormatDetection.h](/mnt/pikachu/STAR-suite/core/features/vbem/source/LibFormatDetection.h)
- [LibFormatDetection.cpp](/mnt/pikachu/STAR-suite/core/features/vbem/source/LibFormatDetection.cpp)
- [TranscriptQuantEC.cpp](/mnt/pikachu/STAR-suite/core/features/vbem/source/TranscriptQuantEC.cpp)
- [alignment_filter.h](/mnt/pikachu/STAR-suite/core/features/vbem/source/libem/alignment_filter.h)
- [ec_builder.cpp](/mnt/pikachu/STAR-suite/core/features/vbem/source/libem/ec_builder.cpp)

Clean rebuild was done before validation:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

Validated fixed run:

- output dir:
  `tests/transcriptvb_auto_fix_output_20260309_235417/`
- key lines:
  - [Log.out](/mnt/pikachu/STAR-suite/tests/transcriptvb_auto_fix_output_20260309_235417/Log.out:392)
    `Library format votes: OSF(3)=32 ISF(5)=463 ISR(13)=308`
  - [Log.out](/mnt/pikachu/STAR-suite/tests/transcriptvb_auto_fix_output_20260309_235417/Log.out:395)
    `Detected library format: IU (formatID=37) from 1000 reads`
  - [Log.out](/mnt/pikachu/STAR-suite/tests/transcriptvb_auto_fix_output_20260309_235417/Log.out:432)
    `dropped_incompat: 0`
  - [Log.out](/mnt/pikachu/STAR-suite/tests/transcriptvb_auto_fix_output_20260309_235417/Log.out:433)
    `dropped_missing_mate_fields: 0`
  - [Log.out](/mnt/pikachu/STAR-suite/tests/transcriptvb_auto_fix_output_20260309_235417/Log.out:434)
    `dropped_unknown_obs_fmt: 0`
  - [Log.final.out](/mnt/pikachu/STAR-suite/tests/transcriptvb_auto_fix_output_20260309_235417/Log.final.out:21)
    `Pairs kept | 98671`

Equivalence check:

- fixed `A` output is essentially identical to the explicit `IU` control
  (`/tmp/pe_bulk_benchmark_real_v4/downsampled/integrated_IU/quant.sf`)
- totals matched at `84157`
- Pearson all transcripts was effectively `~0.99998`

## Remaining Divergence To Investigate

The big `A` autodetect regression is fixed. Remaining work is to inspect the
residual integrated-vs-stepwise divergence on the downsampled PE dataset.

Known residuals:

- small FASTQ-only Y/noY routing difference of `15` fragments between the
  integrated and external stepwise paths
- broader TranscriptVB vs Salmon count-model differences still exist, but those
  are not the autodetect bug

Do not move on to large benchmarking until the downsampled comparison is
understood.

## Artifact Locations

Safe local debug artifacts used in this investigation:

- fixed auto-detect validation:
  `tests/transcriptvb_auto_fix_output_20260309_235417/`
- prior broken integrated run:
  `/tmp/pe_bulk_benchmark_real_v4/downsampled/integrated/`
- explicit IU comparison run:
  `/tmp/pe_bulk_benchmark_real_v4/downsampled/integrated_IU/`

Restricted Y-related artifacts remain local-only. Do not copy them elsewhere.

## Working Tree State

Modified tracked files:

- `core/features/vbem/source/LibFormatDetection.cpp`
- `core/features/vbem/source/LibFormatDetection.h`
- `core/features/vbem/source/TranscriptQuantEC.cpp`
- `core/features/vbem/source/libem/alignment_filter.h`
- `core/features/vbem/source/libem/ec_builder.cpp`
- `core/legacy/source/ReadAlignChunk.cpp`
- `core/legacy/source/ReadAlignChunk.h`
- `core/legacy/source/ReadAlignChunk_processChunks.cpp`

Untracked items already present:

- `latex/`
- `scripts/paper/`
- `tests/transcriptvb/__pycache__/`

Do not clean or revert unrelated local changes.

## Checkout Warning

`AGENTS.local.md` says:

- `/mnt/pikachu/STAR-suite` is a dirty local-work checkout
- `/mnt/pikachu/STAR-suite-benchmark-pe` is the clean clone intended for PE
  benchmark work

This debugging work was done in `/mnt/pikachu/STAR-suite`. Before broader PE
benchmarking or paper scripts, the next agent should either:

1. port the validated code changes into `/mnt/pikachu/STAR-suite-benchmark-pe`,
   or
2. confirm with the user that continuing from this dirty checkout is acceptable

## Recommended Next Steps

1. Stay on the **downsampled PE** fixture only.
2. Compare integrated vs stepwise output after the STAR run and isolate the
   remaining divergence.
3. Add a targeted regression script for the fixed PE `A -> IU` autodetect case,
   preferably using the real downsampled PE fixture.
4. Keep all Y-derived per-read artifacts local and unshared.
