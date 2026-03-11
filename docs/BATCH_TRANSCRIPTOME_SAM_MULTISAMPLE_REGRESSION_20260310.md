# Batch TranscriptomeSAM Multisample Regression

## Summary

On March 10, 2026, a new PE multisample smoke exposed a batch-mode regression in
the bulk quantification path:

- `--batchMode 1`
- PE input
- internal adapter trimming (`--trimCutadapt Yes`)
- Y removal (`--emitNoYBAM yes --emitYNoYFastq yes`)
- `--quantMode TranscriptomeSAM TranscriptVB`

Sample 1 completed, but sample 2 crashed when STAR tried to write the
transcriptome BAM header.

## Root Cause

Batch reset was reopening the main per-sample BAM outputs, but it was not
rebinding and reopening the transcriptome BAM writer (`outQuantBAMfile`) for
the next sample.

The second sample therefore entered `samHeaders()` with a null BGZF handle for
`Aligned.toTranscriptome.out.bam`.

Failure signature:

- crash in `bgzf_write(fp=0x0, ...)`
- call chain through `outBAMwriteHeader()` and `samHeaders()`

Fix location:

- [Parameters_batchReset.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/Parameters_batchReset.cpp)

The fix updates the per-sample transcriptome BAM path during batch reset and
reopens the BGZF writer before the next sample starts.

## Regression Coverage

Added smoke:

- [run_bulk_pe_multisample_equivalence.sh](/mnt/pikachu/STAR-suite/tests/run_bulk_pe_multisample_equivalence.sh)

Helper:

- [split_fastq_consecutive_pe.py](/mnt/pikachu/STAR-suite/tests/split_fastq_consecutive_pe.py)

What the smoke checks:

- one batch STAR run over two consecutive pseudo-samples from the PE downsample
- internal adapter trimming
- Y/noY BAM emission
- Y/noY FASTQ emission
- `TranscriptVB` autodetect resolves `IU`
- batch vs sequential outputs match for:
  - `quant.sf`
  - `quant.genes.sf`
  - coordinate BAM
  - transcriptome BAM
  - Y/noY BAMs
  - Y/noY FASTQs
  - `SJ.out.tab`

The smoke also requires both pseudo-samples to contain at least one emitted
`Y` read, so it exercises the Y split rather than only the batch reset path.

## Determinism Note

The checked-in smoke defaults to:

- `BULK_MULTISAMPLE_THREADS=1`

Reason:

- the batch reset bug itself is independent of thread count
- exact byte-for-byte `TranscriptVB` parity is reproducible at one thread
- multithreaded VBEM introduces small run-to-run floating-point drift, which is
  not useful as a regression signal for this bug

This smoke is therefore designed to be a strict correctness guard, not a
throughput benchmark.

## Current Status

This issue is currently considered closed for active benchmark follow-up.

- the regression is fixed
- the strict batch-vs-sequential smoke remains checked in
- future broader PE multisample benchmarking can be revisited later without
  blocking current paper/debugging work

## Validation Artifact

Clean passing run:

- [bulk_pe_multisample_equivalence_output_final_20260310_075405](/mnt/pikachu/STAR-suite/tests/bulk_pe_multisample_equivalence_output_final_20260310_075405)

Key outputs:

- [SUMMARY.txt](/mnt/pikachu/STAR-suite/tests/bulk_pe_multisample_equivalence_output_final_20260310_075405/SUMMARY.txt)
- [comparison.tsv](/mnt/pikachu/STAR-suite/tests/bulk_pe_multisample_equivalence_output_final_20260310_075405/comparison.tsv)
