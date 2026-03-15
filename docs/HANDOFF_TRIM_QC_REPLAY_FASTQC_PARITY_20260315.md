# Trim QC Replay and FastQC Parity

## Summary

This change adds a standalone alignment replay tool for the existing STAR
FastQC-like trim-QC output path:

- `core/legacy/source/trim_qc_replay.cpp`
- `scripts/run_trim_qc_replay.sh`
- `tests/run_trim_qc_replay_smoke.sh`

The tool reads SAM/BAM/CRAM, reconstructs sequence and quality strings from the
alignment records, and writes the same `*.trim_qc.json` / `*.trim_qc.html`
schema used by the in-STAR `--trimQcReport` path.

## Why this was added

Two use cases motivated it:

1. Post hoc QC on existing alignment files without rerunning STAR.
2. A cheap regression surface for trim-QC output changes, independent of a full
   alignment run.

## FastQC comparison

The replay tool was compared against `FastQC` on:

- `tests/cbub_regress_unsorted/Aligned.out.bam`

Key finding:

- The replay tool can match FastQC closely on BAM input, but only when the BAM
  replay semantics are aligned with FastQC.

Specifically:

1. FastQC includes all BAM records by default, including secondary alignments.
2. FastQC reverse-complements negative-strand BAM records back to original read
   orientation before computing per-base content and quality.

With `trim_qc_replay --all-records`, the comparison matched on:

- total sequences
- sequence length
- GC content
- first-base composition
- first-base mean quality

Scratch comparison artifacts:

- `/tmp/fastqc_bam_compare_V9OdQu/Aligned.out_fastqc.html`
- `/tmp/fastqc_bam_compare_V9OdQu/replay_all.trim_qc.html`
- `/tmp/fastqc_bam_compare_V9OdQu/replay_all.trim_qc.json`

## Root cause of the parity bug

An initial replay implementation bug caused per-base composition and quality to
diverge from FastQC on reverse-strand BAM records.

Root cause:

- `trim_qc_replay.cpp` incorrectly treated `bam_seqi()` values as one-hot A/C/G/T
  masks.
- HTSlib actually uses BAM's nt16 alphabet encoding.
- That made reverse-complement handling incorrect even though the total-count
  surface was otherwise reasonable.

Fix:

- decode BAM bases via `seq_nt16_str`
- then reverse-complement negative-strand records before collecting QC

After this fix, BAM replay matched FastQC on the tested BAM.

## Controlled fixture behavior

The committed smoke fixture is:

- `tests/fixtures/trim_qc_replay_tiny.sam`

Important semantic difference demonstrated by the tiny fixture:

- FastQC on SAM/BAM counts all records unless filtered for another reason.
- The replay tool defaults to primary-only behavior for safer post hoc analysis.
- Use `--all-records` when the goal is explicit FastQC parity on alignment files.

## Commands

Build:

```bash
make -C core/legacy/source trim_qc_replay
```

Replay:

```bash
scripts/run_trim_qc_replay.sh \
  --input tests/cbub_regress_unsorted/Aligned.out.bam \
  --report /tmp/replay \
  --stage bam_replay \
  --all-records
```

Smoke:

```bash
tests/run_trim_qc_replay_smoke.sh
```

## Related files

- `core/legacy/source/STAR.cpp`
- `core/legacy/source/parametersDefault`
- `scripts/lib/star_trim_qc.sh`
- `scripts/run_ucsf_perturb_yremove_batch.sh`
- `scripts/run_ucsf_full_compat_forward_rescue_guides.sh`
