# STAR Suite v1.4.3 Release Notes

Date: 2026-07-01

`v1.4.3` is a patch release for byte-identical Y/no-Y FASTQ emission in
integrated bulk PE workflows. The release artifact version is `v1.4.3`, and
`STAR --version` reports `1.4.3`. The upstream STAR base remains available
through `STAR --upstream-version` (`2.7.11b`), and the genome index
compatibility string remains available through `STAR --genome-compat-version`
(`2.7.4a`).

## Fixes

- Preserve Illumina FASTQ header suffixes such as
  `1:N:0:ATGAGGCC+CAATTAAC` through the FastxInputModule path and the legacy
  FASTQ chunk path before Y/no-Y FASTQ re-emission.
- Emit Y/no-Y FASTQ sequence and quality payloads using the post-trim
  `readLength` value instead of treating `Read0`/`Qual0` as null-terminated
  strings. This prevents adapter/quality-trimmed tails left in the read buffer
  from appearing in sidecar FASTQs.
- Keep the historical `gzip`/`awk`/`samtools` Y/no-Y FASTQ artifact contract
  byte-equivalent to STAR-suite integrated Y/no-Y FASTQ emission.

## Validation

- Clean rebuilt STAR from `core/legacy/source`.
- Passed PPARG 100K paired-end byte-equivalence gate against the historical
  Trim Galore + awk/samtools/gzip path with `--runThreadN 1`.
- Passed the same PPARG 100K paired-end byte-equivalence gate with
  `--runThreadN 32` and Trim Galore `--cores 8`.
- In both gates, decompressed FASTQ MD5s matched exactly for `R1_Y`, `R2_Y`,
  `R1_noY`, and `R2_noY`.
