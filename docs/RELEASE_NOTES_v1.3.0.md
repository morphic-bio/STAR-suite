# STAR Suite v1.3.0 Release Notes

Date: 2026-06-01

`v1.3.0` is a feature release for parallel CBQ integration. The release
artifact version is `1.3.0`, and `STAR --version` reports `1.3.0`. The upstream
STAR base remains available through `STAR --upstream-version` (`2.7.11b`), and
the genome index compatibility string remains available through
`STAR --genome-compat-version` (`2.7.4a`).

## Feature Scope

- Enable native CBQ range-reader paths for STAR core and Flex workflows.
- Extend CBQ integration through process_features and the native libchromap
  multiome ATAC input path.
- Support PBMC multiome GEX+ATAC CBQ production input without gzip ingestion in
  the STAR/Chromap run itself.
- Keep ATAC Y/noY BAM output available on the native libchromap path.

## PBMC Multiome CBQ Benchmark

- Dataset: PBMC 3k 10x Multiome full fixture.
- Native CBQ input: GEX `--readFilesType Binseq PE` plus ATAC
  `--chromapAtacInputFormat cbq`.
- STAR/Chromap wall time: `18:46.41`, mapping speed `582.23M reads/hr`, peak
  RSS `51,048,592 KB`.
- Completed FASTQ.gz comparator used external `zcat`: `20:40.05`, mapping
  speed `520.89M reads/hr`, peak RSS `46,267,956 KB`.
- CBQ saved `1:53.64` in the shared STAR/Chromap phase (`9.2%`, `1.10x`) while
  also emitting ATAC Y/noY BAM outputs.
- Sidecar peak/MEX materialization from the CBQ run took `1:27.53`; wrapper
  total was `20:13.94`.
- Validation: GEX BAM, ATAC BAM, ATAC noY BAM, and ATAC Y BAM pass
  `samtools quickcheck`; noY has zero `chrY` alignments, Y has zero non-Y
  alignments; peak count is `50,274`.

## Validation

- Clean rebuilt Chromap-suite `libchromap.a` and STAR Suite core.
- Passed STAR libchromap CBQ contract smoke.
- Passed Chromap CBQ adapter smoke.
- Passed PBMC 100K multiome FASTQ smoke.
- Passed full PBMC multiome native CBQ smoke with BAM writing, ATAC Y/noY BAM
  output, sidecar peak/MEX materialization, and output inventory validation.
