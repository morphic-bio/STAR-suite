# STAR Suite v1.7.0 Release Notes

Date: 2026-08-08

`v1.7.0` is a feature release for distributed bulk-RNA processing, faster
reference and TranscriptVB setup, the RapidMACS naming transition, and a
long-standing STAR genome-generation memory-safety fix.

The release artifact version is `v1.7.0`, Debian packages use `1.7.0-1`, and
`STAR --version` reports `1.7.0`. The upstream STAR base remains `2.7.11b`;
genome-index compatibility remains `2.7.4a`. Official binaries embed the exact
40-character release commit reported by `STAR --source-revision`.

## TranscriptVB scatter/gather

- Add `--quantVBSidecarOut`, `--quantVBSidecarOnly`, and shard identity fields
  for writing compact TranscriptVB evidence from independent alignment shards.
- Add the packaged `transcriptvb_finalize` executable to validate and gather
  sidecars into `quant.sf` and gene-level abundance without loading Genome, SA,
  or SAindex.
- Preserve in-process quantification semantics by normalizing shard EC weights
  and reconciling transcript-table order with transcriptome-FASTA order before
  finalization.
- Support range emission, range workers, deterministic assembly, sparse
  gating, and an opt-in parallel GC-background pass. The no-GC path remains the
  byte-parity reference for scatter/gather validation.

## Distributed trim QC

- Feed STAR's trim-QC collector from actual reads; previous reports could be
  structurally valid while incorrectly reporting zero reads.
- Add `--trimQcShardOut` for portable additive per-shard counter dumps.
- Package `trim_qc_merge` to combine shard dumps deterministically into the
  same JSON and HTML histograms produced by a single run.
- Package `trim_qc_fastq` for direct FASTQ/FASTQ.gz collection and compatible
  shard output, avoiding BAM reconstruction in distributed bulk pipelines.

## Parallel setup and I/O

- Add configurable concurrent genome-index loading through
  `--genomeLoadThreads`.
- Parallelize the optional TranscriptVB GC-background pass while retaining a
  serial parity mode.
- Merge trim-QC summaries across alignment shards instead of requiring a
  monolithic BAM or read replay.

## RapidMACS and integrated ATAC naming

- Rename the in-process MACS-compatible peak-calling dependency to RapidMACS,
  its library to `librapidmacs`, and the executable to `rapidmacs` throughout
  the STAR/Chromap integration contract.
- Preserve the existing Chromap multiome feature behavior and defaults. Hosted
  STAR Suite tarballs, Debian packages, and Docker builds remain portable
  no-Chromap artifacts because the integrated path depends on the sibling
  Chromap Suite and RapidMACS source trees.

## Genome-generation memory safety

- Fix `genomeGenerate --sjdbGTFfile` when the annotation is valid but produces
  zero splice junctions, such as a single-exon synthetic or real annotation.
- The apparent teardown error (`double free or corruption (!prev)`) was heap
  corruption from an earlier packed suffix-array write at record `nSA` after
  junction insertion had returned without rebuilding the array.
- Skip that padding write unless junctions were actually inserted. Real
  nonzero-junction insertion retains the established padding behavior.

## Release artifact coverage

- Ship `transcriptvb_finalize`, `trim_qc_fastq`, and `trim_qc_merge` in both
  Linux compatibility tarballs, the automatic installer bundle, Debian
  packages, and the multi-architecture Docker image.
- Validate companion-tool inventory during tarball, installer, and Debian
  install/uninstall checks.
- Run installed-binary synthetic gates for zero-junction GTF generation,
  adapter trimming, TranscriptVB scatter/gather parity, standalone trim-QC
  merging, and STAR in-process versus sharded trim-QC identity.
- Continue validating Ubuntu 22.04 (`glibc234`) and Ubuntu 24.04 (`glibc239`)
  compatibility assets, checksums, exact source identity, and installer
  selection.

## Compatibility

- Existing genome indexes remain compatible; no index rebuild is required.
- Existing STAR, STARsolo, Flex, SLAM, CR-compatible, and spatial defaults are
  unchanged unless a new scatter/gather, trim-QC, parallel-GC, or parallel-load
  option is selected.
- The integrated RapidMACS path remains an opt-in source-build integration;
  portable official STAR Suite binaries continue to use `WITH_CHROMAP=0`.

Detailed methodology and contracts are documented in
`docs/STANDALONE_VBEM_FINALIZER_20260807.md`,
`docs/TRANSCRIPTVB_SHARD_SETUP_STAGE_20260807.md`,
`docs/TRIM_QC_SHARD_MERGE_20260807.md`, and
`docs/LIBCHROMAP_CONTRACT.md`.
