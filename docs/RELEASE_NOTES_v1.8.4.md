# STAR Suite v1.8.4 Release Notes

Date: 2026-09-05

`v1.8.4` adds a deliberately narrow BGZF input bridge for fused STAR-Flex
runs that request coordinate-sorted BAM output. Paired BGZF FASTQ can now use
the parallel range reader while feeding the existing STAR mapping and BAM
writing machinery. Other STAR and STARsolo input paths retain their previous
behavior.

The release artifact version is `v1.8.4`, the Debian source package uses
`1.8.4-1`, and `STAR --version` reports `1.8.4`. GitHub Releases provides
binary packages built and validated on Ubuntu 22.04 and Ubuntu 24.04 as
`1.8.4-1~ubuntu22.04.1` and `1.8.4-1~ubuntu24.04.1`. The upstream STAR base
remains `2.7.11b`; genome-index compatibility remains `2.7.4a`. Existing
indexes do not need to be rebuilt.

## Flex coordinate-sorted BAM input

- Feed leased BGZF records into the established STAR mapping chunks and
  coordinate-sorted BAM writer without an intermediate decompression file.
- Discover independent BGZF members from their inline `BC/BSIZE` fields. No
  pre-scan, external index, or sidecar is required.
- Preserve exact mate ownership by record ordinal and retain the Illumina
  pass-filter flag used by STAR's established read-filter behavior.
- Drain one lane at a time and distribute reader workers between mates using
  their compressed sizes.
- Add `--out-samtype bam-sorted` to the Flex recipe wrapper as an explicit
  opt-in; its existing no-BAM default is unchanged.

## Validation

- BGZF `off` and `range` modes produced byte-identical coordinate-sorted BAM
  bodies and byte-identical `Barcodes.stats` on the synthetic Flex gate.
- The resulting BAM passed `samtools quickcheck`.
- Negative gates confirmed that plain gzip, unsorted BAM, and classic
  non-Flex STAR do not select the new bridge.
- The complete BGZF T1-T10 suite passed, including fused alignment, pairing,
  truncation, mixed-lane, gzip, and CBQ regressions.
- A separate Flex CBQ smoke produced byte-identical BAM bodies and barcode
  statistics relative to FASTQ input.

## Compatibility boundaries

- The bridge is selected only for paired Fastx input when Flex inline hashing
  and coordinate-sorted BAM output are active and every input mate is BGZF.
- Plain gzip continues through the established serial zlib reader.
- Classic STAR, STARsolo, bulk, non-Flex, unsorted BAM, two-pass and splice-
  junction staging, transcriptome BAM, spatial R1 tap, and ordered paired-input
  modes remain on their established paths.
- No htslib code is changed, no external dependency is added, and BGZF CRC and
  reader-mode defaults are unchanged.
