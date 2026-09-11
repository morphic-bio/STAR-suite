# STAR Suite v1.9.2 Release Notes

Date: 2026-09-11

STAR Suite 1.9.2 corrects scRNA-seq cell calling. In 1.9.1 a single UMI floor
decided both which barcodes EmptyDrops tests and which OrdMag primary cells are
kept, and the scRNA-seq path set that floor to 100. EmptyDrops therefore tested
barcodes down to 100 UMIs instead of the declared 500, and on low-depth
libraries every primary cell below 100 UMIs was removed. 1.9.2 separates the two
floors. It includes the merged changes through
`fa93df77d79c3e8c20dcdd09240367e04abc4c44`.

`STAR --version` reports `1.9.2`. Debian source packaging uses `1.9.2-1`;
Ubuntu packages use `1.9.2-1~ubuntu22.04.1` and
`1.9.2-1~ubuntu24.04.1`. Upstream STAR remains `2.7.11b`, genome-index
compatibility remains `2.7.4a`, and legacy compatibility remains `2.7.1a`.
Existing indexes do not need rebuilding.

## Cell calling

- **scRNA-seq (EmptyDrops_CR, non-Flex).** EmptyDrops candidates now start at
  the declared `umiMin`, the sixth `EmptyDrops_CR` value (500 by default).
  OrdMag primary cells are no longer trimmed by it.
- **OCM.** Candidates use the declared `umiMin` (500 by default). The primary
  trim stays at 100.
- **Flex.** Unchanged. Flex outputs are byte-identical to 1.9.1.
- libscrna has a separate OrdMag primary floor (`primary_umi_min`). Zero, the
  default, keeps the 1.9.1 behaviour of trimming primaries at `umi_min`, which
  is what Flex, the tag-aware caller and `scrna_simpleed` use. Both floors are
  written to the log and to the OrdMag summary.
- scRNA-seq and OCM cell calls change relative to 1.9.1.

## soloCellFiltering

- `--runMode soloCellFiltering` with `EmptyDrops_CR` exited in 1.9.1 with
  "MT rank annotation does not match this feature matrix". The `-` default of
  `--soloCellFilterMitochondrialGenes` was cleared only after that mode had
  returned. It is now cleared first. Alignment runs are unaffected.

## Benchmark tooling

- `scripts/paper/run_pe_bulk_feature_benchmark.sh` again runs Trim Galore with
  FastQC as the external trimmer, awk/samtools Y removal after alignment, and
  `--trimQcReport` on the integrated arm. These options were used for the June
  2026 bulk controls but had never been committed. The default trim-QC limit is
  250,000 reads per alignment-thread collector (`--trim-qc-max-reads`).

## Validation

Cells were re-called from saved raw matrices with `--runMode soloCellFiltering`
at 32 threads. The comparison build is 1.9.1 with only the soloCellFiltering
fix; it reproduces the in-run 1.9.1 PBMC callset exactly. These are correctness
checks, not timed benchmarks.

| Dataset | 1.9.1 | 1.9.2 |
| --- | --- | --- |
| PBMC 10K v3 vs Cell Ranger 9.0.1 (11,806 cells) | 12,146 cells, Jaccard 0.971358 | 11,863 cells, Jaccard 0.995027 |
| UCSF iPSC2, low depth, vs Cell Ranger (7,325 cells) | 6,946 cells, Jaccard 0.948259 | 7,275 cells, Jaccard 0.989101 |
| Flex GSE268426, GSE325982 and JAX CBQ reruns | reference | byte-identical: 19, 52 and 151 output files |

- PBMC 10K: the 308 STAR-only cells between 100 and 499 UMIs are removed. The
  remaining 58 STAR-only cells and 1 Cell Ranger-only cell all have at least
  500 UMIs.
- iPSC2: Cell Ranger calls 379 cells below 100 UMIs, and 1.9.1 missed all of
  them. 1.9.2 misses 65 and adds 15, all below 100 UMIs, matching the result
  from before the shared floor was introduced (Jaccard 0.989786).
- OCM: the tiny-fixture test passes. With the primary trim at 100, a native
  re-call differs from 1.9.1 only in the recorded candidate floor.
- `tests/test_ordmag_primary_floor.cpp` checks both floors on synthetic
  counts (10 of 10 pass).

The PBMC, iPSC2 and Flex checks used the build before the OCM primary floor
was set back to 100; that setting affects only OCM, which was re-tested with
it. No Cell Ranger source was inspected for these changes. New paper benchmark
executions should pin `v1.9.2` and record its resolved source commit and
binary hash.

## Distribution

Release packaging provides amd64 glibc234/glibc239 tarballs, a compatibility
installer bundle, Ubuntu 22.04/24.04 Debian binaries, Debian source packaging,
runtime manifests and SHA-256 checksums. Hosted tarball/Debian builds retain
the established portable no-Chromap configuration; local production source
builds retain the Chromap-enabled default. `STAR --source-revision` reports
the exact tagged source embedded in each official binary.
