# STAR Suite v1.10.0 Release Notes

Date: 2026-09-29

STAR Suite 1.10.0 removes every dependency on other suites. STAR no longer
links Chromap Suite or RapidMACS, and no longer hosts the joint RNA + ATAC
(multiome) run. That integration moves to Multiomics Suite, which builds the
multiome binary from released STAR Suite, Chromap Suite and RapidMACS
versions. In its place STAR gains a small, generic host interface that knows
nothing about Chromap or ATAC. The release also adds paired hashtag
demultiplexing and per-library ambient-FDR feature calling.

`STAR --version` reports `1.10.0`. Debian source packaging uses `1.10.0-1`;
Ubuntu packages use `1.10.0-1~ubuntu22.04.1` and `1.10.0-1~ubuntu24.04.1`.
Upstream STAR remains `2.7.11b`, genome-index compatibility remains `2.7.4a`,
and legacy compatibility remains `2.7.1a`. Existing indexes do not need
rebuilding.

## Known limitations

- The bundled official recipe catalog remains the pinned 1.9.5 snapshot. Its
  multiome recipes have not migrated to the Multiomics Suite executable and are
  not supported standalone STAR 1.10 workflows. The release validator checks
  snapshot integrity, not runtime compatibility.
- The canonical downstream fixes validated with this release (seeded
  scDblFinder, fail-closed CellBender) are on `morphic-recipes` branch
  `fix/v1100-downstream-validation`, commit `bbf8b54`, and are not yet on its
  default branch.
- Multi-thread TranscriptVB quantification can vary between runs, as it does
  in 1.9.5.a. Cross-version TranscriptVB validation therefore used one thread.

## Removed: the Chromap integration

- **The multiome run moves to Multiomics Suite.** Removed from STAR: the
  Chromap orchestration (asynchronous ATAC worker, permit-telemetry sampler,
  drain-time and saturation controllers), the libchromap contract library and
  its runners, the inline peak/matrix step, libscrna's ATAC evidence reader
  and its five multiome cell-calling tools, the CBQ-to-Chromap FASTQ adapter
  with its harness and smoke test, the multiome tests and scripts, the
  `morphic_multiome` MCP workflow and the integration runbooks.
- **57 parameters are gone:** the 42 `--chromapAtac*`, the 11
  `--multiomeAtac*`, and the ATAC-only permit parameters
  `--dynamicThreadAtacFloor`, `--dynamicThreadAtacController`,
  `--dynamicThreadAtacWorkEstimate` and `--dynamicThreadTelemetryIntervalSec`.
  STAR now rejects them as unknown parameters. Multiomics Suite accepts the
  same names, types and defaults.
- **Nothing is lost.** Every removed file, code block, symbol and parameter is
  listed with its last STAR commit in
  [the hand-over document](https://github.com/morphic-bio/STAR-suite/blob/v1.10.0/docs/HANDOVER_MULTIOMICS_1.10.md),
  so that Multiomics Suite can copy it unchanged. The CBQ adapter moves rather
  than being deleted.
- **What stays in STAR:** the CAT-ATAC guide arm (a feature-barcode arm that
  reads its barcode in the ATAC whitelist namespace), the shared thread-permit
  pool and its saturation controller, and libscrna's EmptyDrops, OrdMag and
  occupancy code.

## Build

- **`make core` no longer needs Chromap Suite.** It builds STAR with its
  bundled HTSlib, which is what `make core-portable` built before;
  `core-portable` remains as an alias. `WITH_CHROMAP` and `CHROMAP_SUITE_DIR`
  are no longer used.
- **Published artifacts were already Chromap-free.** Release tarballs were
  built with `core-static`, Debian packages with `core-portable`, and the
  Docker image with `STAR_WITH_CHROMAP=0`. Only a local `make core` linked
  Chromap. Users who built that way and write BAM will see different
  compressed BAM bytes with identical records, because STAR's bundled HTSlib
  replaces the system HTSlib.
- **`HTSLIB=external`** compiles and links STAR against an installed HTSlib
  (`pkg-config htslib`, or `HTSLIB_CFLAGS`/`HTSLIB_LIBS`), for programs that
  embed STAR alongside other HTSlib users.
- **Reproducible library builds.** The PCG helper's unused build-time seed,
  derived from `__DATE__` and `__TIME__`, is now a fixed literal, and libem no
  longer adds `-march=native`: it uses the compiler's target baseline, and
  `CXXFLAGS` can still select a target. STAR's samplers keep their explicit
  seeds. A regression test in partial-build CI checks both. The unlinked older
  Flex copy of the PCG helper carries the same fix.

## Host interface (new)

A program can now link STAR Suite as a library and run its own work beside
STAR in the same process, sharing STAR's thread permits. See the
[host interface reference](https://github.com/morphic-bio/STAR-suite/blob/v1.10.0/docs/HOST_API.md).

- `make star-host-lib` builds `libstar_suite.a` and `libstar_suite.link`
  (the libraries to link after it). The headers are in
  `core/legacy/source/host/`. Host API version 1.
- `star::host::runMain(argc, argv, hooks)` runs STAR. The `STAR` executable is
  `runMain(argc, argv, nullptr)`, so standalone behaviour is unchanged.
- Parameters STAR does not know go to the host, which accepts or rejects them.
  STAR records accepted ones in `Log.out` and in the command lines it writes to
  BAM headers, including those given in `--parametersFiles`.
- Callbacks run before read mapping (`preflight`, `start`) and after all of
  STAR's own work (`finish`), where STAR used to call Chromap. Further
  callbacks set the extra permit threads and the initial floors, refuse routes
  that cannot share the pool, and require a full pool at exit.
- The third permit domain, formerly `ATAC`, is now a neutral `EXTERNAL`
  domain that STAR lends to the host. Standalone permit logs name it
  `external`; a host supplies its own label (Multiomics Suite uses `atac`).
  The `atacController=` field of the "Dynamic thread interface enabled" line is
  gone.
- The saturation permit controller is now a public header,
  `host/SaturationPermitController.h`, in namespace `star::permits`.

## Feature calling

- **Paired hashtag demultiplexing.** `star_hash_demux_method=pair` in a
  pf-multi library assigns each cell to a sample from a
  `sample, hash_a, hash_b` table (`star_hash_sample_table`), with
  `star_hash_min_pair_ratio` (default 2.0). Cells whose top two tags form an
  unknown pair are reported as `unknown_pair`. The standalone
  `assignBarcodes` accepts `--hash-demux-method pair`. The default `ratio`
  method and its output are unchanged.
- **Per-library ambient-FDR calls.** `star_feature_caller=ambient-fdr` now
  calls a named non-GEX library on its own, with `star_feature_call_fdr`
  (default 0.01) and `star_feature_call_min_umi` (default 1), within the
  EmptyDrops cell set. Output is written to
  `outs/feature_analysis/<library_id>/ambient_fdr/`. CellTag libraries are
  selected by `feature_type=CellTag` and are not mixed with guides. The
  automatic CRISPR path and its `guide_*` outputs are unchanged.

## Correctness and regression coverage

- **Inline CB correction with Velocyto, without BAM.** A memory-saving
  allocation guard incorrectly skipped the per-read CB/UMI storage requested
  by Velocyto's gene-like source. Honor that existing requirement, including
  native OCM composite barcodes. No cell-calling or velocity-counting rule is
  changed. A new synthetic regression checks exact GeneFull and velocity
  layers in plain, inline-corrected and OCM modes.
- Update the feature-permit smoke to validate batched acquisitions and
  concurrent GEX work instead of expecting one permit per record. Add nine
  validator unit tests, including invalid and missing telemetry.
- Explicitly select the legacy caller in the legacy Flex smoke. Its remaining
  hash-on/legacy parity failure remains a historical diagnostic, not the modern
  release gate. The release gate now uses the established H0/H1X2 half-probe
  fixture with matching model/filtered gene axes and preserved 1.9.5.a output.
- Version the standalone March decision replay with an explicit legacy negative
  policy, leaving the default policy and STAR runtime classification unchanged.
- Seed scDblFinder in the canonical downstream recipe's R container and match
  the STAR mirror. Fail on doublet errors rather than silently assigning singlets.
  Require CUDA and successful CellBender output; reject stale/partial failures.

## Validation

The validated source is `v1.10.0-rc2` (`d885bce`); the `v1.10.0` tag adds only
these release notes and the distribution entry. Comparisons are against
1.9.5.a built in the same container image (Ubuntu 22.04, g++ 12.3.0).

- **Build.** Clean build of the exact rc2 source; 11/11 partial Make targets
  pass, including the host library and host API tests. The reproducible-library
  regression passes. Official snapshot integrity passes (11 recipes, 10
  evidence records).
- **Regression (G-S1).** 25/25 production regression cases, 15/15 Tier A
  tests and 9/9 additional feature and build checks pass. All 323 retained
  output checks match 1.9.5.a; only source-directory provenance paths are
  normalized. Two initial test invocations failed because an external decoder
  wrapper passed a second thread argument; run with the test's own argument,
  both passed. STAR was not changed.
- **TranscriptVB, one thread.** Four 100K cases (treated and no-4sU, single-
  and paired-end) give 12 quantification tables; 11 are byte-identical to
  1.9.5.a. In the treated single-end `star_quant.genes.sf`, gene
  `ENSG00000099917` reports `Length` 641.878 against 641.877. Its read count,
  TPM and effective length are identical, as are every other field and the row
  order. The author accepted this rounding difference as an effect of the
  portable build: libem is no longer compiled with `-march=native`.
- **Performance (G-S3).** Three repetitions per version and workload, run
  under the shared host lock with alternating version order. All 18 counted
  attempts match the 1.9.5.a outputs and have clean host-load verdicts. Limits
  are +3% median wall time and +2% median peak RSS.

| Workload | 1.9.5.a median wall (s) | 1.10.0 median wall (s) | Wall change | Peak RSS change |
|---|---:|---:|---:|---:|
| scRNA-seq 100K regression | 56.83 | 56.11 | -1.27% | +0.001% |
| scRNA-seq 10M reads | 40.06 | 39.68 | -0.95% | -0.003% |
| Flex 100K half-probe | 15.50 | 15.46 | -0.26% | -0.092% |

In the first timing series, five of the six Flex attempts, which take about
15 s, were flagged because about 0.3 MB of disk I/O was attributed to other
processes; most of it is the filesystem journal recording the benchmark's own
output. With the author's approval those five were re-run until each had a
clean verdict (130 runs; an earlier series of 60 was set aside because
it ran under a different locale, which reorders the Flex test's checksum
listing). Every attempt is kept in the evidence; Flex re-run wall times ranged
from 15.2 to 15.6 s.

- **Multiomics Suite integration.** Multiomics Suite built its single binary
  against the rc2 host interface with its strict reproducible build. G-M2
  passed: six alternating 20M-read, 32-thread runs gave a median wall change of
  -0.46% and a median RSS change of +0.07% (limits 3% and 2%). In G-M1, three
  of five fixtures matched exactly. The author accepted the differences in the
  other two: the antibody `feature_per_cell.csv` files (raw and filtered)
  contain the same rows in a different order, and the CAT-ATAC ATAC BAM is
  2,647,891 rather than 2,648,549 bytes with byte-identical decompressed
  contents (145,114 records). G-M3, a repeated clean rebuild, was not run and
  does not block STAR.

Earlier candidate evidence, including the rc1 gates, is in the
[validation follow-up](VALIDATION_STAR_1_10_0_20260929.md), the
[runbook](runbooks/RUNBOOK_STAR_1_10_0_HOST_API_20260928.md) and the
[handoff](handoffs/HANDOFF_STAR_1_10_0_HOST_API_20260928.md). The medians
meet the regression limits; they are not a claim of a speedup.
