# STAR Suite

STAR Suite updates the original STAR aligner by integrating four modules — STAR-core (with integrated TranscriptVB quantification), STAR-perturb, STAR-Flex, and STAR-SLAM — to provide complete internal C/C++ pipelines for bulk RNA-seq, scRNA-seq, Perturb-seq, 10x Flex, and SLAM-seq. The integration results in **substantial speedups** (STAR Suite 1.9.4: **1.6–4.2x for bulk RNA-seq** against an external stepwise pipeline; against Cell Ranger 9.0.1, **3.9x for scRNA-seq**, **4.8–5.8x for Perturb-seq**, and **17–23x for Flex** from the delivered FASTQ, up to **42x** from CBQ input) and a simplified toolchain that can be **installed through pre-compiled binaries** for researchers and agents. **No new external dependencies** are required; the suite is built entirely with the existing STAR toolchain and vendored components. **This is a drop-in replacement for the STAR aligner.**

Current production release: **[STAR Suite v1.9.5](https://github.com/morphic-bio/STAR-suite/releases/tag/v1.9.5)**. The suite release tag and
source-packaging version are `v1.9.5` / `1.9.5-1`; `STAR --version` reports
`1.9.5`. The benchmark results below retain their measured 1.9.4 identity.
GitHub Releases also provides Ubuntu 22.04- and 24.04-built `.deb`
packages from the same source revision.
Use `STAR --upstream-version` for the underlying upstream STAR base
(`2.7.11b`) and `STAR --genome-compat-version` for the genome index
compatibility string (`2.7.4a`). Use `STAR --source-revision` to report the
exact full commit embedded in an official binary.

Advanced-user previews live on `dev-release` or version-scoped
`dev-release-vX.Y.Z` branches. Immutable prerelease tags use `vX.Y.Z-rcN`;
stable releases are tagged from `master`.

STAR Suite supports partial compilation: build only the module/tool targets you need instead of building the full suite every time.

Agent quickstart: see `AGENTS.md` for repo-specific guardrails, tests, and recent changes.

## Core Additions over STAR 2.7.11b

- **Speedup** (STAR Suite 1.9.4; see [Benchmarks](#benchmarks)): bulk RNA-seq **1.6x** faster than an external Trim Galore + STAR + Salmon pipeline and **4.2x** with Y-chromosome removal; scRNA-seq **3.9x** faster than Cell Ranger 9.0.1 on the public 10x PBMC 10K library; Perturb-seq **4.8x** (A375) and **5.8x** (MSK 30-KO, three feature libraries in one pass against two Cell Ranger runs); Flex **16.8–22.8x** from the delivered FASTQ and **20.2–41.8x** from CBQ input. Gene-level Spearman and Pearson **0.98–1.0** against the reference outputs.
- **Batch Mode** (`--batchMode 1`): Processes multiple FASTQs in one STAR invocation while reusing the loaded genome. Removes the need for `--genomeLoad` keep-in-memory workflows. Single-pass only (no `--twopassMode`); not supported with Solo (`--soloType`). Use `--outFileNamePrefixAuto 1` for per-sample subdirectories.
- **TranscriptVB Quantification** (`--quantMode TranscriptVB`): Variational Bayes and EM quantification for transcript-level abundance, with parity-oriented behavior against Salmon alignment-mode. Gene-level summarization via `--quantVBgenesMode Tximport`.
- **TranscriptVB Scatter/Gather** (`--quantVBSidecarOnly 1`): Writes compact,
  mergeable per-shard evidence; the packaged `transcriptvb_finalize` executable
  gathers shards into transcript- and gene-level abundance without loading the
  genome index.
- **Sharded Trim QC** (`--trimQcShardOut`): Records additive per-shard QC
  counters during alignment. Packaged `trim_qc_merge` and `trim_qc_fastq`
  executables merge distributed counters or collect the same report directly
  from FASTQ without reconstructing BAM files.
- **Parallel Reference I/O** (`--genomeLoadThreads`): Loads genome-index files
  concurrently, while TranscriptVB can opt into a parallel GC-background pass.
- **Transcriptome Output** (`--quantTranscriptomeSAMoutput`): Replaces the former `--quantTranscriptomeBan` with more explicit control (e.g., `BanSingleEnd_ExtendSoftclip`).
- **Reference Automation** (`--autoIndex Yes`): Automated reference download/build with `--cellrangerStyleIndex Yes` formatting and `--genomeGenerateTranscriptome Yes` for transcript-level quant workflows.
- **Native Gzip FASTQ Handling**: Automatic detection of `.gz` FASTQ inputs with internal zlib streaming — no `--readFilesCommand zcat` needed for correctness. FLEX FASTQ production recipes use this path by default; legacy external helper mode remains available via `--readFilesLegacyZcat Yes`.
- **CBQ/BINSEQ Input** (`--readFilesType Binseq PE|SE`): Native C++ CBQ reader plus an order-preserving FASTQ/FASTQ.gz-to-CBQ encoder for STAR mapper, STARsolo, OCM, Flex, SLAM, and process_features adapter workflows. Exact FASTQ-vs-CBQ parity smokes are registered in the production regression manifest. See [`docs/CBQ_FORMAT_AND_IMPLEMENTATION.md`](docs/CBQ_FORMAT_AND_IMPLEMENTATION.md) for the format and adapter reference. Multiome ATAC now also supports the native libchromap CBQ path with paired-read and barcode CBQs.
- **Cutadapt-Compatible Trimming** (`--trimCutadapt Yes`): Native cutadapt-style trimming for bulk/PE workflows. Compatibility mode: `--trimCutadaptCompat Cutadapt3`.
- **Poly-G Trimming** (`--clip3pPolyG yes|no|auto`): Trims poly-G artifacts common on NovaSeq/NextSeq platforms. Default `auto` activates in CellRanger4 mode. Without this, poly-G reads can inflate specific genes (e.g., LINC00486) and degrade gene-level correlations.
- **Samtools-style BAM Sorting** (`--outBAMsortMethod samtools`): Spill-to-disk sort to reduce peak RAM pressure. Works with all modes including Flex.
- **Y/NoY Separation** (`--emitNoYBAM yes`, `--emitYNoYFastq yes`): Split BAM and FASTQ outputs by chrY alignment. Works with bulk and single-cell runs; for Flex it needs a BAM, which is a legacy route (`--flexLegacy yes`).
- **EmptyDrops_CR Integration**: CR-compatible EmptyDrops path (including libscrna-backed behavior in scRNA/perturb flows).
- **Molecule-first Barcode Assignment**: The opt-in
  `molecule_first_resolver` preserves finite read candidates through
  deterministic read-clique formation and candidate-specific UMI collapse,
  then emits strict, soft expected-count, hard, and gated-hard products before
  cell or spatial-bin calling. Existing STARsolo and Flex decoding remains the
  default compatibility path.
- **Integrated Visium HD 3' GEX** (`--soloSpatialGexIntegrated yes`): Fuses R1
  spatial decoding with post-alignment GeneFull evidence, retains ambiguous
  coordinate candidates through molecule-first resolution, and writes 2, 8,
  and 16 micrometer strict/soft/hard/gated-hard matrices. Bounded downstream
  spooling supports complete slides on smaller-memory hosts. The feature is
  opt-in; ordinary bulk, STARsolo, CR-compatible scRNA, and Flex runs keep
  their existing paths and defaults.
- **Native Visium HD Flex** (`--soloSpatialFlexIntegrated yes`): From 1.9.5,
  resolves probes with the shared half-probe classifier and retains every
  spatial candidate coordinate through molecule resolution. The fused FASTQ
  path loads no genome and emits all four policies at 2, 8, and 16 micrometers
  through the existing memory or bounded spill engine. Ordinary gzip and
  native BGZF are supported. The retired spatial alignment route is rejected;
  see [spatial half-probe processing](docs/SPATIAL_FLEX_HALF_PROBE.md).
- **Flex gDNA QC** (`--soloFlexGdna auto|yes|no`): Computes the 10x-style gDNA
  diagnostic from final filtered barcode/gene/UMI families and writes
  per-sample JSON plus an audit TSV. The diagnostic is matrix-inert and inert
  outside Flex.
- **Solo Features**: `sF` BAM tag for feature type, `--soloCBtype String` for arbitrary barcode strings, `--soloCellReadStats Standard` for improved cell filtering.
- **CR-compat GEX** (`--soloCrGexFeature auto|gene|genefull`): Controls which GEX source is merged in CR-compat mode.
- **Native Velocyto MEX Packaging**: Current production binaries write raw and
  filtered Velocyto MEX under `outs/` internally. `prepare_velocyto_mex.py` is a
  legacy repair/backfill helper for old STAR outputs, not the normal production
  path.
- **OCM Composite Barcode Mode** (`--ocmMultiEnable yes`,
  `--ocmMultiBarcodeMode flex`): 10x OCM runs can promote the barcode to an
  effective `CB16+OCM_TAG8` before correction, UMI collapse, and Velocyto,
  run per-sample CR-compatible EmptyDrops after the OCM split, then emit Cell
  Ranger multi-compatible `outs/multi`,
  `outs/per_sample_outs`, and per-sample downstream mirrors.
- **CB/UB Tag Pairing** (`--soloCbUbRequireTogether yes|no`): Enforce CB/UB tag pairing for tag injection (default `yes`).

## Folder Structure

```
core/
  legacy/                        # Canonical STAR core tree (upstream layout preserved; not deprecated)
  features/                      # Shared overlays and feature tooling
    process_features/            # Perturb feature extraction/calling implementation
    feature_barcodes/            # Standalone barcode tools (assignBarcodes, demux)
    libscrna/                    # EmptyDrops/OrdMag/Occupancy shared library
flex/                    # Flex-specific code + tools
slam/                    # SLAM-seq code + tools
build/                   # Modular make fragments
scripts/                 # Suite-level helper scripts (see scripts/README.md)
docs/                    # Suite-level docs
tests/                   # Suite-level tests (see tests/ARTIFACTS.md for artifact locations)
tools/                   # Suite-level scripts/utilities
mcp_server/              # MCP server for scripted discovery/preflight/run workflows
```

## Modules

- **STAR-core** (`core/`): Legacy STAR (indexing, bulk, Solo) plus shared utilities.
  Build: `make core` for the Chromap-enabled multiome-capable binary at
  `core/legacy/source/STAR`; use `make core-portable` for an explicit
  no-Chromap compatibility build.
- **STAR-perturb** (`core/legacy/` + `core/features/process_features/`): CR-compatible perturb-seq path with integrated feature extraction/calling (`process_features` + `call_features`) and `crispr_analysis/` outputs in CR-compat mode.
  Primary run path: `STAR --pfMultiConfig ... --defaultCrCompat yes` (see STAR-perturb section below).
- **STAR-OCM scRNA-seq** (`core/legacy/`): GEM-X OCM support on the CR-compatible
  GEX path. Production uses `GeneFull Velocyto`, per-sample `EmptyDrops_CR`
  after OCM split, native `CB16+OCM_TAG8` effective barcodes, native per-sample
  MEX/Velocyto materialization, and optional Y/noY side outputs.
- **STAR-Flex** (`flex/`): FlexFilter pipeline and Flex-specific integrations.
  Build tools: `make flex` or `make flex-tools`.
- **STAR-SLAM** (`slam/`): SLAM-seq quantification, SNP masking, trimming/QC.
  Build tools: `make slam` or `make slam-tools`.
- **Feature Barcodes** (`core/features/feature_barcodes/`): Standalone barcode tools (`assignBarcodes`, `demux_bam`, `demux_fastq`) for perturb-seq testing.
  Build tools: `make feature-barcodes-tools`.
- **Process Features** (`core/features/process_features/`): Full feature extraction/calling pipeline (`assignBarcodes`, `call_features`, `demux_bam`, `demux_fastq`) and standalone tool (`star_feature_call`).
  Build tools: `make process-features-tools`, `make star-feature-call`.
- **Shared Feature Toolchains** (`core/features/`): Reusable tool layers used across modules, including `vbem` (TranscriptVB helpers), `yremove_*` (Y/noY splitting), `bamsort`, and `libscrna`.
  Build tools: `make vbem-tools`, `make yremove-tools`, plus in-core integrations.
- **MCP Server (tooling)** (`mcp_server/`): Agent automation service for dataset/test discovery and controlled execution (`list_datasets`, `list_test_suites`, `preflight`, `run_script`, `collect_outputs`), plus **STAR Launchpad** (`/launchpad/`) in the browser for workflow recipes (defaults to **`star_*`** CLI recipes; optional full list). See [`mcp_server/README.md`](mcp_server/README.md).
- **Helper Scripts** (`scripts/`): Standalone Python and Bash tools for FASTQ preflight, QC, parity benchmarking, downstream h5ad processing, and fixture management. These are not compiled into STAR; they run independently. Highlights include `preflight_library_pairing.py` (chemistry detection and library pairing for mislabeled Perturb-seq), `report_additional_parity_metrics.py` (STAR vs CR parity), and `build_gene_full_velocyto_h5ad.py` (Velocyto h5ad packaging). See [`scripts/README.md`](scripts/README.md) for the full catalogue.

## Benchmarks

These are the STAR Suite **1.9.4** results reported in the manuscript. The
protocol, comparator commands and parity definitions are in
[docs/PAPER_BENCHMARK_METHODOLOGY.md](docs/PAPER_BENCHMARK_METHODOLOGY.md); the
benchmark entry points are listed in
[publications/benchmarks/README.md](publications/benchmarks/README.md). Records
from earlier releases remain under `comparisons/` and `docs/` for history and are
not the current numbers.

**Protocol.** Every run used 32 threads with every input, output and scratch
file on local NVMe, one run at a time on a quiet machine with the page cache
dropped first; wall time is the whole process from `/usr/bin/time -v`. All rows
were measured on one server (Intel i9-13900KF, 24 cores / 32 threads, 126 GiB
RAM) except the 320k Flex rows, which ran on a rented AWS instance limited to
the same size (m6id.16xlarge at 32 vCPUs and 126 GiB, local NVMe) because that
Cell Ranger run needs about 3 TB of scratch. The reference is Cell Ranger 9.0.1
with BAM output and secondary analysis disabled (`--create-bam=false`,
`--nosecondary`), timed to completion; the MSK 30-KO reference is the sum of the
two Cell Ranger runs needed for the same outputs. The bulk RNA-seq reference is
an external stepwise pipeline: Trim Galore 0.6.10 with FastQC, upstream STAR
2.7.11b (`--quantMode TranscriptomeSAM`), and Salmon 1.10.3 in alignment mode
with automatic library-type detection (`-l A --gcBias`); the with-Y-removal
reference adds an `awk`/`samtools`/`gzip` Y split. CBQ rows read binary input
converted once from the FASTQ; the conversion is excluded.

| Assay | Benchmark | STAR Suite 1.9.4 | Reference | Speedup |
|---|---|---:|---:|---:|
| Bulk RNA-seq | PPARG 35.1M PE, no Y-removal | 8 min 2 s | 12 min 49 s | 1.6x |
| Bulk RNA-seq | PPARG 35.1M PE, with Y-removal | 9 min 40 s | 40 min 47 s | 4.2x |
| scRNA-seq | 10x PBMC 10K, 3' v3 (638.9M pairs) | 16 min 8 s | 63 min 0 s | 3.9x |
| Perturb-seq | A375, gene expression + CRISPR guides | 2 min 21 s | 11 min 25 s | 4.8x |
| Perturb-seq | MSK 30-KO ES, gene expression + guides + LARRY, one pass | 28 min 39 s | 167 min 2 s | 5.8x |
| 10x Flex | JAX SC2300771 (2.01B pairs), delivered BGZF FASTQ | 3 min 19 s | 57 min 6 s | 17.2x |
| 10x Flex | JAX SC2300771, CBQ | 2 min 50 s | 57 min 6 s | 20.2x |
| 10x Flex | GSE325982 (1.12B pairs), plain gzip via rapidgzip | 4 min 28 s | 74 min 45 s | 16.8x |
| 10x Flex | GSE325982, CBQ | 1 min 47 s | 74 min 45 s | 41.8x |
| 10x Flex | 10x 320k scFFPE (7.30B pairs), delivered BGZF FASTQ (cloud) | 17 min 25 s | 6 h 37 min 38 s | 22.8x |
| 10x Flex | 10x 320k scFFPE, CBQ (cloud) | 10 min 45 s | 6 h 37 min 38 s | 37.0x |

| Assay | Benchmark (reference) | Concordance |
|---|---|---|
| Bulk RNA-seq | PPARG (Salmon) | gene Spearman 0.9992, Pearson 1.0000; transcript read-count Spearman 0.985, Pearson 0.99998 |
| scRNA-seq | 10x PBMC 10K (Cell Ranger 9.0.1) | 11,863 cells against 11,806; barcode Jaccard 0.995; cell Pearson 0.99999; gene Spearman 0.990, Pearson 0.9998 |
| Perturb-seq | A375 (Cell Ranger 9.0.1) | CRISPR calls 100% (1,076 shared cells); feature-UMI Pearson 0.99999; barcode Jaccard 0.992; cell Pearson 0.99995; gene Spearman 0.988, Pearson 0.980 |
| Perturb-seq | MSK 30-KO ES (Cell Ranger 9.0.1) | CRISPR calls 99.0%; feature-UMI Pearson 0.9994; cell Pearson 0.99998; gene Spearman 0.993, Pearson 0.9994 |
| 10x Flex | JAX SC2300771 (Cell Ranger 9.0.1) | barcode Jaccard 0.988 (pooled; 20,657 cells against 20,419); cell Pearson 0.9986; gene Spearman 0.99996, Pearson 0.999999 |
| 10x Flex | GSE325982 (submitters' Cell Ranger 9.0.1) | barcode Jaccard 0.994 (pooled; 38,634 cells against 38,444); cell Pearson 0.9952; gene Spearman 0.99991, Pearson 0.999999 |
| 10x Flex | 10x 320k scFFPE (Cell Ranger 9.0.1) | barcode Jaccard 0.974 (333,411 cells against 325,410); gene Spearman 0.99998, Pearson 0.999999 |
| SLAM-seq | GRAND-SLAM 100K human fixture (GRAND-SLAM) | NTR Pearson 0.9989 / 0.9961 / 0.9944 at >= 20 / 50 / 100 reads (30 s run) |

Gene Spearman and Pearson are computed on raw per-gene totals over the cells
both tools called, over every gene in both annotations (bulk: estimated read
counts against Salmon). Barcode Jaccard is on the called-cell sets; for
scRNA-seq and Perturb-seq, cell Pearson is on per-barcode UMI totals, and for
Flex it is the mean per-cell Pearson of log counts over shared cells, averaged
over samples. A Flex cell is its 16-base barcode together with its sample tag.

On the same PBMC 10K data, upstream STAR 2.7.11b with the community-optimized
CellGENI STARsolo parameters took 30 min 43 s and gave gene Spearman 0.950,
Pearson 0.993 and barcode Jaccard 0.994 against Cell Ranger 9.0.1; STAR Suite's
poly-G clipping and Cell Ranger-matched multimapper policy account for the
difference in gene-level counts. STAR-Flex loads no genome index: peak memory
was 23 GiB on JAX and 72 GiB on 320k with cell-barcode buckets held in memory.
On spinning disk (matched cloud instances, reported in the manuscript's
supplement) STAR-Flex remains 12-18x faster than Cell Ranger from the delivered
files.

## Building & Installing

### From source

```bash
# Core STAR binary
make core

# Module-focused builds
make flex           # core + Flex tools
make slam           # core + SLAM tools

# Individual tool targets
make feature-barcodes-tools    # assignBarcodes/demux (standalone)
make process-features-tools    # full process_features pipeline
make star-feature-call         # standalone feature caller
make vbem-tools                # TranscriptVB helpers
make yremove-tools             # Y/noY splitting tools

# Default build (core + common tools)
make                           # or: make default

# Build everything
make all
```

Selective filtering:

```bash
make default INCLUDE="core slam-tools"
make default EXCLUDE="flex-tools"
```

Run `make help` to see the full target list and descriptions.

### From release artifacts

```bash
# Ubuntu package from a local artifact
sudo apt install ./star-suite_<version>_<arch>.deb

# Installer tarball (auto-detects host glibc level)
tar -xzf STAR-suite-<version>-linux-<arch>-installer.tar.gz
cd STAR-suite-<version>-linux-<arch>-installer
./install.sh

# Manual compatibility tarball
tar -xzf STAR-suite-<version>-linux-<arch>-glibc234.tar.gz
cd STAR-suite-<version>-linux-<arch>-glibc234
./install.sh
```

Release tarballs are validated in clean Ubuntu 22.04 and 24.04 Docker containers before publication. The installer bundle auto-detects the host glibc level and chooses the right bundled binary.

Packaging/release details and artifact policy:
- `docs/Star-binary-distribution.md`
- `docs/Github-actions.md`

Compilation details (module-by-module, clean rebuilds, and clean Ubuntu 24.04 validation):
- `docs/compile_instructions.md`

## Docker

A multi-stage Docker setup (Ubuntu 24.04) provides a clean build environment and separate runtime/test images.

**Builder stage**: Compiles STAR Suite from source with no host leakage. Validates `make core`, `flex`, `slam`, `feature-barcodes-tools`, `default`, and `all`.

**Suite base runtime (`suite-base`)**: Minimal executable image with suite binaries (e.g. `STAR`) and no Python/test-only helpers.

**Test images** (built from `suite-base`):
- `test-tier-a`: self-contained smoke helpers.
- `test-tier-b`: fixture-backed helper stack (e.g. `python3`, `bc`, `samtools`).

### Quickstart

```bash
# Build suite base image (default tag: biodepot/star-suite:latest)
./scripts/docker/build_image.sh

# Override tag or parallel jobs
IMAGE_TAG=myorg/star-suite:v1 MAKE_JOBS=8 ./scripts/docker/build_image.sh

# Reproducibility check: force a clean rebuild (no cache)
docker build --no-cache --target suite-base -f docker/Dockerfile -t biodepot/star-suite:latest --build-arg MAKE_JOBS=8 .

# Run STAR from suite base image
docker run --rm biodepot/star-suite:latest

# Run Tier A smoke tests (builds/uses test-tier-a image)
./scripts/docker/run_smokes_tier_a.sh

# Run Tier B smoke tests (builds/uses test-tier-b image; requires fixtures)
./scripts/docker/run_smokes_tier_b.sh
```

### Fixture mount for Tier B

Tier B tests require data under `/storage`. Mount your fixture root:

```bash
docker run --rm -v /path/to/your/data:/storage biodepot/star-suite:test-tier-b bash -c "tests/run_cbub_regression_test.sh"
```

By default, `./scripts/docker/run_smokes_tier_b.sh` uses `STORAGE=/storage`.
Set `STORAGE=/path` to override (script uses it for the `-v` mount).

Expected layout: `/storage/A375`, `/storage/flex_filtered_reference`, etc. See `plans/docker_plan.md` for full fixture roots.

### STAR_BIN override

Smoke tests honor `STAR_BIN` to decouple from source-relative paths. Docker smoke wrappers set `STAR_BIN=/usr/local/bin/STAR` automatically.

### Validation

See [docs/docker_validation.md](docs/docker_validation.md) for the latest portability check results.

## Module Reference

This section documents the key features and flags for each module. For standard STAR flags not listed here, see `core/legacy/README.md`. Core additions are listed above in [Core Additions over STAR 2.7.11b](#core-additions-over-star-2711b).

### Flex

See [flex/README_flex.md](flex/README_flex.md) for the full pipeline reference.

From STAR Suite 1.9.4, STAR-Flex assigns each read to its probe from a half-probe (H1X2) hash cache and aligns nothing. Exact probe matches, and reads with at most one mismatch in each 25-base half of the probe, are resolved from the cache; halves that point to different probes are rejected as ambiguous. `--flex yes` selects this route by default and requires the cache, which is built once per probe set with `--runMode hashCacheGenerate --hashCacheTiers H0,H1X2` (see `hashCacheTiers` in the parameter reference).

The earlier alignment-based routes are **legacy**: probes embedded as pseudo-chromosomes in a hybrid reference and resolved by STAR's aligner, the alignment-validated H0/H1 cache, and BAM output with CB/UB tags or Y-chromosome splitting. They remain available with `--flexLegacy yes`, only to reproduce results from earlier releases.

Key flags:
- `--flex yes`: Enable Flex pipeline (half-probe route by default).
- `--soloHashScreenFile`: Half-probe (H1X2) cache; required with `--flex yes`.
- `--flexLegacy yes`: Permit a legacy alignment-based Flex route (reproducing earlier releases only).
- `--soloSampleWhitelist`: TSV mapping sample tags to labels.
- `--soloFlexAllowedTags`: TSV of the sample tags to call (the manuscript runs pass the sample whitelist here too).
- `--soloFlexCellCaller tag-aware` (default): calls tags that share a sample label as one sample, keeps CB16+TAG8 cell identities, and estimates cell counts itself; it rejects legacy expected-cell options such as `--soloFlexExpectedCellsPerTag`, which only `legacy` uses.
- `--soloProbeList`: Probe gene list (auto-detected from index if omitted).
- `--soloSampleProbes`: 10x probe barcode sequences file.

Features:
- Sample tag detection, 1MM pseudocount correction for CBs, clique-based UMI deduplication, and occupancy filtering.
- Y-chromosome splitting was tested on the legacy BAM route (`tests/TEST_REPORT_Y_SPLIT_FLEX.md`); the default 1.9.4 route writes no BAM.

#### Flex Parity: CR9-Projected Leiden UMAP

Using a fixed CR9 embedding removes the visual ambiguity from independently fit UMAPs. When full-align and no-align are both projected into the same CR9 PCA/UMAP space, they use the same 13 CR9 Leiden clusters and agree almost perfectly on shared cells: projected-label ARI `0.9979`, NMI `0.9967` on `20,315` shared cells.

| CR9 Reference | STAR-Flex Full Projected To CR9 | STAR-Flex No-Align Projected To CR9 |
|---|---|---|
| ![CR9 reference Leiden UMAP](docs/images/flex_parity/umap_sc2300771_cr9_reference.png) | ![STAR-Flex full projected to CR9 Leiden UMAP](docs/images/flex_parity/umap_sc2300771_fullalign_projected_to_cr9.png) | ![STAR-Flex no-align projected to CR9 Leiden UMAP](docs/images/flex_parity/umap_sc2300771_noalign_projected_to_cr9.png) |

### SLAM

See [slam/docs/SLAM_COMPATIBILITY_MODE.md](slam/docs/SLAM_COMPATIBILITY_MODE.md) and [slam/docs/SLAM_seq.md](slam/docs/SLAM_seq.md).

Integrated SLAM-seq quantification with paired-end support, GRAND-SLAM parity,
count-binomial output, and tximport-ready TranscriptVB gene counts:

Key flags:
- `--slamQuantMode 1`: Enable SLAM quantification.
- `--slamGrandSlamOut 1`: Generate GRAND-SLAM compatible output.
- `--slamCbOut 1 --slamCbFormat star|ezbakr`: Generate model-ready count-binomial output.
- `--slamMinCallableLength 30`: Require a minimum callable post-trim/overlap-consensus evidence length for SLAM transition statistics.
- `--slamCompatMode gedi`: Enable GEDI compatibility (intronic classification, lenient overlap, overlap weighting).
- `--slamCompatIntronic`, `--slamCompatLenientOverlap`: Fine-grained compat control.
- `--autoTrim variance`: Variance-based detection of artifact-prone read ends.
- `--slamTrim5p`, `--slamTrim3p`, `--slamCompatTrim5pMate1`, `--slamCompatTrim3pMate1`, `--slamCompatTrim5pMate2`, `--slamCompatTrim3pMate2`: Manual SE/PE trim guards.
- `--slamErrorRateFromBlank 1`: Seed error rate from a blank (e.g. no4sU) sample.
- `--outFileNamePrefixAuto 1`: Derive sample name from first FASTQ and route outputs into subdirs.
- `--slamDumpBinary 1 --slamDumpWeights 1`: Emit binary dumps for offline re-quantification with `slam_requant`.

Features:
- Full gene-level NTR estimation (Binomial/EM models).
- Paired-end transition coordinate handling with overlap consensus before counting.
- Fixed 2026-05 PE smoke trims for the NW panel: SE R1 `8/12`; PE R1 `8/13`, R2 `19/14`.
- Auto-trimming: variance-based detection of artifact-prone read ends.
- QC: comprehensive interactive HTML reports for T->C rates and error modeling.
- Batch layout organizes outputs into `alignments/`, `counts/`, `qc/`, `y_separated/`.
- Binary dump format documented in `slam/docs/SLAM_DUMP_FORMAT.md`.
- Reproducible PE smoke and production runbooks: `docs/RUNBOOK_SLAM_PE_100K_SMOKE.md`,
  `docs/RUNBOOK_SLAM_PE_DESEQ2_COUNT_SURFACES.md`, and
  `docs/RUNBOOK_SLAM_PE_PRODUCTION.md`.

### STAR-perturb / CR-Compat

See [docs/feature_barcodes.md](docs/feature_barcodes.md) and [docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md](docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md).

CR-compatible Solo behavior with integrated CRISPR feature calling:

Key flags:
- `--pfMultiConfig`: Cell Ranger-style multi processing with feature libraries.
- `--defaultCrCompat yes`: Apply the CR-compat perturb defaults bundle.
- `--dynamicThreadInterface 1`: Enable STAR/PF permit coordination.
- `--dynamicThreadConstMapPermits 32`: Start with full map-side permit budget.
- `--crAssignConsumerThreads 32`: Provision PF worker pool to full host budget.
- `--crAssignSearchThreads 1`: Per-consumer search-thread mode.
- `--crMinUmi`: Minimum UMI threshold for CRISPR feature calling (default `3`; `--defaultCrCompat yes` sets `10`; the A375 benchmark used `10` and the MSK LARRY library `2`).
- `--soloStrand`: Must match the gene-expression library: `Forward` for 10x 3', `Reverse` for 10x 5' libraries sequenced from read 2 only (Cell Ranger chemistry `SC5P-R2*`, e.g. A375).
- `--soloCrGexFeature`: Control merged GEX source (`auto`, `gene`, `genefull`).
- `--soloCrMode CR`: Enable CR-compatible single-cell behavior.
- `--crChemistry`: Barcode chemistry (`auto`, `NXT`, `TRU`). Default `auto` enables per-library auto-detection. Mixed NXT/TRU experiments are handled automatically; per-library overrides via the `star_chemistry` column in `--pfMultiConfig`.

Recommended execution profile (32-thread host):

```bash
--runThreadN 32 --dynamicThreadInterface 1 --dynamicThreadConstMapPermits 32 \
--dynamicThreadTelemetry 1 --crAssignConsumerThreads 32 --crAssignSearchThreads 1
```

Standalone tool (`star_feature_call`):
- `--compat-perturb`: CR9-compatible output layout (`crispr_analysis/`).
- `--feature-ref`, `--whitelist`, `--fastq-dir`, `--output-dir`: FASTQ -> MEX -> calls.
- `--call-only --mex-dir`: call_features-only pass on existing MEX.
- `--emptydrops-use-fdr`, `--min-umi`, `--ratio-test`: calling controls.

### OCM scRNA-seq

See [docs/RUNBOOK_SCRNA_OCM_CR_COMPAT.md](docs/RUNBOOK_SCRNA_OCM_CR_COMPAT.md),
[docs/RUNBOOK_SCRNA_OCM_MULTI_MEX_MATERIALIZER_IMPLEMENTATION_20260519.md](docs/RUNBOOK_SCRNA_OCM_MULTI_MEX_MATERIALIZER_IMPLEMENTATION_20260519.md),
and [docs/RUNBOOK_JAX_SCRNASEQ02_OCM_20260518.md](docs/RUNBOOK_JAX_SCRNASEQ02_OCM_20260518.md).

OCM is treated as a CR-compatible GEX run with an effective sample-aware cell
barcode, not as a guide-feature library. In production, STAR derives an OCM tag
from bases 8-9 of the raw 16 bp barcode, appends the internal TAG8 suffix before
barcode correction/counting, runs per-sample CR-compatible EmptyDrops after the
OCM split, and later strips that suffix from Cell Ranger-compatible output
labels.

Key flags:
- `--ocmMultiEnable yes`: emit OCM multi-compatible outputs.
- `--ocmMultiConfig <config.csv>`: Cell Ranger multi-style config with
  `[samples]` and `ocm_barcode_ids`.
- `--ocmMultiBarcodeMode flex`: production mode; count on `CB16+OCM_TAG8`.
- `--ocmMultiOutputCompat cellranger`: writes `outs/multi`,
  `outs/per_sample_outs`, and downstream `samples/<sample>/run/outs` mirrors.
- `--soloFeatures GeneFull Velocyto`: expression and velocity surface used by
  the downstream h5ad/CellBender path.
- `--soloCellFilter None`: current split-before-ED OCM production mode; the
  native OCM materializer applies CR-compatible EmptyDrops separately per OCM
  sample.

OCM production should also use the dataset-specific whitelist family, the
MSK/UCSF GRCh38 2024-A STAR reference, and Y-removal for KOLF2-derived JAX
samples.

### QC Outputs

- **SLAM QC** (`--slamQcReport <prefix>`): Interactive HTML report (`.html`) and JSON metrics (`.json`) for T->C conversion rates, variance analysis, and trimming overlays.
- **FlexFilter QC** (`flexfilter_summary.tsv`): Cell calling statistics (EmptyDrops/OrdMag), cell counts, UMI thresholds, and filtering rates per sample.

## Sample Commands

**Core alignment:**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq.gz \
  --outFileNamePrefix out/ \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM MD
```

**Batch mode (bulk, single-pass, SE):**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn A_R1.fq.gz,B_R1.fq.gz \
  --outFileNamePrefix /path/to/out_root/ \
  --outFileNamePrefixAuto 1 \
  --batchMode 1 \
  --outSAMtype BAM SortedByCoordinate
```

**Batch mode (bulk, single-pass, PE):**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn A_R1.fq.gz,B_R1.fq.gz A_R2.fq.gz,B_R2.fq.gz \
  --outFileNamePrefix /path/to/out_root/ \
  --outFileNamePrefixAuto 1 \
  --batchMode 1 \
  --outSAMtype BAM SortedByCoordinate
```

**Flex Mode (10x Fixed RNA Profiling):**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/flex_index \
  --readFilesIn reads_R2.fq.gz reads_R1.fq.gz \
  --flex yes \
  --soloHashScreenFile /path/to/flex_h01x2_cache.half.khash \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /path/to/737K-fixed-rna-profiling.txt \
  --soloSampleWhitelist sample_whitelist.tsv \
  --outFileNamePrefix output/
```

**SLAM Mode (Standard):**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq.gz \
  --outFileNamePrefix out/ \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM MD \
  --slamQuantMode 1 \
  --slamSnpBed /path/to/snps.bed
```

**SLAM Mode (GEDI Compatibility):**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn reads.fq.gz \
  --slamQuantMode 1 \
  --slamCompatMode gedi \
  --autoTrim variance \
  --outFileNamePrefix output/
```

**SLAM PE 100K smoke (R1-only SE vs R1/R2 PE):**

```bash
bash scripts/run_slam_100k_se_pe_smoke.sh \
  --sample ARID1A-no4su_S50 \
  --sample ARID1A-6h-1_S43 \
  --threads 16
```

**SLAM PE production panel (safe dry-run default):**

```bash
bash scripts/run_slam_prod_set.sh \
  --pilot \
  --dry-run \
  --globus-dst-endpoint 61fb8b9a-9b52-456e-928c-30c0fb0140bf \
  --globus-dst-root SLAM-seq-PE-results
```

**SLAM Batch Mode (blank-first, SE/PE):**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn blank_R1.fq.gz,0h_R1.fq.gz,6h_R1.fq.gz,24h_R1.fq.gz \
  --outFileNamePrefix /path/to/out_root/ \
  --outFileNamePrefixAuto 1 \
  --slamQuantMode 1 \
  --slamBatchMode 1 \
  --slamErrorRateFromBlank 1 \
  --slamSnpBed /path/to/snps.bed
```

For paired-end, pass **two comma-separated mate lists**:
`--readFilesIn blank_R1.fq.gz,0h_R1.fq.gz,... blank_R2.fq.gz,0h_R2.fq.gz,...`

**STAR-perturb (integrated CR-compat mode):**

The gene-expression FASTQs go in `--readFilesIn` (cDNA read first); the feature
libraries are listed in `--pfMultiConfig`. The option set below is the one used
for the manuscript's Perturb-seq benchmarks
([docs/PAPER_BENCHMARK_METHODOLOGY.md](docs/PAPER_BENCHMARK_METHODOLOGY.md)
Section 1.6); explicit options take precedence over the `--defaultCrCompat`
bundle, which on its own would select `Rescue` multimappers and
`Gene GeneFull` and sets no strand.

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --runThreadN 32 \
  --genomeDir /path/to/index \
  --readFilesIn gex_R2.fastq.gz gex_R1.fastq.gz \
  --readFilesBgzfMode auto \
  --pfMultiConfig /path/to/multi_config.csv \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist /path/to/gex_whitelist.txt \
  --soloStrand Forward \
  --clipAdapterType CellRanger4 --clip3pPolyG yes \
  --alignEndsType Local --chimSegmentMin 1000000 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique --soloCrMultimapRescue yes \
  --soloCellFilter EmptyDrops_CR --soloCbUbRequireTogether no \
  --soloFeatures GeneFull --soloCrGexFeature genefull \
  --crMinUmi 10 \
  --dynamicThreadInterface 1 \
  --dynamicThreadConstMapPermits 32 \
  --crAssignSearchThreads 1 \
  --outSAMtype None \
  --outFileNamePrefix /path/to/outs/
```

Use `--soloStrand Reverse` (and omit `--clip3pPolyG` if the data are not
two-colour NovaSeq/NextSeq) for 10x 5' libraries sequenced from read 2 only,
such as the A375 benchmark.

**OCM scRNA-seq (native composite barcode mode):**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --runThreadN 16 \
  --genomeDir /storage/autoindex_110_44/bulk_index \
  --readFilesIn "${R2_FILES}" "${R1_FILES}" \
  --readFilesCommand zcat \
  --outFileNamePrefix /path/to/library/run/ \
  --outSAMtype BAM Unsorted \
  --emitNoYBAM yes \
  --emitYNoYFastq yes \
  --clipAdapterType CellRanger4 \
  --clip3pPolyG yes \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloCBwhitelist /storage/scRNAseq_output/whitelists/3M-3pgex-may-2023_TRU.txt \
  --soloInlineCBCorrection yes \
  --soloCellFilter None \
  --soloFeatures GeneFull Velocyto \
  --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes \
  --ocmMultiEnable auto \
  --ocmMultiConfig /path/to/cellranger_multi_config.csv \
  --ocmMultiBarcodeMode flex \
  --ocmMultiOutputCompat cellranger
```

**STAR-perturb (standalone feature pipeline):**

```bash
core/legacy/source/star_feature_call \
  --compat-perturb \
  --feature-ref /path/to/feature_reference.csv \
  --whitelist /path/to/whitelist.txt \
  --fastq-dir /path/to/feature_fastqs \
  --filtered-barcodes /path/to/filtered_barcodes.tsv \
  --output-dir /path/to/feature_out \
  --emptydrops-use-fdr \
  --min-umi 10
```

## STAR Launchpad (Recipe Builder)

If you prefer not to assemble command lines by hand, STAR Launchpad is a
browser-based recipe builder served from the STAR-suite MCP server. Select a
workflow, fill in parameters through a guided form, and get a validated,
copy-pasteable shell command.

![STAR Launchpad screenshot](docs/images/launchpad/launch.png)

### Quick start

```bash
# Install dependencies (once)
pip install -r mcp_server/requirements.txt

# Start Launchpad + MCP on one port
bash scripts/launchpad_server.sh up

# Open in your browser
# http://localhost:8765/launchpad/
```

### MCP endpoints

When the server is running on port `8765`, the main endpoints are:

- `http://127.0.0.1:8765/launchpad/` — Launchpad UI
- `http://127.0.0.1:8765/` — MCP streamable-HTTP endpoint
- `http://127.0.0.1:8765/sse` — MCP SSE endpoint

### MCP client setup

Use the endpoint that matches your client. If the server is running on another
host or port, replace the URL accordingly.

#### VS Code / GitHub Copilot

VS Code and GitHub Copilot use MCP over HTTP. Add a workspace config at
`.vscode/mcp.json`:

```json
{
  "servers": {
    "starSuite": {
      "type": "http",
      "url": "http://127.0.0.1:8765/"
    }
  }
}
```

Then reload the VS Code window or use the MCP commands in the Command Palette
to restart the server listing.

#### Cursor

Cursor uses the SSE MCP endpoint. Add this to your Cursor MCP config:

```json
{
  "mcpServers": {
    "star-suite": {
      "url": "http://127.0.0.1:8765/sse"
    }
  }
}
```

If Cursor is already open, reload the window after updating the config.

#### Claude

Claude can register the server over streamable HTTP:

```bash
claude mcp add --transport http star-suite http://127.0.0.1:8765/
```

Useful follow-up commands:

```bash
claude mcp list
claude mcp get star-suite
```

If you are connecting from another machine, use the remote host or public URL
instead of `127.0.0.1`.

### What it does

1. **Pick a recipe** -- choose a workflow (e.g. A375 CR-compatible alignment).
2. **Fill the form** -- parameters are grouped by category with defaults
   pre-filled, descriptions as help text, and constraints shown inline.
3. **Validate** -- the server checks required fields, types, and constraint
   rules. Errors and warnings appear inline.
4. **Generate** -- get the full bash command, environment variable overrides,
   and a checklist of required input files.
5. **Copy and run** -- paste the command into your terminal.

Launchpad does not execute anything. It generates commands; you run them.

Design details: [`plans/star_launchpad_v1_runbook.md`](plans/star_launchpad_v1_runbook.md)

## Codespaces Walkthroughs

STAR Suite includes GitHub Codespaces walkthroughs for the main module entry points.

Start here:
- [Codespaces walkthrough summary](docs/CODESPACES_DEMO_WALKTHROUGHS_20260312.md)
- [Codespaces overview](docs/codespaces/00_overview.md)

Ready now:
- [Optional setup: build the small demo reference](docs/codespaces/01_setup_reference.md)
- [Bulk demo](docs/codespaces/02_bulk.md)
- [SLAM demo](docs/codespaces/03_slam.md)
- [Single-cell fixture builder](docs/codespaces/04_single_cell_fixture.md)

Work in progress:
- [Perturb demo](docs/codespaces/05_perturb.md)
- [Flex demo](docs/codespaces/06_flex.md)

Helpful follow-up guides:
- [If you already use STAR or Cell Ranger](docs/codespaces/07_star_cellranger_users.md)
- [Using your own data](docs/codespaces/08_using_your_own_data.md)

## More Detail

- Core usage: [core/legacy/README.md](core/legacy/README.md)
- Flex pipeline: [flex/README_flex.md](flex/README_flex.md)
- SLAM compatibility: [slam/docs/SLAM_COMPATIBILITY_MODE.md](slam/docs/SLAM_COMPATIBILITY_MODE.md)
- SLAM methodology: [slam/docs/SLAM_seq.md](slam/docs/SLAM_seq.md)
- STAR-perturb feature docs: [docs/feature_barcodes.md](docs/feature_barcodes.md)
- OCM scRNA-seq runbook: [docs/RUNBOOK_SCRNA_OCM_CR_COMPAT.md](docs/RUNBOOK_SCRNA_OCM_CR_COMPAT.md)
- Velocyto CR-compat policy runbook: [docs/RUNBOOK_VELOCYTO_CR_COMPAT_POLICY_20260519.md](docs/RUNBOOK_VELOCYTO_CR_COMPAT_POLICY_20260519.md)
- STAR Suite binary distribution: [docs/Star-binary-distribution.md](docs/Star-binary-distribution.md)
- STAR-perturb A375 parity report: [tests/crispr_feature_calling_comparison_report.md](tests/crispr_feature_calling_comparison_report.md)
- Cell Ranger multi smoke tool: [docs/cr_multi.md](docs/cr_multi.md)
- Docker validation: [docs/docker_validation.md](docs/docker_validation.md)
