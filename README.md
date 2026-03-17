# STAR-suite

STAR-suite reorganizes STAR into module-focused directories while keeping a single source of truth for shared code. Build outputs remain compatible with existing STAR workflows, and the new top-level `Makefile` exposes module targets.

**No new external dependencies are required for the suite modules in this repo.**
The current integrations (including STAR-perturb, STAR-SLAM, and QC updates) are
built with the existing toolchain and vendored components.

STAR-suite supports partial compilation: build only the module/tool targets you need
instead of building the full suite every time.

Agent quickstart: see `AGENTS.md` for repo-specific guardrails, tests, and recent changes.

## Folder Structure

```
core/
  legacy/                        # Upstream STAR layout (single source of truth)
  features/                      # Shared overlays and feature tooling
    process_features/            # Perturb feature extraction/calling implementation
    feature_barcodes/            # Standalone barcode tools (assignBarcodes, demux)
    libscrna/                    # EmptyDrops/OrdMag/Occupancy shared library
flex/                    # Flex-specific code + tools
slam/                    # SLAM-seq code + tools
build/                   # Modular make fragments
docs/                    # Suite-level docs
tests/                   # Suite-level tests (see tests/ARTIFACTS.md for artifact locations)
tools/                   # Suite-level scripts/utilities
mcp_server/              # MCP server for scripted discovery/preflight/run workflows
```

## Modules

- **STAR-core** (`core/`): Legacy STAR (indexing, bulk, Solo) plus shared utilities.
  Build: `make core` (binary at `core/legacy/source/STAR`).
- **STAR-perturb** (`core/legacy/` + `core/features/process_features/`): CR-compatible perturb-seq path with integrated feature extraction/calling (`process_features` + `call_features`) and `crispr_analysis/` outputs in CR-compat mode.
  Primary run path: `STAR --pfMultiConfig ... --defaultCrCompat yes` (see STAR-perturb section below).
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
- **MCP Server (tooling)** (`mcp_server/`): Agent automation service for dataset/test discovery and controlled execution (`list_datasets`, `list_test_suites`, `preflight`, `run_script`, `collect_outputs`). This is repo tooling, not an analysis module.

## Benchmarks

### Perturb-seq (STAR vs Cell Ranger 9)

| Dataset | Cells (STAR / CR) | Jaccard | Gene Pearson | Cell Pearson | CRISPR match | Speedup |
|---|---|---|---|---|---|---|
| UCSF EBs2_2 (full, NXT) | 13,721 / 13,760 | 0.98 | 0.995 | 1.000 | 98.9% (11,032/11,157) | 2.5x |
| A375 1k CRISPR 5' (GeneFull) | 1,186 / 1,163 | 0.98 | 0.933 | 0.999 | 96.8% (1,048/1,083) | 3.9x |
| MSK 30polyKO (full, NXT) | 30,497 / 32,256 | 0.94 | 0.993 | 1.000 | 98.5% (22,200/22,531) | 4.0x |

- UCSF EBs2_2: Gene Pearson on 24,343 filtered genes; Cell Pearson on 13,571 common barcodes; CRISPR top-guide match on 11,157 common evaluated barcodes; speedup = 20 min vs 52 min (32 threads, no BAM, dynamic permits, bootstrap tiered hash).
- A375: Gene Pearson on 17,253 filtered genes; Cell Pearson on 1,163 common barcodes; CRISPR exact-match at min-UMI 10; speedup = 4 min vs 15 min (32 threads, no BAM, dynamic permits, bootstrap tiered hash).
- MSK: Gene Pearson on 17,448 filtered genes; Cell Pearson on 30,417 common barcodes; CRISPR set-equivalent calls on 22,531 evaluated rows (30 guides, min-UMI 2); speedup = 42 min vs 168 min (32 threads, with BAM, dynamic permits). CR requires two separate runs (GEX+gRNA 58 min + GEX+LARRY 110 min); STAR handles all three libraries in a single pass.

### Flex (STAR vs Cell Ranger 7.2)

| Dataset | Cells (STAR / CR) | Jaccard | Gene Pearson | Cell Pearson | Speedup |
|---|---|---|---|---|---|
| JAX SC2300771 full (4 samples) | 20,291 / 20,444 | 0.98 | 0.998 | 1.000 | pending |

- Gene Pearson on 18,021 common genes; Cell Pearson on 20,173 shared barcodes.
- Speedup not yet measured (optimization pending).

### PE Bulk (Integrated STAR-suite vs External Stepwise Pipeline)

"External stepwise" = Trim Galore + STAR align + remove\_y\_reads + Salmon quant (sequential).

| Dataset | Transcript Pearson | Gene Pearson | Speedup |
|---|---|---|---|
| JAX PE (full, 32 threads) | 0.995 | 0.997 | 2.1x |

- TranscriptVB vs Salmon (alignment-mode VB) on expressed transcripts.
- Integrated: 61 s vs external stepwise: 125 s (32 threads).

### SLAM-seq (STAR-SLAM vs GrandSLAM/GEDI)

**NTR parity (compat mode, GEDI is reference):**

| Dataset | Sample | NTR Pearson | NTR Spearman |
|---|---|---|---|
| NW-5-21 ARID1A 1M (compat, no trim) | 0h | 0.978 | 0.990 |
| NW-5-21 ARID1A 1M (compat, no trim) | 6h | 0.972 | 0.986 |
| NW-5-21 ARID1A 1M (compat, no trim) | 24h | 0.967 | 0.985 |
| 100K fixture (SNP BED, ≥20 reads) | -- | 0.999 | 0.981 |

- Comparison uses SNP-masked BAMs; GEDI is reference.
- `slam_requant` replay: Pearson/Spearman 1.0 (exact parity with STAR output).
- Compat mode (`--slamCompatMode gedi`) adds negligible overhead (<0.1% wall time, <1% memory).
- Direct speedup comparison is not reported because GRAND-SLAM depends on alignment being completed first (it operates on pre-aligned BAMs), whereas STAR-SLAM performs alignment and quantification in a single pass. On the ARID1A time-course (167M reads, 4 samples), GEDI quantification alone adds ~14% to the alignment time (~5.5 min on top of ~40 min alignment).

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

## Codespaces Walkthroughs

STAR-suite includes GitHub Codespaces walkthroughs for the main module entry points.

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

## Docker

A multi-stage Docker setup (Ubuntu 24.04) provides a clean build environment and separate runtime/test images.

**Builder stage**: Compiles STAR-suite from source with no host leakage. Validates `make core`, `flex`, `slam`, `feature-barcodes-tools`, `default`, and `all`.

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

This section documents the key features and flags for each module. For standard STAR flags not listed here, see `core/legacy/README.md`.

### Core

STAR-suite is built on top of STAR 2.7.11b. The following features are STAR-suite additions to the core:

- **Batch Mode** (`--batchMode 1`): Processes multiple FASTQs in one STAR invocation while reusing the loaded genome. Removes the need for `--genomeLoad` keep-in-memory workflows. Single-pass only (no `--twopassMode`); not supported with Solo (`--soloType`). Use `--outFileNamePrefixAuto 1` for per-sample subdirectories.
- **TranscriptVB Quantification** (`--quantMode TranscriptVB`): Variational Bayes and EM quantification for transcript-level abundance, with parity-oriented behavior against Salmon alignment-mode. Gene-level summarization via `--quantVBgenesMode Tximport`.
- **Transcriptome Output** (`--quantTranscriptomeSAMoutput`): Replaces the former `--quantTranscriptomeBan` with more explicit control (e.g., `BanSingleEnd_ExtendSoftclip`).
- **Reference Automation** (`--autoIndex Yes`): Automated reference download/build with `--cellrangerStyleIndex Yes` formatting and `--genomeGenerateTranscriptome Yes` for transcript-level quant workflows.
- **Cutadapt-Compatible Trimming** (`--trimCutadapt Yes`): Native cutadapt-style trimming for bulk/PE workflows. Compatibility mode: `--trimCutadaptCompat Cutadapt3`.
- **Poly-G Trimming** (`--clip3pPolyG yes|no|auto`): Trims poly-G artifacts common on NovaSeq/NextSeq platforms. Default `auto` activates in CellRanger4 mode. Without this, poly-G reads can inflate specific genes (e.g., LINC00486) and degrade gene-level correlations.
- **Samtools-style BAM Sorting** (`--outBAMsortMethod samtools`): Spill-to-disk sort to reduce peak RAM pressure. Works with all modes including Flex.
- **Y/NoY Separation** (`--emitNoYBAM yes`, `--emitYNoYFastq yes`): Split BAM and FASTQ outputs by chrY alignment. Works with bulk, single-cell, and Flex.
- **EmptyDrops_CR Integration**: CR-compatible EmptyDrops path (including libscrna-backed behavior in scRNA/perturb flows).
- **Solo Features**: `sF` BAM tag for feature type, `--soloCBtype String` for arbitrary barcode strings, `--soloCellReadStats Standard` for improved cell filtering.
- **CR-compat GEX** (`--soloCrGexFeature auto|gene|genefull`): Controls which GEX source is merged in CR-compat mode.
- **CB/UB Tag Pairing** (`--soloCbUbRequireTogether yes|no`): Enforce CB/UB tag pairing for tag injection (default `yes`).

### Flex

See [flex/README_flex.md](flex/README_flex.md) for the full pipeline reference.

STAR-Flex uses a pseudo-chromosome alignment approach: probe sequences are embedded as pseudo-chromosomes in a hybrid reference genome, and STAR's native alignment machinery handles gene assignment. Core features (trimming, spill-to-disk sorting, Y-chromosome splitting, TranscriptVB) all work with Flex.

Key flags:
- `--flex yes`: Enable Flex pipeline.
- `--soloFlexExpectedCellsPerTag`: Expected cells per sample tag.
- `--soloSampleWhitelist`: TSV mapping sample tags to labels.
- `--soloProbeList`: Probe gene list (auto-detected from index if omitted).
- `--soloSampleProbes`: 10x probe barcode sequences file.

Features:
- Sample tag detection, 1MM pseudocount correction for CBs, clique-based UMI deduplication, and occupancy filtering.
- Y-chromosome splitting tested and validated (`tests/TEST_REPORT_Y_SPLIT_FLEX.md`).

### SLAM

See [slam/docs/SLAM_COMPATIBILITY_MODE.md](slam/docs/SLAM_COMPATIBILITY_MODE.md) and [slam/docs/SLAM_seq.md](slam/docs/SLAM_seq.md).

Integrated SLAM-seq quantification with GRAND-SLAM parity:

Key flags:
- `--slamQuantMode 1`: Enable SLAM quantification.
- `--slamGrandSlamOut 1`: Generate GRAND-SLAM compatible output.
- `--slamCompatMode gedi`: Enable GEDI compatibility (intronic classification, lenient overlap, overlap weighting).
- `--slamCompatIntronic`, `--slamCompatLenientOverlap`: Fine-grained compat control.
- `--autoTrim variance`: Variance-based detection of artifact-prone read ends.
- `--slamTrim5p`, `--slamTrim3p`: Manual trim guards.
- `--slamErrorRateFromBlank 1`: Seed error rate from a blank (e.g. no4sU) sample.
- `--outFileNamePrefixAuto 1`: Derive sample name from first FASTQ and route outputs into subdirs.
- `--slamDumpBinary 1 --slamDumpWeights 1`: Emit binary dumps for offline re-quantification with `slam_requant`.

Features:
- Full gene-level NTR estimation (Binomial/EM models).
- Auto-trimming: variance-based detection of artifact-prone read ends.
- QC: comprehensive interactive HTML reports for T->C rates and error modeling.
- Batch layout organizes outputs into `alignments/`, `counts/`, `qc/`, `y_separated/`.
- Binary dump format documented in `slam/docs/SLAM_DUMP_FORMAT.md`.

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
- `--crMinUmi`: Minimum UMI threshold for CRISPR feature calling (default `10`; lower to `2-3` for lineage barcodes).
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
  --readFilesCommand zcat \
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
  --readFilesCommand zcat \
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
  --readFilesCommand zcat \
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
  --readFilesCommand zcat \
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

**SLAM Batch Mode (blank-first, SE/PE):**

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --genomeDir /path/to/genome_index \
  --readFilesIn blank_R1.fq.gz,0h_R1.fq.gz,6h_R1.fq.gz,24h_R1.fq.gz \
  --readFilesCommand zcat \
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

```bash
core/legacy/source/STAR \
  --runMode alignReads \
  --runThreadN 32 \
  --genomeDir /path/to/index \
  --pfMultiConfig /path/to/multi_config.csv \
  --dynamicThreadInterface 1 \
  --dynamicThreadConstMapPermits 32 \
  --dynamicThreadTelemetry 1 \
  --crAssignConsumerThreads 32 \
  --crAssignSearchThreads 1 \
  --defaultCrCompat yes \
  --crChemistry auto \
  --outFileNamePrefix /path/to/outs/
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

## More Detail

- Core usage: [core/legacy/README.md](core/legacy/README.md)
- Flex pipeline: [flex/README_flex.md](flex/README_flex.md)
- SLAM compatibility: [slam/docs/SLAM_COMPATIBILITY_MODE.md](slam/docs/SLAM_COMPATIBILITY_MODE.md)
- SLAM methodology: [slam/docs/SLAM_seq.md](slam/docs/SLAM_seq.md)
- STAR-perturb feature docs: [docs/feature_barcodes.md](docs/feature_barcodes.md)
- STAR-perturb A375 parity report: [tests/crispr_feature_calling_comparison_report.md](tests/crispr_feature_calling_comparison_report.md)
- Cell Ranger multi smoke tool: [docs/cr_multi.md](docs/cr_multi.md)
- Docker validation: [docs/docker_validation.md](docs/docker_validation.md)
