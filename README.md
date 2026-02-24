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
    feature_barcodes/            # assignBarcodes/demux tooling
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
- **Feature Barcodes** (`core/features/feature_barcodes/`): Vendored `process_features` tools for perturb-seq testing (`assignBarcodes`, `demux_bam`, `demux_fastq`).
  Build tools: `make feature-barcodes-tools`.
- **Shared Feature Toolchains** (`core/features/`): Reusable tool layers used across modules, including `vbem` (TranscriptVB helpers), `yremove_*` (Y/noY splitting), `bamsort`, and `libscrna`.
  Build tools: `make vbem-tools`, `make yremove-tools`, plus in-core integrations.
- **MCP Server (tooling)** (`mcp_server/`): Agent automation service for dataset/test discovery and controlled execution (`list_datasets`, `list_test_suites`, `preflight`, `run_script`, `collect_outputs`). This is repo tooling, not an analysis module.

## Installation

Build from repo root:

```bash
# Core STAR binary
make core

# Module-focused builds
make flex
make slam
make feature-barcodes-tools

# Build everything
make all
```

Selective default build:

```bash
make default INCLUDE="core slam-tools"
make default EXCLUDE="flex-tools"
```

Install alternatives (release artifacts):

```bash
# Ubuntu package from a local artifact
sudo apt install ./star-suite_<version>_<arch>.deb

# Optional PPA path (once published)
sudo add-apt-repository ppa:<org>/star-suite
sudo apt update
sudo apt install star-suite
```

```bash
# Static tarball install
tar -xzf STAR-static-<version>-linux-<arch>.tar.gz
sudo install -m 0755 STAR-static-<version>-linux-<arch>/bin/STAR /usr/local/bin/STAR
```

Packaging/release details and artifact policy:
- `docs/Star-binary-distribution.md`
- `docs/Github-actions.md`

Compilation details (module-by-module, clean rebuilds, and clean Ubuntu 24.04 validation):
- `docs/compile_instructions.md`

## Build Targets

The top-level `Makefile` supports a default build, full build, and conditional include/exclude filters.

- **Default build**: `make` (same as `make default`)
  - Builds the “usual culprits” (core + common tools).
  - Optional filters:
    - `make default INCLUDE="core flex-tools"`
    - `make default EXCLUDE="slam-tools yremove-tools"`
- **Build everything**: `make all`
  - Includes everything in the suite (core + all tools).

Run `make help` to see the full target list and descriptions.

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

### Validation Status

Recent portability checks on this branch (Ubuntu 24.04 Docker):

- `docker build --no-cache --target suite-base -f docker/Dockerfile --build-arg MAKE_JOBS=8 .` passed.
- `scripts/docker/run_smokes_tier_a.sh` passed.
- `scripts/docker/run_smokes_tier_b.sh` with `STORAGE=/storage` produced:
  - pass: `run_flex_smoke`, `run_cbub_regression_test`
  - skip: `run_cr_compat_integration_smoke` (ASAN-only path), `run_cr_multi_smoke` (external script dependency), `run_unsorted_cbub_smoke_test` (missing probes fixture)
  - fail: `test_cr_compat_crispr_calling` (missing expected `outs/*` files in this containerized fixture run)
- Runtime split verified:
  - `suite-base` and `test-tier-a` have no Python runtime
  - `test-tier-b` includes Python runtime for fixture-backed tests
- `suite-base` contains: `STAR`, `star_feature_call`, `slam_requant`, `pileup_snp`, `flexfilter` (`run_flexfilter_mex` alias), and `remove_y_reads`.

## Technical Updates

### Core Updates
Recent updates to the Core module (STAR 2.7.11b and prior) include:
- **Batch Mode (single-pass, non-Solo)**: `--batchMode 1` processes multiple FASTQs in one STAR invocation while reusing the loaded genome. This removes the need for `--genomeLoad` keep-in-memory workflows that are often brittle in containerized and HPC job environments. It is also important when analyses require shared static inputs across many samples (for example SLAM SNP masks and blank-derived background/error settings), so each sample is processed under the same fixed context.
  - **Limits**: batch mode is single-pass only (no `--twopassMode`) and not supported with Solo (`--soloType`).
  - **Output routing**: use `--outFileNamePrefixAuto 1` for per-sample subdirectories under one output root.
- **Transcriptome Output**: Replaced `--quantTranscriptomeBan` with `--quantTranscriptomeSAMoutput` for more explicit control (e.g., `BanSingleEnd_ExtendSoftclip`).
- **TranscriptVB Quantification**: Variational Bayes and EM quantification for transcript-level abundance (`--quantMode TranscriptVB`), with parity-oriented behavior against Salmon alignment-mode.
- **Reference Automation**: Automated reference download/build (`--autoIndex`, `--autoCksumUpdate`, `--cellrangerLegacyGtfFilter`) plus automatic `transcriptome.fa` generation during indexing for transcript-level quant workflows. CellRanger parity is guaranteed for formatted reference inputs (`cellranger_ref/genome.fa`, `cellranger_ref/genes.gtf`) only.
- **Cutadapt-Compatible Trimming**: Native cutadapt-style trimming path (`--trimCutadapt Yes`) for bulk/PE workflows.
- **Samtools-style BAM Sorting**: Spill-to-disk sort (`--outBAMsortMethod samtools`) to reduce peak RAM pressure versus in-memory bin sorting.
- **Y/NoY Separation**: Split BAM and FASTQ outputs by chrY alignment (`--emitNoYBAM`, `--emitYNoYFastq`).
- **EmptyDrops_CR Integration**: CR-compatible EmptyDrops path (including libscrna-backed behavior in scRNA/perturb flows).
- **Solo Features**:
  - `sF` BAM tag for feature type and gene counts.
  - `--soloCBtype String` for arbitrary barcode strings.
  - Improved cell filtering and statistics with `--soloCellReadStats Standard`.

### Flex Updates
STAR-Flex extends STAR-core with Flex-specific behavior:
- **Flex Pipeline**: Inline hash-based processing for 10x Genomics Flex (Fixed RNA Profiling). Includes sample tag detection, 1MM pseudocount correction for CBs, clique-based UMI deduplication, and occupancy filtering.

### SLAM Updates
Integrated SLAM-seq quantification with GRAND-SLAM parity:
- **Quantification**: Full gene-level NTR estimation (Binomial/EM models).
- **Compatibility Mode**: `--slamCompatMode gedi` enables GEDI-compatible behaviors (intronic classification, lenient overlap, overlap weighting) for parity testing.
- **Auto-Trimming**: Variance-based detection of artifact-prone read ends (`--autoTrim variance`).
- **QC**: Comprehensive reports for T->C rates and error modeling.
- **Batch Layout + Blank-First**: `--outFileNamePrefixAuto 1` organizes SLAM outputs into `alignments/`, `counts/`, `qc/`, `y_separated/` under a single root, and `--slamErrorRateFromBlank 1` can seed the background error rate from a blank (e.g. no4sU).
- **Binary Dump + Requant**: `--slamDumpBinary 1 --slamDumpWeights 1` emits `<sample>_slam_dump.bin` and `<sample>_slam_weights.bin` in `alignments/` (batch + auto prefix layout). The `slam_requant` tool can re‑quantify these dumps with **exact parity** to `SlamQuant.out` (Pearson/Spearman 1.0 in the 1M parity check).
- **Binary Dump Format**: bitwise header + record layout is documented in `slam/docs/SLAM_DUMP_FORMAT.md`.

### STAR-perturb (CR-Compat + process_features/call_features)
STAR-suite includes a perturb-seq path that combines CR-compatible Solo behavior
with integrated CRISPR feature calling. This is the path used for STAR-perturb
work and CR compatibility comparisons.

- **Default STAR-perturb execution profile (recommended):**
  - `--runThreadN 32`
  - `--dynamicThreadInterface 1`
  - `--dynamicThreadConstMapPermits 32`
  - `--dynamicThreadTelemetry 1`
  - `--crAssignConsumerThreads 32`
  - `--crAssignSearchThreads 1`
  - This profile keeps both map and feature sides provisioned to use full CPU
    when unconstrained, with permit gating preventing oversubscription.
- **UCSF full-sample benchmark (2026-02-24, 32-thread host):**
  - STAR dynamic (`32x32`): `20:26.93`
  - STAR full fixture baseline: `30:20`
  - Cell Ranger 9 full run (`32 cores`): `33:57.18`
  - Speedup: ~`1.66x` throughput vs CR full run (~`39.8%` wall-time reduction).

- **Integrated CR-compat in STAR** (GEX + feature merge + CRISPR calling):
  - Use `--pfMultiConfig <multi_config.csv>`
  - Recommended bundle: `--defaultCrCompat yes`
  - Key controls:
    - `--crMinUmi 10` (default; lower to `2-3` for lineage-barcode style assays)
    - `--soloCrGexFeature GeneFull` (or `Gene` when explicitly required)
- **Standalone feature pipeline tool** (`core/legacy/source/star_feature_call`):
  - Full pipeline: FASTQ -> MEX -> calls
  - Call-only mode: MEX -> calls
  - `--compat-perturb` writes CR9-style `crispr_analysis/` outputs.
- **A375 small-set parity result**:
  - On the A375 1k CRISPR 5' small set, STAR CRISPR calling matched Cell Ranger at
    `1083/1083` common barcodes (**100.0% exact-match**) when using min-UMI `10`.
  - Reference report: `tests/crispr_feature_calling_comparison_report.md`.

### QC Outputs
STAR-Flex and STAR-SLAM now generate detailed QC reports:
- **SLAM QC** (`--slamQcReport <prefix>`): Generates an interactive HTML report (`.html`) and JSON metrics (`.json`) visualizing:
  - T->C conversion rates per read position.
  - Variance analysis for auto-trimming (Stdev curves, segmented regression fits).
  - Trimming overlays showing chosen 5'/3' cut sites.
- **FlexFilter QC** (`flexfilter_summary.tsv`):
  - Cell calling statistics (EmptyDrops/OrdMag results).
  - Cell counts, UMI thresholds, and filtering rates per sample.


## Summary of Flags

### Core (Legacy)
Standard STAR flags apply. See `core/legacy/README.md`.
- `--runMode`: `alignReads`, `genomeGenerate`, `soloCellFiltering`
- `--genomeDir`: Path to genome index
- `--readFilesIn`: Input read files
- `--outSAMtype`: Output SAM/BAM format (e.g., `BAM SortedByCoordinate`)
- `--batchMode`: Batch multiple FASTQs in one run (bulk, single-pass only; no Solo or 2-pass)
- `--soloType`: Single-cell mode (e.g., `CB_UMI_Simple`, `SmartSeq`)
- `--soloCbUbRequireTogether`: Enforce CB/UB tag pairing for tag injection (`yes`/`no`, default `yes`)
- `--soloCrGexFeature`: CR-compat merged GEX source (`auto`, `gene`, `genefull`)

### Flex
See `flex/README_flex.md` for full reference.
- **Pipeline**:
  - `--flex yes`: Enable Flex pipeline.
  - `--soloFlexExpectedCellsPerTag`: Expected cells per sample tag.
  - `--soloSampleWhitelist`: TSV mapping sample tags to labels.
- **Trimming**:
  - `--trimCutadapt Yes`: Enable cutadapt-style trimming.
  - `--trimCutadaptCompat`: Compatibility mode (e.g., `Cutadapt3`).
- **Quantification**:
  - `--quantMode TranscriptVB`: Enable VB/EM quantification.
- **Y-Split**:
  - `--emitNoYBAM yes`: Emit `_Y.bam` and `_noY.bam`.
  - `--emitYNoYFastq yes`: Emit split FASTQ files.
- **Reference**:
  - `--autoIndex Yes`: Enable automated reference download/build.
  - `--cellrangerStyleIndex Yes`: Use CellRanger-style reference formatting.
  - `--cellrangerLegacyGtfFilter Auto|Yes|No`: Select legacy vs updated CellRanger-style filtering mode.
- **Sorting**:
  - `--outBAMsortMethod samtools`: Enable spill-to-disk sorting.

### SLAM
See `slam/docs/SLAM_COMPATIBILITY_MODE.md` and `slam/docs/SLAM_seq.md`.
- **Quantification**:
  - `--slamQuantMode 1`: Enable SLAM quantification.
  - `--slamGrandSlamOut 1`: Generate GRAND-SLAM compatible output.
  - `--slamErrorRateFromBlank 1`: Seed error rate from the detection pass (useful when a blank is first).
- **Compatibility**:
  - `--slamCompatMode gedi`: Enable GEDI compatibility.
  - `--slamCompatIntronic`, `--slamCompatLenientOverlap`: Fine-grained control.
- **Trimming**:
  - `--autoTrim variance`: Enable variance-based auto-trimming.
  - `--slamTrim5p`, `--slamTrim3p`: Manual trim guards.
 - **Batch Layout**:
   - `--outFileNamePrefixAuto 1`: Derive sample name from first FASTQ and route outputs into subdirs under `--outFileNamePrefix`.

### STAR-perturb / CR-Compat
See `docs/feature_barcodes.md` and `docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md`.
- `--pfMultiConfig`: Enable Cell Ranger-style multi processing with feature libraries.
- `--defaultCrCompat yes`: Apply the CR-compat perturb defaults bundle.
- `--dynamicThreadInterface 1`: Enable STAR/PF permit coordination.
- `--dynamicThreadConstMapPermits 32`: Start with full map-side permit budget.
- `--crAssignConsumerThreads 32`: Provision PF worker pool to full host budget.
- `--crAssignSearchThreads 1`: Per-consumer search-thread mode used in UCSF full benchmark.
- `--crMinUmi`: Minimum UMI threshold for CRISPR feature calling (default `10`).
- `--soloCrGexFeature`: Control merged GEX source (`auto`, `gene`, `genefull`).
- `--soloCrMode CR`: Enable CR-compatible single-cell behavior.

Standalone tool (`star_feature_call`) key flags:
- `--compat-perturb`: CR9-compatible output layout (`crispr_analysis/`).
- `--feature-ref`, `--whitelist`, `--fastq-dir`, `--output-dir`: FASTQ -> MEX -> calls.
- `--call-only --mex-dir`: call_features-only pass on existing MEX.
- `--emptydrops-use-fdr`, `--min-umi`, `--ratio-test`: calling controls.

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
Notes:
- Batch mode is **single-pass only** (not compatible with `--twopassMode`).
- Batch mode is **not supported with Solo** (`--soloType`).

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
