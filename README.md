# STAR-suite

STAR-suite reorganizes STAR into module-focused directories while keeping a single source of truth for shared code. Build outputs remain compatible with existing STAR workflows, and the new top-level `Makefile` exposes module targets.

Agent quickstart: see `AGENTS.md` for repo-specific guardrails, tests, and recent changes.

## Folder Structure

```
core/
  legacy/                # Upstream STAR layout (single source of truth)
  features/              # Shared feature overlays (vbem, yremove, bamsort)
flex/                    # Flex-specific code + tools
slam/                    # SLAM-seq code + tools
build/                   # Modular make fragments
docs/                    # Suite-level docs
  tests/                   # Suite-level tests (see tests/ARTIFACTS.md for artifact locations)
tools/                   # Suite-level scripts/utilities
```

## Modules

- **STAR-core** (`core/`): Legacy STAR (indexing, bulk, Solo) plus shared utilities.
  Build: `make core` (binary at `core/legacy/source/STAR`).
- **STAR-Flex** (`flex/`): FlexFilter pipeline and Flex-specific integrations.
  Build tools: `make flex` or `make flex-tools`.
- **STAR-SLAM** (`slam/`): SLAM-seq quantification, SNP masking, trimming/QC.
  Build tools: `make slam` or `make slam-tools`.
- **Feature Barcodes** (`core/features/feature_barcodes/`): Vendored `process_features` tools for perturb-seq testing (`assignBarcodes`, `demux_bam`, `demux_fastq`).
  Build tools: `make feature-barcodes-tools`.

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

## Technical Updates

### Core Updates
Recent updates to the Core module (STAR 2.7.11b and prior) include:
- **Transcriptome Output**: Replaced `--quantTranscriptomeBan` with `--quantTranscriptomeSAMoutput` for more explicit control (e.g., `BanSingleEnd_ExtendSoftclip`).
- **Solo Features**:
  - `sF` BAM tag for feature type and gene counts.
  - `--soloCBtype String` for arbitrary barcode strings.
  - Improved cell filtering and statistics with `--soloCellReadStats Standard`.

### Flex Updates
STAR-Flex extends upstream STAR with features for 10x Flex, bulk RNA-seq, and reference management:
- **Flex Pipeline**: Inline hash-based processing for 10x Genomics Flex (Fixed RNA Profiling). Includes sample tag detection, 1MM pseudocount correction for CBs, clique-based UMI deduplication, and occupancy filtering.
- **Cutadapt Parity**: Native implementation of cutadapt v5.1 trimming algorithm (`--trimCutadapt Yes`) for bulk RNA-seq.
- **TranscriptVB**: Variational Bayes and EM quantification for transcript-level abundance, offering parity with Salmon alignment-mode.
- **Reference Building**:
  - **AutoIndex**: Automated reference download (`--autoIndex`), checksum verification (`--autoCksumUpdate`), and CellRanger-style formatting (`--cellrangerStyleIndex`).
  - **Transcriptome FASTA**: Automatic generation of `transcriptome.fa` during indexing for Salmon/TranscriptVB compatibility.
- **BAM Sorting**:
  - **Samtools-style Sort**: Spill-to-disk sorting (`--outBAMsortMethod samtools`) reduces memory usage compared to the default bin-based sorter. Configurable with `--limitBAMsortRAM`.
- **Y-Chromosome Splitting**: Split BAMs and FASTQs by Y-chromosome alignment (`--emitNoYBAM`, `--emitYNoYFastq`) for sex-specific analysis.

### SLAM Updates
Integrated SLAM-seq quantification with GRAND-SLAM parity:
- **Quantification**: Full gene-level NTR estimation (Binomial/EM models).
- **Compatibility Mode**: `--slamCompatMode gedi` enables GEDI-compatible behaviors (intronic classification, lenient overlap, overlap weighting) for parity testing.
- **Auto-Trimming**: Variance-based detection of artifact-prone read ends (`--autoTrim variance`).
- **QC**: Comprehensive reports for T->C rates and error modeling.
- **Batch Layout + Blank-First**: `--outFileNamePrefixAuto 1` organizes SLAM outputs into `alignments/`, `counts/`, `qc/`, `y_separated/` under a single root, and `--slamErrorRateFromBlank 1` can seed the background error rate from a blank (e.g. no4sU).
- **Batch Mode (Core)**: The multi-FASTQ “batch mode” originated in SLAM-seq workflows with a blank sample, and has been generalized to bulk single-pass runs. Use `--batchMode 1` to process multiple FASTQs in one invocation while reusing the loaded genome. **Exclusions**: not supported with Solo (`--soloType`) or 2-pass (`--twopassMode`). Output is routed into per-sample subdirectories via `--outFileNamePrefixAuto 1` (forced in batch mode).
- **Binary Dump + Requant**: `--slamDumpBinary 1 --slamDumpWeights 1` emits `<sample>_slam_dump.bin` and `<sample>_slam_weights.bin` in `alignments/` (batch + auto prefix layout). The `slam_requant` tool can re‑quantify these dumps with **exact parity** to `SlamQuant.out` (Pearson/Spearman 1.0 in the 1M parity check).
- **Binary Dump Format**: bitwise header + record layout is documented in `slam/docs/SLAM_DUMP_FORMAT.md`.

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

## More Detail

- Core usage: [core/legacy/README.md](core/legacy/README.md)
- Flex pipeline: [flex/README_flex.md](flex/README_flex.md)
- SLAM compatibility: [slam/docs/SLAM_COMPATIBILITY_MODE.md](slam/docs/SLAM_COMPATIBILITY_MODE.md)
- SLAM methodology: [slam/docs/SLAM_seq.md](slam/docs/SLAM_seq.md)
- Cell Ranger multi smoke tool: [docs/cr_multi.md](docs/cr_multi.md)
