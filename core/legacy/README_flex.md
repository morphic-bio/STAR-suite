# STAR-Flex: Extended Features for 10x Flex and Trimming

This document describes the STAR-Flex fork, which extends upstream STAR with additional features.

## Overview

STAR-Flex extends STAR with additional features:

1. **Cutadapt-style trimming** (`--trimCutadapt Yes`) for **bulk RNA-seq** with perfect parity to Trim Galore/cutadapt v5.1. For legacy datasets processed with Trim Galore + cutadapt 3.2, enable `--trimCutadaptCompat Cutadapt3` (see [docs/cutadapt_3.2_parity_report.md](docs/cutadapt_3.2_parity_report.md)). This is a general-purpose trimming feature usable with any STAR workflow (bulk RNA-seq, single-cell, etc.). See [docs/trimming.md](docs/trimming.md) for details.

2. **Inline hash-based pipeline for 10x Genomics Flex** (Fixed RNA Profiling) samples using probes for transcript detection and RTL tags for multiplexing. We generate a hybrid reference with the regular genome and with synthetic chromosomes for each of the probes. This allows to use the STAR alignment routines to quantify probe alignment and use the genomic hits to confirm the match and detect off-probe noise. However, the rest of the workflow diverges from the standard STAR solo workflow, largely due to the presence of RTL tags for multiplexing samples. Because these are on the same mate as the probe and not the cell barcode, STAR's barcode and UMI correction, and UMI deduping routines could not be used. Furthermore, the noise characteristics of Flex are different that the native STAR's multimapping ad emptyDrops functions could not be used. A fast inline path was created to handle Flex processing after STAR alignment.

The Flex pipeline includes:
- **Sample tag detection** during alignment identifies multiplexed sample barcodes
- **Inline hash capture** stores CB/UMI/gene tuples directly in memory
- **Cell Barcode (CB) correction** applies 1MM pseudocount-based correction (Cell Ranger compatible)
- **UMI correction** uses clique-based 1MM deduplication
- **Cell filtering** via OrdMag (simple EmptyDrops) or full EmptyDrops per sample
- **Tag occupancy filtering** via Monte Carlo estimation of the expected distribution of samples per cell barcode
- **MEX output** produces raw and per-sample filtered matrices

When `--flex no` (default), STAR behavior is identical to upstream.

## STAR-Flex Extras

This fork adds several features beyond upstream STAR:

### Bulk RNA-seq Features

- **[Cutadapt-Style Trimming](docs/trimming.md)**: Perfect parity with Trim Galore/cutadapt v5.1 for quality and adapter trimming. For legacy datasets trimmed with Trim Galore + cutadapt 3.2, enable `--trimCutadaptCompat Cutadapt3` (see [docs/cutadapt_3.2_parity_report.md](docs/cutadapt_3.2_parity_report.md)).

- **TranscriptVB (VB/EM) + Salmon parity workflow**: STAR-Flex can quantify transcripts with variational Bayes (VB, default) or EM, and can run Salmon in alignment mode for cross-tool parity checks.
  - **STAR TranscriptVB** (VB by default; set `--quantVBem 1` for EM):
    ```bash
    STAR \
      --runMode alignReads \
      --genomeDir /path/to/star_index \
      --readFilesIn R2.fastq.gz R1.fastq.gz \
      --readFilesCommand zcat \
      --quantMode TranscriptVB \
      --quantVBLibType A \
      --quantVBAutoDetectWindow 1000 \
      --quantVBem 0 \
      --outFileNamePrefix out/
    ```
  - **Generate transcriptome BAM for Salmon (alignment-mode quant)**:
    ```bash
    STAR \
      --runMode alignReads \
      --genomeDir /path/to/star_index \
      --readFilesIn R2.fastq.gz R1.fastq.gz \
      --readFilesCommand zcat \
      --quantMode TranscriptomeSAM \
      --outSAMtype BAM Unsorted \
      --outFileNamePrefix out_trbam/
    ```
    This produces `out_trbam/Aligned.toTranscriptome.out.bam`.
  - **Salmon VB in alignment mode** (requires `transcriptome.fa`):
    ```bash
    salmon quant \
      -t /path/to/transcriptome.fa \
      -l A \
      -a out_trbam/Aligned.toTranscriptome.out.bam \
      --useVBOpt \
      -o salmon_out
    ```
    For gene-level Salmon output, also pass a gene map: `-g tx2gene.tsv` (format `transcript_id<TAB>gene_id`).
  - **tximport-style gene counts in STAR** (DESeq2-friendly):
    ```bash
    STAR \
      --runMode alignReads \
      --genomeDir /path/to/star_index \
      --readFilesIn R2.fastq.gz R1.fastq.gz \
      --readFilesCommand zcat \
      --quantMode TranscriptVB \
      --quantVBgenesMode Tximport \
      --outFileNamePrefix out/
    ```
    This writes `out/quant.genes.tximport.sf` in the standard `quant.sf`-like schema.
  - **Testing and expected correlations**:
    - TranscriptVB vs Salmon (alignment-mode VB) parity tests: `tests/transcriptvb/salmon_parity_test.sh`
    - tximport parity tests (STAR vs R tximport, and CLI tool): `tests/tximport/star_tximport_e2e_test.sh`
    - On the JAX PE validation runs, gene-level correlations between STAR and Salmon tximport-style summaries are typically very high (Spearman/Pearson ~0.99+), and base-vs-auto (trimming) parity is expected to be ~1.0.

- **[SLAM-seq (STAR-SLAM)](docs/SLAM_seq.md)**: Integrated SLAM-seq quantification with GRAND-SLAM parity.
  - **Features**:
    - **Gene Quantification**: Full gene-level NTR estimation (Binomial/EM models).
    - **GRAND-SLAM Output**: Generates `<prefix>SlamQuant.grandslam.tsv` for seamless tool compatibility (`--slamGrandSlamOut 1`).
    - **Auto-Trimming**: Robust variance-based detection of artifact-prone read ends (`--autoTrim variance`).
    - **QC Reports**: Comprehensive HTML reports for T->C rates, error modeling, and trimming (`--slamQcReport`).
  - **Usage**:
    ```bash
    STAR \
      --runMode alignReads \
      --genomeDir /path/to/star_index \
      --readFilesIn reads.fastq.gz \
      --slamQuantMode 1 \
      --slamGrandSlamOut 1 \
      --autoTrim variance \
      --slamQcReport output/qc_report \
      --outFileNamePrefix output/
    ```
  - **Documentation**: See **[docs/SLAM_seq.md](docs/SLAM_seq.md)** for detailed theory, usage guide, and parameter reference.
  - **Benchmarks**: See [docs/SLAM_METHOD_TECHNICAL_NOTE.md](docs/SLAM_METHOD_TECHNICAL_NOTE.md) for parity benchmarks.

- **Samtools-style spill-to-disk BAM sorting** (`--outBAMsortMethod samtools`): Optional coordinate-sorting backend that uses a spill-to-disk strategy (bounded by `--limitBAMsortRAM`) to reduce temporary disk usage compared to STAR’s legacy bin-based sorter.
  - Rationale: the legacy STAR sorter partitions alignments into many genomic bins and can create large temporary files; the spill-to-disk sorter keeps in-memory buffers up to the RAM cap and only writes spill files as needed.
  - Usage:
    ```bash
    STAR \
      --runMode alignReads \
      --genomeDir /path/to/star_index \
      --readFilesIn R2.fastq.gz R1.fastq.gz \
      --readFilesCommand zcat \
      --outSAMtype BAM SortedByCoordinate \
      --outBAMsortMethod samtools \
      --outBAMsortingThreadN 8 \
      --limitBAMsortRAM 24000000000 \
      --outFileNamePrefix out/
    ```
  - Notes:
    - Default remains `--outBAMsortMethod star`.
    - Works with `--emitNoYBAM yes` / `--keepBAM yes`.

- **[Y-Chromosome BAM Split](docs/Y_CHROMOSOME_BAM_SPLIT.md)**: Split BAMs by Y-chromosome alignment. Developed for **Morphic requirements for KOLF cell lines** (not connected to Flex pipeline). Works with both bulk and single-cell workflows.

### Index-Time Features

- **[AutoIndex + CellRanger-style references](docs/autoindex_cellranger.md)**: Optional reference download + integrity verification, CellRanger-style FASTA/GTF formatting, and automatic index creation in `--genomeDir` (`--autoIndex`, `--forceIndex`, `--forceAllIndex`).

- **[Transcriptome FASTA Generation](#transcriptome-fasta-generation)**: Generate `transcriptome.fa` during index creation for Salmon quantification parity and TranscriptVB error modeling. Eliminates the need to run gffread/rsem-prepare-reference separately.

### Flex-Specific Features

- **[Flex Pipeline](docs/flex_methodology.md)**: Inline hash pipeline for 10x Genomics Flex (Fixed RNA Profiling) samples.

For complete parameter reference, see [docs/flex_parameters.md](docs/flex_parameters.md) (STAR-Flex-only flags) and upstream `README.md` (all other parameters).

For detailed technical documentation of the flex data flow and algorithms, see [docs/flex_methodology.md](docs/flex_methodology.md).

## Quick Start

```bash
STAR \
  --genomeDir /path/to/flex_reference \
  --readFilesIn R2.fastq.gz R1.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /path/to/737K-fixed-rna-profiling.txt \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist sample_whitelist.tsv \
  --soloProbeList probe_list.txt \
  --soloSampleProbes probe-barcodes-fixed-rna-profiling-rna.txt \
  --soloSampleProbeOffset 68 \
  --soloFlexOutputPrefix output/per_sample \
  --soloMultiMappers Rescue \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloFeatures Gene \
  --outFileNamePrefix output/
```

### Example: Y-Chromosome BAM Split

To split BAM output into Y and noY files:

```bash
STAR \
  --genomeDir /path/to/flex_reference \
  --readFilesIn R2.fastq.gz R1.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /path/to/737K-fixed-rna-profiling.txt \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist sample_whitelist.tsv \
  --soloProbeList probe_list.txt \
  --soloSampleProbes probe-barcodes-fixed-rna-profiling-rna.txt \
  --soloSampleProbeOffset 68 \
  --soloFlexOutputPrefix output/per_sample \
  --soloMultiMappers Rescue \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloFeatures Gene \
  --outSAMtype BAM SortedByCoordinate \
  --emitNoYBAM yes \
  --outFileNamePrefix output/
```

This produces:
- `output/Aligned.sortedByCoord.out_Y.bam` - Reads with any Y-chromosome alignment
- `output/Aligned.sortedByCoord.out_noY.bam` - Reads with no Y-chromosome alignments
- Primary BAM (`output/Aligned.sortedByCoord.out.bam`) is suppressed by default

To emit a read-name list for FASTQ filtering (with or without Y/noY BAMs):

```bash
STAR \
  ... \
  --emitYReadNames yes \
  --outFileNamePrefix output/
```

This writes `output/Aligned.out_Y.names.txt` by default (override with `--YReadNamesOutput`).

To emit Y/noY FASTQ files directly during alignment:

```bash
STAR \
  ... \
  --emitYNoYFastq yes \
  --emitYNoYFastqCompression gz \
  --outFileNamePrefix output/
```

This creates FASTQs named after the input files, with `_Y` / `_noY` inserted before the last `_R1` or `_R2`.
For example, `Sample_R1_001.fastq.gz` becomes `Sample_Y_R1_001.fastq.gz` and `Sample_noY_R1_001.fastq.gz`
(output written under the `--outFileNamePrefix` directory).
If no `_R1`/`_R2` token is found, STAR falls back to `Y_reads.mateN.fastq(.gz)` and `noY_reads.mateN.fastq(.gz)` under the output prefix.
You can override names explicitly with `--YFastqOutputPrefix` and `--noYFastqOutputPrefix`.
If a separate barcode read is present (e.g., scRNA-seq R3), only the true mates (R1/R2) are emitted.

Edge cases to be aware of:
- If the reference has no Y contigs, the Y FASTQs are empty and a warning is logged.
- FASTA inputs produce `.fa(.gz)` outputs with `>` headers and no `+`/quality lines.
- Multiple input files per mate derive output names from the first file for each mate.
- `--emitYNoYFastqCompression none` writes uncompressed `.fastq`/`.fa` outputs.
- Unmapped reads are routed to noY.

You can use `--emitYNoYFastq yes` with `--outSAMtype None` to emit FASTQ files without BAM output.

To keep the primary BAM alongside the split files:

```bash
STAR \
  ... \
  --emitNoYBAM yes \
  --keepBAM yes \
  --outFileNamePrefix output/
```

**Note**: The Y/noY split is a general-purpose feature developed for **Morphic requirements for KOLF cell lines** (not Flex-specific). It works with both single-cell RNA-seq (Flex mode) and bulk RNA-seq modes. In single-cell mode, R1/R2 are not traditional paired-end mates, so routing is based on each read's own alignments. In bulk paired-end mode, if either mate has a Y-chromosome alignment, both mates route to `_Y.bam`.

## Required Inputs

| Input | Description |
|-------|-------------|
| Flex reference genome | Hybrid genome with probe pseudo-chromosomes (see [Building References](#building-references)) |
| CB whitelist | 10x barcode whitelist (e.g., `737K-fixed-rna-profiling.txt`) |
| Sample whitelist | TSV mapping sample tag sequences to labels |
| Probe list | Gene list from probe set |
| Sample probe barcodes | 10x probe barcode sequences file |

## Parameters

### Master Switch

| Flag | Default | Description |
|------|---------|-------------|
| `--flex` | `no` | Enable flex pipeline (`yes`/`no`) |

### Y-Chromosome BAM Split

| Flag | Default | Description |
|------|---------|-------------|
| `--emitNoYBAM` | `no` | Enable Y-chromosome BAM splitting (`yes`/`no`). When enabled, emits two additional BAM files: `<out>_noY.bam` (reads with no Y-chromosome alignments) and `<out>_Y.bam` (reads with any Y-chromosome alignment). Primary BAM is suppressed by default unless `--keepBAM yes` is specified. |
| `--emitYReadNames` | `no` | Emit list of read names with any Y-chromosome alignment (one per line). Can be used with or without Y/noY BAMs. |
| `--emitYNoYFastq` | `no` | Emit Y/noY FASTQ files directly during alignment (`yes`/`no`). |
| `--emitYNoYFastqCompression` | `gz` | Compression for Y/noY FASTQ output (`gz`/`none`). |
| `--YFastqOutputPrefix` | - | Optional: override output prefix for Y FASTQ files (default: derived from input name; falls back to `Y_reads.mateN`). |
| `--noYFastqOutputPrefix` | - | Optional: override output prefix for noY FASTQ files (default: derived from input name; falls back to `noY_reads.mateN`). |
| `--keepBAM` | `no` | Keep primary BAM output when `--emitNoYBAM yes` is enabled (`yes`/`no`) |
| `--noYOutput` | - | Optional: override default path for noY BAM output (default: `<out>_noY.bam`) |
| `--YOutput` | - | Optional: override default path for Y BAM output (default: `<out>_Y.bam`) |
| `--YReadNamesOutput` | - | Optional: override output path for Y read names list (default: `<out>Aligned.out_Y.names.txt`) |

### Sample Detection

| Flag | Default | Description |
|------|---------|-------------|
| `--soloSampleWhitelist` | - | Path to sample tag whitelist TSV |
| `--soloProbeList` | auto | Path to probe gene list (auto-detects from genome index if not specified) |
| `--soloSampleProbes` | - | Path to 10x sample probe barcodes |
| `--soloSampleProbeOffset` | 0 | Offset in read for sample probe sequence |
| `--soloSampleSearchNearby` | `yes` | Search nearby positions for sample tag |
| `--soloSampleStrictMatch` | `no` | Require strict match for sample tag |

### FlexFilter (Cell Calling)

| Flag | Default | Description |
|------|---------|-------------|
| `--soloFlexExpectedCellsPerTag` | 0 | Expected cells per sample tag |
| `--soloFlexExpectedCellsTotal` | 0 | Total expected cells (alternative to per-tag) |
| `--soloFlexAllowedTags` | - | Optional: restrict to specific sample tags |
| `--soloFlexOutputPrefix` | - | Output prefix for per-sample MEX |

### EmptyDrops Parameters (Advanced)

| Flag | Default | Description |
|------|---------|-------------|
| `--soloFlexEdNiters` | 10000 | Monte Carlo simulation iterations |
| `--soloFlexEdFdrThreshold` | 0 (disabled) | FDR threshold for cell calling; if set (>0), FDR gate is used |
| `--soloFlexEdPvalueThreshold` | 0.05 | Raw p-value threshold when FDR gate is disabled (default behavior) |
| `--soloFlexEdLower` | 100 | Lower UMI bound for ambient profile |

## Output Structure

```
output/
├── Solo.out/Gene/raw/          # Raw MEX (all barcodes)
│   ├── barcodes.tsv
│   ├── features.tsv
│   └── matrix.mtx
├── per_sample/                  # Per-sample filtered MEX (labels from whitelist)
│   ├── SampleA/Gene/filtered/
│   ├── SampleB/Gene/filtered/
│   └── flexfilter_summary.tsv   # Cell calling statistics
├── Aligned.sortedByCoord.out_Y.bam      # Y-chromosome reads (if --emitNoYBAM yes)
└── Aligned.sortedByCoord.out_noY.bam    # Non-Y reads (if --emitNoYBAM yes)
```

When `--emitNoYBAM yes` is enabled:
- `_Y.bam`: Contains all reads where any alignment (primary, secondary, or supplementary) touches a Y-chromosome contig
- `_noY.bam`: Contains all reads with no Y-chromosome alignments
- Primary BAM (`Aligned.sortedByCoord.out.bam` or `Aligned.out.bam`) is suppressed by default unless `--keepBAM yes` is specified
- Works with both `BAM Unsorted` and `BAM SortedByCoordinate` output types

## Building References

The flex pipeline requires a hybrid reference genome that includes pseudo-chromosomes for probe sequences. Scripts are provided in `scripts/` to build these references:

### Integrated Index Generation (Recommended)

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/flex_index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/genes.gtf \
  --sjdbOverhang 100 \
  --flexGeneProbeSet /path/to/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv \
  --runThreadN 8
```

#### Required Inputs

| Input | Description |
|-------|-------------|
| `--genomeFastaFiles` | Base genome FASTA file |
| `--sjdbGTFfile` | Gene annotation GTF file (can be gzipped) |
| `--flexGeneProbeSet` | 10x Flex probe CSV file (50bp gene probes) |

#### Flex Index Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--flexGeneProbeSet` | - | Path to 50bp gene probe CSV file |
| `--flexGeneProbeLength` | 50 | Expected probe length (fails if mismatch) |

#### Output Structure

```
flex_index/
├── probe_gene_list.txt           # Unique gene IDs with probes (auto-detected for --soloProbeList)
├── flex_probe_artifacts/         # Probe processing artifacts
│   ├── filtered_probe_set.csv    # Probes matching GTF genes
│   ├── probes_only.fa            # Probe-only FASTA
│   ├── probes_only.gtf           # Probe-only GTF entries
│   ├── genome.filtered.fa        # Hybrid FASTA (used for indexing)
│   ├── genes.filtered.gtf        # Hybrid GTF (used for indexing)
│   ├── probe_genes_exons.bed     # Probe coordinates
│   ├── probe_list.txt            # Unique gene IDs
│   └── metadata/
│       └── reference_manifest.json
├── Genome                        # Standard STAR index files
├── SA
├── SAindex
└── ... (other STAR index files)
```

#### Probe Filtering Rules

The integrated preprocessor applies these filters:
1. **50bp A/C/G/T only** - Fails if any probe has invalid length or characters
2. **Skip DEPRECATED** - Excludes probes marked as deprecated
3. **Gene match** - Keeps only probes whose gene_id exists in the target GTF
4. **Deterministic ordering** - Stable sort by gene_id then probe_id

### Alternative: Shell Scripts

For custom workflows or debugging, standalone shell scripts are available:

```bash
# Filter probes and build hybrid reference
./scripts/filter_probes_to_gtf.sh \
  --probe-set /path/to/probes.csv \
  --gtf /path/to/genes.gtf.gz \
  --base-fasta /path/to/genome.fa \
  --output-dir ./probe_artifacts
```

The legacy `build_filtered_reference.sh` and `make_filtered_star_index.sh` scripts are also available. See [scripts/README.md](scripts/README.md) for details.

### Using the Flex Index

After building, use the index with the probe gene list:

```bash
STAR \
  --genomeDir /path/to/flex_index \
  --flex yes \
  ... # other flex parameters
  # --soloProbeList is auto-detected from probe_gene_list.txt in the index directory
```

## AutoIndex + CellRanger-Style References

STAR-Flex includes an index-time workflow to reproduce the “CellRanger-style” reference preparation (download → integrity checks → format FASTA/GTF → genomeGenerate).

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --autoIndex Yes \
  --cellrangerStyleIndex Yes \
  --autoCksumUpdate Yes \
  --sjdbOverhang 100 \
  --runThreadN 16
```

Key outputs and paths:
- Formatted inputs: `${genomeDir}/cellranger_ref/genome.fa`, `${genomeDir}/cellranger_ref/genes.gtf`
- Download cache (default): `${genomeDir}/cellranger_ref_cache` (override with `--cellrangerStyleCacheDir`)
- Rebuild controls: `--forceIndex Yes` (re-index), `--forceAllIndex Yes` (re-download + re-index)

See [docs/autoindex_cellranger.md](docs/autoindex_cellranger.md) for URL selection (`--cellrangerRefRelease` / `--faUrl` / `--gtfUrl`), checksum flags, and parity test scripts.

## Transcriptome FASTA Generation

STAR-Flex can generate `transcriptome.fa` during index creation, eliminating the need for separate gffread/rsem-prepare-reference runs. This is required for:
- **Salmon quantification** (identical output for parity)
- **TranscriptVB error modeling** (fragment length distribution estimation)

### Basic Usage

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/genes.gtf \
  --sjdbOverhang 100 \
  --genomeGenerateTranscriptome Yes \
  --runThreadN 8
```

This produces `${genomeDir}/transcriptome.fa` alongside the standard index files.

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--genomeGenerateTranscriptome` | `No` | Enable transcriptome FASTA generation (`Yes`/`No`) |
| `--genomeGenerateTranscriptomeFasta` | `-` | Custom output path (default: `${genomeDir}/transcriptome.fa`) |
| `--genomeGenerateTranscriptomeOverwrite` | `No` | Overwrite existing file (`Yes`/`No`) |

### CellRanger-Style Index

When `--cellrangerStyleIndex Yes`, STAR-Flex formats the annotation inputs into `${genomeDir}/cellranger_ref/`:
- `${genomeDir}/cellranger_ref/genome.fa`
- `${genomeDir}/cellranger_ref/genes.gtf`

When combined with `--genomeGenerateTranscriptome Yes`, the transcriptome is written to both:
- `${genomeDir}/transcriptome.fa` (standard path)
- `${genomeDir}/cellranger_ref/transcriptome.fa` (CellRanger-compatible path)

```bash
STAR --runMode genomeGenerate \
  --genomeDir /path/to/index \
  --genomeFastaFiles /path/to/genome.fa \
  --sjdbGTFfile /path/to/genes.gtf \
  --sjdbOverhang 100 \
  --genomeGenerateTranscriptome Yes \
  --cellrangerStyleIndex Yes \
  --runThreadN 8
```

### Output Format

The transcriptome FASTA follows Salmon conventions:
- **Headers**: Transcript IDs without version suffixes (e.g., `>ENST00000456328` not `>ENST00000456328.2`)
- **Line width**: 70 characters
- **Ordering**: Matches `transcriptInfo.tab` for Salmon parity
- **Negative strand**: Exons concatenated in genomic order, then reverse-complemented

### Verification

Test with the included chr21+chr22 subset:

```bash
./test/run_transcriptome_generation.sh --all
```

This runs:
1. **Synthetic tests**: Basic transcriptome generation with small fixtures
2. **Default path tests**: Validates `${genomeDir}/transcriptome.fa` output
3. **CellRanger tests**: Real GENCODE chr21+chr22 with CellRanger filtering

## Standalone FlexFilter Tool

A standalone tool `run_flexfilter_mex` is available for offline MEX processing. This allows re-running the OrdMag/EmptyDrops cell calling pipeline on existing composite MEX files without re-running STAR alignment.

**Use cases:**
- Parameter tuning (adjust expected cells, EmptyDrops thresholds)
- Reprocessing with different filtering settings
- Integration with non-STAR pipelines (any tool producing composite CB+TAG MEX)
- Batch reprocessing of archived STAR outputs

### Building

The tool is optional and not built by the default `make STAR` target:

```bash
cd source
make flexfilter
```

This produces `tools/flexfilter/run_flexfilter_mex`.

### Input Requirements

The tool expects a composite MEX directory containing:
- `matrix.mtx` - Matrix Market sparse matrix (or `InlineHashDedup_matrix.mtx`)
- `barcodes.tsv` - Composite barcodes in CB16+TAG8 format (24 characters)
- `features.tsv` - Gene IDs (tab-separated)

The composite barcode format concatenates the 16bp cell barcode with the 8bp sample tag:
```
AAACCCAAGAAACACTACGTACGT  # CB16 (AAACCCAAGAAACACT) + TAG8 (ACGTACGT)
```

### Basic Usage

```bash
./tools/flexfilter/run_flexfilter_mex \
  --mex-dir /path/to/Solo.out/Gene/raw \
  --total-expected 12000 \
  --output-prefix /path/to/filtered_output
```

### Key Parameters

| Parameter | Description |
|-----------|-------------|
| `--mex-dir` | Path to composite MEX directory (required) |
| `--total-expected` | Total expected cells across all samples (required) |
| `--output-prefix` | Output directory prefix (required) |
| `--sample-whitelist` | TSV file mapping sample names to tag sequences |
| `--ed-lower-bound` | Lower UMI bound for EmptyDrops (default: 500) |
| `--ed-fdr` | FDR threshold for EmptyDrops (default: 0.01) |
| `--disable-occupancy` | Skip occupancy post-filter (for testing) |

### Output Structure

```
output_prefix/
├── SampleA/Gene/filtered/
│   ├── matrix.mtx
│   ├── barcodes.tsv
│   ├── features.tsv
│   └── EmptyDrops/
│       └── emptydrops_results.tsv
├── SampleB/Gene/filtered/
│   └── ...
└── flexfilter_summary.tsv
```

### Example Workflow: Reprocess with Different Expected Cells

```bash
# Original STAR run produced Solo.out/Gene/raw/
# Reprocess with higher cell expectation
./tools/flexfilter/run_flexfilter_mex \
  --mex-dir /storage/run1/Solo.out/Gene/raw \
  --total-expected 20000 \
  --output-prefix /storage/run1/refiltered_20k

# Or with explicit sample whitelist
./tools/flexfilter/run_flexfilter_mex \
  --mex-dir /storage/run1/Solo.out/Gene/raw \
  --sample-whitelist samples.tsv \
  --total-expected 15000 \
  --output-prefix /storage/run1/refiltered_explicit
```

Sample whitelist format (`samples.tsv`):
```
Sample_A	ACGTACGT
Sample_B	TGCATGCA
Sample_C	GGCCGGCC
```
Labels in the first column are used verbatim for per-sample directories, and the order in the whitelist is preserved.

### Testing

```bash
# Requires tests/gold_standard/ fixtures
./tools/flexfilter/test_smoke.sh

# Validate output format
./tools/flexfilter/validate_output.py /path/to/output
```

See [tools/flexfilter/README.md](tools/flexfilter/README.md) for complete CLI reference and advanced options.

## remove_y_reads - FASTQ Splitter

A standalone C tool that splits FASTQ files based on a Y-only BAM produced by STAR's `--emitNoYBAM` feature. Given the `_Y.bam` output from STAR's Y-chromosome split, this tool partitions original FASTQ files into Y/noY sets while preserving read order.

**Use cases:**
- Split FASTQ files after STAR alignment with Y/noY BAM output
- Prepare separate inputs for sex-specific analyses
- Filter out Y-chromosome reads from FASTQ files
- Downstream analysis requiring separate Y/non-Y FASTQs

**Key features:**
- Uses htslib for BAM reading and kseq.h for robust FASTQ parsing
- Dual-hash collision protection (FNV-1a + djb2) for read name lookup
- File-level threading with semaphore-bounded concurrency
- Preserves original read order in outputs
- Handles gzipped and uncompressed FASTQs

### Building

The tool is optional and not built by the default `make STAR` target:

```bash
cd source
make remove_y_reads
```

This produces `tools/remove_y_reads/remove_y_reads`.

Alternatively, build directly:

```bash
cd tools/remove_y_reads
make
```

### Basic Usage

```bash
./tools/remove_y_reads/remove_y_reads \
    -y Aligned.sortedByCoord.out_Y.bam \
    --threads 4 \
    --gzip-level 6 \
    -o output_dir \
    sample_R1.fastq.gz sample_R2.fastq.gz
```

**Output:** For each input FASTQ, produces `<stem>_Y.fastq.gz` and `<stem>_noY.fastq.gz`.

### Key Parameters

| Flag | Description |
|------|-------------|
| `-y, --ybam` | Y-only BAM file (required) |
| `-o, --outdir` | Output directory (default: alongside input) |
| `-t, --threads` | Number of parallel workers (default: 1) |
| `-z, --gzip-level` | Compression level 1-9 (default: 6) |
| `-h, --help` | Show help message |

### Features

- **Read order preservation**: Outputs maintain the same read order as input FASTQs
- **Name normalization**: Automatically handles FASTQ name formats (strips `@`, `/1`, `/2`, comments)
- **Collision detection**: Uses hash + length to guard against rare hash collisions
- **Multi-threaded**: Process multiple FASTQ files in parallel (file-level parallelism)
- **Gzip support**: Handles both compressed and uncompressed FASTQ files
- **Dynamic parsing**: Uses kseq.h for robust parsing of arbitrarily long reads

### Example Workflow

```bash
# Step 1: Run STAR with Y/noY split
STAR \
  --genomeDir /path/to/reference \
  --readFilesIn R1.fastq.gz R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --emitNoYBAM yes \
  --outFileNamePrefix output/

# Step 2: Split original FASTQs based on Y BAM
./tools/remove_y_reads/remove_y_reads \
    -y output/Aligned.sortedByCoord.out_Y.bam \
    --threads 4 \
    -o output/fastq_split \
    R1.fastq.gz R2.fastq.gz

# Result: output/fastq_split/R1_Y.fastq.gz, R1_noY.fastq.gz, etc.
```

### Testing

```bash
# Basic self-contained test
./tests/run_remove_y_reads_test.sh

# Comprehensive test (single-threaded, multithreaded, multiple files)
./tests/run_y_removal_comprehensive_test.sh
```

Test report generated at `tests/TEST_REPORT_REMOVE_Y_FASTQ.md`.

For detailed technical documentation, see [docs/Y_CHROMOSOME_BAM_SPLIT.md](docs/Y_CHROMOSOME_BAM_SPLIT.md).

## Building STAR-Flex

Standard STAR build process:

```bash
cd source
make -j8
```

The flex objects are automatically included in the build.

## Testing

See [docs/TESTING_flex.md](docs/TESTING_flex.md) for detailed testing instructions.

Quick test:

```bash
./tests/run_flex_multisample_test.sh
```

Gold standard comparison files are bundled in `tests/gold_standard/`.

## Code Organization

```
source/
├── libflex/                    # Core flex filtering library
│   ├── FlexFilter.cpp/h        # Main filter orchestration
│   ├── EmptyDropsMultinomial.cpp/h  # Full EmptyDrops
│   ├── OrdMagStage.cpp/h       # Simple EmptyDrops (OrdMag)
│   └── OccupancyGuard.cpp/h    # Occupancy-based filtering
├── solo/
│   └── CbCorrector.cpp/h       # CB correction with pseudocounts
├── SampleDetector.cpp/h        # Sample tag detection
├── InlineCBCorrection.cpp/h    # Inline hash CB correction
├── UMICorrector.cpp/h          # Clique-based UMI correction
├── MexWriter.cpp/h             # MEX matrix output
├── GeneResolver.cpp/h          # Probe-to-gene mapping
├── SoloFeature_flexfilter.cpp  # FlexFilter integration
├── SoloFeature_writeMexFromInlineHashDedup.cpp
└── UmiCodec.h                  # UMI encoding/decoding helpers

tools/flexfilter/               # Standalone FlexFilter CLI
├── run_flexfilter_mex.cpp      # Main CLI wrapper
├── Makefile                    # Build configuration
├── README.md                   # CLI documentation
├── test_smoke.sh               # Smoke test script
└── validate_output.py          # Output validation script

tools/remove_y_reads/           # Standalone FASTQ Y-splitter CLI
├── remove_y_reads.c            # Main implementation (C + htslib)
└── Makefile                    # Build configuration
```

## Compatibility

- Baseline: STAR 2.7.11b
- When `--flex no` (default), behavior is identical to upstream STAR
- Upstream `README.md` and `CHANGES.md` are not modified

## Code Statistics

STAR-Flex adds approximately **33,800 lines of C/C++ code** on top of upstream STAR 2.7.11b.

### Summary (excluding third-party libraries)

| Repository | Files | Lines |
|------------|-------|-------|
| Upstream STAR (2.7.11b) | 250 | 28,228 |
| STAR-Flex | 356 | 62,016 |
| **Difference** | **+106** | **+33,788** |

### Bulk/PE vs Flex (heuristic file classification)

| Category | Files | Lines |
|----------|-------|-------|
| Flex-specific | 25 | 9,670 |
| Bulk/PE-specific | 44 | 10,993 |
| Shared/core (includes overlap) | 287 | 41,353 |

Bulk/PE-specific includes (among others) `source/libtrim/` (trimming) and `source/libem/` (TranscriptVB/EM + tximport-style gene summaries).

### Key Bulk/PE Modules (selected)

| Path | Files | Lines |
|------|-------|-------|
| `source/libtrim/` | 4 | 962 |
| `source/libem/` | 26 | 4,873 |
| `source/TranscriptomeFasta.*` | 2 | 231 |
| `source/SamtoolsSorter.*` | 2 | 584 |
| `source/CellRangerFormatter.*` | 2 | 1,863 |

### Third-party Libraries (included in repo)

| Library | Lines | Files | Notes |
|---------|-------|-------|-------|
| PCG random (pcg_*.hpp) | 3,623 | 3 | New in STAR-Flex |
| klib/khash.h | 617 | 1 | New in STAR-Flex |
| SimpleGoodTuring | 300 | 1 | Present in upstream |

*Line counts include only compiled source files (`.c`, `.cpp`, `.h`, `.hpp`) under `source/`. Excludes `source/htslib/`, `source/opal/`, reference data, scripts, and test fixtures. Generated by `scripts/count_cpp_lines.sh --upstream-ref 2.7.11b`.*
