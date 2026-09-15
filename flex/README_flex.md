# STAR-Flex: 10x Fixed RNA Profiling Pipeline

This document describes STAR-Flex, the Flex-specific module in STAR Suite.

> **STAR Suite 1.9.4:** `--flex yes` assigns reads to probes from a half-probe (H1X2) hash cache and aligns nothing; no genome index is loaded. The cache is required; build it once per probe set with `--runMode hashCacheGenerate --hashCacheTiers H0,H1X2`. The alignment-based routes described in parts of this document (pseudo-chromosome alignment, the alignment-validated H0/H1 cache, BAM output with CB/UB tags, Y-chromosome splitting) are **legacy** and need `--flexLegacy yes`; use them only to reproduce results from an earlier release. Current speed and concordance numbers are in the top-level [README Benchmarks](../README.md#benchmarks).

## Overview

STAR-Flex adds a **probe-cache pipeline for 10x Genomics Flex** (Fixed RNA Profiling) samples using probes for transcript detection and RTL tags for multiplexing. Each read is assigned to a probe from a precomputed half-probe (H1X2) cache. An exact 50-base match is resolved directly. Otherwise each 25-base half of the probe window is looked up in exact and single-mismatch tables, the other half is scored by fast Hamming distance, and the read is accepted when the whole probe has at most 10 mismatches. A read anchored by more than one probe, or otherwise ambiguous, is rejected. No read is aligned. The cache is built once from a Flex reference index (the genome plus one pseudo-chromosome per probe), and cache entries that could come from more than one probe are marked ambiguous when it is built. Aligning reads against that hybrid reference is the legacy route. The rest of the workflow diverges from the standard STAR Solo workflow because RTL tags are on the same mate as the probe (not the cell barcode), so STAR's barcode/UMI correction and deduplication routines cannot be used. A fast inline path handles Flex-specific processing.

The Flex pipeline includes:
- **Sample tag detection** identifies multiplexed sample barcodes
- **Inline hash capture** stores CB/UMI/gene tuples directly in memory, in cell-barcode buckets
- **Cell Barcode (CB) correction** applies 1MM pseudocount-based correction (Cell Ranger compatible)
- **UMI correction** uses clique-based 1MM deduplication
- **Cell calling** with the tag-aware caller (`--soloFlexCellCaller tag-aware`, the default): OrdMag and EmptyDrops run once per sample, and tags given the same sample label are called together as one group. A Flex cell is its 16-base barcode together with its sample tag (CB16+TAG8)
- **GEM occupancy filtering** after all samples are called: the number of distinct sample tags seen per 16-base barcode is fitted with a zero-truncated Poisson, and barcodes above the 0.999 quantile are removed
- **MEX output** produces raw and per-sample filtered matrices

When `--flex no` (default), STAR behavior is identical to upstream.

## Core Features Available in Flex Mode

The following features were originally developed in the STAR-Flex fork and are now
part of STAR-core. They work with all STAR modes (bulk, single-cell, Flex). See the
main suite `README.md` for full documentation and flags.

- **Cutadapt-style trimming** (`--trimCutadapt Yes`): See [trimming docs](../core/features/vbem/docs/trimming.md).
- **TranscriptVB quantification** (`--quantMode TranscriptVB`): VB/EM transcript-level quantification with Salmon parity.
- **SLAM-seq** (`--slamQuantMode 1`): See [slam/docs/SLAM_seq.md](../slam/docs/SLAM_seq.md).
- **Spill-to-disk BAM sorting** (`--outBAMsortMethod samtools`): Bounded-RAM coordinate sorting. Works with the legacy Flex alignment route (`--flexLegacy yes`); the default Flex route writes no BAM.
- **Y-chromosome BAM/FASTQ splitting** (`--emitNoYBAM yes`, `--emitYNoYFastq yes`): Split reads by chrY alignment. Developed for MorPHiC KOLF cell lines. Tested and validated with Flex in both sorted and unsorted modes (see `tests/TEST_REPORT_Y_SPLIT_FLEX.md`); for Flex this needs alignment and therefore the legacy route (`--flexLegacy yes`). See [Y-chromosome BAM split docs](../core/features/yremove_fastq/docs/Y_CHROMOSOME_BAM_SPLIT.md).

## STAR-Flex Extras

### Index-Time Features

- **[AutoIndex + CellRanger-style references](../core/legacy/docs/autoindex_cellranger.md)**: Optional reference download + integrity verification, CellRanger-style FASTA/GTF formatting, and automatic index creation in `--genomeDir` (`--autoIndex`, `--forceIndex`, `--forceAllIndex`).

- **[Transcriptome FASTA Generation](#transcriptome-fasta-generation)**: Generate `transcriptome.fa` during index creation for Salmon quantification parity and TranscriptVB error modeling. Eliminates the need to run gffread/rsem-prepare-reference separately.

### Flex-Specific Features

- **[Flex Pipeline](docs/flex_methodology.md)**: Inline hash pipeline for 10x Genomics Flex (Fixed RNA Profiling) samples.

For complete parameter reference, see [flex parameter docs](../core/legacy/docs/flex_parameters.md) (STAR-Flex-only flags) and upstream `README.md` (all other parameters).

For detailed technical documentation of the flex data flow and algorithms, see [docs/flex_methodology.md](docs/flex_methodology.md).

## Quick Start

Build the half-probe cache once per probe set, from a Flex reference index (see [Building References](#building-references)):

```bash
STAR --runMode hashCacheGenerate \
  --runThreadN 16 \
  --genomeDir /path/to/flex_index \
  --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
  --soloCBwhitelist /path/to/737K-fixed-rna-profiling.txt \
  --flex yes --soloFeatures Gene \
  --soloProbeList probe_list.txt \
  --soloSampleWhitelist sample_whitelist.tsv \
  --soloSampleProbes probe-barcodes-fixed-rna-profiling-rna.txt \
  --soloSampleProbeOffset 68 \
  --hashCacheTiers H0,H1X2 \
  --hashCacheOutput /path/to/flex_h01x2_sequence_cache.bin \
  --outSAMtype None \
  --outFileNamePrefix cache_run/
```

Cache generation is memory-hungry: for the human probe set it peaked at 74.5 GiB with 16 threads (see `docs/RELEASE_NOTES_v1.9.4.md`). The FH01SEQ1 `.bin` it writes is accepted directly; `flex/tools/hash_screen_replay/flex_hash_cache_pack --half` converts it to the compact `.half.khash` form (see [the cache format notes](../docs/FLEX_KHASH_CACHE.md)).

Then run Flex. This is the command used for the STAR Suite 1.9.4 benchmarks (see `docs/PAPER_BENCHMARK_METHODOLOGY.md`, Section 1.6):

```bash
STAR \
  --runThreadN 32 --flex yes --flexPipeline yes --flexNoAlign 1 \
  --readFilesIn R2.fastq.gz R1.fastq.gz \
  --soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17 --soloCBlen 16 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 --soloCBwhitelist /path/to/737K-fixed-rna-profiling.txt \
  --soloSampleWhitelist sample_whitelist.tsv --soloFlexAllowedTags sample_whitelist.tsv \
  --soloSampleProbes probe-barcodes-fixed-rna-profiling-rna.txt --soloSampleProbeOffset 68 \
  --soloProbeList probe_list.txt \
  --soloHashScreenFile /path/to/flex_h01x2_sequence_cache.bin \
  --soloFeatures Gene --soloCellFilter None --soloMultiMappers Rescue \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloStrand Unstranded --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 0 \
  --soloKeysCompat cr --soloBucketMode ram --soloBucketCount 256 \
  --soloRunFlexFilter yes --soloFlexCellCaller tag-aware --soloCellFilterBootstrapThreads 32 \
  --soloFlexEdFdrThreshold 0.01 --outSAMtype None --outSJtype None --dynamicThreadInterface 1 \
  --outFileNamePrefix output/
```

Choose the input reader to match the files:

- Sequencer-delivered BGZF FASTQ: `--readFilesBgzfMode range --bgzfReaderThreads 32 --bgzfCrcCheck 1` (parallel in-process reader).
- Plain single-stream gzip: `--readFilesBgzfMode off --readFilesCommand rapidgzip -d -c -P 8` (rapidgzip is an external program).
- CBQ (experimental): `--readFilesType Binseq PE --readFilesCbqRangeMode range`.

No genome index is loaded when the command is count-only: `--outSAMtype None`, `--outSJtype None`, `--chimSegmentMin 0`, `--soloFeatures Gene` and no other genome-backed output. Otherwise `Log.out` reports `Flex count-only no-genome: not active` with the reason, and STAR loads the index given by `--genomeDir`.

`--soloBucketMode ram` keeps the cell-barcode buckets in memory; the default `auto` switches to spill files once `--soloBucketMemGB` (32 GiB) is crossed. With buckets in memory, peak memory in the 1.9.4 benchmarks was 23 GiB on JAX SC2300771 and 72 GiB on the 10x 320k scFFPE dataset. EmptyDrops draws are split across bootstrap streams, so keep `--soloCellFilterBootstrapThreads` fixed to reproduce a call set exactly. Legacy expected-cell options such as `--soloFlexExpectedCellsPerTag` are rejected by the tag-aware caller.

The sample whitelist has one `label<TAB>TAG8` line per tag. The cell caller takes its sample labels from the `--soloFlexAllowedTags` file (the same format), so tags pooled into one sample are given the same label and are called together as one group; for example, the 10x 320k scFFPE benchmark pooled its 16 tags in pairs into 8 samples.

### Example: Y-Chromosome BAM Split (legacy route)

Flex BAM output needs genomic alignment, so it runs only on the legacy route (`--flexLegacy yes`) with a Flex reference index. To split BAM output into Y and noY files:

```bash
STAR \
  --genomeDir /path/to/flex_reference \
  --readFilesIn R2.fastq.gz R1.fastq.gz \
  --readFilesCommand zcat \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist /path/to/737K-fixed-rna-profiling.txt \
  --flex yes \
  --flexLegacy yes \
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

**Note**: The Y/noY split is a general-purpose core feature developed for **MorPHiC requirements for KOLF cell lines**. It works with all modes: Flex (legacy alignment route only), single-cell, and bulk RNA-seq. Validated with Flex in both sorted and unsorted modes (see `tests/TEST_REPORT_Y_SPLIT_FLEX.md`). In single-cell mode, R1/R2 are not traditional paired-end mates, so routing is based on each read's own alignments. In bulk paired-end mode, if either mate has a Y-chromosome alignment, both mates route to `_Y.bam`.

## Required Inputs

| Input | Description |
|-------|-------------|
| Half-probe cache | H0,H1X2 probe cache passed with `--soloHashScreenFile` (see [Quick Start](#quick-start)); required by `--flex yes` |
| Flex reference index | Hybrid genome with probe pseudo-chromosomes (see [Building References](#building-references)); needed once to build the cache, and at run time only for the legacy alignment route |
| CB whitelist | 10x barcode whitelist (e.g., `737K-fixed-rna-profiling.txt`) |
| Sample whitelist | TSV with one `label<TAB>TAG8` line per sample tag |
| Probe list | Gene list from probe set; must be the list the cache was built with |
| Sample probe barcodes | 10x probe barcode sequences file |

## Parameters

### Master Switch

| Flag | Default | Description |
|------|---------|-------------|
| `--flex` | `no` | Enable flex pipeline (`yes`/`no`). `yes` selects the half-probe route: `--flexNoAlign 1`, fully fused pipeline, `--outSAMtype None` unless set explicitly |
| `--soloHashScreenFile` | - | Half-probe (H1X2) cache. If omitted, STAR uses `STAR_FLEX_HASH_SCREEN_CACHE`, then looks next to the probe list for `flex_h01x2_cache.half.khash` and then `flex_h01x2_sequence_cache.bin` |
| `--flexLegacy` | `no` | `yes` permits the legacy routes (H0/H1 cache without H1X2, aligning cache misses, `--no-hash-screen yes`, SAM/BAM output, spatial Flex), kept only to reproduce earlier releases |

### Y-Chromosome BAM Split

For Flex these need alignment and therefore `--flexLegacy yes`.

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
| `--soloSampleProbeOffset` | 68 | The one fixed offset used for the 8-base sample tag |
| `--soloSampleTagMismatch` | 1 | After an exact-table miss, query a separately constructed Hamming-1 cache; accept only unique-owner keys |
| `--soloSampleSearchNearby` | `no` | Compatibility flag; `yes` is rejected because neighboring offsets are not searched |
| `--soloSampleStrictMatch` | `no` | Set `yes` to disable the H1 tier and require an exact table entry |

See [the fixed sample-tag policy](../docs/FLEX_SAMPLE_TAG_POLICY.md) for the
lookup order and ambiguity rules.

### FlexFilter (Cell Calling)

| Flag | Default | Description |
|------|---------|-------------|
| `--soloRunFlexFilter` | `no` | Run cell calling after the raw MEX is written; `--flex yes` sets `yes` unless given explicitly |
| `--soloFlexCellCaller` | `tag-aware` | `tag-aware`: one OrdMag + EmptyDrops model per sample label (pooled tags together), full CB16+TAG8 identities, joint occupancy fit. `legacy`: the earlier per-tag EmptyDrops + partition-occupancy caller |
| `--soloFlexAllowedTags` | - | Optional `label<TAB>TAG8` (or TAG8-only) file; restricts calling to these tags and supplies the sample labels used for grouping |
| `--soloFlexOutputPrefix` | - | Output prefix for per-sample MEX |
| `--soloCellFilterBootstrapThreads` | 0 | Bootstrap streams (0 = `--runThreadN`); the call set depends on this value |
| `--soloFlexFilteredGeneList` | - | Optional gene-ID allowlist for the per-sample filtered MEX; does not restrict calling |
| `--soloFlexExpectedCellsPerTag` | 0 | Legacy caller only: expected cells per sample tag (rejected by `tag-aware`) |
| `--soloFlexExpectedCellsTotal` | 0 | Legacy caller only: total expected cells (rejected by `tag-aware`) |

### EmptyDrops Parameters (Advanced)

| Flag | Default | Description |
|------|---------|-------------|
| `--soloFlexEdNiters` | 0 | Monte Carlo simulation iterations; 0 = 100,000 for `tag-aware`, 10,000 for `legacy` |
| `--soloFlexEdFdrThreshold` | 0 | FDR threshold for cell calling; 0 = 0.01 for `tag-aware`, 0.001 for `legacy` |
| `--soloFlexEdLower` | 0 | Legacy caller only: start of the ambient rank window (0 = 45,000). The tag-aware caller uses a fixed window of ranks 45,000-90,000 per tag and rejects this option |

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
├── Aligned.sortedByCoord.out_Y.bam      # Y-chromosome reads (legacy route, if --emitNoYBAM yes)
└── Aligned.sortedByCoord.out_noY.bam    # Non-Y reads (legacy route, if --emitNoYBAM yes)
```

With the tag-aware caller, barcodes in the per-sample filtered MEX keep the full CB16+TAG8 identity, so cells from tags pooled into one sample stay distinct.

When `--emitNoYBAM yes` is enabled (legacy route):
- `_Y.bam`: Contains all reads where any alignment (primary, secondary, or supplementary) touches a Y-chromosome contig
- `_noY.bam`: Contains all reads with no Y-chromosome alignments
- Primary BAM (`Aligned.sortedByCoord.out.bam` or `Aligned.out.bam`) is suppressed by default unless `--keepBAM yes` is specified
- Works with both `BAM Unsorted` and `BAM SortedByCoordinate` output types

## Building References

The half-probe cache is built from a hybrid reference index that includes a pseudo-chromosome for each probe: `--runMode hashCacheGenerate` reads the probe sequences from it (see [Quick Start](#quick-start)). A default `--flex yes` run does not load the index. Aligning reads to the hybrid reference is the legacy route (`--flexLegacy yes`).

Scripts are provided in `scripts/` to build these references:

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

After building, use the index with the probe gene list to generate the half-probe cache (see [Quick Start](#quick-start)), or to run a legacy alignment route:

```bash
STAR \
  --genomeDir /path/to/flex_index \
  --flex yes \
  --flexLegacy yes \
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

See [autoindex docs](../core/legacy/docs/autoindex_cellranger.md) for URL selection (`--cellrangerRefRelease` / `--faUrl` / `--gtfUrl`), checksum flags, and parity test scripts.

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

The tool runs the **legacy** per-tag caller (expected-cell allocation, per-tag EmptyDrops, partition-occupancy filter) and by default strips output barcodes to 16 bases. It does not reproduce the tag-aware caller that `--flex yes` uses by default in STAR Suite 1.9.4; to reproduce STAR's calls, run STAR itself.

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
