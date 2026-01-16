# STAR-Flex Parameter Reference

This document lists **STAR-Flex-only parameters** that are not present in upstream STAR. For all other parameters, see upstream `README.md` and `parametersDefault`.

## Bulk RNA-seq Features

### Cutadapt-Style Trimming

**Note**: Trimming is a general-purpose feature for bulk RNA-seq (and other workflows), not specific to the Flex pipeline.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--trimCutadapt` | `-` | Enable cutadapt-style trimming (`Yes`/`-`). When enabled, all `clip*` parameters are ignored and ClipMate clipping is bypassed. |
| `--trimCutadaptQuality` | `20` | Quality threshold for 3' trimming (Phred scale, integer) |
| `--trimCutadaptMinLength` | `20` | Minimum read length after trimming. Pairs with either mate shorter than this are dropped (integer) |
| `--trimCutadaptAdapter` | `-` | Custom adapter sequences for R1 and R2 (space-separated). Default (`-`) uses TruSeq adapters: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` (R1) and `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` (R2) |
| `--trimCutadaptCompat` | `-` | Compatibility mode for adapter matching (`-`/`Off`/`Cutadapt3`). Default (`-`/`Off`) uses cutadapt 5.1 parity. `Cutadapt3` enables cutadapt 3.x compatibility mode for matching Trim Galore/cutadapt 3.2 behavior. |

**Note**: Trimming achieves perfect parity with Trim Galore/cutadapt v5.1 by default. The `Cutadapt3` compatibility mode reproduces cutadapt 3.x behavior for datasets processed with older Trim Galore versions. See [docs/trimming.md](trimming.md) for algorithm details and usage examples.

### TranscriptVB (VB/EM) + tximport-Style Gene Summaries

**Note**: This is a bulk RNA-seq quantification feature (Salmon-style EC-based VB/EM) integrated into STAR-Flex.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--quantMode TranscriptVB` | - | Enable TranscriptVB quantification (value of `--quantMode`) |
| `--quantVBem` | `0` | Use EM instead of VB (0/1) |
| `--quantVBLibType` | `A` | Library format (A=auto-detect; IU/ISF/ISR/U fixed) |
| `--quantVBAutoDetectWindow` | `1000` | Reads used for auto-detect (when `--quantVBLibType A`) |
| `--quantVBgenes` | `1` | Write gene summary `quant.genes.sf` (0/1) |
| `--quantVBgenesMode` | `Legacy` | Gene summarization mode: `Legacy` (simple sum) or `Tximport` (lengthScaledTPM, writes `quant.genes.tximport.sf`) |

See [docs/TranscriptVB_quantification.md](TranscriptVB_quantification.md) for details, warnings, and example commands.

### Samtools-Style BAM Sorting Backend

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outBAMsortMethod` | `star` | Coordinate-sorting backend for BAM output: `star` (legacy) or `samtools` (spill-to-disk) |

### Y-Chromosome BAM Split

**Note**: This feature was developed for **Morphic requirements for KOLF cell lines**. It is **not connected to the Flex pipeline** and is a general-purpose feature usable with any STAR workflow.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--emitNoYBAM` | `no` | Enable Y-chromosome BAM splitting (`yes`/`no`). When enabled, emits `<out>_Y.bam` and `<out>_noY.bam`. Primary BAM is suppressed unless `--keepBAM yes` is specified. |
| `--emitYReadNames` | `no` | Emit list of read names with any Y-chromosome alignment (one per line). Useful for FASTQ filtering without generating Y/noY BAMs. |
| `--emitYNoYFastq` | `no` | Emit Y/noY FASTQ files directly during alignment (`yes`/`no`). Works with both bulk PE and Flex/solo modes. |
| `--emitYNoYFastqCompression` | `gz` | Compression for FASTQ output (`gz`/`none`). Default is gzip-compressed. |
| `--YFastqOutputPrefix` | - | Optional: override output prefix for Y FASTQ files (default: derived from input name, inserting `_Y` before the last `_R1`/`_R2`; if missing, falls back to `Y_reads.mateN`) |
| `--noYFastqOutputPrefix` | - | Optional: override output prefix for noY FASTQ files (default: derived from input name, inserting `_noY` before the last `_R1`/`_R2`; if missing, falls back to `noY_reads.mateN`) |
| `--keepBAM` | `no` | Keep primary BAM output when `--emitNoYBAM yes` is enabled (`yes`/`no`) |
| `--noYOutput` | - | Optional: override default path for noY BAM output |
| `--YOutput` | - | Optional: override default path for Y BAM output |
| `--YReadNamesOutput` | - | Optional: override output path for Y read names list (default: `<outFileNamePrefix>Aligned.out_Y.names.txt`) |

See [docs/Y_CHROMOSOME_BAM_SPLIT.md](Y_CHROMOSOME_BAM_SPLIT.md) for details.

## Index-Time Features

### AutoIndex + CellRanger-Style Reference Formatting

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--autoIndex` | `No` | Build index in `--genomeDir` only if missing |
| `--forceIndex` | `No` | Rebuild index even if it exists (keeps cached inputs) |
| `--forceAllIndex` | `No` | Re-download/re-format/rebuild (wipes cache + formatted files) |
| `--cellrangerStyleIndex` | `No` | Format FASTA/GTF using CellRanger-style rules before indexing |
| `--cellrangerRefRelease` | `2024-A` | Default URL set used by AutoIndex when no inputs provided (`2024-A`/`2020-A`) |
| `--cellrangerStyleCacheDir` | `-` | Cache directory for downloads and checksum cache (default `${genomeDir}/cellranger_ref_cache`) |
| `--cellrangerStyleDownloadOnly` | `No` | Download-only mode (exit after download; skip formatting and indexing) |
| `--faUrl` | `-` | FASTA URL (must be provided together with `--gtfUrl`) |
| `--gtfUrl` | `-` | GTF URL (must be provided together with `--faUrl`) |
| `--autoCksumUpdate` | `No` | Auto-fill missing checksums from CHECKSUMS files for trusted URLs |
| `--replaceUnverifiableFiles` | `No` | Replace existing final files that cannot be verified (checksum mismatches still error) |
| `--allUntrustedUrl` | `No` | Allow untrusted URLs without checksums (disables integrity checking for those URLs) |

See [docs/autoindex_cellranger.md](autoindex_cellranger.md) for usage examples and output layout.

### Transcriptome FASTA Generation

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--genomeGenerateTranscriptome` | `No` | Generate `${genomeDir}/transcriptome.fa` during genomeGenerate |
| `--genomeGenerateTranscriptomeFasta` | `-` | Override transcriptome FASTA output path |
| `--genomeGenerateTranscriptomeOverwrite` | `No` | Overwrite existing transcriptome FASTA (`Yes`/`No`) |

## Flex-Specific Features

### Flex Pipeline

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--flex` | `no` | Enable flex pipeline (`yes`/`no`) |

**Note**: Flex pipeline is specific to 10x Genomics Flex (Fixed RNA Profiling) samples. See [README_flex.md](../README_flex.md) for complete flex parameter documentation.

## Upstream Parameters

For all other parameters (genome generation, alignment, solo, etc.), refer to upstream STAR documentation:
- `README.md` - User guide
- `parametersDefault` - Complete parameter list with defaults
- `doc/STARmanual.pdf` - Full manual
