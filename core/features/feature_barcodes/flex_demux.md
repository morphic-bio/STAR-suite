# Flex Demux - Single Cell Feature Assignment and Demultiplexing Suite

## Overview

Flex Demux is a comprehensive, high-performance toolkit for single-cell targeted sequencing analysis, providing fast barcode assignment, feature matching, and sample demultiplexing capabilities. The suite consists of three main tools optimized for different stages of the single-cell analysis pipeline.

## Core Tools

### 1. assignBarcodes - Feature Barcode Assignment
The primary tool for assigning feature barcodes from FASTQ files to known sequence barcodes with advanced error correction and UMI deduplication.

**Key Features:**
- Exhaustive search with fuzzy matching for both ATAC-seq and RNA-seq
- Advanced error correction for sequence and feature barcodes
- Sophisticated UMI deduplication with connected component analysis
- Multi-level parallelization for high throughput processing
- Interactive QC plots and heatmaps

### 2. demux_fastq - Sample-Level Demultiplexing  
Lightweight helper for grouping raw FASTQ reads into per-sample folders based on 8-base probe barcodes embedded in reads.

**Key Features:**
- Fast probe barcode extraction and sample assignment
- Producer-consumer parallelization model
- Direct 64-bit comparison optimization for small probe sets
- Automatic file organization and compression

### 3. demux_bam - BAM to Probe Matrix Conversion
Processes STAR Solo-aligned BAM files to produce probe × cell-barcode matrices with UMI counting.

**Key Features:**
- Extracts STAR Solo tags (CB, UB, GX/GE)
- Configurable filtering (MAPQ, duplicates, alignment type)
- Memory-efficient counting without heavy deduplication tables
- Matrix Market format output
- Optional CB/UB binary stream ingestion via `--CBUB_file`

## Installation

### Prerequisites
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get -y install build-essential zlib1g-dev libglib2.0-dev libcairo2-dev

# Optional: Disable heatmap generation to skip cairo dependency
make NO_HEATMAP=1
```

### Build
```bash
cd process_features
make
sudo cp assignBarcodes demux_fastq demux_bam /usr/local/bin
```

### Docker
```bash
docker pull biodepot/process_features:latest
# or build from Dockerfile
docker build -t flex_demux .
```

## Quick Start Examples

### Feature Assignment (assignBarcodes)
```bash
./assignBarcodes \
    -d ./output_dir/ \
    -w /path/to/10x_whitelist.txt \
    -f /path/to/features.csv \
    -t 8 -S 4 -c 2 \
    -m 2 -u 12 -o 26 \
    --limit_search 5 \
    /path/to/sample_fastqs/
```

### Sample Demultiplexing (demux_fastq)
```bash
./demux_fastq \
    --probe_barcodes tables/probe-barcodes.txt \
    --sample_map tables/probe-to-sample-map.txt \
    --outdir demux_out \
    --probe_read R2 --probe_offset 68 \
    --threads 4 \
    /path/to/fastq_dir/
```

### BAM Processing (demux_bam)
```bash
./demux_bam \
    --bam input.bam \
    --outdir probe_matrix_out \
    --sample_probes tables/probe-barcodes.txt \
    --probe_offset 68 \
    --search_nearby \
    -S 4 -t 1
```

#### Using the STAR cb/ub stream
When STARsolo produces `Aligned.out.cb_ub.bin`, enable the faster CBUB path:

```bash
./demux_bam \
    --bam Aligned.out.bam \
    --CBUB_file Aligned.out.cb_ub.bin \
    --whitelist whitelists/737K-fixed-rna-profiling.txt \
    --sample_probes tables/probe-barcodes.txt \
    --probe_offset 68 \
    -S 4 -t 1
```

`demux_bam` synchronizes each BAM alignment with one binary record, so truncated or oversized streams are flagged immediately. The whitelist is mandatory in this mode because STAR encodes barcodes as 1-based indices.

##### CBUB binary layout

| Component | Description |
|-----------|-------------|
| Header | Four little-endian `uint64_t`: `status_bits`, `cb_bits`, `umi_bits`, `record_count`. Only `status_bits = 1` is accepted; `umi_bits` must be even. |
| Record | Packed `status_bits + cb_bits + umi_bits` bits (LSB-first, byte-aligned). `status=0` or `cb_index=0` are treated as missing tags. UMIs decode 2-bit bases (`00 A`, `01 C`, `10 G`, `11 T`). |

The BAM still provides gene tags (`GX`/`GE`), so downstream filtering and Matrix Market emission remain unchanged.

## Performance Characteristics

- **assignBarcodes**: ~1M reads/second on standard laptop for targeted sequencing
- **demux_fastq**: ~200 MB/s decompression, ~80k reads/s per consumer thread
- **demux_bam**: Memory usage O(#active_CB × #probes × 4 bytes)

## Output Formats

All tools produce standardized outputs:
- **Matrix Market format** (.mtx) for count matrices
- **TSV files** for barcodes and features
- **Interactive HTML plots** for quality control
- **PNG heatmaps** for visual QC assessment
- **Comprehensive statistics** files

## Quality Control Features

- Interactive UMI count histograms with Plotly visualization
- Feature counts heatmaps showing count distributions
- Feature richness heatmaps for multiplet detection
- Detailed run statistics and matched sequence reports
- Configurable filtering thresholds for QC metrics

## Architecture

The suite uses a modular C codebase with:
- **Memory pools** for efficient allocation
- **Producer-consumer** threading models
- **OpenMP parallelization** for compute-intensive operations
- **GLib data structures** for hash tables and arrays
- **HTSlib integration** for BAM file processing

This toolkit provides a complete solution for single-cell targeted sequencing analysis from raw FASTQ files through final count matrices, with emphasis on performance, accuracy, and comprehensive quality control.
