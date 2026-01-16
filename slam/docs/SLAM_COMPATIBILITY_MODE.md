# SLAM Compatibility Mode

This document describes the optional GEDI-compatibility mode for SLAM quantification in STAR-Flex. This mode mirrors specific GEDI behaviors for parity testing and optional alignment with GEDI results.

## Overview

The SLAM compatibility mode provides an optional set of behaviors that approximate how GEDI (Genome-wide Dynamics of Intronic and Exonic transcription) processes SLAM-seq data. This enables:

1. **Parity testing**: Compare STAR-Flex SLAM results with GEDI outputs
2. **Migration**: Users transitioning from GEDI can maintain consistent behavior
3. **Research**: Study the effect of different algorithmic choices on NTR estimates

**Important**: This is a clean-room re-implementation based on observed behaviors and public descriptions. No GEDI source code was used.

## Quick Start

Enable GEDI-compatible mode with a single flag:

```bash
STAR --quantMode Slam \
     --slamCompatMode gedi \
     # ... other parameters
```

Or enable individual features granularly:

```bash
STAR --quantMode Slam \
     --slamCompatIntronic 1 \
     --slamCompatLenientOverlap 1 \
     --slamCompatOverlapWeight 1 \
     # ... other parameters
```

## Parameters

### Mode Selection

| Parameter | Values | Default | Description |
|-----------|--------|---------|-------------|
| `--slamCompatMode` | `none`, `gedi` | `none` | Compatibility mode preset |

When `--slamCompatMode gedi` is set:
- `slamCompatIntronic` → enabled
- `slamCompatLenientOverlap` → enabled
- `slamCompatOverlapWeight` → enabled
- `slamCompatIgnoreOverlap` → **NOT** enabled (see [Known Divergences](#known-divergences))

### Individual Feature Flags

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--slamCompatIntronic` | 0/1 | not set | GEDI-style intronic classification |
| `--slamCompatLenientOverlap` | 0/1 | not set | 50% overlap + SJ concordance acceptance |
| `--slamCompatOverlapWeight` | 0/1 | not set | Divide weight by gene count |
| `--slamCompatIgnoreOverlap` | 0/1 | not set | Skip all PE overlap positions (both mates) |
| `--slamCompatTrim5p` | int | 0 | 5' trim guard (bases) |
| `--slamCompatTrim3p` | int | 0 | 3' trim guard (bases) |

**Note**: Granular flags override mode defaults, including explicit `0` to disable a feature. For example, to use GEDI mode but disable overlap weighting:

```bash
--slamCompatMode gedi --slamCompatOverlapWeight 0
```

This works because flags use a sentinel value internally (`-1` = not set). When you explicitly provide `0`, it overrides the mode default.

## Feature Details

### 1. Intronic Classification (`slamCompatIntronic`)

When enabled, reads that would otherwise have no gene assignment can be classified as intronic if they meet specific criteria:

**Requirements for intronic classification:**
1. No exonic gene assignment from standard `classifyAlign`
2. Read is single-part (unspliced, `nExons == 1`)
3. Read intersects introns of **>1 transcripts** for the candidate gene

This approximates GEDI's "location count" check for intronic reads.

**Candidate gene source**: Uses `GeneFull_ExonOverIntron` overlap set (genes where read overlaps gene body but wasn't assigned exonically).

### 2. Lenient Overlap Acceptance (`slamCompatLenientOverlap`)

Relaxes the strict exonic concordance requirement for gene assignment:

**Standard behavior**: Read must be fully exon-concordant with transcript (all blocks inside exons, splice junctions match).

**Lenient behavior**: Accepts reads meeting:
1. ≥50% of mapped length overlaps transcript exons
2. Splice junction consistency (via `alignToTranscriptMinOverlap`)
3. Status must be `Concordant` or `ExonIntron` (NOT purely `Intron`)

This catches reads that slightly span exon boundaries but are otherwise consistent.

### 3. Overlap-Gene Weighting (`slamCompatOverlapWeight`)

Adjusts alignment weights based on ambiguous gene assignments:

**Standard behavior**: Weight = 1/nTr (divided by number of alignments)

**Compat behavior**: Weight = (1/nTr) / geneCount

Where `geneCount` is the **read-level** union of all genes assigned across all alignments (`annFeat.fSet.size()`). This matches GEDI's `gset` semantics.

**Semantics**: Uses read-level union (not per-alignment) to match GEDI behavior.

### 4. Position Filtering (`slamCompatIgnoreOverlap`, `slamCompatTrim5p/3p`)

Controls which positions are counted for T→C conversions:

**Overlap filtering** (`slamCompatIgnoreOverlap`):
- When enabled, **all** positions where PE mates overlap are skipped (both mate 1 and mate 2 positions)
- This differs from STAR's default behavior which only skips mismatches on mate 2 in overlap regions
- **This is a deliberate divergence from GEDI** (see below)

**Trim guards** (`slamCompatTrim5p`, `slamCompatTrim3p`):
- Skip positions within N bases of mate 5' or 3' ends
- Applied per-mate (not concatenated read position)
- Handles edge effects and sequencing errors at read boundaries

**Coordinate conversion**: STAR uses concatenated read positions internally. The compat layer converts to per-mate local coordinates:
```
Mate 1: mateLocalPos = readPos
Mate 2: mateLocalPos = readPos - readLength[0]
```

## Known Divergences from GEDI

### 1. `slamCompatIgnoreOverlap` (NOT default in gedi mode)

GEDI treats PE overlap positions specially but still counts them. STAR's `slamCompatIgnoreOverlap` skips **all** overlap positions entirely (both mates). This is available as an explicit override but is **not enabled** by `--slamCompatMode gedi`.

To enable complete overlap skipping:
```bash
--slamCompatMode gedi --slamCompatIgnoreOverlap 1
```

Note: STAR's default behavior (without compat mode) skips only mate 2 mismatches in overlap regions. The `slamCompatIgnoreOverlap` option is more aggressive, skipping all positions in overlap regions.

### 2. Intronic Classification Threshold

GEDI requires the read to overlap introns of >1 transcripts. STAR uses the same threshold. However, single-transcript genes can never be classified as intronic under this rule.

### 3. Splice Junction Concordance Implementation

STAR uses `alignToTranscriptMinOverlap` with `minOverlap=0` for lenient acceptance. This may differ slightly from GEDI's internal SJ concordance check.

## Diagnostics and Monitoring

### Log Output

When compat mode is enabled, a summary line is written to `Log.out`:

```
SLAM compat(gedi): alignsIntronic=1234 alignsLenient=567 alignsWeightAdj=890 posSkipOvlp=456 posSkipTrim=123
```

### Counter Semantics

| Counter | Level | Description |
|---------|-------|-------------|
| `compatAlignsReclassifiedIntronic` | Per-alignment | Alignments reclassified as intronic |
| `compatAlignsLenientAccepted` | Per-alignment | Alignments accepted via lenient overlap (at most once per alignment, even if multiple transcripts matched) |
| `compatAlignsOverlapWeightApplied` | Per-alignment | Alignments where weight was adjusted |
| `compatPositionsSkippedOverlap` | Per-position | Positions skipped due to PE overlap |
| `compatPositionsSkippedTrim` | Per-position | Positions skipped due to trim guards |

**Note**: `compatAlignsLenientAccepted` counts alignments, not transcript matches. If one alignment matches 3 transcripts via lenient acceptance, it contributes 1 to the counter.

### Diagnostics File

When SLAM diagnostics are enabled, compat counters appear in the diagnostics file under `CompatModeCounters:` section.

## Architecture

### Module Structure

```
source/
├── SlamCompat.h         # Config struct and helper class
├── SlamCompat.cpp       # Stateless helpers with cached transcript data
├── Parameters.h         # CLI parameter definitions (quant.slam.compat*)
└── Integration points:
    ├── ReadAlignChunk.cpp           # Owns SlamCompat instance
    ├── ReadAlign_outputAlignments.cpp  # Weighting, intronic classification
    ├── ReadAlign_slamQuant.cpp      # Position filtering
    └── Transcriptome_classifyAlign.cpp  # Lenient overlap acceptance
```

### Ownership Model

- `ReadAlignChunk` owns `SlamCompat` instance (created per thread)
- `ReadAlign` holds raw pointer for runtime access
- `Transcriptome::classifyAlign` receives `SlamCompat*` as optional parameter
- All counters live in `SlamDiagnostics` (single source of truth)

### Stateless Design

`SlamCompat` helpers are stateless functions. The class only caches:
- `geneToTranscripts_`: Gene ID → list of transcript IDs
- `transcriptIntrons_`: Transcript ID → list of intron intervals (genomic coords)

This allows safe parallel execution across threads.

## Testing

### Unit Tests

Run the SlamCompat unit tests:

```bash
./tests/run_slam_compat_test.sh
```

Tests cover:
- `compatOverlapWeight`: Gene count weighting
- `compatShouldCountPos`: Trim guard logic with mate boundaries
- `compatIsIntronic`: Single-part gate and multi-transcript requirement
- `computeExonOverlap`: Exon overlap calculation

### Integration Testing

Lenient overlap acceptance flows through `Transcriptome::classifyAlign` since `alignToTranscriptMinOverlap` is a static function. Currently, unit tests cover the helper functions (`computeExonOverlap`, etc.) but full integration testing through `classifyAlign` requires constructing a complete `Transcriptome` object with transcript data. This is an integration path that exercises the full lenient acceptance logic but is not currently covered by standalone unit tests.

## Example Workflow

### Comparing STAR-Flex vs GEDI

```bash
# Run STAR with GEDI compat mode
STAR --runMode alignReads \
     --quantMode Slam \
     --slamCompatMode gedi \
     --genomeDir /path/to/genome \
     --readFilesIn reads_R1.fastq reads_R2.fastq \
     --outFileNamePrefix slam_compat_

# Compare NTR correlations between outputs
python compare_ntr.py gedi_output.tsv slam_compat_SlamQuant.out
```

### Debugging Discrepancies

Enable debug mode to trace individual reads:

```bash
STAR --quantMode Slam \
     --slamCompatMode gedi \
     --slamDebugGeneList genes_of_interest.txt \
     --slamDebugReadList reads_of_interest.txt \
     --slamDebugOutPrefix debug_
```

This produces per-read logs with gene assignments, weights, and drop reasons. Note that debug output currently logs overall read processing status but does not explicitly annotate which compat-specific features (lenient overlap, intronic classification, etc.) affected each read decision. Compat feature effects can be inferred from the aggregate counters in the diagnostics file.

## Validation Results

### Fixture Parity Tests

Tests were run on the SLAM fixture dataset (100K reads) comparing STAR-Flex with GEDI-compat mode against GRAND-SLAM reference.

#### Run A: GEDI Compat Mode Defaults

**Configuration:**
```bash
--slamQuantMode 1 --slamSnpDetect 1 --slamCompatMode gedi
```

**Correlation Metrics:**

| Threshold | Filter | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|--------|---------|-------------|--------------|--------------|---------------|
| >=20 | readcount | 384 | 0.982929 | 0.963909 | 0.987470 | 0.986444 |
| >=50 | readcount | 76 | 0.996176 | 0.979539 | 0.997971 | 0.993322 |
| >=100 | readcount | 23 | 0.993175 | 0.988131 | 0.995595 | 0.987154 |

#### Run B: GEDI Compat Mode with Overlap Weight Override

**Configuration:**
```bash
--slamQuantMode 1 --slamSnpDetect 1 --slamCompatMode gedi --slamCompatOverlapWeight 0
```

**Correlation Metrics:**

| Threshold | Filter | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |
|-----------|--------|---------|-------------|--------------|--------------|---------------|
| >=20 | readcount | 384 | 0.989526 | 0.967828 | 0.991410 | 0.990479 |
| >=50 | readcount | 76 | 0.996190 | 0.979457 | 0.997982 | 0.993322 |
| >=100 | readcount | 23 | 0.993175 | 0.988131 | 0.995595 | 0.987154 |

#### Observations

1. **Override semantics work correctly**: Run B shows different correlations than Run A, confirming that `--slamCompatOverlapWeight 0` successfully overrides the gedi mode default.

2. **Correlation changes**: Disabling overlap weighting improves correlations slightly:
   - NTR Pearson: 0.982929 → 0.989526 (+0.0066)
   - k/nT Pearson: 0.987470 → 0.991410 (+0.0039)

3. **High correlations maintained**: Both configurations show excellent correlations (>0.98) with GRAND-SLAM reference, indicating the compat mode is functioning as intended.

4. **Threshold stability**: Correlations remain high across thresholds (>=20, >=50, >=100), with best performance at >=50 threshold.

## Auto-Trim and Multi-File Processing

### Auto-Trim Mode

STAR-Flex supports automatic detection of optimal 5' and 3' trim values based on variance analysis of the first N reads:

```bash
STAR --quantMode Slam \
     --autoTrim variance \
     --autoTrimDetectionReads 100000 \
     --autoTrimMinReads 1000 \
     # ... other parameters
```

The variance-based auto-trim analyzes per-position quality scores and T→C conversion rates to identify problematic regions at read ends, then computes optimal trim values using knee detection.

### Multi-File Trim Scope

When processing multiple input files in a single run, the `--trimScope` parameter controls how auto-trim is applied:

| Value | Behavior |
|-------|----------|
| `first` (default) | Compute trim from first file, apply to all files |
| `per-file` | Compute and apply separate trim values for each file |

**Example with per-file trim:**

```bash
STAR --quantMode Slam \
     --autoTrim variance \
     --trimScope per-file \
     --readFilesIn file1.fq,file2.fq,file3.fq \
     # ... other parameters
```

### Important: FILE Marker Requirement

The `trimScope=per-file` mode requires FILE markers in the read stream to detect file boundaries. These markers are automatically generated when:

1. Using comma-separated `--readFilesIn` (e.g., `file1.fq,file2.fq,file3.fq`)
2. Using `--readFilesCommand` with multiple input files

**Warning**: If FILE markers are not present (e.g., when using a single concatenated input), `trimScope=per-file` cannot reliably detect file boundaries and will emit a warning.

### Read ID Handling

In multi-file mode with `trimScope=per-file`, read IDs are kept globally unique across all files. This ensures:
- Consistent behavior with downstream tools expecting unique read IDs
- Correct BAM record ordering for coordinate-sorted output
- Proper handling by samtools and other BAM processing tools

### Timing and Statistics

Statistics in `Log.final.out` reflect the overall run:
- Mapping speed is calculated from the overall start time (not per-file)
- Read counts accumulate across all files
- Detection pass reads are excluded from final statistics

## Changelog

- **v1.1.0**: Auto-trim and multi-file support
  - Variance-based auto-trim with knee detection
  - Per-file trim scope for multi-file runs
  - Global read ID uniqueness across files
  - Detection pass output suppression
  - FILE marker requirement documentation

- **v1.0.0**: Initial implementation
  - Intronic classification with single-part gate
  - Lenient overlap acceptance (50% threshold + SJ concordance)
  - Read-level overlap-gene weighting
  - Position filtering with trim guards
  - Full diagnostics integration
  - Override semantics with sentinel values (-1 = not set, 0/1 = explicit override)
  - Per-alignment lenient acceptance counter (not per-transcript)