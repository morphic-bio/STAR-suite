# Cutadapt-Style Trimming in STAR-Flex

**Note**: This is a **STAR-Flex-only feature** for **bulk RNA-seq** (and other workflows), not present in upstream STAR. This feature is **not specific to the Flex pipeline** and can be used with any STAR workflow.

STAR-Flex includes an optional cutadapt-parity trimming feature that matches Trim Galore defaults for quality and adapter trimming. It is designed for general-purpose use in bulk RNA-seq pipelines.

## Quick Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--trimCutadapt` | `-` | Enable trimming (`Yes`/`-`) |
| `--trimCutadaptQuality` | `20` | Quality threshold (Phred scale) |
| `--trimCutadaptMinLength` | `20` | Minimum read length after trimming |
| `--trimCutadaptAdapter` | `-` | Custom R1/R2 adapters (space-separated) or `-` for TruSeq defaults |
| `--trimCutadaptCompat` | `-` | Compatibility mode (`-`/`Off`/`Cutadapt3`). Default uses cutadapt 5.1 parity. `Cutadapt3` enables cutadapt 3.x compatibility mode. |

**Important**: When `--trimCutadapt Yes` is enabled, all `clip*` parameters are ignored and ClipMate clipping is bypassed.

For complete parameter reference, see [docs/flex_parameters.md](flex_parameters.md).

## Overview

When enabled, this feature:
- Performs 3' quality trimming using cutadapt's modified Mott algorithm
- Removes Illumina TruSeq adapters using Hamming distance matching
- Drops read pairs where either mate falls below the minimum length threshold
- Provides detailed trimming statistics in `Log.final.out`

## Parameters

### `--trimCutadapt`

Enable cutadapt-style trimming.

- **Default**: `-` (disabled)
- **Options**: 
  - `-`: Trimming disabled (default)
  - `Yes`: Enable with Trim Galore defaults (Q20, length 20, TruSeq adapters)

**Note**: When `trimCutadapt` is enabled, existing `clip*` parameters are **ignored**. ClipMate clipping is bypassed; use `trimCutadapt*` parameters instead.

### `--trimCutadaptQuality`

Quality threshold for 3' trimming (Phred scale).

- **Default**: `20`
- **Type**: Integer

### `--trimCutadaptMinLength`

Minimum read length after trimming. Pairs with either mate shorter than this are dropped.

- **Default**: `20`
- **Type**: Integer

### `--trimCutadaptAdapter`

Custom adapter sequences for R1 and R2 (space-separated).

- **Default**: `-` (uses default TruSeq adapters)
- **Format**: Two adapter sequences separated by space
- **Example**: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`

### `--trimCutadaptCompat`

Compatibility mode for adapter matching algorithm.

- **Default**: `-` or `Off` (cutadapt 5.1 parity)
- **Options**:
  - `-` or `Off`: Use default cutadapt 5.1 matching algorithm (default)
  - `Cutadapt3`: Use cutadapt 3.x compatibility mode for matching Trim Galore/cutadapt 3.2 behavior

**When to use**: Enable `Cutadapt3` mode when you need to match trimming results from datasets processed with Trim Galore using cutadapt 3.2 (or earlier 3.x versions). This mode adjusts adapter matching tie-breaking to prefer earlier matches (shorter trims) and relaxes acceptance criteria for short overlaps, reproducing the behavior differences between cutadapt 3.x and 5.x.

**Note**: Compatibility mode only affects adapter matching behavior, not quality trimming or minimum length logic.

## Default Adapter Sequences

When using default adapters (when `trimCutadaptAdapter` is `-`):

- **R1**: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`
- **R2**: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`

These match Trim Galore's default TruSeq adapter detection.

## Trimming Algorithm

### Quality Trimming

Uses cutadapt's modified Mott algorithm:
- Scans from 3' to 5'
- Accumulates score = (quality - cutoff) for each base
- Finds position where cumulative sum is maximized
- Resets negative sums (Mott modification)

### Adapter Matching

Default mode (`--trimCutadaptCompat -/Off`) uses cutadapt's semiglobal alignment algorithm (exact port from cutadapt v5.1). Legacy compatibility mode (`--trimCutadaptCompat Cutadapt3`) uses a faithful port of cutadapt v3.2 adapter matching (see [docs/cutadapt_3.2_parity_report.md](cutadapt_3.2_parity_report.md)).

- **Algorithm**: Semiglobal dynamic programming with edit distance (Levenshtein)
- **Source**: Ported verbatim from `cutadapt v5.1` (`_align.pyx`, `Aligner.locate()` method); Cutadapt3 mode ports cutadapt v3.2 matching/tie-breaking
- **Minimum overlap**: 1 base (Trim Galore default `--stringency 1`)
- **Maximum error rate**: 0.1 (10% errors allowed)
- **Error threshold**: `floor(aligned_adapter_length × 0.1)` where `aligned_adapter_length` is the adapter portion actually matched
- **Scoring**: `MATCH=+1, MISMATCH=-1, INSERTION=-2, DELETION=-2`
- **Flags**: `QUERY_START | QUERY_STOP | REFERENCE_END` (allows adapter to start anywhere in read, end early, and read suffix can be skipped)

**Note on adapter overlap**: The default minimum overlap of 1 base matches Trim Galore's `--stringency 1` setting. This ensures that even single-base adapter matches are trimmed, providing maximum parity with Trim Galore output. Increasing the overlap threshold (e.g., to 3 bases) would be less aggressive and may retain additional bases in edge cases where only 1-2 bases match the adapter sequence.

**Note**: Adapter matching is implemented in `source/libtrim/adapter_trim.cpp` and supports both default cutadapt v5.1 parity and optional cutadapt v3.2 compatibility. If trimming logic is changed, re-run parity tests (`make test_trim_parity`) and (for Cutadapt3 mode) the cutadapt 3.2 comparison described in [docs/cutadapt_3.2_parity_report.md](cutadapt_3.2_parity_report.md).

## Statistics Output

When trimming is enabled, `Log.final.out` includes a new section:

```
                        TRIMMING (cutadapt-style) |
              Reads processed for trimming |   1000000
                             Reads trimmed |   850000
     Reads dropped (below minimum length) |   5000
                    Bases quality-trimmed |   12500000
                    Bases adapter-trimmed |   8750000
```

**Note**: Statistics count individual reads (not pairs). For paired-end data, each pair contributes 2 to the read counts. For example, 500,000 pairs = 1,000,000 reads processed.

## Usage Example

```bash
STAR \
    --runThreadN 8 \
    --genomeDir /path/to/genome \
    --readFilesIn sample_R1.fastq sample_R2.fastq \
    --trimCutadapt Yes \
    --trimCutadaptQuality 20 \
    --trimCutadaptMinLength 20 \
    --outFileNamePrefix output/
```

## Single-End Support

Single-end trimming is supported. When `trimCutadapt` is enabled for single-end data:
- Uses R1 adapter sequence (or first adapter if custom adapters specified)
- Applies same quality and adapter trimming logic
- Drops reads below minimum length threshold

## Parity Validation

The trimming implementation achieves **perfect parity** with Trim Galore/cutadapt v5.1 (default mode):
- **FASTQ-level parity**: 0 diff lines in integration tests (all 9 tests pass)
- **Alignment parity**: Guaranteed by FASTQ parity (identical FASTQs → identical alignments)

**Test Suite**: Run `make test_trim_parity` from `source/` directory to validate parity.

For legacy datasets trimmed with Trim Galore + cutadapt 3.2, enable `--trimCutadaptCompat Cutadapt3` and validate against the cutadapt 3.2 parity procedure in [docs/cutadapt_3.2_parity_report.md](cutadapt_3.2_parity_report.md).

**End-to-End Test**: An optional end-to-end alignment test script is available at `tools/trimvalidate/run_end_to_end_parity.sh` for additional validation (requires reference genome).

## Comparison with ClipMate

| Feature | ClipMate | trimCutadapt |
|---------|----------|--------------|
| Quality trimming | No | Yes (Mott algorithm) |
| Adapter matching | Hamming/Opal | Semiglobal DP (cutadapt v5.1; optional v3.2 compat) |
| 5' adapter | Supported | Not supported (future) |
| Single-end | Supported | Supported |
| Statistics | Limited | Detailed (counts reads, not pairs) |
| Parity with Trim Galore | No | Yes (perfect parity) |

**Important**: When `trimCutadapt` is enabled, all `clip*` parameters are ignored. STAR will log a warning at startup.

## Future Enhancements

Planned (but not yet implemented):
- 5' quality trimming (optional)
- 5' adapter trimming
- Multiple adapter support
- Per-mate minimum length thresholds
