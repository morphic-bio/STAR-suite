# Auto-Trim and Multi-File Implementation Summary

## Overview

This document summarizes the implementation of auto-trim detection with multi-file support for STAR-SLAM. The feature enables automatic detection of optimal 5' and 3' trim values based on per-position T→C conversion rate variance analysis using segmented regression.

## Key Features Implemented

### 1. Variance-Based Auto-Trim Detection (Segmented Regression)

**Algorithm:** Segmented regression with **two breakpoints** on the smoothed T→C rate standard deviation curve.

**Status:** ✅ Implemented (replaced kneedle algorithm)

**Signal Used:** Standard deviation of T→C conversion rate at each read position
- High stdev = inconsistent conversion = artifact (trim this region)
- Low stdev = consistent signal = keep this region

**Methodology:**
1. Build per-position T→C stdev curve from variance analysis
2. Smooth with median window (default 5 bp)
3. Fit a 3-segment piecewise linear model with 2 breakpoints (b1, b2) by minimizing total SSE
4. Use prefix sums for efficient O(n²) segment fitting
5. Enforce minimum segment length (default 3) and max-trim (default 15) constraints
6. Set `trim5p = b1`, `trim3p = read_len - 1 - b2`

**Data Collection:**
- Per-position statistics: T count, T→C count, sum of rates, sum of squares
- Computes `stddevTcRate()` for each position
- Linear interpolation for missing positions before smoothing

**Files:**
- `source/SlamVarianceAnalysis.h` - Data structures, `SegmentFit` struct, API
- `source/SlamVarianceAnalysis.cpp` - Segmented regression implementation

### 2. Multi-File Trim Scope

**Two modes via `--trimScope` parameter:**

| Mode | Behavior |
|------|----------|
| `first` (default) | Compute trim from first file, apply to all files |
| `per-file` | Compute and apply separate trim values for each file |

**Per-file processing flow:**
1. For each file: run single-threaded detection pass
2. Compute trim values from segmented regression
3. Rewind input and run multi-threaded mapping with trims applied
4. Repeat for next file

### 3. Detection Pass Output Suppression

Detection pass does NOT produce outputs:
- BAM output guards in `unsortedOneAlign()`, `coordOneAlign()`
- SAM output guard in `outputTranscriptSAM()`
- Flag: `P.quant.slam.autoTrimDetectionPass`

### 4. Statistics Handling

**Timing:** Overall `timeStartMap` preserved for accurate speed reporting in `Log.final.out`

**Read IDs:** Monotonic across files via `cumulativeReadOffset` tracking

**Stats accumulation:** Detection reads excluded from final totals; mapping stats accumulate across files

### 5. FILE Marker Requirement

`trimScope=per-file` requires FILE markers in read stream:
- Auto-generated with comma-separated `--readFilesIn`
- Safety warnings if markers not found after 100M lines

### 6. QC Output

**QC JSON:** Per-position variance stats and segmented regression info
- Smoothed stdev curve
- Breakpoints b1, b2
- Segment fit parameters (slope, intercept, SSE)

**QC HTML:** Interactive Plotly visualization
- T→C stdev plot with segment fits
- Trim position markers
- Quality statistics plot

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--autoTrim` | `-` | Auto-trim mode: `variance` or `-` (disabled) |
| `--trimScope` | `first` | `first` or `per-file` |
| `--trimSource` | `-` | Path to file for computing shared trims (overrides first file) |
| `--autoTrimDetectionReads` | `100000` | Max reads for detection pass |
| `--autoTrimMinReads` | `1000` | Min reads required for trim computation |
| `--autoTrimMaxReads` | `100000` | Max reads for variance analysis |
| `--autoTrimSmoothWindow` | `5` | Median smoothing window for stdev curve |
| `--autoTrimSegMinLen` | `3` | Minimum segment length for regression |
| `--autoTrimMaxTrim` | `15` | Maximum trim at either end |

### Trim Source Override (`--trimSource`)

The `--trimSource` parameter allows specifying a different file for computing shared auto-trim values:

- When set, auto-trim is computed from the specified file instead of the first input file
- The computed trims are then applied to all input files in the run
- Useful when the first file is not representative (e.g., low depth or outlier sample)
- Ignored when `--trimScope per-file` is set (per-file mode always computes trims independently)
- The file must exist and be readable; STAR exits with an error otherwise

**Example:**
```bash
# Use file2.fq.gz for trim computation, apply trims to all files
STAR --readFilesIn file1.fq.gz,file2.fq.gz,file3.fq.gz \
     --autoTrim variance \
     --trimSource file2.fq.gz \
     ...
```

## Test Results (SE50 SLAM-seq data)

**Single file test with `--autoTrim variance`:**
```
SLAM auto-trim (segmented regression) computed:
    trim5p=10 trim3p=4 mode=auto_segmented
    breakpoints: b1=10 b2=45 total_sse=0.0085467
    reads_analyzed=13376 detection_reads_processed=100000 min_reads=5000
    smooth_window=5 min_seg_len=3 max_trim=15 scope=first
```

**Interpretation:**
- Breakpoint b1=10: First 10 positions have elevated T→C stdev (5' artifact)
- Breakpoint b2=45: Last 5 positions (50-45=5, clamped to trim3p=4) have elevated stdev (3' artifact)
- Middle segment (positions 10-45) has consistent, low stdev = true signal
- Total SSE = 0.0085 indicates good fit

**QC outputs generated:**
- `slam_qc.json` - Full variance stats and regression info
- `slam_qc.html` - Interactive visualization

## Files Modified

### Core Implementation
- `source/SlamVarianceAnalysis.h` - Variance stats struct with `tcRateSumSq`, `varianceTcRate()`, `stddevTcRate()`, `SegmentFit`, segmented regression API
- `source/SlamVarianceAnalysis.cpp` - Segmented regression implementation with prefix sums, median smoothing, linear interpolation
- `source/SlamQuant.h/.cpp` - Extended `enableVarianceAnalysis()` with smoothing/segment parameters
- `source/STAR.cpp` - Per-file processing loop, stats accumulation, timing preservation, QC output
- `source/ReadAlignChunk.cpp` - Pass new parameters to variance analyzer
- `source/ReadAlignChunk_processChunks.cpp` - FILE marker skip logic, PE support, boundary detection

### Output Guards
- `source/BAMoutput.cpp` - Detection pass guards
- `source/ReadAlign_outputTranscriptSAM.cpp` - SAM output guard

### QC Output
- `source/SlamQcOutput.h/.cpp` - JSON and HTML output with segmented regression info

### Parameters
- `source/Parameters.h` - New parameters: `autoTrimSmoothWindow`, `autoTrimSegMinLen`, `autoTrimMaxTrim`
- `source/Parameters.cpp` - Parameter registration
- `source/parametersDefault` - Parameter descriptions

### Documentation
- `docs/SLAM_COMPATIBILITY_MODE.md` - Auto-trim and multi-file section
- `plans/autotrim_multifile_implementation_summary.md` - This document

## Commits

1. `ed2a1e9` - Fix auto-trim detection pass output and stats handling
2. `7fb8e02` - Fix per-file auto-trim timing, read IDs, and FILE marker handling
3. `250f492` - Add documentation for auto-trim multi-file processing
4. `c219839` - Fix auto-trim to use T→C rate variance/stdev instead of mean rate
5. (pending) - Implement segmented regression, add QC output, new parameters

## Known Limitations

1. **Coordinate-sorted BAM with per-file mode:** May have issues due to multiple mapping phases. Use unsorted BAM or `trimScope=first` for coordinate-sorted output.

2. **FILE marker requirement:** `trimScope=per-file` only works with comma-separated `--readFilesIn` or `--readFilesCommand` that generates FILE markers.

3. **Single-threaded detection:** Detection pass runs single-threaded for simplicity. Main mapping is multi-threaded.

## Algorithm Details: Segmented Regression

The segmented regression algorithm finds the optimal breakpoints by:

1. **Smoothing:** Apply median filter with window size W (default 5)
2. **Interpolation:** Fill missing positions with linear interpolation
3. **Prefix sums:** Precompute cumulative sums for efficient segment fitting:
   - `prefixX[i]` = Σx (x = position index)
   - `prefixXX[i]` = Σx²
   - `prefixY[i]` = Σy (y = smoothed stdev)
   - `prefixXY[i]` = Σxy
   - `prefixYY[i]` = Σy²

4. **Grid search:** For each candidate (b1, b2) pair:
   - Segment 1: [0, b1-1] - 5' artifact region
   - Segment 2: [b1, b2] - middle signal region
   - Segment 3: [b2+1, n-1] - 3' artifact region
   - Fit linear model to each segment using prefix sums
   - Compute SSE for each segment
   - Track best (b1, b2) with minimum total SSE

5. **Constraints:**
   - Each segment must have at least `minSegLen` positions
   - b1 ≤ maxTrim (5' trim cannot exceed max)
   - n-1-b2 ≤ maxTrim (3' trim cannot exceed max)

6. **Output:** trim5p = b1, trim3p = n-1-b2

This matches the methodology in `analyze_phred_by_position.py` for consistency.
