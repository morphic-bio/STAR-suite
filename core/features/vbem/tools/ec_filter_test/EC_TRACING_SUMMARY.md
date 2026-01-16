# EC-Level Read Tracing Implementation Summary

## Completed

### 1. Analysis Script Created
**File**: `tools/ec_filter_test/analyze_ec_reads.py`

- Parses trace files from both tools
- Groups reads by EC label
- Identifies ECs with count differences
- Outputs detailed report with read names

**Usage**:
```bash
python3 analyze_ec_reads.py \
    --our-trace /tmp/our_ec_trace.txt \
    --salmon-ec-file /tmp/salmon_eq_classes.txt \
    --our-ec-file /tmp/our_eq_classes.txt \
    --output report.txt
```

### 2. Read Characteristic Analysis Script
**File**: `tools/ec_filter_test/analyze_differing_reads.py`

- Extracts BAM characteristics for reads in differing ECs
- Analyzes flags, MAPQ, NH, AS tags, etc.
- Helps identify patterns in filtered reads

**Usage**:
```bash
python3 analyze_differing_reads.py \
    --analysis-file report.txt \
    --bam-file input.bam \
    --ec-label "93" \
    --output read_analysis.txt
```

### 3. Trace Files Generated
- **Our CLI**: `/tmp/our_ec_trace.txt` (30,130 reads traced)
- **Salmon**: Trace file generation not working (needs investigation)

## Findings

### ECs with Largest Differences
1. EC (92,): 162 reads difference (13,350 vs 13,512)
2. EC (93,): 65 reads difference (5,400 vs 5,465)
3. EC (81,): 42 reads difference (2,494 vs 2,536)
4. EC (38,): 39 reads difference (2,127 vs 2,166)
5. EC (38,39): 36 reads difference (2,585 vs 2,621)

### Read Characteristics Analyzed
- **Sample reads from EC (93,)**: All proper pairs, MAPQ 255, NH:1, **no AS tags**
- **AS tag distribution**: 3,588 reads with AS, 26,542 without AS
- **Conclusion**: Missing AS tags is NOT the cause (too many reads without AS are included)

## Remaining Work

### 1. Fix Salmon Trace File Generation
**Issue**: `SALMON_TRACE_FILE` environment variable not generating trace file

**Possible causes**:
- Salmon binary not built with trace instrumentation
- Trace file path issue
- Trace code not executing

**Next steps**:
- Check if Salmon trace code is compiled in
- Verify trace file path/write permissions
- Rebuild Salmon if needed

### 2. Complete Read-by-Read Comparison
Once Salmon trace is working:
- Compare read sets for each differing EC
- Identify exact 404 reads that differ
- Categorize filtering reasons

### 3. Root Cause Analysis
After identifying the 404 reads:
- Check if they share common characteristics
- Verify if they're filtered by Salmon's checks we're missing
- Implement fixes or document as expected behavior

## Current Status

- ✅ Analysis infrastructure in place
- ✅ Our trace file working
- ⏳ Salmon trace file generation needs fixing
- ⏳ Read-by-read comparison pending

## Next Steps

1. **Investigate Salmon trace file generation**
   - Check if trace code exists in SalmonQuantifyAlignments.cpp
   - Verify compilation and execution
   - Fix or rebuild if needed

2. **Run full comparison** once both trace files available
3. **Identify root cause** of 404-read difference
4. **Implement fix** or document as expected

## Analysis Results

### Read Characteristics Summary
- **Single-transcript ECs** (92, 93, 81, 38, etc.):
  - All proper pairs (flags 83/163 or 99/147)
  - MAPQ 255 (perfect mapping quality)
  - NH:1 (single alignment location)
  - No AS tags (STAR doesn't always output AS)
  
- **Multi-transcript ECs** (38,39, 81,82):
  - Lower MAPQ (3.0 average)
  - NH:2 (multiple alignment locations)
  - More complex flag distributions
  - No AS tags

### Key Finding
**No obvious filtering pattern identified**. The differing reads appear normal and should be included by Salmon. This suggests:
1. The difference is due to subtle filtering logic we haven't identified
2. Salmon may be filtering based on transcript mass or other factors we don't check
3. The difference may be in how reads are grouped/aggregated rather than filtered

## Scripts Created

1. **analyze_ec_reads.py**: Main EC-level comparison script
2. **analyze_differing_reads.py**: Detailed BAM analysis for specific ECs
3. **identify_differing_reads.py**: Pattern identification across differing ECs

All scripts are ready to use once Salmon trace file generation is fixed.

## ROOT CAUSE IDENTIFIED (2024-12-19)

### The 404-Read Difference Explained

**Cause**: Salmon's default `--incompatPrior` behavior filters strand-incompatible alignments.

**Details**:
- Default `incompatPrior` is very small (~0)
- When `incompatPrior < 1e-100`, Salmon sets `ignoreIncompat = true`
- This causes Salmon to skip alignments where strand orientation doesn't match the library format
- Our CLI included all alignments regardless of strand compatibility

**Evidence**:
```
# Default Salmon (ignores incompatible):
Counted 29,726 total reads in the equivalence classes

# With --incompatPrior 1.0:
Counted 30,130 total reads in the equivalence classes
```

### Resolution

**Option 1**: Run Salmon with `--incompatPrior 1.0` to match our CLI behavior:
```bash
salmon quant ... --incompatPrior 1.0
```
Result: **100% parity** (30,130 reads, 102 ECs)

**Option 2**: Add strand compatibility filtering to our CLI to match Salmon's default:
- Check `isCompatible(libFormat, expectedFormat, ...)` for each alignment
- Skip alignments where `!isCompat` and `incompatPrior` is below threshold

### Current Parity Status (with --incompatPrior 1.0)

| Metric | Value |
|--------|-------|
| Reads processed | 30,130 (both) |
| ECs created | 102 (both) |
| Label parity | 100% |
| Weight parity | 100% |
| Count parity | 100% |
