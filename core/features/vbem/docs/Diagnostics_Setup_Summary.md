# Count Gap Diagnostics - Setup Summary

**Date:** January 1, 2026  
**Status:** Tests configured and ready to run

---

## What Was Created

### 1. Automated Test Script
**File:** `/mnt/pikachu/STAR-Flex/test_count_gap_diagnostics.sh`

Runs three diagnostic tests:
- **Test 1:** Fixed libType (no auto-detect) - isolates detection window effects
- **Test 2:** Score filtering disabled - checks if filtering drops valid reads
- **Test 3:** Trace analysis of STAR-only transcripts - identifies drop reasons

### 2. Test Matrix Document
**File:** `/mnt/pikachu/STAR-Flex/docs/Count_Gap_Diagnostics_Test_Matrix.md`

Complete documentation including:
- Test hypotheses
- Expected outputs
- Analysis scripts
- Interpretation guide

---

## Quick Start

### Run All Tests
```bash
bash /mnt/pikachu/STAR-Flex/test_count_gap_diagnostics.sh
```

**Expected runtime:** ~15-20 minutes per test (3 tests total)

**Output location:** `/storage/production/bulk_vb_diagnostics/`

### Run Individual Tests

**Test 1 only:**
```bash
mkdir -p /storage/production/bulk_vb_diagnostics/test1_fixed_libtype
# Edit test script to run only Test 1 section
```

**Test 2 only:**
```bash
mkdir -p /storage/production/bulk_vb_diagnostics/test2_no_score_filter
# Edit test script to run only Test 2 section
```

---

## Test 1: Fixed libType (No Auto-Detect)

### Purpose
Determine if auto-detection window differences contribute to the count gap.

### What It Tests
- Baseline: Auto-detect with window=12,866
- Test: Fixed libType=ISF (no detection)

### Expected Result
If gap remains ~3,132: Auto-detect is NOT the issue  
If gap changes significantly: Auto-detect affects counts

### Key Metrics
- Read gap (should be ~3,132 if unchanged)
- EC count (should be ~20,025 if unchanged)
- Dropped incompat (should be ~2,226 if unchanged)

---

## Test 2: Score Filtering Disabled

### Purpose
Determine if STAR's score filtering drops alignments that Salmon keeps.

### What It Tests
- Baseline: Default score filtering (min_score_fraction=0.65)
- Test: Score filtering disabled (min_score_fraction=0)

### Expected Result
If EC count increases: Score filtering was reducing ECs  
If STAR-only transcripts decrease: Score filtering was creating discrepancies

### Key Metrics
- EC count (should increase if filtering was dropping reads)
- STAR-only transcripts (should decrease if filtering was the cause)
- Read gap (may decrease if filtering was dropping valid reads)

---

## Test 3: Trace STAR-Only Transcripts

### Purpose
Identify WHY 3,007 transcripts receive reads only in STAR.

### What It Tests
- Samples reads from top 100 STAR-only transcripts
- Traces these reads through STAR's EC building
- Identifies drop reasons: incompat, score filter, missing mate fields

### Expected Result
Most STAR-only reads should be due to:
- Incompatibility detection (if STAR is more aggressive)
- Score filtering (if Test 2 shows improvement)
- Missing mate fields (if STAR handles these differently)

### Key Metrics
- % dropped as incompatible
- % dropped by score filter
- % dropped due to missing mate fields
- % that are actually compatible (should be 0 if logic is correct)

---

## Analysis Scripts

### Compare Test Results
```python
import pandas as pd

# Load baseline
baseline_star = pd.read_csv("/storage/production/bulk_vb_matched_window/quant.sf", sep="\t")
baseline_salmon = pd.read_csv("/storage/production/bulk_vb_deterministic/salmon_deterministic/quant.sf", sep="\t")

# Load test results
test_star = pd.read_csv("/storage/production/bulk_vb_diagnostics/test1_fixed_libtype/star_quant.sf", sep="\t")
test_salmon = pd.read_csv("/storage/production/bulk_vb_diagnostics/test1_fixed_libtype/salmon/quant.sf", sep="\t")

# Compare gaps
baseline_gap = baseline_salmon['NumReads'].sum() - baseline_star['NumReads'].sum()
test_gap = test_salmon['NumReads'].sum() - test_star['NumReads'].sum()

print(f"Baseline gap: {baseline_gap:.0f}")
print(f"Test gap: {test_gap:.0f}")
print(f"Difference: {abs(test_gap - baseline_gap):.0f}")
```

### Analyze Trace for Drop Reasons
```python
# Parse trace file for sample reads
sample_reads = set(open('sample_reads.txt').readlines())

drop_reasons = {
    'incompat': 0,
    'compatible': 0,
    'unknown': 0
}

with open('test3/star_quant.ec_trace.tsv') as f:
    for line in f:
        qname = line.split('\t')[0]
        if qname in sample_reads:
            if 'isCompat=0' in line and 'isCompat=1' not in line:
                drop_reasons['incompat'] += 1
            elif 'isCompat=1' in line:
                drop_reasons['compatible'] += 1
            else:
                drop_reasons['unknown'] += 1

print("Drop reason distribution:")
for reason, count in drop_reasons.items():
    print(f"  {reason}: {count}")
```

---

## Next Steps After Tests Complete

1. **If Test 1 gap unchanged:**
   - Auto-detect is NOT the issue
   - Focus on filtering/EC building

2. **If Test 2 shows EC increase:**
   - Score filtering was dropping valid reads
   - Consider adjusting filtering threshold

3. **If Test 3 shows >70% incompat drops:**
   - Incompatibility detection is main issue
   - Need to match Salmon's filtering logic exactly

4. **If all tests show similar gaps:**
   - Issue is deeper (FLD training, EC building, EM resolution)
   - May need to investigate FLD divergence next

---

## Files Created

| File | Purpose |
|------|---------|
| `test_count_gap_diagnostics.sh` | Automated test execution |
| `Count_Gap_Diagnostics_Test_Matrix.md` | Complete test documentation |
| `Diagnostics_Setup_Summary.md` | This file |

---

## Status

✅ Test scripts created  
✅ Test matrix documented  
⏳ Tests running (Test 1 in progress)  
⏳ Results pending  

Check test status:
```bash
ls -la /storage/production/bulk_vb_diagnostics/test*/star_Log.final.out
```

