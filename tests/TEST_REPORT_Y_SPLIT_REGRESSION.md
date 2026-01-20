# Y/NoY Split Regression Test Report

Generated: Tue Jan 20 10:09:11 PM UTC 2026

## Overview

This report validates that STAR-Flex baseline behavior matches upstream STAR Solo, and that Y/noY split outputs are stable across test runs.

---

## Test 1: Bulk RNA-Seq Baseline Regression

### Configuration

- STAR-Flex Binary: `/mnt/pikachu/STAR-suite/core/legacy/source/STAR.release`
- STAR Solo Binary: `/usr/local/bin/STAR`
- Dataset: `/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R1_001.fastq.gz`
- Reference: `/storage/flex_filtered_reference/star_index`

### Test Procedure

1. Run STAR-Flex baseline (no --emitNoYBAM) - unsorted and sorted
2. Run STAR Solo baseline (no split) - unsorted and sorted
3. Compare total mapped reads and chrY counts


### Results

#### Unsorted BAM

| Metric | STAR-Flex | STAR Solo | Match |
|--------|-----------|-----------|-------|
| Total reads | 270089 | 270089 | ✓ |
| chrY reads | 570 | 570 | ✓ |

#### Sorted BAM

| Metric | STAR-Flex | STAR Solo | Match |
|--------|-----------|-----------|-------|
| Total reads | 270089 | 270089 | ✓ |
| chrY reads | 570 | 570 | ✓ |

### Conclusion

- Unsorted total reads: PASS
- Unsorted chrY reads: PASS
- Sorted total reads: PASS
- Sorted chrY reads: PASS


---

## Test 2: Flex Regression - BASELINE CREATED

Previous artifacts not found. Created new baseline at `/tmp/ychrom_flex_test_current`:

- Y reads: 499
- noY reads: 150790

Future runs will compare against this baseline.


---

## Summary

- Bulk RNA-seq regression checks passed: `4`
- Bulk RNA-seq regression checks failed: `0`
- Flex regression checks passed: `0`
- Flex regression checks failed: `0`
- **Total passed: `4`**
- **Total failed: `0`**

