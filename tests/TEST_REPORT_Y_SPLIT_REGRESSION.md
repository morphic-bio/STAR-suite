# Y/NoY Split Regression Test Report

Generated: Wed Dec 24 06:32:30 AM UTC 2025

## Overview

This report validates that STAR-Flex baseline behavior matches upstream STAR Solo, and that Y/noY split outputs are stable across test runs.

---

## Test 1: Bulk RNA-Seq Baseline Regression

### Configuration

- STAR-Flex Binary: `/mnt/pikachu/STAR-Flex/tests/../source/STAR`
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

## Test 2: Flex Regression (Previous vs Current)

### Configuration

- Previous artifacts: `/tmp/ychrom_flex_test`
- Current run: `/tmp/ychrom_flex_test_current`
- Dataset: `/storage/downsampled/SC2300771`

### Results

| Metric | Previous | Current | Match |
|--------|----------|---------|-------|
| Y reads | 499 | 499 | ✓ |
| noY reads | 150790 | 150790 | ✓ |
| Y idxstats | - | - | ✓ |
| noY idxstats | - | - | ✓ |

### Conclusion

- Y read counts: PASS
- noY read counts: PASS
- Y idxstats: PASS
- noY idxstats: PASS


---

## Summary

- Bulk RNA-seq regression checks passed: `4`
- Bulk RNA-seq regression checks failed: `0`
- Flex regression checks passed: `4`
- Flex regression checks failed: `0`
- **Total passed: `8`**
- **Total failed: `0`**

