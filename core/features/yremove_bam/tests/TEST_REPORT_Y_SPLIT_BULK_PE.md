# Bulk Paired-End Y/NoY Split Validation Test Report

Generated: Wed Dec 24 06:32:30 AM UTC 2025

## Test Configuration

- STAR Binary: `/mnt/pikachu/STAR-Flex/tests/../source/STAR`
- Dataset: `/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R1_001.fastq.gz`, `/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R2_001.fastq.gz`
- Reference: `/storage/flex_filtered_reference/star_index`
- Base Directory: `/tmp/ychrom_bulk_pe_test`
- Mode: Bulk Paired-End RNA-seq (NOT single-cell/Flex)

## Baseline Statistics

- Total reads: `270089`
- chrY reads: `570`

## Split Run Results

### Unsorted Mode

- Y reads: `5639`
- noY reads: `264450`
- Total: `270089`

### Sorted Mode

- Y reads: `5639`
- noY reads: `264450`
- Total: `270089`

## Validation Results

- Passed: `9`
- Failed: `0`

### Checks Performed

1. Count consistency: Y + noY = baseline total
2. noY exclusivity: Zero chrY reads in `_noY.bam`
3. Y exclusivity: All reads in `_Y.bam` have chrY alignments (direct or via mate)
4. Mate consistency: No read names appear in both Y and noY files
5. Mate-checking: Mates route together (if one mate is on Y, both go to _Y.bam)

## Files

- Baseline BAM: `/tmp/ychrom_bulk_pe_test/baseline/Aligned.out.bam`
- Unsorted Y: `/tmp/ychrom_bulk_pe_test/unsorted_split/Aligned.out_Y.bam`
- Unsorted noY: `/tmp/ychrom_bulk_pe_test/unsorted_split/Aligned.out_noY.bam`
- Sorted Y: `/tmp/ychrom_bulk_pe_test/sorted_split/Aligned.sortedByCoord.out_Y.bam`
- Sorted noY: `/tmp/ychrom_bulk_pe_test/sorted_split/Aligned.sortedByCoord.out_noY.bam`

