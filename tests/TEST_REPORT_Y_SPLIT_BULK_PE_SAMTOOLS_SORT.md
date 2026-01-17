# Bulk Paired-End Y/NoY Split Test Report (samtools sorter)

Generated: Tue Dec 23 07:02:54 PM UTC 2025

## Test Configuration

- STAR Binary: `/mnt/pikachu/STAR-Flex/tests/../source/STAR`
- Dataset: `/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R1_001.fastq.gz`, `/storage/PE/downsampled/21033-09-01-13-01_S1_L007_R2_001.fastq.gz`
- Reference: `/storage/flex_filtered_reference/star_index`
- Base Directory: `/tmp/ychrom_bulk_pe_samtools_sort_test`
- Mode: Bulk Paired-End RNA-seq with samtools sorter
- Sort Method: `--outBAMsortMethod samtools`

## Baseline Statistics

- Total reads: `270089`
- chrY reads: `570`

## Split Run Results

### Sorted Mode (samtools sorter)

- Y reads: `5639`
- noY reads: `264450`
- Total: `270089`

## Validation Results

- Passed: `5`
- Failed: `0`

### Checks Performed

1. Count consistency: Y + noY = baseline total
2. noY exclusivity: Zero chrY reads in `_noY.bam`
3. Y exclusivity: All reads in `_Y.bam` have chrY alignments (direct or via mate)
4. Mate consistency: No read names appear in both Y and noY files
5. Mate-checking: Mates route together (if one mate is on Y, both go to _Y.bam)

## Files

- Baseline BAM: `/tmp/ychrom_bulk_pe_samtools_sort_test/baseline/Aligned.sortedByCoord.out.bam`
- Sorted Y: `/tmp/ychrom_bulk_pe_samtools_sort_test/sorted_split/Aligned.sortedByCoord.out_Y.bam`
- Sorted noY: `/tmp/ychrom_bulk_pe_samtools_sort_test/sorted_split/Aligned.sortedByCoord.out_noY.bam`

