# remove_y_reads FASTQ Splitter Test Report

Generated: Sun Apr  5 08:36:43 AM UTC 2026

## Test Configuration

- Tool: `/mnt/pikachu/STAR-suite/tests/../core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads`
- Test Directory: `/tmp/y_removal_comprehensive_test_681787`


## Test Data

- Y BAM: `/tmp/y_removal_comprehensive_test_681787/y_reads.bam` (3 unique read names)
- Test FASTQ 1: `/tmp/y_removal_comprehensive_test_681787/test1.fastq.gz` (5 reads)
- Test FASTQ 2: `/tmp/y_removal_comprehensive_test_681787/test2.fastq.gz` (3 reads)
- Expected Y reads: read1, read3, read7
- Expected noY reads: read2, read4, read5, read6, read8

---

## Test Results


### single_threaded Mode (--threads 1)

| Check | Result |
|-------|--------|
| Output files exist | ✓ PASS |
| Count consistency | ✓ PASS |
| Y hash membership | ✓ PASS |
| noY hash membership | ✓ PASS |
| Order preservation | ✓ PASS |
| Overall counts | Y=2, noY=3, original=5 |


### multithreaded_single Mode (--threads 4)

| Check | Result |
|-------|--------|
| Output files exist | ✓ PASS |
| Count consistency | ✓ PASS |
| Y hash membership | ✓ PASS |
| noY hash membership | ✓ PASS |
| Order preservation | ✓ PASS |
| Overall counts | Y=2, noY=3, original=5 |


### multithreaded_multi Mode (--threads 2)

| Check | Result |
|-------|--------|
| Output files exist | ✓ PASS |
| Count consistency | ✓ PASS |
| Y hash membership | ✓ PASS |
| noY hash membership | ✓ PASS |
| Order preservation | ✓ PASS |
| Overall counts | Y=3, noY=5, original=8 |


---

## Summary

- **Total checks passed**: 22
- **Total checks failed**: 0
- **Status**: ✓ ALL TESTS PASSED

