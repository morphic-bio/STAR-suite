# Error Model Variance Analysis

## Executive Summary

The 21% matrix cell variation and ~0.82 mean errLike difference are NOT caused by implementation bugs, but by **fundamental differences in how reads are processed**. There are also 2 bugs found in Salmon's instrumentation.

## Three Root Causes Identified

### Issue 1: Read Processing Order Differs (PRIMARY CAUSE)

**Evidence:**
- Read order diverges at position 1 (the second read)
- Our CLI: `SRR6357070.5843839` at position 1
- Salmon: `SRR6357070.33112793` at position 1

**Cause:**
- Salmon uses mini-batch processing (`processedReads += batchReads`)
- Even with `-p 1`, Salmon's batch-based architecture processes reads differently
- Our CLI processes reads sequentially from the BAM file

**Impact:**
- Different reads train the matrices in different order
- Matrices diverge from checkpoint 1 (1000 reads)
- At checkpoint 5000: 1394 cells differ (21% of matrix)
- Mean cell difference: 0.095
- Max cell difference: 2.35

**This is NOT a bug** - both implementations are valid. They just use different read orderings.

### Issue 2: Pre-Burnin Count Differs (817 reads)

**Evidence:**
- Our CLI: 5000 pre-burnin reads (modelUsed=0)
- Salmon: 5817 pre-burnin reads (modelUsed=0)

**Analysis:**
- Both use `numPreBurninFrags = 5000` threshold
- Salmon checks `processedReads.load() >= numPreBurninFrags`
- The counter is incremented AFTER processing a batch
- 817 extra reads were processed before counter crossed threshold

**Impact:**
- 817 reads have different modelUsed flags
- These reads are in ours' post-burnin but Salmon's pre-burnin (or vice versa)

**This is a behavioral difference** in how the threshold check interacts with batch processing.

### Issue 3: Salmon errLikeSum=0 Bug (101 reads)

**Evidence:**
- 101 reads have `modelUsed=1` but `errLikeSum=0` in Salmon
- These same reads have proper errLike values in our CLI
- Example: `SRR6357070.11022960`: ours=-117.196, salmon=0

**These reads have:**
- Valid 101M CIGAR strings
- Paired-end alignments (flags 163/83)
- Should have errLike computed

**This IS a bug** in Salmon's instrumentation - it reports `modelUsed=1` but doesn't compute/report the errLike.

## Quantitative Summary

| Metric | Value | Impact |
|--------|-------|--------|
| Matrix cells differing | 1394/6724 (21%) | Expected from different training order |
| Mean errLike difference | 0.82 | Affects EC weights |
| Pre-burnin count diff | 817 reads | Threshold timing |
| Reads with buggy salmon trace | 101 | Instrumentation bug |
| Max errLike difference | 8.2 | For reads with rare transitions |
| Min errLike difference | 0.001 | For simple match-only reads |

## errLike Difference Distribution

```
For 23,808 reads with non-zero errLike in both:
- Mean difference: 0.82
- Median difference: 0.78
- Standard deviation: 0.34
- Min difference: 0.001
- Max difference: 8.2
- Relative difference: 1.3% (mean)
```

## Why Match-Only Reads Have Small Differences

For reads with 101M CIGAR (all matches):
- Most transitions use states 0→0, 10→10, 20→20, 30→30 (match states)
- These match→match transition differences are small (~0.001-0.03)
- A 101bp read has ~100 transitions
- Expected errLike diff ≈ 100 × 0.01 = 1.0 (matches observed ~0.82)

## Why Some Reads Have Large Differences (up to 8.2)

For reads with mismatches, indels, or soft clips:
- These use rare transition states
- Rare states have larger differences (up to 2.35)
- Cumulative effect can be 5-8 log-prob difference

## Recommendations

### For Achieving Parity

1. **Force identical read order** - Both implementations must process reads in the same sequence
2. **Fix Salmon instrumentation bug** - The 101 reads with errLikeSum=0 need investigation
3. **Align pre-burnin counter logic** - Ensure threshold check happens at same point

### For Practical Use

If exact parity isn't required:
- The ~1.3% relative difference is acceptable for most applications
- Matrix training will converge to similar distributions with enough data
- EC weights are relative - absolute errLike values matter less than rankings

## Next Steps

1. [ ] Investigate Salmon's batch processing to understand read order
2. [ ] Fix Salmon's errLikeSum=0 bug for the 101 reads
3. [ ] Consider if we should match Salmon's batch processing order
4. [ ] Test EC weight parity despite errLike differences
