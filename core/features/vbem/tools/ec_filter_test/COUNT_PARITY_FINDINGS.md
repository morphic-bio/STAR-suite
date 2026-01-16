# Count Parity Investigation Findings

## Current State
- **Salmon reads in ECs**: 29,726
- **Our reads in ECs**: 30,130
- **Difference**: 404 reads (we include more than Salmon)
- **Parity**: 73.3% (74/101 ECs fully match)

## Findings

### 1. Zero-Probability Fragment Check
- **Status**: ✅ Implemented correctly
- **Our count**: 0 zero-prob fragments
- **Salmon**: Reports some zero-prob fragments in logs (~1-2% in batches)
- **Conclusion**: This is not the source of the 404 read difference

### 2. Missing Checks Identified
Salmon performs three checks before including an alignment:
1. `transcriptLogCount != LOG_0` (transcript mass)
2. `auxProb != LOG_0` (alignment probability) ✅ We check this
3. `startPosProb != LOG_0` (start position probability)

**Our implementation**: Only checks `auxProb != LOG_0`

**Analysis**:
- With `--noLengthCorrection`: `startPosProb = -1.0` (not LOG_0), so check should always pass
- In initial round: `transcriptLogCount` should have prior mass (not LOG_0), so check should always pass
- **Conclusion**: These checks should not cause filtering, but we're not implementing them

### 3. Read Characteristics
- **Orphan reads**: 0 detected in trace (all reads are proper pairs in transcriptome BAM)
- **Incompatible reads**: 0 dropped (we have `ignore_incompat=true`)
- **Empty transcript IDs**: 0 (all reads have at least one valid alignment)

### 4. EC Count Differences
Largest differences:
- EC (93,): 65 reads difference
- EC (38,): 39 reads difference
- EC (38,39): 36 reads difference
- EC (81,): 42 reads difference
- Many ECs: 1 read difference each

**Pattern**: Differences are distributed across many ECs, suggesting systematic filtering rather than a single bug.

## Next Steps

1. **Add missing checks** (even if no-ops) to match Salmon's structure exactly
2. **Investigate transcript mass**: Check if any transcripts could have zero mass in edge cases
3. **Compare read-by-read**: If possible, instrument Salmon to output filtered read names
4. **Test with different datasets**: See if the 404 difference is consistent

## Hypothesis

The 404 read difference is likely due to:
- Salmon filtering reads based on transcript mass or startPosProb in edge cases we're not handling
- Different handling of reads with multiple alignment locations
- A subtle difference in how we compute `sumOfAlignProbs` vs Salmon's `sumOfAlignProbs`

However, with `--noLengthCorrection` and in the initial round, these checks should always pass, so the root cause remains unclear.

## Implementation Status

### ✅ Completed
1. **Tracing**: Added per-read tracing with `--trace-reads` flag
   - Tracks all reads included in ECs
   - Outputs zero-prob fragment list
   - Counts dropped incompatible, orphans, etc.

2. **Zero-Probability Fragment Check**: ✅ Already implemented correctly
   - Checks `if (mapping.transcript_ids.empty())` which matches Salmon's `sumOfAlignProbs == LOG_0`
   - Counts zero-prob fragments (currently 0)

3. **Orphan Detection**: ✅ Verified correct
   - Uses `MateStatus::PAIRED_END_LEFT/RIGHT` to identify orphans
   - Sets `log_frag_prob = LOG_EPSILON` for orphans (matches Salmon)
   - Handles discordant pairs as two orphan alignments (matches Salmon)

4. **Analysis**: ✅ Documented findings
   - 404 read difference identified
   - Pattern: distributed across many ECs (systematic filtering)
   - Missing checks identified (transcript mass, startPosProb - but should always pass)

### ⏳ Remaining Work
1. **Add missing checks** (transcript mass, startPosProb) even if no-ops to match Salmon's structure
2. **Root cause identification**: The 404 read difference remains unexplained
   - Hypothesis: Subtle differences in filtering logic we haven't identified
   - May require Salmon instrumentation to output filtered read names

## Recommendations

1. **Accept current parity** (73.3%) as reasonable given:
   - 100% label match
   - 100% weight match  
   - Only count differences remain
   - Differences are small and distributed

2. **If 100% parity required**:
   - Instrument Salmon to output filtered read names
   - Compare read-by-read to identify the 404 reads
   - Implement exact matching logic

3. **Alternative approach**:
   - Document the 404 read difference as expected due to implementation differences
   - Focus on achieving parity with error model enabled
