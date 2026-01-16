# STAR Agent Next Steps: Error Model Comparison

## Status: ✅ Salmon Instrumentation Complete

The Salmon agent has successfully fixed the instrumentation issue. Trace files and matrix dumps are now being generated correctly.

## Generated Files (Ready for Comparison)

### Trace File
- **Location**: `/tmp/salmon_debug_trace.txt`
- **Status**: ✅ Generated
- **Lines**: 29,726 READ entries (one per read)
- **Format**: Matches specification

### Matrix Files
- **Location**: `/tmp/salmon_debug_matrices_*.tsv`
- **Status**: ✅ Generated
- **Count**: 72 files total
- **Checkpoints**: Includes checkpoint dumps at 1000, 2000, 5000, 10000, 20000, 50000 reads, plus final

## Our CLI Outputs (Already Generated)

### Trace File
- **Location**: `/tmp/error_model_test/our_trace.txt`
- **Lines**: 29,726 READ entries
- **Status**: ✅ Ready for comparison

### Matrix Files
- **Location**: `/tmp/error_model_test/our_matrices_*.tsv`
- **Count**: 72 files (checkpoints at 1000, 2000, 3000, 4000, 5000, plus final)
- **Status**: ✅ Ready for comparison

## Next Steps: Run Comparisons

### Step 1: Compare Trace Files

```bash
cd /mnt/pikachu/STAR-Flex/tools/ec_filter_test

# Compare trace files
python3 compare_error_model.py \
  --ours /tmp/error_model_test/our_trace.txt \
  --salmon /tmp/salmon_debug_trace.txt \
  --tolerance 1e-6 \
  --report /tmp/trace_comparison_report.txt
```

**Expected Output**:
- Comparison statistics
- Matching reads count
- Diverging reads count
- Delta statistics for errLike values

### Step 2: Compare Matrix Files

Compare matrices at each checkpoint:

```bash
# Compare initial matrices (if available)
python3 compare_matrices.py \
  --ours /tmp/error_model_test/our_matrices_left_bin0_pre1000.tsv \
  --salmon /tmp/salmon_debug_matrices_left_bin0_checkpoint1000.tsv \
  --tolerance 1e-10 \
  --output /tmp/matrix_comparison_1000.txt

# Compare at 5000 reads checkpoint
python3 compare_matrices.py \
  --ours /tmp/error_model_test/our_matrices_left_bin0_pre5000.tsv \
  --salmon /tmp/salmon_debug_matrices_left_bin0_checkpoint5000.tsv \
  --tolerance 1e-10 \
  --output /tmp/matrix_comparison_5000.txt

# Compare final matrices
python3 compare_matrices.py \
  --ours /tmp/error_model_test/our_matrices_left_bin0_final.tsv \
  --salmon /tmp/salmon_debug_matrices_left_bin0_final.tsv \
  --tolerance 1e-10 \
  --output /tmp/matrix_comparison_final.txt
```

**Note**: Salmon uses checkpoint names like `checkpoint1000`, `checkpoint5000`, etc., while our CLI uses `pre1000`, `pre5000`, etc.

### Step 3: Analyze Divergences

If comparisons show differences:

1. **Trace Divergences**:
   - Check which reads have different errLike values
   - Verify modelUsed/modelUpdated flags match
   - Check if differences are within tolerance

2. **Matrix Divergences**:
   - Identify which cells differ
   - Check if differences are systematic or random
   - Verify if differences grow over time (training divergence)

### Step 4: Identify Root Causes

Common areas to check if divergences found:

1. **CIGAR Parsing**: Verify both implementations parse CIGAR strings identically
2. **State Indexing**: Verify state calculation matches (ref_base * 9 + read_base)
3. **Transition Probabilities**: Verify probability calculations match
4. **Log-space Operations**: Verify logAdd and normalization match
5. **Read Position Binning**: Verify bin calculation matches
6. **Left/Right Determination**: Verify read orientation logic matches

## Regeneration Commands (If Needed)

### Regenerate Salmon Outputs

```bash
export SALMON_ERROR_MODEL_TRACE=/tmp/salmon_debug_trace.txt
export SALMON_TRACE_LEVEL=1
export SALMON_DUMP_MATRICES=/tmp/salmon_debug_matrices

/mnt/pikachu/salmon/build/src/salmon quant \
  -t /tmp/nfcore_ec_parity_test/transcriptome.fasta \
  -l A \
  -a /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam \
  --dumpEqWeights \
  --noLengthCorrection \
  --noEffectiveLengthCorrection \
  --noFragLengthDist \
  -p 1 \
  -o /tmp/salmon_debug_out
```

### Regenerate Our CLI Outputs

```bash
cd /mnt/pikachu/STAR-Flex/tools/ec_filter_test

./ec_filter_cli \
  --input /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam \
  --transcripts /tmp/nfcore_ec_parity_test/transcriptome.fasta \
  --use-error-model \
  --error-model-trace /tmp/error_model_test/our_trace.txt \
  --trace-level 2 \
  --dump-matrices /tmp/error_model_test/our_matrices \
  -o /tmp/error_model_test/our_ec.txt
```

## File Naming Differences

### Salmon Matrix Files
- Format: `salmon_debug_matrices_<left/right>_bin<N>_checkpoint<M>.tsv`
- Example: `salmon_debug_matrices_left_bin0_checkpoint1000.tsv`

### Our CLI Matrix Files
- Format: `our_matrices_<left/right>_bin<N>_pre<M>.tsv`
- Example: `our_matrices_left_bin0_pre1000.tsv`

**Note**: When comparing, match:
- `checkpoint1000` ↔ `pre1000`
- `checkpoint5000` ↔ `pre5000`
- `final` ↔ `final`

## Expected Comparison Results

### Trace Comparison (Ideal)
- **Total reads**: 29,726 (should match)
- **Common reads**: 29,726 (all reads should match)
- **Matching reads**: 29,726 (within tolerance)
- **Diverging reads**: 0

### Matrix Comparison (Ideal)
- **Mean difference**: < 1e-10
- **Max difference**: < 1e-10
- **Matching cells**: 100% (6724 cells per matrix = 82×82)

## Troubleshooting

### If Trace Files Don't Match

1. **Check read order**: Reads may be processed in different order
   - Solution: Sort by qname before comparison
2. **Check errLike values**: May differ due to floating-point precision
   - Solution: Use tolerance (1e-6 recommended)
3. **Check modelUsed flags**: May differ if burn-in thresholds differ
   - Solution: Verify both use same thresholds (5000 pre-burnin)

### If Matrix Files Don't Match

1. **Check numObserved**: May differ if checkpoint timing differs
   - Solution: Compare matrices at same numObserved values
2. **Check initialization**: Initial matrices should match exactly
   - Solution: Compare initial matrices first
3. **Check training**: Differences may accumulate during training
   - Solution: Compare at early checkpoints to find divergence point

## Success Criteria

✅ **Trace Comparison**: 
- All reads match within tolerance
- errLike values match within 1e-6
- modelUsed/modelUpdated flags match

✅ **Matrix Comparison**:
- Initial matrices match exactly (< 1e-10)
- Checkpoint matrices match within tolerance
- Final matrices match within tolerance

## Files Reference

### Comparison Scripts
- `/mnt/pikachu/STAR-Flex/tools/ec_filter_test/compare_error_model.py` - Trace comparison
- `/mnt/pikachu/STAR-Flex/tools/ec_filter_test/compare_matrices.py` - Matrix comparison

### Output Directories
- `/tmp/error_model_test/` - Our CLI outputs
- `/tmp/salmon_debug_*` - Salmon outputs

### Documentation
- `SALMON_INSTRUMENTATION_COMPLETE.md` - Salmon instrumentation status
- `ERROR_MODEL_TESTING_GUIDE.md` - Testing guide
- `TEST_RESULTS_SUMMARY.md` - Previous test results

## Next Actions

1. ✅ **Run trace comparison** - Compare `/tmp/error_model_test/our_trace.txt` vs `/tmp/salmon_debug_trace.txt`
2. ✅ **Run matrix comparisons** - Compare matrices at each checkpoint
3. ✅ **Analyze results** - Identify any divergences
4. ✅ **Fix issues** - Harmonize any differences found
5. ✅ **Re-test** - Verify parity after fixes
