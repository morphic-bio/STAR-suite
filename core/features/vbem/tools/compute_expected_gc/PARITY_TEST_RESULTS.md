# GC Bias Expected Distribution Parity Test Results

## Test Date
2025-12-20

## Test Setup
- **Salmon**: Run with `--gcBias` flag
- **Our Tool**: `compute_expected_gc --transcriptome`
- **Test Data**: `/tmp/nfcore_ec_parity_test/transcriptome.fasta` (125 transcripts)

## Results

### Salmon Output
- **Format**: Binary GCFragModel (gzipped)
- **Dimensions**: 1 conditional bin × 25 GC bins
- **Distribution**: **Uniform** (0.04 probability per bin)
- **Mean GC%**: 48.00%
- **Status**: ⚠️ Suspicious - uniform distribution suggests expected GC may not have been computed properly

### Our Tool Output
- **Format**: TSV (101 bins, 0-100% GC)
- **Distribution**: **Non-uniform** with peak around 40% GC
- **Mean GC%**: 41.68%
- **Non-zero bins**: 62 out of 101
- **Status**: ✓ Appears to be correctly computed

## Comparison

| Metric | Salmon | Our Tool | Difference |
|--------|--------|----------|------------|
| **Bins** | 25 | 101 | Different resolution |
| **Distribution** | Uniform (0.04/bin) | Non-uniform (peak ~40%) | **Significant** |
| **Mean GC%** | 48.00% | 41.68% | 6.32% difference |
| **Max diff** | - | - | 0.055 (after interpolation) |

## Analysis

### Why Salmon's Distribution is Uniform

Salmon's uniform distribution (0.04 = 1/25 per bin) suggests:

1. **Uniform Prior Not Updated**: Salmon may initialize expected GC with a uniform prior and not update it properly
2. **Missing FLD**: Expected GC computation requires fragment length distribution - if FLD wasn't available, it might fall back to uniform
3. **Computation Not Triggered**: The expected GC computation might not have run despite `--gcBias` flag

### Our Tool's Distribution

Our tool shows a realistic distribution:
- Peak around 40% GC (typical for transcriptomes)
- Non-zero values in reasonable range (16-68% GC)
- Proper normalization (sums to 1.0)

## Potential Issues

1. **Fragment Length Distribution**:
   - Salmon uses empirical FLD from sample
   - Our tool uses default normal-like FLD (mean=200, sd=80)
   - **Impact**: Different FLDs will produce different expected GC distributions

2. **Bin Resolution**:
   - Salmon: 25 bins (4% per bin)
   - Our tool: 101 bins (1% per bin)
   - **Impact**: Higher resolution should be fine, but need to ensure compatibility

3. **Salmon's Uniform Output**:
   - May indicate a bug or missing computation
   - Need to verify Salmon actually computed expected GC
   - May need to check if FLD was available

## Recommendations

1. **Investigate Salmon's Uniform Output**:
   - Check if FLD was properly loaded
   - Verify expected GC computation actually ran
   - Check Salmon logs for warnings/errors

2. **Use Salmon's FLD**:
   - Extract FLD from Salmon output (`fld.gz`)
   - Use it as input to our tool: `--fld salmon_fld.txt`
   - Re-run comparison

3. **Verify Algorithm**:
   - Review Salmon's expected GC computation code
   - Ensure our algorithm matches Salmon's approach
   - Check for differences in normalization or prior handling

4. **Test with Larger Dataset**:
   - Current test uses only 125 transcripts
   - May need more data for stable expected GC computation
   - Test with full transcriptome

## Next Steps

1. Extract and use Salmon's FLD for our tool
2. Re-run comparison with matching FLD
3. Investigate why Salmon produced uniform distribution
4. If Salmon's output is correct (uniform prior), document this as expected behavior
5. If our tool is correct, verify algorithm matches Salmon's intended behavior

## Files

- Salmon output: `/tmp/salmon_gc_test/aux_info/exp_gc.gz`
- Our output: `/tmp/our_gc_test/expected_gc.tsv`
- Comparison script: `compare_with_salmon.py`
- Test script: `test_gc_parity.sh`

## Conclusion

**Parity check INCONCLUSIVE** - Salmon's uniform distribution suggests the expected GC may not have been computed properly, or Salmon uses a uniform prior that wasn't updated. Our tool produces a realistic non-uniform distribution. Further investigation needed to determine which is correct.
