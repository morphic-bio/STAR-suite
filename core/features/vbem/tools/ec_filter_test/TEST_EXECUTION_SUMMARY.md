# Error Model Test Execution Summary

## Test Date
December 19, 2025

## Test Status

### ✅ Our CLI - Working
- **Trace Generation**: ✅ Working (29,726 trace lines generated)
- **Matrix Dumps**: ✅ Working (final matrices dumped)
- **Trace Format**: ✅ Correct format
- **Model Usage**: ✅ `modelUsed=1` entries present with errLike values

### ⚠️ Salmon - Instrumentation Not Working
- **Trace Generation**: ❌ No trace file created
- **Matrix Dumps**: ❌ No matrix files created
- **Environment Variables**: ⚠️ May not be read correctly

## Test Results

### Our CLI Trace File
- **File**: `/tmp/error_model_test/our_trace.txt`
- **Lines**: 29,726 (one per read)
- **Format**: `READ <qname> numAlns=<N> errLikeSum=<X> modelUsed=<0/1> modelUpdated=<0/1>`
- **Model Usage**:
  - `modelUsed=0`: Pre-burnin reads (first ~5000 reads)
  - `modelUsed=1`: Post pre-burnin reads (remaining ~24,726 reads)
- **errLike Values**: 
  - `errLikeSum=0` for pre-burnin reads
  - Non-zero errLike values (e.g., -58.8786, -77.6697) for post pre-burnin reads

### Matrix Files Generated
- **Our CLI**: `/tmp/error_model_test/our_matrices_final_*.tsv`
  - Files for left/right reads
  - Files for bins 0-5 (6 bins total)
  - Format: TSV with transition probabilities

### Salmon Trace File
- **Status**: ❌ Not generated
- **Expected**: Should match our format
- **Issue**: Environment variables may not be read, or instrumentation not enabled

## Issues Identified

1. **Salmon Instrumentation Not Active**
   - Environment variables (`SALMON_ERROR_MODEL_TRACE`, `SALMON_TRACE_LEVEL`) set but no trace file created
   - May need to verify Salmon binary was rebuilt with instrumentation
   - May need to check if environment variables are read at runtime

2. **Matrix Dump Checkpoints**
   - Our CLI: Checkpoint dumps at 1000, 2000, 3000, 4000, 5000 reads
   - Salmon: May need checkpoint implementation (see `INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`)

## Next Steps

1. **Verify Salmon Instrumentation**
   - Check if Salmon binary includes instrumentation code
   - Verify environment variable reading in `AlignmentModel` constructor
   - Test with a minimal example

2. **Compare When Salmon Works**
   - Run trace comparison once Salmon generates trace files
   - Compare errLike values for matching reads
   - Identify divergences

3. **Matrix Comparison**
   - Compare initial matrices (should be identical)
   - Compare checkpoint matrices (should match at each checkpoint)
   - Compare final matrices (may show differences if error models diverge)

4. **Error Model Harmonization**
   - Fix any identified differences in error model computation
   - Re-run tests to verify parity

## Files Generated

### Our CLI Outputs
- `/tmp/error_model_test/our_trace.txt` - Full trace (29,726 lines)
- `/tmp/error_model_test/our_matrices_final_*.tsv` - Final transition matrices
- `/tmp/error_model_test/our_ec.txt` - EC output (116 equivalence classes)

### Salmon Outputs
- `/tmp/error_model_test/salmon_out/` - Salmon output directory
- ⚠️ No trace files generated
- ⚠️ No matrix dump files generated

## Sample Trace Entries

### Pre-Burnin (modelUsed=0)
```
READ SRR6357070.6710070 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.5843839 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
```

### Post Pre-Burnin (modelUsed=1)
```
READ SRR6357070.16033236 numAlns=1 errLikeSum=-58.8786 modelUsed=1 modelUpdated=1
READ SRR6357070.12271739 numAlns=1 errLikeSum=-77.6697 modelUsed=1 modelUpdated=1
READ SRR6357070.23037324 numAlns=1 errLikeSum=-57.4266 modelUsed=1 modelUpdated=1
```

## Test Commands Used

### Our CLI
```bash
./ec_filter_cli \
  --input /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam \
  --transcripts /tmp/nfcore_ec_parity_test/transcriptome.fasta \
  --use-error-model \
  --error-model-trace /tmp/error_model_test/our_trace.txt \
  --trace-level 2 \
  --dump-matrices /tmp/error_model_test/our_matrices \
  -o /tmp/error_model_test/our_ec.txt
```

### Salmon (Expected)
```bash
SALMON_ERROR_MODEL_TRACE=/tmp/error_model_test/salmon_trace.txt \
SALMON_TRACE_LEVEL=2 \
SALMON_DUMP_MATRICES=/tmp/error_model_test/salmon_matrices \
salmon quant \
  -t /tmp/nfcore_ec_parity_test/transcriptome.fasta \
  -l A \
  -a /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam \
  --dumpEqWeights \
  --noLengthCorrection \
  --noEffectiveLengthCorrection \
  --noFragLengthDist \
  -p 1 \
  -o /tmp/error_model_test/salmon_out
```

## Notes

- All tests run single-threaded (`-p 1`) for determinism
- Test data: 29,726 reads from nf-core test dataset
- Error model parameters: `alpha=1.0`, `readBins=6`
- Pre-burnin threshold: 5000 reads (model starts being used after this)
- Full burnin threshold: 5M reads (not reached in test data, so `modelUpdated=1` for all reads)
