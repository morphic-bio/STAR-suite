# Error Model Test Results Summary

## Test Execution Date
December 19, 2025

## Test Status: ✅ Our CLI Working, ⚠️ Salmon Instrumentation Not Active

### Our CLI Status: ✅ FULLY FUNCTIONAL

#### Trace Generation
- **Status**: ✅ Working perfectly
- **File**: `/tmp/error_model_test/our_trace.txt`
- **Lines**: 29,726 (one per read)
- **Format**: Correct (`READ <qname> numAlns=<N> errLikeSum=<X> modelUsed=<0/1> modelUpdated=<0/1>`)
- **Model Usage Breakdown**:
  - `modelUsed=0`: 5,000 reads (pre-burnin phase)
  - `modelUsed=1`: 24,726 reads (post pre-burnin phase)
- **errLike Values**:
  - Pre-burnin: `errLikeSum=0` (model not used)
  - Post pre-burnin: Non-zero values (e.g., -58.8786, -77.6697, -59.3901)

#### Matrix Dumps
- **Status**: ✅ Working perfectly
- **Files Generated**: 72 matrix files total
  - 6 bins (0-5) × 2 orientations (left/right) × 6 checkpoints (pre1000, pre2000, pre3000, pre4000, pre5000, final)
- **File Naming**: `our_matrices_<left/right>_bin<N>_<checkpoint>.tsv`
- **Checkpoints**: Successfully dumped at 1000, 2000, 3000, 4000, 5000 reads, and final
- **File Size**: ~132KB per matrix file (82×82 transition matrix)

#### Sample Trace Entries

**Pre-Burnin (modelUsed=0)**:
```
READ SRR6357070.6710070 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.5071605 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
```

**Post Pre-Burnin (modelUsed=1)**:
```
READ SRR6357070.16033236 numAlns=1 errLikeSum=-58.8786 modelUsed=1 modelUpdated=1
READ SRR6357070.37604693 numAlns=1 errLikeSum=-59.3901 modelUsed=1 modelUpdated=1
READ SRR6357070.40473550 numAlns=1 errLikeSum=-57.9368 modelUsed=1 modelUpdated=1
```

### Salmon Status: ⚠️ INSTRUMENTATION NOT WORKING

#### Trace Generation
- **Status**: ❌ Not generating trace files
- **Issue**: Environment variables (`SALMON_ERROR_MODEL_TRACE`, `SALMON_TRACE_LEVEL`) set but no trace file created
- **Possible Causes**:
  1. Salmon binary not rebuilt with instrumentation code
  2. Environment variables not read correctly in `AlignmentModel` constructor
  3. Error model not enabled/used in alignment mode

#### Matrix Dumps
- **Status**: ❌ Not generating matrix dump files
- **Issue**: Environment variable (`SALMON_DUMP_MATRICES`) set but no files created
- **Note**: According to `INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`, checkpoint dumps require additional implementation

## Test Data

- **BAM File**: `/tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam`
- **Transcriptome**: `/tmp/nfcore_ec_parity_test/transcriptome.fasta`
- **Total Reads**: 29,726
- **ECs Built**: 116 equivalence classes

## Error Model Parameters

- **Alpha**: 1.0 (Salmon default)
- **Read Bins**: 6 (Salmon default)
- **Pre-Burnin Threshold**: 5,000 reads
- **Full Burnin Threshold**: 5,000,000 reads (not reached in test data)

## Files Generated

### Our CLI Outputs
```
/tmp/error_model_test/
├── our_trace.txt                    # Full trace (29,726 lines)
├── our_matrices_left_bin0_pre1000.tsv
├── our_matrices_left_bin0_pre2000.tsv
├── our_matrices_left_bin0_pre3000.tsv
├── our_matrices_left_bin0_pre4000.tsv
├── our_matrices_left_bin0_pre5000.tsv
├── our_matrices_left_bin0_final.tsv
├── ... (similar for bins 1-5 and right orientation)
└── our_ec.txt                       # EC output (116 classes)
```

### Salmon Outputs
```
/tmp/error_model_test/
└── salmon_out/                      # Standard Salmon output
    └── (no trace or matrix files)
```

## Next Steps

### Immediate Actions Required

1. **Verify Salmon Instrumentation**
   - Check if `/mnt/pikachu/salmon/build/src/salmon` includes instrumentation code
   - Verify `AlignmentModel` constructor reads environment variables
   - Test with minimal example to isolate issue

2. **Once Salmon Works**
   - Run trace comparison: `compare_error_model.py`
   - Run matrix comparison: `compare_matrices.py`
   - Identify divergences in errLike computation
   - Identify divergences in transition matrices

3. **Harmonization**
   - Fix identified differences
   - Re-run tests to verify parity
   - Document any intentional differences

## Test Commands

### Our CLI (Working)
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

### Salmon (Expected - Not Working Yet)
```bash
SALMON_ERROR_MODEL_TRACE=/tmp/error_model_test/salmon_trace.txt \
SALMON_TRACE_LEVEL=2 \
SALMON_DUMP_MATRICES=/tmp/error_model_test/salmon_matrices \
/mnt/pikachu/salmon/build/src/salmon quant \
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

## Observations

1. **Our Implementation**: Fully functional with correct trace and matrix dump generation
2. **Pre-Burnin Behavior**: Correctly shows `modelUsed=0` and `errLikeSum=0` for first 5,000 reads
3. **Post Pre-Burnin Behavior**: Correctly shows `modelUsed=1` and non-zero errLike values
4. **Matrix Evolution**: Matrices dumped at each checkpoint show model training progression
5. **Salmon Blocking Issue**: Instrumentation not active, preventing comparison

## Conclusion

Our CLI error model instrumentation is **fully functional and ready for comparison**. Once Salmon's instrumentation is verified and working, we can proceed with detailed comparison and harmonization.
