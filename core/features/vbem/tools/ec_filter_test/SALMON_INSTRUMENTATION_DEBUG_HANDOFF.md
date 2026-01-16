# Salmon Error Model Instrumentation - Debug Handoff

## Problem Summary

The error model instrumentation was implemented in Salmon (per `INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`), but **trace files and matrix dumps are not being generated** when running Salmon with the expected environment variables.

## Expected Behavior

When running Salmon with environment variables set, it should:
1. **Generate trace files** containing per-read error model information
2. **Generate matrix dump files** containing transition matrices at checkpoints

## Current Behavior

- ❌ **No trace files generated** - `SALMON_ERROR_MODEL_TRACE` environment variable set but no file created
- ❌ **No matrix dump files generated** - `SALMON_DUMP_MATRICES` environment variable set but no files created
- ✅ Salmon runs successfully and produces normal output (quantification works)

## Test Case

### Command Used
```bash
SALMON_ERROR_MODEL_TRACE=/tmp/test_salmon_trace.txt \
SALMON_TRACE_LEVEL=2 \
SALMON_DUMP_MATRICES=/tmp/test_salmon_matrices \
/mnt/pikachu/salmon/build/src/salmon quant \
  -t /tmp/nfcore_ec_parity_test/transcriptome.fasta \
  -l A \
  -a /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam \
  --dumpEqWeights \
  --noLengthCorrection \
  --noEffectiveLengthCorrection \
  --noFragLengthDist \
  -p 1 \
  -o /tmp/test_salmon_out
```

### Expected Outputs
- `/tmp/test_salmon_trace.txt` - Trace file with READ and ALN level details
- `/tmp/test_salmon_matrices_*_bin*.tsv` - Matrix dump files

### Actual Outputs
- ❌ `/tmp/test_salmon_trace.txt` - **File not created**
- ❌ `/tmp/test_salmon_matrices_*.tsv` - **No files created**
- ✅ `/tmp/test_salmon_out/` - Normal Salmon output directory (quantification works)

## Implementation Reference

According to `INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`, the instrumentation should:

1. **Read environment variables in `AlignmentModel` constructor**:
   - `SALMON_ERROR_MODEL_TRACE` - Path to trace file
   - `SALMON_TRACE_LEVEL` - Trace level (1-3, default: 1)
   - `SALMON_DUMP_MATRICES` - Prefix for matrix dump files

2. **Generate trace output** via `traceRead()`, `traceAlignment()`, `traceTransition()` methods

3. **Generate matrix dumps** via `dumpMatrices()` method

## Debugging Checklist

### 1. Verify Binary Includes Instrumentation Code

**Check**: Does the compiled binary include the instrumentation code?

```bash
# Check if trace methods exist in binary
nm /mnt/pikachu/salmon/build/src/salmon | grep -i "trace\|dumpMatrices"

# Expected: Should see symbols like:
# - setTraceOutput
# - traceRead
# - traceAlignment
# - dumpMatrices
```

**If missing**: Binary needs to be rebuilt with instrumentation code.

### 2. Verify Environment Variable Reading

**Check**: Does `AlignmentModel` constructor read environment variables?

**Location**: `/mnt/pikachu/salmon/src/AlignmentModel.cpp` constructor

**Expected code** (from `INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`):
```cpp
// Check environment variables on construction
const char* traceFile = std::getenv("SALMON_ERROR_MODEL_TRACE");
if (traceFile) {
    // Open trace file and set trace level
}

const char* dumpPrefix = std::getenv("SALMON_DUMP_MATRICES");
if (dumpPrefix) {
    matrixDumpPrefix_ = dumpPrefix;
}
```

**Action**: Verify this code exists and is executed.

### 3. Verify Error Model is Enabled

**Check**: Is the error model actually being used in alignment mode?

**Location**: `/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp`

**Expected**: Error model tracing should be enabled when:
- `useErrorModel` is true
- `useAuxParams` is true  
- `useASWithoutCIGAR` is false

**Check code** (from implementation summary):
```cpp
bool emTraceEnabled = (salmonOpts.useErrorModel && useAuxParams && !useASWithoutCIGAR);
```

**Action**: Verify `emTraceEnabled` is true during execution.

### 4. Verify Trace Context is Set

**Check**: Is `setTraceContext()` being called before computing errLike?

**Location**: `/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp` (around line ~552)

**Expected code**:
```cpp
if (emTraceEnabled) {
    std::string qname = std::string(aln->getName());
    uint32_t tid = static_cast<uint32_t>(transcriptID);
    alnMod.setTraceContext(qname, tid);
}
```

**Action**: Verify this code exists and is executed.

### 5. Verify Trace Methods are Called

**Check**: Are `traceRead()` and `traceAlignment()` methods actually being called?

**Location**: 
- `traceRead()`: `/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp` (around line ~775)
- `traceAlignment()`: `/mnt/pikachu/salmon/src/AlignmentModel.cpp` (in `logLikelihood()` method)

**Action**: Add debug output or breakpoints to verify these methods are called.

### 6. Verify File Permissions

**Check**: Can Salmon write to the specified trace file path?

```bash
# Test write permissions
touch /tmp/test_salmon_trace.txt
echo "test" > /tmp/test_salmon_trace.txt
rm /tmp/test_salmon_trace.txt
```

**Action**: Verify Salmon has write permissions to trace file directory.

### 7. Verify Error Model is Active

**Check**: Is the error model actually being used (not just enabled)?

**Note**: According to implementation summary, error model tracing is independent of the existing `traceFile` mechanism and uses a separate `emTrace_errLike` vector.

**Action**: Verify error model is computing errLike values (check logs for error model activity).

## Minimal Test Case

Create a minimal test to isolate the issue:

```bash
# 1. Create minimal BAM with 1 read
samtools view -h /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam | head -4 | samtools view -Sb - > /tmp/minimal_test.bam

# 2. Run Salmon with instrumentation
SALMON_ERROR_MODEL_TRACE=/tmp/minimal_trace.txt \
SALMON_TRACE_LEVEL=1 \
/mnt/pikachu/salmon/build/src/salmon quant \
  -t /tmp/nfcore_ec_parity_test/transcriptome.fasta \
  -l A \
  -a /tmp/minimal_test.bam \
  --dumpEqWeights \
  --noLengthCorrection \
  --noEffectiveLengthCorrection \
  --noFragLengthDist \
  -p 1 \
  -o /tmp/minimal_out

# 3. Check if trace file exists
ls -lh /tmp/minimal_trace.txt
cat /tmp/minimal_trace.txt
```

## Expected Trace Format

If instrumentation works, trace file should contain lines like:

```
READ SRR6357070.6710070 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.5843839 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
```

For trace level 2, also expect:
```
ALN SRR6357070.6710070 tid=0 fg=-123.456 bg=-125.789 errLike=2.333 bin=2
```

## Comparison with Working Implementation

Our CLI implementation (which works correctly) shows:

**Trace file**: `/tmp/error_model_test/our_trace.txt`
- ✅ 29,726 lines generated
- ✅ Correct format
- ✅ Shows model usage correctly

**Matrix files**: `/tmp/error_model_test/our_matrices_*.tsv`
- ✅ 72 files generated (6 bins × 2 orientations × 6 checkpoints)
- ✅ Correct TSV format with headers

## Files to Check

1. **`/mnt/pikachu/salmon/src/AlignmentModel.cpp`**
   - Constructor: Environment variable reading
   - `logLikelihood()`: `traceAlignment()` calls
   - `traceRead()`, `traceAlignment()`, `traceTransition()`: Implementation

2. **`/mnt/pikachu/salmon/include/AlignmentModel.hpp`**
   - Public methods: `setTraceOutput()`, `setTraceContext()`, `traceRead()`
   - Private members: `traceStream_`, `traceLevel_`, `currentQname_`, etc.

3. **`/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp`**
   - `emTraceEnabled` flag setup
   - `setTraceContext()` calls
   - `traceRead()` calls
   - `emTrace_errLike` vector population

## Known Limitations (from Implementation Summary)

1. **numObserved tracking**: Matrix dumps show `numObserved=0` because `AlignmentModel` doesn't track this
2. **Checkpoint dumps**: Matrix dumps at checkpoints (1000, 2000, etc. reads) require calling `dumpMatrices()` from `SalmonQuantifyAlignments.cpp` when `processedReads` reaches those values - **this is not yet implemented**

## Action Items for Salmon Agent

1. ✅ **Verify binary includes instrumentation** - Check if code was compiled
2. ✅ **Verify environment variable reading** - Check constructor code
3. ✅ **Verify error model is enabled** - Check `emTraceEnabled` flag
4. ✅ **Verify trace methods are called** - Add debug output
5. ✅ **Fix any issues found** - Implement missing pieces
6. ✅ **Test with minimal case** - Verify trace file generation
7. ✅ **Test with full case** - Verify full trace and matrix dumps

## Success Criteria

Instrumentation is working when:
- ✅ Trace file is created at path specified by `SALMON_ERROR_MODEL_TRACE`
- ✅ Trace file contains READ-level entries (trace level 1+)
- ✅ Trace file contains ALN-level entries (trace level 2+)
- ✅ Matrix dump files are created (if `SALMON_DUMP_MATRICES` is set)
- ✅ Trace format matches specification (see `SALMON_HANDOFF.md`)

## Reference Documents

- **Implementation Summary**: `/mnt/pikachu/salmon/INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`
- **Original Handoff**: `/mnt/pikachu/STAR-Flex/tools/ec_filter_test/SALMON_HANDOFF.md`
- **Instrumentation Guide**: `/mnt/pikachu/STAR-Flex/tools/ec_filter_test/SALMON_INSTRUMENTATION_GUIDE.md`
- **Test Results**: `/mnt/pikachu/STAR-Flex/tools/ec_filter_test/TEST_RESULTS_SUMMARY.md`

## Contact

For questions about expected behavior or format, refer to:
- Our working implementation: `/mnt/pikachu/STAR-Flex/tools/ec_filter_test/ec_filter_cli.cpp`
- Our trace output: `/tmp/error_model_test/our_trace.txt` (29,726 lines)
- Our matrix files: `/tmp/error_model_test/our_matrices_*.tsv`

## Quick Verification Command

Once fixed, this should generate a trace file:

```bash
SALMON_ERROR_MODEL_TRACE=/tmp/verify_trace.txt \
SALMON_TRACE_LEVEL=1 \
/mnt/pikachu/salmon/build/src/salmon quant \
  -t /tmp/nfcore_ec_parity_test/transcriptome.fasta \
  -l A \
  -a /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam \
  --dumpEqWeights \
  --noLengthCorrection \
  --noEffectiveLengthCorrection \
  --noFragLengthDist \
  -p 1 \
  -o /tmp/verify_out

# Should create:
ls -lh /tmp/verify_trace.txt
# Should contain trace entries:
head -10 /tmp/verify_trace.txt
```

Expected output format:
```
READ <qname> numAlns=<N> errLikeSum=<X> modelUsed=<0/1> modelUpdated=<0/1>
```
