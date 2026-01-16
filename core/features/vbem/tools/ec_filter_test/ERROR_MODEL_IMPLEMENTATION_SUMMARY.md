# Error Model Harmonization Implementation Summary

## Status: Implementation Complete (Our Side)

All instrumentation and comparison tools for our CLI have been implemented. Salmon instrumentation guide has been created but requires manual implementation.

## Completed Components

### 1. Our CLI - Multi-Level Trace System ✅

**Files Modified**:
- `source/libem/alignment_model.h` - Added trace support, ErrorModelTraceLevel enum
- `source/libem/alignment_model.cpp` - Implemented trace methods, matrix dump methods
- `tools/ec_filter_test/ec_filter_cli.cpp` - Added CLI arguments and integration

**Features**:
- `--error-model-trace <file>` - Enable error model tracing
- `--trace-level <1|2|3>` - Control verbosity:
  - Level 1: Per-read summary (`READ qname numAlns=X errLikeSum=Y modelUsed=0/1 modelUpdated=0/1`)
  - Level 2: Per-alignment details (`ALN qname tid=X fg=Y bg=Z errLike=W bin=B`)
  - Level 3: Per-transition details (`TRANS qname tid=X readPos=P bin=B prevState=S1 curState=S2 logProb=X`)
- `--dump-matrices <prefix>` - Dump transition matrices at checkpoints

**Trace Integration Points**:
- `setTraceContext()` called before computing errLike for each read
- `traceAlignment()` called after computing log-likelihood
- `traceTransition()` called for each state transition (level 3)
- `traceRead()` called after processing each read group
- `clearTraceContext()` called after each read

**Matrix Dump Checkpoints**:
- Automatically dumps at 1000, 2000, 3000, 4000, 5000 reads (pre-burnin)
- Dumps final matrices at end
- Format: `<prefix>_<checkpoint>_<left/right>_bin<N>.tsv`

### 2. Comparison Scripts ✅

**Created Files**:
- `tools/ec_filter_test/compare_error_model.py` - Compares trace files
  - Parses READ and ALN lines
  - Identifies first divergence
  - Reports errLike delta statistics
  - Generates comparison report

- `tools/ec_filter_test/compare_matrices.py` - Compares matrix dumps
  - Loads TSV matrix files
  - Computes per-cell differences
  - Identifies cells with largest divergence
  - Reports match percentage

### 3. Constants Verification ✅

**Created File**:
- `tools/ec_filter_test/verify_constants.md` - Documents constant verification

**Verified Matches**:
- ✅ State constants (ALN_A=0, ALN_C=1, etc.) - All match
- ✅ Start state index (81) - Matches
- ✅ Two-bit encoding array - Matches exactly
- ✅ Number of states (82 total, 9 base states) - Matches

**Potential Differences**:
- ⚠️ Left/right read determination: We use `<=`, Salmon uses `<` when positions equal
- ⚠️ Alpha default: Need to verify Salmon's default
- ⚠️ Read bins: Need to verify if Salmon uses 6 by default

### 4. Documentation ✅

**Created Files**:
- `tools/ec_filter_test/SALMON_INSTRUMENTATION_GUIDE.md` - Complete guide for adding matching instrumentation to Salmon
- `tools/ec_filter_test/ERROR_MODEL_TESTING_GUIDE.md` - Step-by-step testing workflow
- `tools/ec_filter_test/verify_constants.md` - Constants verification results

## Salmon Instrumentation Status

**Status**: Guide created, requires manual implementation

The `SALMON_INSTRUMENTATION_GUIDE.md` contains complete instructions for:
1. Adding trace support to `AlignmentModel.hpp` and `AlignmentModel.cpp`
2. Environment variable support (`SALMON_ERROR_MODEL_TRACE`, `SALMON_TRACE_LEVEL`, `SALMON_DUMP_MATRICES`)
3. Integration points in `SalmonQuantifyAlignments.cpp`

**Note**: Salmon files are outside the workspace, so these changes need to be applied manually.

## Usage Examples

### Our CLI with Full Instrumentation

```bash
./ec_filter_cli --input test.bam --transcripts transcriptome.fa \
    --use-error-model \
    --error-model-trace /tmp/our_trace.txt \
    --trace-level 2 \
    --dump-matrices /tmp/our_matrices \
    -o /tmp/our_ec.txt
```

### Comparison Workflow

```bash
# Compare traces
python3 compare_error_model.py \
    --ours /tmp/our_trace.txt \
    --salmon /tmp/salmon_trace.txt \
    --tolerance 1e-6 \
    --report /tmp/comparison.txt

# Compare matrices
python3 compare_matrices.py \
    --ours /tmp/our_matrices_final_left_bin0.tsv \
    --salmon /tmp/salmon_matrices_final_left_bin0.tsv \
    --tolerance 1e-10
```

## Next Steps

1. **Apply Salmon Instrumentation**: Follow `SALMON_INSTRUMENTATION_GUIDE.md` to add matching trace support
2. **Run Baseline Tests**: Execute Stage 1-2 tests from `ERROR_MODEL_TESTING_GUIDE.md`
3. **Identify Divergence**: Use Stage 3-5 tests to pinpoint where matrices/errLike diverge
4. **Fix Issues**: Address any differences found (likely left/right determination, alpha defaults, etc.)
5. **Verify Parity**: Confirm 99%+ parity with error model enabled

## Key Implementation Details

### Trace Context Management

Our CLI sets trace context before computing errLike:
```cpp
alignment_model->setTraceContext(qname, tid);
// ... compute errLike ...
alignment_model->traceRead(qname, numAlns, errLikeSum, modelUsed, modelUpdated);
alignment_model->clearTraceContext();
```

### Matrix Dump Timing

Matrices are dumped:
- After `incrementObserved()` so counts are accurate
- At checkpoints: 1000, 2000, 3000, 4000, 5000 reads
- At end: final matrices

### Trace Format

All trace formats match exactly between our CLI and Salmon (when Salmon instrumentation is added):
- `READ` lines: Per-read summary
- `ALN` lines: Per-alignment errLike components
- `TRANS` lines: Per-transition probabilities

This enables direct line-by-line comparison using `compare_error_model.py`.

## Files Summary

**Modified**:
- `source/libem/alignment_model.h` - Trace support
- `source/libem/alignment_model.cpp` - Trace implementation
- `tools/ec_filter_test/ec_filter_cli.cpp` - CLI integration

**Created**:
- `tools/ec_filter_test/compare_error_model.py` - Trace comparison
- `tools/ec_filter_test/compare_matrices.py` - Matrix comparison
- `tools/ec_filter_test/SALMON_INSTRUMENTATION_GUIDE.md` - Salmon guide
- `tools/ec_filter_test/ERROR_MODEL_TESTING_GUIDE.md` - Testing guide
- `tools/ec_filter_test/verify_constants.md` - Constants verification
- `tools/ec_filter_test/ERROR_MODEL_IMPLEMENTATION_SUMMARY.md` - This file

## Testing Readiness

All components are ready for testing once Salmon instrumentation is added. The testing guide provides step-by-step instructions for:
1. Baseline verification
2. Single-read test
3. Pre-burnin phase verification
4. Post pre-burnin errLike comparison
5. Full run comparison

Each stage builds on the previous one to systematically identify any divergence points.
