# Prompt for Salmon Agent: Debug Error Model Instrumentation

## Task

Debug why Salmon's error model instrumentation is not generating trace files when environment variables are set. The instrumentation code exists (per `INSTRUMENTATION_IMPLEMENTATION_SUMMARY.md`), but trace files are not being created.

## Test Command to Run

```bash
# Set up test environment
export SALMON_ERROR_MODEL_TRACE=/tmp/salmon_debug_trace.txt
export SALMON_TRACE_LEVEL=1
export SALMON_DUMP_MATRICES=/tmp/salmon_debug_matrices

# Run Salmon with test data
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

## Expected Outcome (Success)

After running the command above, you should see:

### 1. Trace File Created
```bash
ls -lh /tmp/salmon_debug_trace.txt
# Expected: File exists with non-zero size
```

### 2. Trace File Contains Entries
```bash
head -10 /tmp/salmon_debug_trace.txt
```

**Expected output format**:
```
READ SRR6357070.6710070 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.5843839 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.15476987 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.39609240 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.12656563 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
```

**Format specification**:
- Line format: `READ <qname> numAlns=<N> errLikeSum=<X> modelUsed=<0/1> modelUpdated=<0/1>`
- Values are space-separated
- `qname` is the read name from BAM
- `numAlns` is number of alignments for this read
- `errLikeSum` is sum of errLike values (log space, may be 0 for pre-burnin)
- `modelUsed` is 1 if error model is being used for errLike computation, 0 otherwise
- `modelUpdated` is 1 if error model can be updated, 0 otherwise

### 3. Trace File Has Many Entries
```bash
wc -l /tmp/salmon_debug_trace.txt
# Expected: Should have ~29,726 lines (one per read in test data)
```

### 4. Matrix Files Created (if SALMON_DUMP_MATRICES set)
```bash
ls -lh /tmp/salmon_debug_matrices_*.tsv
# Expected: Multiple TSV files (one per bin/orientation)
```

## Current Behavior (Failure)

Currently, when running the test command:
- ❌ `/tmp/salmon_debug_trace.txt` - **File not created**
- ❌ `/tmp/salmon_debug_matrices_*.tsv` - **No files created**
- ✅ `/tmp/salmon_debug_out/` - Normal Salmon output (quantification works)

## Debugging Steps

### Step 1: Verify Binary Includes Instrumentation

```bash
# Check if trace methods exist in binary
nm /mnt/pikachu/salmon/build/src/salmon | grep -i "trace\|dumpMatrices\|isTraceEnabled"

# Expected symbols:
# - setTraceOutput
# - traceRead
# - traceAlignment
# - dumpMatrices
# - isTraceEnabled (if this method exists)
```

**If missing**: Binary needs to be rebuilt with instrumentation code.

### Step 2: Verify Environment Variables Are Read

**Check**: `/mnt/pikachu/salmon/src/AlignmentModel.cpp` constructor

**Expected**: Should read environment variables:
```cpp
const char* traceFile = std::getenv("SALMON_ERROR_MODEL_TRACE");
if (traceFile) {
    // Open trace file and call setTraceOutput()
}

const char* traceLevel = std::getenv("SALMON_TRACE_LEVEL");
// Parse and set trace level

const char* dumpPrefix = std::getenv("SALMON_DUMP_MATRICES");
if (dumpPrefix) {
    matrixDumpPrefix_ = dumpPrefix;
}
```

**Action**: Add debug output to verify environment variables are read:
```cpp
std::cerr << "DEBUG: SALMON_ERROR_MODEL_TRACE = " 
          << (traceFile ? traceFile : "NULL") << std::endl;
```

### Step 3: Verify isTraceEnabled() Method

**Check**: Does `AlignmentModel` have an `isTraceEnabled()` method?

**Location**: `/mnt/pikachu/salmon/include/AlignmentModel.hpp`

**Expected**: Public method that returns whether tracing is enabled:
```cpp
bool isTraceEnabled() const { return traceStream_ != nullptr; }
```

**If missing**: Add this method to check if trace stream is open.

### Step 4: Verify Trace Context is Set

**Check**: `/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp`

**Expected**: `setTraceContext()` should be called before computing errLike:
```cpp
if (alnMod.isTraceEnabled()) {
    std::string qname = std::string(aln->getName());
    uint32_t tid = static_cast<uint32_t>(transcriptID);
    alnMod.setTraceContext(qname, tid);
}
```

**Action**: Add debug output to verify this is called:
```cpp
if (alnMod.isTraceEnabled()) {
    std::cerr << "DEBUG: Setting trace context for " << qname << std::endl;
    alnMod.setTraceContext(qname, tid);
}
```

### Step 5: Verify traceRead() is Called

**Check**: `/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp` (around line ~775)

**Expected**: `traceRead()` should be called after processing read group:
```cpp
if (alnMod.isTraceEnabled()) {
    // ... compute errLikeSum ...
    alnMod.traceRead(qname, numAlns, errLikeSum, modelUsed, modelUpdated);
    alnMod.clearTraceContext();
}
```

**Action**: Add debug output to verify this is called:
```cpp
if (alnMod.isTraceEnabled()) {
    std::cerr << "DEBUG: Calling traceRead for " << qname << std::endl;
    alnMod.traceRead(qname, numAlns, errLikeSum, modelUsed, modelUpdated);
}
```

### Step 6: Verify traceRead() Implementation

**Check**: `/mnt/pikachu/salmon/src/AlignmentModel.cpp`

**Expected**: `traceRead()` should write to trace stream:
```cpp
void AlignmentModel::traceRead(...) {
    if (traceStream_ && traceLevel_ >= ErrorModelTraceLevel::READ) {
        *traceStream_ << "READ " << qname 
                      << " numAlns=" << numAlns
                      << " errLikeSum=" << errLikeSum
                      << " modelUsed=" << (modelUsed ? 1 : 0)
                      << " modelUpdated=" << (modelUpdated ? 1 : 0)
                      << "\n";
        traceStream_->flush();
    }
}
```

**Action**: Add debug output to verify method is called:
```cpp
std::cerr << "DEBUG: traceRead called, traceStream_ = " 
          << (traceStream_ ? "valid" : "null") << std::endl;
```

### Step 7: Test with Minimal Case

Create a minimal test to isolate the issue:

```bash
# Create minimal BAM with 1 read
samtools view -h /tmp/nfcore_ec_parity_test/star_align/Aligned.toTranscriptome.out.bam | head -4 | samtools view -Sb - > /tmp/minimal_test.bam

# Run Salmon
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

# Check result
ls -lh /tmp/minimal_trace.txt
cat /tmp/minimal_trace.txt
```

**Expected**: Should see at least 1 READ line in trace file.

## Reference: Working Implementation

For comparison, our CLI implementation works correctly:

**Trace file**: `/tmp/error_model_test/our_trace.txt`
- ✅ 29,726 lines generated
- ✅ Correct format
- ✅ Shows model usage correctly

**Sample entries**:
```
READ SRR6357070.6710070 numAlns=1 errLikeSum=0 modelUsed=0 modelUpdated=1
READ SRR6357070.16033236 numAlns=1 errLikeSum=-58.8786 modelUsed=1 modelUpdated=1
```

## Success Criteria

The instrumentation is working when:

1. ✅ Trace file is created at path specified by `SALMON_ERROR_MODEL_TRACE`
2. ✅ Trace file contains READ-level entries (trace level 1+)
3. ✅ Trace file has correct format matching specification
4. ✅ Trace file has ~29,726 lines (one per read)
5. ✅ Matrix dump files are created (if `SALMON_DUMP_MATRICES` is set)

## Files to Check

1. **`/mnt/pikachu/salmon/src/AlignmentModel.cpp`**
   - Constructor: Environment variable reading
   - `traceRead()`: Implementation
   - `isTraceEnabled()`: Implementation (if exists)

2. **`/mnt/pikachu/salmon/include/AlignmentModel.hpp`**
   - Public methods: `setTraceOutput()`, `setTraceContext()`, `traceRead()`, `isTraceEnabled()`
   - Private members: `traceStream_`, `traceLevel_`, etc.

3. **`/mnt/pikachu/salmon/src/SalmonQuantifyAlignments.cpp`**
   - `isTraceEnabled()` checks
   - `setTraceContext()` calls
   - `traceRead()` calls

## Notes

- The recent change to use `alnMod.isTraceEnabled()` instead of `salmonOpts.useErrorModel` is correct - tracing should work independently of whether error model is enabled
- Make sure `isTraceEnabled()` method exists and returns true when trace stream is open
- Environment variables must be read in `AlignmentModel` constructor
- Trace stream must be opened and flushed after each write

## Expected Debug Output

When debugging, you should see output like:
```
DEBUG: SALMON_ERROR_MODEL_TRACE = /tmp/salmon_debug_trace.txt
DEBUG: Trace file opened successfully
DEBUG: Setting trace context for SRR6357070.6710070
DEBUG: Calling traceRead for SRR6357070.6710070
DEBUG: traceRead called, traceStream_ = valid
```

If you don't see these, that's where the issue is.
