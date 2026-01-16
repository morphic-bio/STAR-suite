# Salmon Error Model Instrumentation - Implementation Handoff

## Purpose

Add trace and matrix dump support to Salmon's `AlignmentModel` to enable side-by-side comparison with our EC filter CLI for error model harmonization. This instrumentation will help identify where the error models diverge during training and inference.

## Context

We're harmonizing our CIGAR-based error model implementation with Salmon's. To identify divergence points, we need matching instrumentation in both tools that outputs:
1. **Trace files**: Per-read, per-alignment, and per-transition details
2. **Matrix dumps**: Transition probability matrices at checkpoints

Our CLI already has this instrumentation. This document describes the exact changes needed in Salmon to match.

## Files to Modify

### 1. `include/AlignmentModel.hpp`

**Location**: Add includes at the top of the file (after existing includes)

```cpp
#include <fstream>
#include <string>
#include <cstdint>
#include <memory>
```

**Location**: Add public methods before the `private:` section (around line 46, after `void normalize();`)

```cpp
  void normalize();
  
  // Trace support (for error model harmonization)
  void setTraceOutput(const std::string& traceFile, int traceLevel);
  void setTraceContext(const std::string& qname, uint32_t tid);
  void clearTraceContext();
  void setMatrixDumpPrefix(const std::string& prefix);
  bool dumpMatrices(const std::string& prefix, const std::string& suffix) const;
```

**Location**: Add private members before the closing brace (around line 83, before `};`)

```cpp
  std::mutex outputMutex_;
  
  // Trace support (private helpers)
  void traceRead(const std::string& qname, size_t numAlns, double errLikeSum, bool modelUsed, bool modelUpdated);
  void traceAlignment(const std::string& qname, uint32_t tid, double fg, double bg, double errLike, uint32_t bin);
  void traceTransition(const std::string& qname, uint32_t tid, uint32_t readPos, uint32_t bin, uint32_t prevState, uint32_t curState, double logProb);
  
  std::unique_ptr<std::ofstream> traceStream_;
  int traceLevel_;  // 0=none, 1=read, 2=alignment, 3=transition
  std::string currentQname_;
  uint32_t currentTid_;
  std::string matrixDumpPrefix_;
  
  void dumpMatrixToFile(const AtomicMatrix<double>& matrix, const std::string& filename, bool isLeft, uint32_t bin, size_t numObserved) const;
```

### 2. `src/AlignmentModel.cpp`

**Location**: Add includes at the top (after existing includes, around line 6)

```cpp
#include <cstdlib>
#include <iomanip>
```

**Location**: Modify constructor (around line 15)

**Find**:
```cpp
AlignmentModel::AlignmentModel(double alpha, uint32_t readBins)
    : transitionProbsLeft_(readBins), transitionProbsRight_(readBins),
      isEnabled_(true), readBins_(readBins) {

  for (size_t i = 0; i < readBins; ++i) {
    transitionProbsLeft_[i] = std::move(AtomicMatrix<double>(
        numAlignmentStates(), numAlignmentStates(), alpha));
    transitionProbsRight_[i] = std::move(AtomicMatrix<double>(
        numAlignmentStates(), numAlignmentStates(), alpha));
  }
}
```

**Replace with**:
```cpp
AlignmentModel::AlignmentModel(double alpha, uint32_t readBins)
    : transitionProbsLeft_(readBins), transitionProbsRight_(readBins),
      isEnabled_(true), readBins_(readBins),
      traceStream_(nullptr), traceLevel_(0), currentQname_(""), currentTid_(0), matrixDumpPrefix_("") {

  for (size_t i = 0; i < readBins; ++i) {
    transitionProbsLeft_[i] = std::move(AtomicMatrix<double>(
        numAlignmentStates(), numAlignmentStates(), alpha));
    transitionProbsRight_[i] = std::move(AtomicMatrix<double>(
        numAlignmentStates(), numAlignmentStates(), alpha));
  }
  
  // Check for trace environment variable
  const char* traceFile = std::getenv("SALMON_ERROR_MODEL_TRACE");
  const char* traceLevelStr = std::getenv("SALMON_TRACE_LEVEL");
  if (traceFile && traceFile[0] != '\0') {
    int level = 1;  // Default level 1
    if (traceLevelStr) {
      level = std::atoi(traceLevelStr);
      if (level < 0 || level > 3) level = 1;
    }
    setTraceOutput(std::string(traceFile), level);
  }
  
  // Check for matrix dump environment variable
  const char* matrixDump = std::getenv("SALMON_DUMP_MATRICES");
  if (matrixDump && matrixDump[0] != '\0') {
    setMatrixDumpPrefix(std::string(matrixDump));
  }
}
```

**Location**: Add trace call in `logLikelihood()` method (around line 201-205)

**Find**:
```cpp
      curStateIdx = curRefBase * numStates + curReadBase;
      double tp = transitionProbs[readPosBin](prevStateIdx, curStateIdx);
      logLike += tp;
      bgLogLike += transitionProbs[readPosBin](0, 0);
      prevStateIdx = curStateIdx;
```

**Replace with**:
```cpp
      curStateIdx = curRefBase * numStates + curReadBase;
      double tp = transitionProbs[readPosBin](prevStateIdx, curStateIdx);
      logLike += tp;
      bgLogLike += transitionProbs[readPosBin](0, 0);
      
      // Trace transition if enabled
      if (traceStream_ && traceLevel_ >= 3) {
        traceTransition(currentQname_, currentTid_, static_cast<uint32_t>(readIdx), readPosBin, prevStateIdx, curStateIdx, tp);
      }
      
      prevStateIdx = curStateIdx;
```

**Location**: Add trace call before return in `logLikelihood()` (around line 227)

**Find**:
```cpp
  return {logLike, bgLogLike};
}

double AlignmentModel::logLikelihood(const ReadPair& hit, const ReadPair& primary, Transcript& ref) {
```

**Replace with**:
```cpp
  // Trace alignment if enabled
  if (traceStream_ && traceLevel_ >= 2) {
    uint32_t bin = 0;
    if (readLen > 0) {
      double invLen = static_cast<double>(readBins_) / readLen;
      bin = static_cast<uint32_t>(readLen / 2 * invLen);
      if (bin >= readBins_) bin = readBins_ - 1;
    }
    traceAlignment(currentQname_, currentTid_, logLike, bgLogLike, logLike - bgLogLike, bin);
  }
  
  return {logLike, bgLogLike};
}

double AlignmentModel::logLikelihood(const ReadPair& hit, const ReadPair& primary, Transcript& ref) {
```

**Location**: Add implementation at end of file (before closing brace, around line 533)

**Add**:
```cpp
//=============================================================================
// Trace and dump support
//=============================================================================

void AlignmentModel::setTraceOutput(const std::string& traceFile, int traceLevel) {
  traceStream_ = std::make_unique<std::ofstream>(traceFile);
  if (!traceStream_->is_open()) {
    traceStream_.reset();
    traceLevel_ = 0;
    return;
  }
  traceLevel_ = traceLevel;
}

void AlignmentModel::setTraceContext(const std::string& qname, uint32_t tid) {
  currentQname_ = qname;
  currentTid_ = tid;
}

void AlignmentModel::clearTraceContext() {
  currentQname_.clear();
  currentTid_ = 0;
}

void AlignmentModel::traceRead(const std::string& qname, size_t numAlns, double errLikeSum, bool modelUsed, bool modelUpdated) {
  if (traceStream_ && traceLevel_ >= 1) {
    *traceStream_ << "READ " << qname 
                  << " numAlns=" << numAlns
                  << " errLikeSum=" << errLikeSum
                  << " modelUsed=" << (modelUsed ? 1 : 0)
                  << " modelUpdated=" << (modelUpdated ? 1 : 0)
                  << "\n";
    traceStream_->flush();
  }
}

void AlignmentModel::traceAlignment(const std::string& qname, uint32_t tid, double fg, double bg, double errLike, uint32_t bin) {
  if (traceStream_ && traceLevel_ >= 2) {
    *traceStream_ << "ALN " << qname
                  << " tid=" << tid
                  << " fg=" << fg
                  << " bg=" << bg
                  << " errLike=" << errLike
                  << " bin=" << bin
                  << "\n";
    traceStream_->flush();
  }
}

void AlignmentModel::traceTransition(const std::string& qname, uint32_t tid, uint32_t readPos, uint32_t bin, uint32_t prevState, uint32_t curState, double logProb) {
  if (traceStream_ && traceLevel_ >= 3) {
    *traceStream_ << "TRANS " << qname
                  << " tid=" << tid
                  << " readPos=" << readPos
                  << " bin=" << bin
                  << " prevState=" << prevState
                  << " curState=" << curState
                  << " logProb=" << logProb
                  << "\n";
    traceStream_->flush();
  }
}

void AlignmentModel::setMatrixDumpPrefix(const std::string& prefix) {
  matrixDumpPrefix_ = prefix;
}

bool AlignmentModel::dumpMatrices(const std::string& prefix, const std::string& suffix) const {
  if (prefix.empty()) {
    return false;
  }
  
  // Note: numObserved is not tracked in AlignmentModel
  // This would need to be passed from SalmonQuantifyAlignments.cpp or tracked here
  // For now, use 0 as placeholder - this is a known limitation
  size_t numObserved = 0;
  
  // Dump left matrices
  for (uint32_t bin = 0; bin < readBins_; ++bin) {
    std::string filename = prefix + "_left_bin" + std::to_string(bin) + suffix + ".tsv";
    dumpMatrixToFile(transitionProbsLeft_[bin], filename, true, bin, numObserved);
  }
  
  // Dump right matrices
  for (uint32_t bin = 0; bin < readBins_; ++bin) {
    std::string filename = prefix + "_right_bin" + std::to_string(bin) + suffix + ".tsv";
    dumpMatrixToFile(transitionProbsRight_[bin], filename, false, bin, numObserved);
  }
  
  return true;
}

void AlignmentModel::dumpMatrixToFile(const AtomicMatrix<double>& matrix, const std::string& filename, bool isLeft, uint32_t bin, size_t numObserved) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    return;
  }
  
  out << "# Transition matrix: " << (isLeft ? "left" : "right") 
      << ", bin=" << bin 
      << ", numObserved=" << numObserved 
      << "\n";
  out << "# Rows: prevState (0-81), Cols: curState (0-81)\n";
  
  for (uint32_t row = 0; row < numAlignmentStates(); ++row) {
    for (uint32_t col = 0; col < numAlignmentStates(); ++col) {
      if (col > 0) out << "\t";
      double logProb = matrix(row, col);
      out << std::setprecision(17) << logProb;
    }
    out << "\n";
  }
  
  out.close();
}
```

### 3. Integration in `src/SalmonQuantifyAlignments.cpp` (Optional but Recommended)

To enable full tracing, you'll need to call `setTraceContext()` before computing errLike for each read group. This requires access to the `AlignmentModel` instance and read names.

**Location**: In the main quantification loop where errLike is computed (around line 545-552)

**Find the section where errLike is computed**:
```cpp
double errLike = salmon::math::LOG_1;
if (useASWithoutCIGAR) {
  errLike = -salmonOpts.scoreExp * (bestAS - alnScore);
} else if (useAuxParams and salmonOpts.useErrorModel) {
  errLike = alnMod.logLikelihood(*aln, *alnPrimary, transcript);
  ++sidx;
}
```

**Before this section**, add:
```cpp
// Set trace context for this alignment
if (salmonOpts.useErrorModel && alnMod.traceStream_) {
  const char* qname = bam_name(aln->read1 ? aln->read1 : aln->read2);
  uint32_t tid = aln->transcriptID();
  alnMod.setTraceContext(std::string(qname), tid);
}
```

**After processing all alignments in a read group**, add:
```cpp
// Trace read summary if enabled
if (salmonOpts.useErrorModel && alnMod.traceStream_ && alnMod.traceLevel_ >= 1) {
  // Count alignments and sum errLike
  size_t numAlns = alnGroup->alignments().size();
  double errLikeSum = 0.0;
  bool modelUsed = (processedReads >= salmonOpts.numPreBurninFrags);
  bool modelUpdated = (processedReads < salmonOpts.numBurninFrags);
  
  for (const auto& aln : alnGroup->alignments()) {
    // errLike was computed above, would need to store it or recompute
    // For now, this is a placeholder - full implementation would track errLike per alignment
  }
  
  const char* qname = bam_name(alnGroup->alignments().front()->read1 ? 
                                alnGroup->alignments().front()->read1 : 
                                alnGroup->alignments().front()->read2);
  alnMod.traceRead(std::string(qname), numAlns, errLikeSum, modelUsed, modelUpdated);
  alnMod.clearTraceContext();
}
```

**Note**: The full integration requires tracking errLike values per alignment, which may require refactoring. The trace context setting is the minimum needed for ALN and TRANS level tracing to work.

## Trace Format Specification

The trace format must match exactly for comparison scripts to work:

### Level 1 (READ) - Per-read summary
```
READ <qname> numAlns=<N> errLikeSum=<X> modelUsed=<0/1> modelUpdated=<0/1>
```

### Level 2 (ALN) - Per-alignment details
```
ALN <qname> tid=<X> fg=<Y> bg=<Z> errLike=<W> bin=<B>
```

### Level 3 (TRANS) - Per-transition details
```
TRANS <qname> tid=<X> readPos=<P> bin=<B> prevState=<S1> curState=<S2> logProb=<X>
```

**Important**: 
- All values are space-separated
- Numeric values use full precision (no scientific notation unless necessary)
- Each line ends with `\n` and is flushed immediately
- `qname` is the read name from BAM
- `tid` is the transcript ID
- `bin` is the read position bin (0 to readBins-1)

## Matrix Dump Format

Matrix files are TSV format:
- Header lines start with `#`
- First header: `# Transition matrix: <left/right>, bin=<N>, numObserved=<M>`
- Second header: `# Rows: prevState (0-81), Cols: curState (0-81)`
- Data: 82 rows × 82 columns, tab-separated
- Values: Log probabilities with 17 decimal precision

**File naming**: `<prefix>_<left/right>_bin<N>_<suffix>.tsv`

Example: `/tmp/salmon_matrices_left_bin0_pre5000.tsv`

## Environment Variables

- `SALMON_ERROR_MODEL_TRACE`: Path to trace output file (enables tracing)
- `SALMON_TRACE_LEVEL`: Trace level 1, 2, or 3 (default: 1 if trace file specified)
- `SALMON_DUMP_MATRICES`: Prefix for matrix dump files (enables matrix dumping)

## Usage Example

```bash
SALMON_ERROR_MODEL_TRACE=/tmp/salmon_trace.txt \
SALMON_TRACE_LEVEL=2 \
SALMON_DUMP_MATRICES=/tmp/salmon_matrices \
salmon quant -t transcriptome.fa -l A -a aligned.bam \
    --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
    --noFragLengthDist -p 1 -o salmon_out
```

This will generate:
- `/tmp/salmon_trace.txt` - Trace file with ALN-level details
- `/tmp/salmon_matrices_*_bin*.tsv` - Matrix dump files (if checkpoints are added)

## Testing

After implementation, verify:

1. **Compilation**: Code compiles without errors
2. **Trace file generation**: Run with `SALMON_ERROR_MODEL_TRACE` set and verify file is created
3. **Trace format**: Check that trace lines match the format specification above
4. **Matrix dump**: Run with `SALMON_DUMP_MATRICES` set and verify TSV files are created
5. **Integration**: Verify trace context is set correctly (qname and tid appear in trace output)

## Known Limitations

1. **numObserved tracking**: Matrix dumps show `numObserved=0` because `AlignmentModel` doesn't track this. This can be addressed by:
   - Adding a `numObserved_` member to `AlignmentModel` and incrementing it
   - Or passing `numObserved` as a parameter to `dumpMatrices()`

2. **Full READ-level tracing**: Requires tracking errLike per alignment in the quantification loop, which may need refactoring. ALN and TRANS level tracing work with just `setTraceContext()`.

3. **Checkpoint dumps**: Matrix dumps at checkpoints (1000, 2000, etc. reads) require calling `dumpMatrices()` from `SalmonQuantifyAlignments.cpp` when `processedReads` reaches those values.

## Questions or Issues

If you encounter issues:
- Check that all includes are present (`<fstream>`, `<string>`, `<iomanip>`, `<memory>`)
- Verify `std::unique_ptr` is available (C++11)
- Ensure `AtomicMatrix` has `operator()` for reading values
- Check that `numAlignmentStates()` is accessible (it's a static method)

## Success Criteria

Implementation is complete when:
- ✅ Code compiles without errors
- ✅ Trace file is generated with correct format when env vars are set
- ✅ Matrix dump files are generated with correct format
- ✅ Trace context (qname, tid) appears correctly in trace output
- ✅ Trace format matches specification exactly (for comparison scripts)
