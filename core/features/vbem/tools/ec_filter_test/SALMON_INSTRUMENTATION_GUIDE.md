# Salmon Error Model Instrumentation Guide

This document describes the changes needed to add matching trace and matrix dump support to Salmon's AlignmentModel for error model harmonization.

## Files to Modify

### 1. `/mnt/pikachu/salmon/include/AlignmentModel.hpp`

Add includes at top:
```cpp
#include <fstream>
#include <string>
#include <cstdint>
#include <memory>
```

Add public methods before `private:`:
```cpp
  // Trace support (for error model harmonization)
  void setTraceOutput(const std::string& traceFile, int traceLevel);
  void setTraceContext(const std::string& qname, uint32_t tid);
  void clearTraceContext();
  void setMatrixDumpPrefix(const std::string& prefix);
  bool dumpMatrices(const std::string& prefix, const std::string& suffix) const;
```

Add private members before closing brace:
```cpp
  // Trace support
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

### 2. `/mnt/pikachu/salmon/src/AlignmentModel.cpp`

Add includes:
```cpp
#include <cstdlib>
#include <iomanip>
```

Modify constructor to initialize trace members and check env vars:
```cpp
AlignmentModel::AlignmentModel(double alpha, uint32_t readBins)
    : transitionProbsLeft_(readBins), transitionProbsRight_(readBins),
      isEnabled_(true), readBins_(readBins),
      traceStream_(nullptr), traceLevel_(0), currentQname_(""), currentTid_(0), matrixDumpPrefix_("") {
  // ... existing initialization ...
  
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

Add trace call in `logLikelihood()` method around line 205:
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

Add trace call before return in `logLikelihood()` around line 227:
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
```

Add implementation at end of file (before closing brace):
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
  // This would need to be passed from SalmonQuantifyAlignments.cpp
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

### 3. Integration in SalmonQuantifyAlignments.cpp

To enable tracing, you'll need to call `setTraceContext()` before computing errLike and `traceRead()` after processing each read group. This requires access to the AlignmentModel instance and read names, which may need to be passed through the quantification pipeline.

## Usage

After implementing, use environment variables:

```bash
SALMON_ERROR_MODEL_TRACE=/tmp/salmon_em_trace.txt \
SALMON_TRACE_LEVEL=2 \
SALMON_DUMP_MATRICES=/tmp/salmon_matrices \
salmon quant -t transcriptome.fa -l A -a aligned.bam \
    --dumpEqWeights --noLengthCorrection --noEffectiveLengthCorrection \
    --noFragLengthDist -p 1 -o salmon_out
```

## Notes

- Matrix dumps will show `numObserved=0` until `numObserved` tracking is added to AlignmentModel or passed as a parameter
- Trace context (qname, tid) needs to be set from SalmonQuantifyAlignments.cpp where read groups are processed
- The trace format matches our CLI exactly for direct comparison
