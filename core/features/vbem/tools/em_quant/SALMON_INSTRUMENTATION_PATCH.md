# Salmon Instrumentation Patch Guide

This guide explains how to add debug tracing to Salmon to compare with em_quant traces and find divergence points.

## Overview

Add instrumentation to Salmon's `VBEMUpdate_` function to log the same values as em_quant:
- Per-iteration: alpha, logNorm, expTheta, expected_count
- Per-EC: ec_id, denom, expTheta, aux, contribution

## Files to Modify

### 1. Add Debug Parameters to SalmonOpts

**File**: `/mnt/pikachu/salmon/include/salmon/SalmonOpts.hpp` (or similar)

Add debug tracing options:
```cpp
// Debug tracing
bool debugTrace = false;
std::string debugTraceFile;
std::vector<std::string> debugTranscripts;
```

### 2. Parse Command-Line Arguments

**File**: `/mnt/pikachu/salmon/src/ProgramOptionsGenerator.cpp` (or similar)

Add parsing for debug flags:
```cpp
po::options_description debugOptions("Debug Options");
debugOptions.add_options()
    ("debugTrace", po::value<std::string>(), "Enable debug trace output to file")
    ("debugTranscripts", po::value<std::string>(), "Comma-separated transcript IDs to trace");
```

### 3. Modify VBEMUpdate_ Function (Single-Threaded)

**File**: `/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp`  
**Location**: Line 104 (single-threaded version)

Add debug parameters and logging:

```cpp
template <typename VecT>
void VBEMUpdate_(std::vector<std::vector<uint32_t>>& txpGroupLabels,
                 std::vector<std::vector<double>>& txpGroupCombinedWeights,
                 const std::vector<uint64_t>& txpGroupCounts,
                 std::vector<double>& priorAlphas, 
                 const VecT& alphaIn, VecT& alphaOut, VecT& expTheta,
                 // ADD: Debug parameters
                 const std::vector<Transcript>& transcripts,
                 uint32_t iterNum,
                 std::ofstream* debugStream,
                 const std::unordered_set<std::string>& debugTranscripts) {

  assert(alphaIn.size() == alphaOut.size());
  size_t M = alphaIn.size();
  size_t numEQClasses = txpGroupLabels.size();
  double alphaSum = {0.0};
  for (size_t i = 0; i < M; ++i) {
    alphaSum += alphaIn[i] + priorAlphas[i];
  }

  double logNorm = boost::math::digamma(alphaSum);

  // Precompute expTheta
  for (size_t i = 0; i < M; ++i) {
    auto ap = alphaIn[i] + priorAlphas[i];
    if (ap > ::digammaMin) {
        expTheta[i] = std::exp(boost::math::digamma(ap) - logNorm);
    } else {
      expTheta[i] = 0.0;
    }
    alphaOut[i] = 0.0;
  }

  // ADD: Log per-iteration values for debug transcripts (BEFORE E-step)
  if (debugStream && debugStream->is_open()) {
    for (size_t i = 0; i < M; ++i) {
      const std::string& txpName = transcripts[i].RefName;
      if (debugTranscripts.find(txpName) != debugTranscripts.end()) {
        *debugStream << iterNum << "\t" << txpName << "\t"
                     << (alphaIn[i] + priorAlphas[i]) << "\t"
                     << logNorm << "\t"
                     << expTheta[i] << "\t0.0\n";
        debugStream->flush();
      }
    }
  }

  // E-step loop
  for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
    uint64_t count = txpGroupCounts[eqID];
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    size_t groupSize = txpGroupCombinedWeights[eqID].size();
    
    if (BOOST_LIKELY(groupSize > 1)) {
      double denom = 0.0;
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        if (expTheta[tid] > 0.0) {
          double v = expTheta[tid] * aux;
          denom += v;
        }
      }
      if (denom <= ::minEQClassWeight) {
        // Skip
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          if (expTheta[tid] > 0.0) {
            double v = expTheta[tid] * aux;
            double contribution = v * invDenom;
            alphaOut[tid] += contribution;
            
            // ADD: Log EC-level details for debug transcripts
            if (debugStream && debugStream->is_open()) {
              const std::string& txpName = transcripts[tid].RefName;
              if (debugTranscripts.find(txpName) != debugTranscripts.end()) {
                *debugStream << "EC\t" << iterNum << "\t" << eqID << "\t"
                            << txpName << "\t" << denom << "\t"
                            << expTheta[tid] << "\t" << aux << "\t"
                            << contribution << "\n";
                debugStream->flush();
              }
            }
          }
        }
      }
    } else {
      // Single-transcript EC
      auto tid = txps.front();
      alphaOut[tid] += count;
      
      // ADD: Log single-transcript EC
      if (debugStream && debugStream->is_open()) {
        const std::string& txpName = transcripts[tid].RefName;
        if (debugTranscripts.find(txpName) != debugTranscripts.end()) {
          *debugStream << "EC\t" << iterNum << "\t" << eqID << "\t"
                      << txpName << "\tSINGLE\t1\t1\t" << count << "\n";
          debugStream->flush();
        }
      }
    }
  }
  
  // ADD: Log per-iteration values AFTER E-step (with expected_count)
  if (debugStream && debugStream->is_open()) {
    for (size_t i = 0; i < M; ++i) {
      const std::string& txpName = transcripts[i].RefName;
      if (debugTranscripts.find(txpName) != debugTranscripts.end()) {
        // Recompute logNorm after M-step (will be done in caller)
        double newAlphaSum = 0.0;
        for (size_t j = 0; j < M; ++j) {
          newAlphaSum += alphaOut[j] + priorAlphas[j];
        }
        double newLogNorm = boost::math::digamma(newAlphaSum);
        double newExpTheta = (alphaOut[i] + priorAlphas[i] > ::digammaMin) 
            ? std::exp(boost::math::digamma(alphaOut[i] + priorAlphas[i]) - newLogNorm)
            : 0.0;
        
        *debugStream << iterNum << "\t" << txpName << "\t"
                     << (alphaOut[i] + priorAlphas[i]) << "\t"
                     << newLogNorm << "\t"
                     << newExpTheta << "\t"
                     << alphaOut[i] << "\n";
        debugStream->flush();
      }
    }
  }
}
```

### 4. Modify Main Optimize Function

**File**: `/mnt/pikachu/salmon/src/CollapsedEMOptimizer.cpp`  
**Location**: Line 738 (`optimize` function)

Add debug setup and pass to VBEMUpdate_:

```cpp
template <typename ExpT>
bool CollapsedEMOptimizer::optimize(ExpT& readExp, SalmonOpts& sopt,
                                    double relDiffTolerance, uint32_t maxIter) {
  // ... existing code ...

  // ADD: Setup debug tracing
  std::ofstream debugStream;
  std::unordered_set<std::string> debugTranscripts;
  if (sopt.debugTrace && !sopt.debugTraceFile.empty()) {
    debugStream.open(sopt.debugTraceFile);
    if (debugStream.is_open()) {
      // Write header
      debugStream << "iter\ttranscript\talpha\tlogNorm\texpTheta\texpected_count\n";
      debugStream << "# EC-level trace: EC\titer\tec_id\ttranscript\tdenom\texpTheta\taux\tcontribution\n";
      
      // Parse debug transcript IDs
      std::istringstream iss(sopt.debugTranscripts);
      std::string txp;
      while (std::getline(iss, txp, ',')) {
        debugTranscripts.insert(txp);
      }
      
      // Log initial state (iter -1)
      double alphaSum = 0.0;
      for (size_t i = 0; i < transcripts.size(); ++i) {
        alphaSum += alphas[i].load() + priorAlphas[i];
      }
      double logNorm = boost::math::digamma(alphaSum);
      
      for (size_t i = 0; i < transcripts.size(); ++i) {
        const std::string& txpName = transcripts[i].RefName;
        if (debugTranscripts.find(txpName) != debugTranscripts.end()) {
          double expTheta = (alphas[i].load() + priorAlphas[i] > ::digammaMin)
              ? std::exp(boost::math::digamma(alphas[i].load() + priorAlphas[i]) - logNorm)
              : 0.0;
          debugStream << "-1\t" << txpName << "\t"
                      << (alphas[i].load() + priorAlphas[i]) << "\t"
                      << logNorm << "\t" << expTheta << "\t0.0\n";
        }
      }
      debugStream.flush();
    }
  }

  // ... existing code ...

  uint32_t itNum = 0;
  while (itNum < minIter or (itNum < maxIter and !converged)) {
    if (useVBEM) {
      // MODIFY: Pass debug parameters
      VBEMUpdate_(arena, eqVec, priorAlphas, alphas, alphasPrime, expTheta,
                  transcripts, itNum, &debugStream, debugTranscripts);
    } else {
      EMUpdate_(arena, eqVec, priorAlphas, alphas, alphasPrime);
    }

    // ... existing convergence check ...

    ++itNum;
  }

  // ADD: Close debug stream
  if (debugStream.is_open()) {
    debugStream.close();
  }

  // ... rest of function ...
}
```

### 5. Update Function Signatures

Update all `VBEMUpdate_` call sites to include the new parameters (or use default arguments).

## Testing

After adding instrumentation:

```bash
salmon quant \
    -i /tmp/salmon_chr22_idx \
    -l A \
    -1 /tmp/sub_R1.fq -2 /tmp/sub_R2.fq \
    --dumpEqWeights \
    --threads 1 \
    --noLengthCorrection \
    --noFragLengthDist \
    --noEffectiveLengthCorrection \
    --debugTrace /tmp/salmon_trace.txt \
    --debugTranscripts ENST00000465752,ENST00000464034 \
    -o /tmp/salmon_out
```

## Comparison

Once both traces are generated:

```bash
cd tools/em_quant
python3 compare_traces.py /tmp/em_trace.txt /tmp/salmon_trace.txt \
    --transcript ENST00000465752
```

This will identify the first iteration/EC where divergence occurs.

## Notes

1. **Thread Safety**: The single-threaded version is easier to instrument. For multi-threaded version, use `#pragma omp critical` or similar for debug output.

2. **Transcript Names**: Use `transcripts[i].RefName` or similar to get transcript names from IDs.

3. **Iteration Number**: Pass iteration number from the main loop to VBEMUpdate_.

4. **Expected Counts**: Log `alphaOut[i]` as expected_count (before adding prior).

5. **Format Matching**: Ensure output format exactly matches em_quant format for easy comparison.
