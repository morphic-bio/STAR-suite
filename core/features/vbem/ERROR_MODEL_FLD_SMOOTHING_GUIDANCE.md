# Error Model Wiring & FLD Smoothing Implementation Guidance

## Status

**Error Model**: Infrastructure added (AlignmentModel and Transcriptome members in TranscriptQuantEC), but computation not yet wired.

**FLD Smoothing**: Not implemented. FLD uses raw counts with no smoothing, uniform-weight updates, and no stochastic acceptance.

## Error Model Wiring

### Current State
- `TranscriptQuantEC` has `alignment_model_` and `transcriptome_` members (optional, can be nullptr)
- `RawAlignment::err_like` is always `0.0` and `has_err_like` is `false`
- Error model is effectively disabled

### Required Steps

#### 1. Load Transcriptome Sequences
**Location**: `source/STAR.cpp` or `source/ReadAlign.cpp`

```cpp
// In STAR.cpp, after transcriptomeMain is loaded:
std::unique_ptr<libem::Transcriptome> libem_transcriptome;
if (P.quant.transcriptVB.errorModel != "off") {
    // Find transcriptome FASTA file
    std::string fasta_path = P.genomeDir + "/transcriptome.fa";  // Adjust path as needed
    libem_transcriptome = std::make_unique<libem::Transcriptome>();
    if (!libem_transcriptome->loadFromFasta(fasta_path)) {
        // Fallback: error or disable error model
        libem_transcriptome.reset();
    }
}
```

#### 2. Create AlignmentModel Instance
**Location**: `source/STAR.cpp` or `source/ReadAlign.cpp`

```cpp
std::unique_ptr<AlignmentModel> alignment_model;
if (P.quant.transcriptVB.errorModel != "off" && libem_transcriptome) {
    alignment_model = std::make_unique<AlignmentModel>(0.001, 4);  // alpha=0.001, readBins=4
}
```

#### 3. Pass to TranscriptQuantEC Constructor
**Location**: `source/ReadAlign.cpp` (line 24)

```cpp
quantEC = new TranscriptQuantEC(chunkTr->nTr, iChunk, 
                                 P.quant.transcriptVB.traceFile,
                                 P.quant.transcriptVB.traceLimit, P,
                                 alignment_model.get(),  // Pass AlignmentModel
                                 libem_transcriptome.get());  // Pass Transcriptome
```

**Note**: `alignment_model` and `libem_transcriptome` need to be stored in `ReadAlign` or passed through `Parameters` to be accessible.

#### 4. Modify addReadAlignments to Accept Read Sequences
**Location**: `source/TranscriptQuantEC.h` and `source/TranscriptQuantEC.cpp`

**Option A**: Add read sequences to signature:
```cpp
void addReadAlignments(Transcript* trMult, uint nTr, 
                       const std::vector<double>& auxProbs, 
                       const char* qname,
                       char** Read0 = nullptr,  // Add read sequences
                       uint* readLengthOriginal = nullptr);  // Add read lengths
```

**Option B**: Store read sequences in `ReadAlign` and access via member (preferred if sequences are already available).

#### 5. Build BAM-Style CIGAR from Transcript
**Location**: `source/TranscriptQuantEC.cpp::transcriptToRawAlignment()`

```cpp
// Convert Transcript::cigar to BAM-packed format
std::vector<uint32_t> packed_cigar;
if (!tr->cigar.empty()) {
    for (const auto& op : tr->cigar) {
        uint32_t opLen = op[1];  // Length
        uint8_t opCode = op[0];   // Operation (BAM_CIGAR_M, etc.)
        packed_cigar.push_back(opLen << BAM_CIGAR_SHIFT | opCode);
    }
}
```

**Alternative**: Use `Transcript::generateCigarP()` and parse the string (less efficient).

#### 6. Pack Read Sequences
**Location**: `source/TranscriptQuantEC.cpp::transcriptToRawAlignment()` or `addReadAlignments()`

```cpp
// Pack read sequences using nuclPackBAM (from SequenceFuns.h)
uint32_t seqLen = readLengthOriginal[mate];
uint32_t packedLen = (seqLen + 1) / 2;
std::vector<uint8_t> packedSeq(packedLen);
nuclPackBAM(Read0[mate], reinterpret_cast<char*>(packedSeq.data()), seqLen);
```

#### 7. Compute err_like
**Location**: `source/TranscriptQuantEC.cpp::transcriptToRawAlignment()` or before `computeAuxProbs()`

```cpp
if (alignment_model_ && transcriptome_) {
    const TranscriptSequence* ref = transcriptome_->getTranscript(transcript_id);
    if (ref && !packed_cigar.empty()) {
        if (tr->readNmates == 2) {
            // Paired-end: use logLikelihoodPairedAndUpdate
            // Need CIGAR and sequences for both mates
            double err_like = alignment_model_->logLikelihoodPairedAndUpdate(
                packed_cigar1.data(), packed_cigar1.size(),
                packedSeq1.data(), readLen1,
                packed_cigar2.data(), packed_cigar2.size(),
                packedSeq2.data(), readLen2,
                ref, refPos1, refPos2,
                0.0,  // logProb (for weighting updates)
                0.0   // logMass (for weighting updates)
            );
            aln.err_like = err_like;
            aln.has_err_like = true;
        } else {
            // Single-end: use logLikelihoodAndUpdate
            double err_like = alignment_model_->logLikelihoodAndUpdate(
                packed_cigar.data(), packed_cigar.size(),
                packedSeq.data(), readLen,
                ref, refPos,
                true,  // isLeft
                0.0,   // logProb
                0.0    // logMass
            );
            aln.err_like = err_like;
            aln.has_err_like = true;
        }
    }
}
```

## FLD Smoothing

### Current State
- FLD uses raw counts (`FLDAccumulator`) with no smoothing
- Updates use uniform weights (`1/nAlignT`) for all alignments
- No stochastic acceptance or forgetting mass

### Required Steps

#### 1. Port FragmentLengthDistribution Class
**Location**: `source/libem/fragment_length_distribution.h/cpp` (new files)

Port Salmon's `FragmentLengthDistribution` class with:
- Gaussian prior (mean=250, sd=25)
- Binomial kernel smoothing
- `maxFragLen=1000` (already fixed in `FLDAccumulator`)
- Methods: `addObservation()`, `getPMF()`, `getLogCMF()`

#### 2. Integrate with FLDAccumulator
**Option A**: Replace `FLDAccumulator` with `FragmentLengthDistribution`
**Option B**: Use `FragmentLengthDistribution` to smooth `FLDAccumulator` output

**Location**: `source/TranscriptQuantEC.cpp::getParams()`

```cpp
// After FLD observations are collected, smooth them
if (observedFLD_.isValid() && fld_smoother_) {
    std::vector<double> smoothed_pmf = fld_smoother_->getPMF();
    // Use smoothed PMF instead of raw PMF
}
```

#### 3. Tie to ForgettingMassCalculator
**Location**: `source/libem/forgetting_mass.h` (already exists)

```cpp
#include "forgetting_mass.h"

// In TranscriptQuantEC:
ForgettingMassCalculator forgetting_mass_;
```

#### 4. Add Stochastic Acceptance
**Location**: `source/TranscriptQuantEC.cpp::addFLDObservation()`

```cpp
void TranscriptQuantEC::addFLDObservation(int32_t frag_len, double weight, double log_prob) {
    // Only update during burn-in
    if (num_processed_fragments_ >= num_burnin_frags_) {
        return;  // Post-burnin: no updates
    }
    
    // Stochastic acceptance: r < exp(logProb)
    if (log_prob != 0.0) {  // If logProb is available
        double r = std::log(static_cast<double>(rand()) / RAND_MAX);
        if (r >= log_prob) {
            return;  // Reject update
        }
    }
    
    // Apply forgetting mass
    double forgetting_mass = forgetting_mass_.getMass(num_processed_fragments_);
    double adjusted_weight = weight * forgetting_mass;
    
    observedFLD_.add(frag_len, adjusted_weight);
    fld_dirty_ = true;
}
```

#### 5. Cache log-CMF
**Location**: `source/libem/fragment_length_distribution.cpp`

```cpp
class FragmentLengthDistribution {
    mutable std::vector<double> cached_log_cmf_;
    mutable bool log_cmf_dirty_ = true;
    
    void updateLogCMF() const {
        if (log_cmf_dirty_) {
            cached_log_cmf_.resize(pmf_.size());
            double cumsum = 0.0;
            for (size_t i = 0; i < pmf_.size(); ++i) {
                cumsum += pmf_[i];
                cached_log_cmf_[i] = std::log(cumsum);
            }
            log_cmf_dirty_ = false;
        }
    }
    
    double getLogCMF(int32_t len) const {
        updateLogCMF();
        if (len < 0 || len >= static_cast<int32_t>(cached_log_cmf_.size())) {
            return LOG_0;
        }
        return cached_log_cmf_[len];
    }
};
```

## Testing

After implementation:
1. Run 20K trace diff to verify error model usage
2. Run 100K EC parity to verify FLD smoothing impact
3. Run EM TPM correlations to verify overall parity improvements
4. Verify error model is active (check `has_err_like=true` in traces)
5. Verify FLD smoothing is applied (check PMF smoothness)

## Notes

- Error model requires transcript sequences (FASTA file), which may not always be available
- FLD smoothing requires additional computational overhead but improves accuracy
- Stochastic acceptance requires `logProb` to be available, which may require refactoring `addReadAlignments()`
