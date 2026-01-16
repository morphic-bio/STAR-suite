#include "alignment_model.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <cmath>

namespace libem {

// Log constants
static constexpr double LOG_0 = -std::numeric_limits<double>::infinity();
static constexpr double LOG_1 = 0.0;
static constexpr double LOG_EPSILON = std::log(0.375e-10);

// Log-sum-exp helper
static inline double logAdd(double a, double b) {
    if (a == LOG_0) return b;
    if (b == LOG_0) return a;
    if (a > b) {
        return a + std::log(1.0 + std::exp(b - a));
    } else {
        return b + std::log(1.0 + std::exp(a - b));
    }
}

//=============================================================================
// TransitionMatrix implementation
//=============================================================================

TransitionMatrix::TransitionMatrix() : alpha_(0.001) {
    storage_.resize(NUM_ALIGNMENT_STATES * NUM_ALIGNMENT_STATES);
    rowsums_.resize(NUM_ALIGNMENT_STATES);
    
    double logAlpha = std::log(alpha_);
    double logRowSum = std::log(NUM_ALIGNMENT_STATES * alpha_);
    
    std::fill(storage_.begin(), storage_.end(), logAlpha);
    std::fill(rowsums_.begin(), rowsums_.end(), logRowSum);
}

TransitionMatrix::TransitionMatrix(double alpha) : alpha_(alpha) {
    storage_.resize(NUM_ALIGNMENT_STATES * NUM_ALIGNMENT_STATES);
    rowsums_.resize(NUM_ALIGNMENT_STATES);
    
    double logAlpha = std::log(alpha_);
    double logRowSum = std::log(NUM_ALIGNMENT_STATES * alpha_);
    
    std::fill(storage_.begin(), storage_.end(), logAlpha);
    std::fill(rowsums_.begin(), rowsums_.end(), logRowSum);
}

TransitionMatrix::TransitionMatrix(const TransitionMatrix& other)
    : storage_(other.storage_)
    , rowsums_(other.rowsums_)
    , alpha_(other.alpha_)
{
}

TransitionMatrix::TransitionMatrix(TransitionMatrix&& other) noexcept
    : storage_(std::move(other.storage_))
    , rowsums_(std::move(other.rowsums_))
    , alpha_(other.alpha_)
{
}

TransitionMatrix& TransitionMatrix::operator=(const TransitionMatrix& other) {
    if (this != &other) {
        storage_ = other.storage_;
        rowsums_ = other.rowsums_;
        alpha_ = other.alpha_;
    }
    return *this;
}

TransitionMatrix& TransitionMatrix::operator=(TransitionMatrix&& other) noexcept {
    if (this != &other) {
        storage_ = std::move(other.storage_);
        rowsums_ = std::move(other.rowsums_);
        alpha_ = other.alpha_;
    }
    return *this;
}

double TransitionMatrix::operator()(uint32_t prevState, uint32_t curState) const {
    size_t k = prevState * NUM_ALIGNMENT_STATES + curState;
    return storage_[k] - rowsums_[prevState];
}

void TransitionMatrix::increment(uint32_t prevState, uint32_t curState, double logMass) {
    size_t k = prevState * NUM_ALIGNMENT_STATES + curState;
    storage_[k] = logAdd(storage_[k], logMass);
    rowsums_[prevState] = logAdd(rowsums_[prevState], logMass);
}

void TransitionMatrix::computeRowSums() {
    for (size_t row = 0; row < NUM_ALIGNMENT_STATES; ++row) {
        double rowSum = LOG_0;
        for (size_t col = 0; col < NUM_ALIGNMENT_STATES; ++col) {
            size_t k = row * NUM_ALIGNMENT_STATES + col;
            rowSum = logAdd(rowSum, storage_[k]);
        }
        rowsums_[row] = rowSum;
    }
}

//=============================================================================
// TranscriptSequence implementation
//=============================================================================

// SAM base encoding: A=1, C=2, G=4, T=8, N=15
static uint8_t charToSamEncode(char c) {
    switch (c) {
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 4;
        case 'T': case 't': return 8;
        default: return 15;  // N
    }
}

TranscriptSequence::TranscriptSequence(const std::string& seq) : length_(seq.size()) {
    // Pack 2 bases per byte (4-bit each)
    sequence_.resize((length_ + 1) / 2);
    for (size_t i = 0; i < length_; ++i) {
        uint8_t base = charToSamEncode(seq[i]);
        size_t byte = i >> 1;
        if (i & 1) {
            // Lower nibble
            sequence_[byte] |= base;
        } else {
            // Upper nibble
            sequence_[byte] = (base << 4);
        }
    }
}

uint8_t TranscriptSequence::baseAt(size_t idx) const {
    if (idx >= length_) return 15;  // N for out of bounds
    size_t byte = idx >> 1;
    size_t shift = (!(idx & 1)) << 2;  // 4 if even index, 0 if odd
    return (sequence_[byte] >> shift) & 0x0F;
}

uint8_t TranscriptSequence::twoBitAt(size_t idx) const {
    uint8_t base = baseAt(idx);
    return samToTwoBit[base];
}

int32_t TranscriptSequence::gcFrac(int32_t start, int32_t end) const {
    if (start < 0 || end < 0 || start >= static_cast<int32_t>(length_) || 
        end >= static_cast<int32_t>(length_) || start > end) {
        return 0;
    }
    
    int gc_count = 0;
    int total_bases = end - start + 1;
    
    for (int32_t i = start; i <= end; ++i) {
        uint8_t base = baseAt(static_cast<size_t>(i));
        // SAM 4-bit encoding: C=2, G=4
        if (base == 2 || base == 4) {
            gc_count++;
        }
    }
    
    return static_cast<int32_t>(std::lrint(100.0 * gc_count / total_bases));
}

//=============================================================================
// Transcriptome implementation
//=============================================================================

bool Transcriptome::loadFromFasta(const std::string& fasta_path) {
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open FASTA file: " << fasta_path << std::endl;
        return false;
    }
    
    std::string line;
    std::string current_name;
    std::string current_seq;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        if (line[0] == '>') {
            // Save previous transcript
            if (!current_name.empty() && !current_seq.empty()) {
                // Extract first word of name (before space)
                size_t space_pos = current_name.find(' ');
                std::string short_name = (space_pos != std::string::npos) 
                    ? current_name.substr(0, space_pos) 
                    : current_name;
                
                name_to_idx_[short_name] = transcripts_.size();
                transcripts_.emplace_back(current_seq);
            }
            
            // Start new transcript
            current_name = line.substr(1);  // Remove '>'
            current_seq.clear();
        } else {
            // Append to sequence
            current_seq += line;
        }
    }
    
    // Save last transcript
    if (!current_name.empty() && !current_seq.empty()) {
        size_t space_pos = current_name.find(' ');
        std::string short_name = (space_pos != std::string::npos) 
            ? current_name.substr(0, space_pos) 
            : current_name;
        
        name_to_idx_[short_name] = transcripts_.size();
        transcripts_.emplace_back(current_seq);
    }
    
    std::cerr << "Loaded " << transcripts_.size() << " transcript sequences" << std::endl;
    return true;
}

const TranscriptSequence* Transcriptome::getTranscript(uint32_t tid) const {
    if (tid >= transcripts_.size()) return nullptr;
    return &transcripts_[tid];
}

const TranscriptSequence* Transcriptome::getTranscriptByName(const std::string& name) const {
    auto it = name_to_idx_.find(name);
    if (it == name_to_idx_.end()) return nullptr;
    return &transcripts_[it->second];
}

bool Transcriptome::reorderByNames(const std::vector<std::string>& names) {
    // Reorder transcripts to match the order in 'names' (BAM header order)
    std::vector<TranscriptSequence> reordered;
    reordered.reserve(names.size());
    
    size_t found = 0;
    size_t not_found = 0;
    
    for (const auto& name : names) {
        auto it = name_to_idx_.find(name);
        if (it != name_to_idx_.end() && it->second < transcripts_.size()) {
            reordered.push_back(std::move(transcripts_[it->second]));
            found++;
        } else {
            // Transcript not found in FASTA - add empty placeholder
            reordered.emplace_back();
            not_found++;
        }
    }
    
    if (not_found > 0) {
        std::cerr << "Warning: " << not_found << " BAM transcripts not found in FASTA, "
                  << found << " found\n";
    }
    
    // Replace with reordered transcripts
    transcripts_ = std::move(reordered);
    
    // Rebuild name_to_idx_ to match new ordering
    name_to_idx_.clear();
    for (size_t i = 0; i < names.size(); ++i) {
        name_to_idx_[names[i]] = static_cast<uint32_t>(i);
    }
    
    return not_found == 0;
}

//=============================================================================
// AlignmentModel implementation
//=============================================================================

AlignmentModel::AlignmentModel(double alpha, uint32_t readBins)
    : transitionProbsLeft_(readBins, TransitionMatrix(alpha))
    , transitionProbsRight_(readBins, TransitionMatrix(alpha))
    , isEnabled_(true)
    , useModel_(false)     // Don't use model until PRE_BURNIN_FRAGS
    , canUpdate_(true)     // Allow updates until BURNIN_FRAGS
    , readBins_(readBins)
    , numObserved_(0)
    , traceStream_(nullptr)
    , traceLevel_(ErrorModelTraceLevel::NONE)
    , currentQname_("")
    , currentTid_(0)
    , matrixDumpPrefix_("")
{
}

void AlignmentModel::setBasesFromCIGAROp(uint8_t op, size_t& curRefBase, size_t& curReadBase) {
    switch (op) {
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            // Both ref and read bases are used as-is
            break;
        case BAM_CINS:
            // Insertion: read base present, ref is gap
            curRefBase = ALN_DASH;
            break;
        case BAM_CDEL:
            // Deletion: ref base present, read is gap
            curReadBase = ALN_DASH;
            break;
        case BAM_CREF_SKIP:
            // Intron/skip: ref base present, read is skip marker
            curReadBase = ALN_REF_SKIP;
            break;
        case BAM_CSOFT_CLIP:
            // Soft clip: read base present, ref is soft clip marker
            curRefBase = ALN_SOFT_CLIP;
            break;
        case BAM_CHARD_CLIP:
            // Hard clip: both are hard clip markers
            curRefBase = ALN_HARD_CLIP;
            curReadBase = ALN_HARD_CLIP;
            break;
        case BAM_CPAD:
            // Padding: both are pad markers
            curRefBase = ALN_PAD;
            curReadBase = ALN_PAD;
            break;
        default:
            // Unknown operation
            break;
    }
}

AlnModelProb AlignmentModel::computeLogLikelihood(
    const uint32_t* cigar,
    uint32_t cigarLen,
    const uint8_t* readSeq,
    int32_t readLen,
    const TranscriptSequence* ref,
    int32_t refPos,
    std::vector<TransitionMatrix>& transitionProbs,
    uint32_t readPosBinStart
) {
    // readPosBinStart is unused but kept for API consistency
    (void)readPosBinStart;
    size_t readIdx = 0;
    int32_t transcriptIdx = refPos;
    size_t transcriptLen = ref ? ref->length() : 0;
    
    // Handle reads starting before reference
    if (transcriptIdx < 0) {
        readIdx = static_cast<size_t>(-transcriptIdx);
        transcriptIdx = 0;
    }
    
    size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);
    
    // Check bounds
    if (uTranscriptIdx >= transcriptLen) {
        return {LOG_0, LOG_1};
    }
    
    if (cigarLen == 0 || !cigar) {
        return {LOG_EPSILON, LOG_1};
    }
    
    // Foreground and background log-likelihoods
    double logLike = LOG_1;
    double bgLogLike = LOG_1;
    
    bool advanceInRead = false;
    bool advanceInReference = false;
    uint32_t readPosBin = 0;
    uint32_t prevStateIdx = START_STATE_IDX;
    uint32_t curStateIdx = 0;
    double invLen = static_cast<double>(readBins_) / readLen;
    
    for (uint32_t cigarIdx = 0; cigarIdx < cigarLen; ++cigarIdx) {
        uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
        uint8_t op = cigar[cigarIdx] & BAM_CIGAR_MASK;
        
        // Get initial bases for this CIGAR operation
        size_t curReadBase = 0;
        size_t curRefBase = 0;
        
        if (cigarConsumesQuery(op) && readIdx < static_cast<size_t>(readLen)) {
            // Get read base (4-bit packed, 2 per byte)
            size_t byteIdx = readIdx >> 1;
            size_t shift = (!(readIdx & 1)) << 2;  // 4 if even, 0 if odd
            uint8_t base4bit = (readSeq[byteIdx] >> shift) & 0x0F;
            curReadBase = samToTwoBit[base4bit];
        }
        
        if (cigarConsumesRef(op) && ref && uTranscriptIdx < transcriptLen) {
            curRefBase = ref->twoBitAt(uTranscriptIdx);
        }
        
        advanceInRead = false;
        advanceInReference = false;
        
        for (size_t i = 0; i < opLen; ++i) {
            // Advance if needed
            if (advanceInRead) {
                if (readIdx >= static_cast<size_t>(readLen)) {
                    // Read index out of bounds - return current likelihood
                    return {logLike, bgLogLike};
                }
                size_t byteIdx = readIdx >> 1;
                size_t shift = (!(readIdx & 1)) << 2;
                uint8_t base4bit = (readSeq[byteIdx] >> shift) & 0x0F;
                curReadBase = samToTwoBit[base4bit];
                readPosBin = static_cast<uint32_t>(readIdx * invLen);
                if (readPosBin >= readBins_) readPosBin = readBins_ - 1;
                advanceInRead = false;
            }
            
            if (advanceInReference) {
                if (uTranscriptIdx >= transcriptLen) {
                    // Reference index out of bounds - return current likelihood
                    return {logLike, bgLogLike};
                }
                if (ref) {
                    curRefBase = ref->twoBitAt(uTranscriptIdx);
                }
                advanceInReference = false;
            }
            
            // Apply CIGAR operation to determine state
            setBasesFromCIGAROp(op, curRefBase, curReadBase);
            
            // Compute state index: refBase * 9 + readBase
            curStateIdx = curRefBase * NUM_BASE_STATES + curReadBase;
            
            // Add transition probability to log-likelihood
            double tp = transitionProbs[readPosBin](prevStateIdx, curStateIdx);
            logLike += tp;
            
            // Background: always use match-match transition (state 0 -> 0)
            bgLogLike += transitionProbs[readPosBin](0, 0);
            
            // Trace transition if enabled
            if (traceStream_ && traceLevel_ >= ErrorModelTraceLevel::TRANSITION) {
                traceTransition(currentQname_, currentTid_, static_cast<uint32_t>(readIdx), readPosBin, prevStateIdx, curStateIdx, tp);
            }
            
            prevStateIdx = curStateIdx;
            
            // Mark what to advance for next iteration
            if (cigarConsumesQuery(op)) {
                ++readIdx;
                advanceInRead = true;
            }
            if (cigarConsumesRef(op)) {
                ++uTranscriptIdx;
                advanceInReference = true;
            }
        }
    }
    
    return {logLike, bgLogLike};
}

double AlignmentModel::logLikelihood(
    const uint32_t* cigar,
    uint32_t cigarLen,
    const uint8_t* readSeq,
    int32_t readLen,
    const TranscriptSequence* ref,
    int32_t refPos,
    bool isLeft
) {
    if (!isEnabled_) {
        return LOG_1;
    }
    
    auto& transitionProbs = isLeft ? transitionProbsLeft_ : transitionProbsRight_;
    auto result = computeLogLikelihood(cigar, cigarLen, readSeq, readLen, ref, refPos, transitionProbs, 0);
    
    // Trace alignment if enabled
    if (traceStream_ && traceLevel_ >= ErrorModelTraceLevel::ALIGNMENT) {
        uint32_t bin = 0;  // Use first bin for trace (could compute average)
        if (readLen > 0) {
            double invLen = static_cast<double>(readBins_) / readLen;
            bin = static_cast<uint32_t>(readLen / 2 * invLen);  // Middle position
            if (bin >= readBins_) bin = readBins_ - 1;
        }
        traceAlignment(currentQname_, currentTid_, result.fg, result.bg, result.fg - result.bg, bin);
    }
    
    // Return foreground - background
    return result.fg - result.bg;
}

double AlignmentModel::logLikelihoodAndUpdate(
    const uint32_t* cigar,
    uint32_t cigarLen,
    const uint8_t* readSeq,
    int32_t readLen,
    const TranscriptSequence* ref,
    int32_t refPos,
    bool isLeft,
    double logProb,
    double logMass
) {
    if (!isEnabled_) {
        return LOG_1;
    }
    
    // Salmon behavior:
    // - Before PRE_BURNIN_FRAGS: return LOG_1, update if canUpdate_
    // - After PRE_BURNIN_FRAGS: compute errLike, update if canUpdate_
    
    // Update model if we're still in the update phase
    if (canUpdate_) {
        update(cigar, cigarLen, readSeq, readLen, ref, refPos, isLeft, logProb, logMass);
    }
    
    // Return LOG_1 if not yet using model (pre-burnin phase)
    if (!useModel_) {
        return LOG_1;
    }
    
    // Compute and return error likelihood
    auto& transitionProbs = isLeft ? transitionProbsLeft_ : transitionProbsRight_;
    auto result = computeLogLikelihood(cigar, cigarLen, readSeq, readLen, ref, refPos, transitionProbs, 0);

    return result.fg - result.bg;
}

double AlignmentModel::logLikelihoodPairedAndUpdate(
    const uint32_t* cigar1, uint32_t cigarLen1,
    const uint8_t* readSeq1, int32_t readLen1,
    const uint32_t* cigar2, uint32_t cigarLen2,
    const uint8_t* readSeq2, int32_t readLen2,
    const TranscriptSequence* ref,
    int32_t refPos1, int32_t refPos2,
    double logProb,
    double logMass
) {
    if (!isEnabled_) {
        return LOG_1;
    }
    
    // Determine which read is left/right based on position
    bool read1IsLeft = (refPos1 <= refPos2);
    
    // Update model if we're still in the update phase
    if (canUpdate_) {
        if (cigar1 && cigarLen1 > 0) {
            update(cigar1, cigarLen1, readSeq1, readLen1, ref, refPos1, read1IsLeft, logProb, logMass);
        }
        if (cigar2 && cigarLen2 > 0) {
            update(cigar2, cigarLen2, readSeq2, readLen2, ref, refPos2, !read1IsLeft, logProb, logMass);
        }
    }
    
    // Return LOG_1 if not yet using model (pre-burnin phase)
    if (!useModel_) {
        return LOG_1;
    }
    
    // Compute and return error likelihood
    double logLike = LOG_1;
    double bg = LOG_1;
    
    if (cigar1 && cigarLen1 > 0) {
        auto& probs = read1IsLeft ? transitionProbsLeft_ : transitionProbsRight_;
        auto result = computeLogLikelihood(cigar1, cigarLen1, readSeq1, readLen1, ref, refPos1, probs, 0);
        logLike += result.fg;
        bg += result.bg;
    }
    
    if (cigar2 && cigarLen2 > 0) {
        auto& probs = read1IsLeft ? transitionProbsRight_ : transitionProbsLeft_;
        auto result = computeLogLikelihood(cigar2, cigarLen2, readSeq2, readLen2, ref, refPos2, probs, 0);
        logLike += result.fg;
        bg += result.bg;
    }
    
    return logLike - bg;
}

double AlignmentModel::logLikelihoodPaired(
    const uint32_t* cigar1, uint32_t cigarLen1,
    const uint8_t* readSeq1, int32_t readLen1,
    const uint32_t* cigar2, uint32_t cigarLen2,
    const uint8_t* readSeq2, int32_t readLen2,
    const TranscriptSequence* ref,
    int32_t refPos1, int32_t refPos2
) {
    if (!isEnabled_) {
        return LOG_1;
    }
    
    double logLike = LOG_1;
    double bg = LOG_1;
    
    // Determine which read is left/right based on position
    bool read1IsLeft = (refPos1 <= refPos2);
    
    if (cigar1 && cigarLen1 > 0) {
        auto& probs = read1IsLeft ? transitionProbsLeft_ : transitionProbsRight_;
        auto result = computeLogLikelihood(cigar1, cigarLen1, readSeq1, readLen1, ref, refPos1, probs, 0);
        logLike += result.fg;
        bg += result.bg;
    }
    
    if (cigar2 && cigarLen2 > 0) {
        auto& probs = read1IsLeft ? transitionProbsRight_ : transitionProbsLeft_;
        auto result = computeLogLikelihood(cigar2, cigarLen2, readSeq2, readLen2, ref, refPos2, probs, 0);
        logLike += result.fg;
        bg += result.bg;
    }
    
    return logLike - bg;
}

void AlignmentModel::update(
    const uint32_t* cigar,
    uint32_t cigarLen,
    const uint8_t* readSeq,
    int32_t readLen,
    const TranscriptSequence* ref,
    int32_t refPos,
    bool isLeft,
    double logProb,
    double logMass
) {
    if (logMass == LOG_0 || !isEnabled_) {
        return;
    }
    
    auto& transitionProbs = isLeft ? transitionProbsLeft_ : transitionProbsRight_;
    
    size_t readIdx = 0;
    int32_t transcriptIdx = refPos;
    size_t transcriptLen = ref ? ref->length() : 0;
    
    if (transcriptIdx < 0) {
        readIdx = static_cast<size_t>(-transcriptIdx);
        transcriptIdx = 0;
    }
    
    size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);
    
    if (cigarLen == 0 || !cigar) {
        return;
    }
    
    bool advanceInRead = false;
    bool advanceInReference = false;
    uint32_t readPosBin = 0;
    uint32_t prevStateIdx = START_STATE_IDX;
    uint32_t curStateIdx = 0;
    double invLen = static_cast<double>(readBins_) / readLen;
    
    for (uint32_t cigarIdx = 0; cigarIdx < cigarLen; ++cigarIdx) {
        uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
        uint8_t op = cigar[cigarIdx] & BAM_CIGAR_MASK;
        
        size_t curReadBase = 0;
        size_t curRefBase = 0;
        
        if (cigarConsumesQuery(op) && readIdx < static_cast<size_t>(readLen)) {
            size_t byteIdx = readIdx >> 1;
            size_t shift = (!(readIdx & 1)) << 2;
            uint8_t base4bit = (readSeq[byteIdx] >> shift) & 0x0F;
            curReadBase = samToTwoBit[base4bit];
        }
        
        if (cigarConsumesRef(op) && ref && uTranscriptIdx < transcriptLen) {
            curRefBase = ref->twoBitAt(uTranscriptIdx);
        }
        
        advanceInRead = false;
        advanceInReference = false;
        
        for (size_t i = 0; i < opLen; ++i) {
            if (advanceInRead) {
                if (readIdx >= static_cast<size_t>(readLen)) {
                    return;
                }
                size_t byteIdx = readIdx >> 1;
                size_t shift = (!(readIdx & 1)) << 2;
                uint8_t base4bit = (readSeq[byteIdx] >> shift) & 0x0F;
                curReadBase = samToTwoBit[base4bit];
                readPosBin = static_cast<uint32_t>(readIdx * invLen);
                if (readPosBin >= readBins_) readPosBin = readBins_ - 1;
                advanceInRead = false;
            }
            
            if (advanceInReference) {
                if (uTranscriptIdx >= transcriptLen) {
                    return;
                }
                if (ref) {
                    curRefBase = ref->twoBitAt(uTranscriptIdx);
                }
                advanceInReference = false;
            }
            
            setBasesFromCIGAROp(op, curRefBase, curReadBase);
            curStateIdx = curRefBase * NUM_BASE_STATES + curReadBase;
            
            // Update transition count
            transitionProbs[readPosBin].increment(prevStateIdx, curStateIdx, logMass + logProb);
            
            prevStateIdx = curStateIdx;
            
            if (cigarConsumesQuery(op)) {
                ++readIdx;
                advanceInRead = true;
            }
            if (cigarConsumesRef(op)) {
                ++uTranscriptIdx;
                advanceInReference = true;
            }
        }
    }
}

//=============================================================================
// Trace and dump support
//=============================================================================

void AlignmentModel::setTraceOutput(std::ostream* traceStream, ErrorModelTraceLevel level) {
    traceStream_ = traceStream;
    traceLevel_ = level;
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

void AlignmentModel::traceAlignment(const std::string& qname, uint32_t tid, double fg, double bg, double errLike, uint32_t bin) {
    if (traceStream_ && traceLevel_ >= ErrorModelTraceLevel::ALIGNMENT) {
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
    if (traceStream_ && traceLevel_ >= ErrorModelTraceLevel::TRANSITION) {
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

bool AlignmentModel::dumpMatrices(const std::string& prefix, const std::string& suffix) const {
    if (prefix.empty()) {
        return false;
    }
    
    // Dump left matrices
    for (uint32_t bin = 0; bin < readBins_; ++bin) {
        std::string filename = prefix + "_left_bin" + std::to_string(bin) + suffix + ".tsv";
        dumpMatrixToFile(transitionProbsLeft_[bin], filename, true, bin);
    }
    
    // Dump right matrices
    for (uint32_t bin = 0; bin < readBins_; ++bin) {
        std::string filename = prefix + "_right_bin" + std::to_string(bin) + suffix + ".tsv";
        dumpMatrixToFile(transitionProbsRight_[bin], filename, false, bin);
    }
    
    return true;
}

void AlignmentModel::dumpMatrixToFile(const TransitionMatrix& matrix, const std::string& filename, bool isLeft, uint32_t bin) const {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Warning: Could not open matrix dump file: " << filename << std::endl;
        return;
    }
    
    out << "# Transition matrix: " << (isLeft ? "left" : "right") 
        << ", bin=" << bin 
        << ", numObserved=" << numObserved_ 
        << "\n";
    out << "# Rows: prevState (0-81), Cols: curState (0-81)\n";
    
    for (uint32_t row = 0; row < NUM_ALIGNMENT_STATES; ++row) {
        for (uint32_t col = 0; col < NUM_ALIGNMENT_STATES; ++col) {
            if (col > 0) out << "\t";
            double logProb = matrix(row, col);
            out << std::setprecision(17) << logProb;
        }
        out << "\n";
    }
    
    out.close();
}

} // namespace libem
