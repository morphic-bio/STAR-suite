#ifndef ALIGNMENT_MODEL_H
#define ALIGNMENT_MODEL_H

#include <vector>
#include <string>
#include <cstdint>
#include <cmath>
#include <mutex>
#include <unordered_map>
#include <memory>
#include <fstream>
#include <functional>

namespace libem {

// Alignment state constants (matching Salmon's AlignmentCommon.hpp)
enum AlignmentState : uint8_t {
    ALN_A = 0,
    ALN_C = 1,
    ALN_G = 2,
    ALN_T = 3,
    ALN_DASH = 4,       // Deletion/insertion gap
    ALN_SOFT_CLIP = 5,
    ALN_HARD_CLIP = 6,
    ALN_PAD = 7,
    ALN_REF_SKIP = 8
};

// State space: 9 base states × 9 = 81 states + 1 start state = 82 total
constexpr uint32_t NUM_ALIGNMENT_STATES = 82;
constexpr uint32_t NUM_BASE_STATES = 9;
constexpr uint32_t START_STATE_IDX = 81;

// SAM base encoding to two-bit (A=0, C=1, G=2, T=3)
// Index by SAM 4-bit encoding: seq_nt16_table
constexpr uint8_t samToTwoBit[] = {
    0, /*A*/ 0, /*C*/ 1, 0, /*G*/ 2, 0, 0, 0, /*T*/ 3, 0, 0, 0, 0, 0, 0, 0
};

// Transition probability matrix (stored in log space)
// Rows: previous state, Columns: current state
// Salmon uses numAlignmentStates() = 82: state = refBase * 9 + readBase (81 states + 1 start state)
// Note: Using non-atomic doubles since our CLI is single-threaded
class TransitionMatrix {
public:
    TransitionMatrix();
    TransitionMatrix(double alpha);  // Initialize with Dirichlet prior alpha
    
    // Copy/move constructors (needed for vector storage)
    TransitionMatrix(const TransitionMatrix& other);
    TransitionMatrix(TransitionMatrix&& other) noexcept;
    TransitionMatrix& operator=(const TransitionMatrix& other);
    TransitionMatrix& operator=(TransitionMatrix&& other) noexcept;
    
    // Get normalized log probability of transition
    double operator()(uint32_t prevState, uint32_t curState) const;
    
    // Increment transition count (log space)
    void increment(uint32_t prevState, uint32_t curState, double logMass);
    
    // Compute row sums for normalization
    void computeRowSums();
    
private:
    std::vector<double> storage_;  // nRow * nCol (non-atomic for simplicity)
    std::vector<double> rowsums_;  // nRow
    double alpha_;
};

// Transcript sequence storage (SAM 4-bit encoding, 2 bases per byte)
class TranscriptSequence {
public:
    TranscriptSequence() : length_(0) {}
    TranscriptSequence(const std::string& seq);
    
    // Get base at position (returns SAM 4-bit encoding)
    uint8_t baseAt(size_t idx) const;
    
    // Get two-bit encoding at position (A=0, C=1, G=2, T=3)
    uint8_t twoBitAt(size_t idx) const;
    
    // Compute GC fraction (0-100) for a fragment [start, end] (inclusive)
    int32_t gcFrac(int32_t start, int32_t end) const;
    
    size_t length() const { return length_; }
    
private:
    std::vector<uint8_t> sequence_;  // Packed 4-bit encoding
    size_t length_;
};

// Transcriptome: collection of transcript sequences
class Transcriptome {
public:
    Transcriptome() {}
    
    // Load from FASTA file
    bool loadFromFasta(const std::string& fasta_path);
    
    // Get transcript sequence by index
    const TranscriptSequence* getTranscript(uint32_t tid) const;
    
    // Get transcript by name
    const TranscriptSequence* getTranscriptByName(const std::string& name) const;
    
    // Number of transcripts
    size_t size() const { return transcripts_.size(); }
    
    // Reorder transcripts to match BAM header order
    // Returns true if reordering was successful, false if names don't match
    bool reorderByNames(const std::vector<std::string>& names);
    
private:
    std::vector<TranscriptSequence> transcripts_;
    std::unordered_map<std::string, uint32_t> name_to_idx_;
};

// Result of log-likelihood computation
struct AlnModelProb {
    double fg;  // Foreground log-likelihood
    double bg;  // Background log-likelihood
    
    AlnModelProb() : fg(0.0), bg(0.0) {}
    AlnModelProb(double fg_, double bg_) : fg(fg_), bg(bg_) {}
};

// Salmon thresholds for error model:
// - Before PRE_BURNIN_FRAGS: return LOG_1, update model
// - Between PRE_BURNIN_FRAGS and BURNIN_FRAGS: compute errLike, update model
// - After BURNIN_FRAGS: compute errLike, don't update (model is "burned in")
constexpr size_t PRE_BURNIN_FRAGS = 5000;      // When to start USING the error model
constexpr size_t BURNIN_FRAGS = 5000000;       // When to STOP updating (5 million)

// Trace levels for error model instrumentation
enum class ErrorModelTraceLevel {
    NONE = 0,      // No tracing
    READ = 1,      // Per-read summary
    ALIGNMENT = 2, // Per-alignment details
    TRANSITION = 3 // Per-transition details
};

// Alignment error model
// Computes log-likelihood of alignment given CIGAR string and reference
class AlignmentModel {
public:
    AlignmentModel(double alpha = 0.001, uint32_t readBins = 4);
    
    // Compute log-likelihood AND update model if still in burn-in
    // This is the main method to use during processing
    // Returns fg - bg (foreground minus background log-likelihood)
    double logLikelihoodAndUpdate(
        const uint32_t* cigar,
        uint32_t cigarLen,
        const uint8_t* readSeq,
        int32_t readLen,
        const TranscriptSequence* ref,
        int32_t refPos,
        bool isLeft,
        double logProb = 0.0,        // Log probability for weighting update
        double logMass = 0.0         // Log mass for weighting update
    );
    
    // Compute log-likelihood of alignment (no update)
    // Returns fg - bg (foreground minus background log-likelihood)
    double logLikelihood(
        const uint32_t* cigar,       // CIGAR array
        uint32_t cigarLen,           // Number of CIGAR operations
        const uint8_t* readSeq,      // Read sequence (4-bit packed)
        int32_t readLen,             // Read length
        const TranscriptSequence* ref,  // Reference transcript
        int32_t refPos,              // Alignment position on reference
        bool isLeft                  // True for left read, false for right
    );
    
    // Compute log-likelihood for paired-end read (with update if in burn-in)
    double logLikelihoodPairedAndUpdate(
        const uint32_t* cigar1, uint32_t cigarLen1,
        const uint8_t* readSeq1, int32_t readLen1,
        const uint32_t* cigar2, uint32_t cigarLen2,
        const uint8_t* readSeq2, int32_t readLen2,
        const TranscriptSequence* ref,
        int32_t refPos1, int32_t refPos2,
        double logProb = 0.0,
        double logMass = 0.0
    );
    
    // Compute log-likelihood for paired-end read (no update)
    double logLikelihoodPaired(
        const uint32_t* cigar1, uint32_t cigarLen1,
        const uint8_t* readSeq1, int32_t readLen1,
        const uint32_t* cigar2, uint32_t cigarLen2,
        const uint8_t* readSeq2, int32_t readLen2,
        const TranscriptSequence* ref,
        int32_t refPos1, int32_t refPos2
    );
    
    // Update model with observed alignment (for learning)
    void update(
        const uint32_t* cigar,
        uint32_t cigarLen,
        const uint8_t* readSeq,
        int32_t readLen,
        const TranscriptSequence* ref,
        int32_t refPos,
        bool isLeft,
        double logProb,
        double logMass
    );
    
    // Increment observed count and update state flags
    void incrementObserved() {
        ++numObserved_;
        // After PRE_BURNIN_FRAGS: start using the model
        if (!useModel_ && numObserved_ >= PRE_BURNIN_FRAGS) {
            useModel_ = true;
        }
        // After BURNIN_FRAGS: stop updating the model
        if (canUpdate_ && numObserved_ >= BURNIN_FRAGS) {
            canUpdate_ = false;
        }
    }
    
    // Get number of observed fragments
    size_t numObserved() const { return numObserved_; }
    
    // Enable/disable the model entirely
    void setEnabled(bool enabled) { isEnabled_ = enabled; }
    bool isEnabled() const { return isEnabled_; }
    
    // Check if model should be used for errLike computation (after PRE_BURNIN_FRAGS)
    bool useModel() const { return useModel_; }
    
    // Check if model should still be updated (before BURNIN_FRAGS)  
    bool canUpdate() const { return canUpdate_; }
    
    // Legacy API - now means "past pre-burnin, using model"
    bool burnedIn() const { return useModel_; }
    
    // Trace and dump support
    void setTraceOutput(std::ostream* traceStream, ErrorModelTraceLevel level);
    void setTraceContext(const std::string& qname, uint32_t tid);
    void clearTraceContext();
    void setMatrixDumpPrefix(const std::string& prefix) { matrixDumpPrefix_ = prefix; }
    
    // Dump transition matrices to files
    bool dumpMatrices(const std::string& prefix, const std::string& suffix) const;
    
    // Trace methods (public for CLI access)
    void traceRead(const std::string& qname, size_t numAlns, double errLikeSum, bool modelUsed, bool modelUpdated);
    
private:
    // Core log-likelihood computation
    AlnModelProb computeLogLikelihood(
        const uint32_t* cigar,
        uint32_t cigarLen,
        const uint8_t* readSeq,
        int32_t readLen,
        const TranscriptSequence* ref,
        int32_t refPos,
        std::vector<TransitionMatrix>& transitionProbs,
        uint32_t readPosBinStart = 0  // Starting bin for trace (default 0)
    );
    
    // Set bases from CIGAR operation
    static void setBasesFromCIGAROp(
        uint8_t op,
        size_t& curRefBase,
        size_t& curReadBase
    );
    
    std::vector<TransitionMatrix> transitionProbsLeft_;
    std::vector<TransitionMatrix> transitionProbsRight_;
    
    bool isEnabled_;
    bool useModel_;    // True after PRE_BURNIN_FRAGS: use model for errLike
    bool canUpdate_;   // True until BURNIN_FRAGS: allow model updates
    uint32_t readBins_;
    size_t numObserved_;  // Number of fragments observed (for burn-in tracking)
    
    std::mutex outputMutex_;
    
        // Trace support (private helpers)
        void traceAlignment(const std::string& qname, uint32_t tid, double fg, double bg, double errLike, uint32_t bin);
        void traceTransition(const std::string& qname, uint32_t tid, uint32_t readPos, uint32_t bin, uint32_t prevState, uint32_t curState, double logProb);
    
    // Trace state
    std::ostream* traceStream_;
    ErrorModelTraceLevel traceLevel_;
    std::string currentQname_;
    uint32_t currentTid_;
    std::string matrixDumpPrefix_;
    
    // Matrix dump support
    void dumpMatrixToFile(const TransitionMatrix& matrix, const std::string& filename, bool isLeft, uint32_t bin) const;
};

// CIGAR operation constants (from htslib)
#ifndef BAM_CIGAR_SHIFT
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  0xf
#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#endif

// Check if CIGAR op consumes query sequence
inline bool cigarConsumesQuery(uint8_t op) {
    return op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP ||
           op == BAM_CEQUAL || op == BAM_CDIFF;
}

// Check if CIGAR op consumes reference
inline bool cigarConsumesRef(uint8_t op) {
    return op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
           op == BAM_CEQUAL || op == BAM_CDIFF;
}

} // namespace libem

#endif // ALIGNMENT_MODEL_H
