#ifndef CODE_CbBayesianResolver
#define CODE_CbBayesianResolver

#include <cstdint>
#include <cstddef>
#include <vector>
#include <unordered_map>
#include <string>

// CB Bayesian Resolver - Phase 2 ambiguity resolution
// Mirrors Cell Ranger's correct_barcode_bayesian(), compute_log_likelihood(), compute_prior()

// CB Context: sequence and quality scores
struct CBContext {
    std::string cbSeq;      // Observed CB sequence (16 bases)
    std::string cbQual;     // Phred quality scores (same length as cbSeq)
    
    CBContext() {}
    CBContext(const std::string &seq, const std::string &qual)
        : cbSeq(seq), cbQual(qual) {}
};

// Candidate: whitelist CB with optional frequency
struct Candidate {
    uint32_t whitelistIdx;  // 1-based whitelist index
    std::string whitelistSeq; // Whitelist CB sequence (for comparison)
    double frequency;       // Prior frequency (0.0 = use uniform prior)
    
    Candidate() : whitelistIdx(0), frequency(0.0) {}
    Candidate(uint32_t idx, const std::string &seq, double freq = 0.0)
        : whitelistIdx(idx), whitelistSeq(seq), frequency(freq) {}
};

struct BayesianResult {
    enum Status {
        Resolved,      // Confident resolution (posterior >= 0.9 and >= 2x runner-up)
        Ambiguous,     // Multiple candidates but no clear winner
        Unresolved     // Low confidence or single candidate
    };
    
    Status status;
    uint32_t bestIdx;       // 1-based whitelist index (only valid if Resolved)
    double posteriorBest;   // Best posterior probability
    double posteriorRunner; // Second-best posterior probability
    
    BayesianResult() : status(Unresolved), bestIdx(0), posteriorBest(0.0), posteriorRunner(0.0) {}
};

class CbBayesianResolver {
public:
    // Constructor: initialize with whitelist size and optional whitelist sequences
    // whitelistSize: number of CBs in whitelist
    // whitelistSeqs: optional pointer to whitelist sequences (for per-base comparison)
    explicit CbBayesianResolver(size_t whitelistSize, 
                                const std::vector<std::string> *whitelistSeqs = nullptr);
    
    // Resolve ambiguous CB using Bayesian inference (matching Cell Ranger's logic)
    // context: observed CB sequence and quality scores
    // candidates: candidate whitelist CBs with optional frequencies
    // umiCounts: UMI24 -> count histogram for this ambiguous CB
    BayesianResult resolve(const CBContext &context,
                          const std::vector<Candidate> &candidates,
                          const std::unordered_map<uint32_t, uint32_t> &umiCounts) const;
    
    // Get whitelist size
    size_t whitelistSize() const { return whitelistSize_; }

private:
    size_t whitelistSize_;
    const std::vector<std::string> *whitelistSeqs_; // Optional whitelist sequences
    
    // Convert Phred quality score to error probability
    // phred: Phred quality score (ASCII character, typically '!' to '~')
    // Returns: error probability (0.0 to 1.0)
    double phredToErrorProb(char phred) const;
    
    // Compute log-likelihood for a candidate CB given observed sequence and quality
    // context: observed CB sequence and quality scores
    // candidateSeq: whitelist CB sequence to compare against
    // Returns: log-likelihood (log P(observed | candidate))
    double computeLogLikelihood(const CBContext &context, const std::string &candidateSeq) const;
    
    // Compute prior probability for a candidate CB
    // candidate: candidate with optional frequency
    // Returns: prior probability (P(candidate))
    double computePrior(const Candidate &candidate) const;
    
    // UMI weighting: weight for a UMI count (matching process_features' weight_table if applicable)
    // For now: uniform weighting (return 1.0)
    double umiWeight(uint32_t umi24, uint32_t count) const;
    
    // Constants for confidence thresholds (matching Cell Ranger)
    static constexpr double MIN_POSTERIOR = 0.9;      // Minimum posterior for resolution
    static constexpr double MIN_RATIO = 2.0;          // Minimum ratio vs runner-up
    static constexpr char PHRED_BASE = 33;            // Phred base offset (ASCII '!')
    static constexpr double MIN_ERROR_PROB = 1e-10;   // Minimum error probability (avoid log(0))
    static constexpr double MAX_ERROR_PROB = 0.75;    // Maximum error probability (cap at 75%)
};

#endif

