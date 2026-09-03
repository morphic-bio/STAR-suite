#ifndef CODE_CbBayesianResolver
#define CODE_CbBayesianResolver

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>

// CB Bayesian Resolver - Phase 2 ambiguity resolution
// Mirrors Cell Ranger's correct_barcode_bayesian(), compute_log_likelihood(), compute_prior()

namespace cb_bayesian {
constexpr char kPhredBase = 33;
constexpr double kMinErrorProb = 1e-10;
constexpr double kMaxErrorProb = 0.75;
constexpr char kDefaultPadQual = 'H';

inline void normalizeCbQual(std::string &cbQual, const std::string &cbSeq)
{
    if (cbQual.length() != cbSeq.length()) {
        if (cbQual.length() < cbSeq.length()) {
            cbQual.append(cbSeq.length() - cbQual.length(), kDefaultPadQual);
        } else {
            cbQual = cbQual.substr(0, cbSeq.length());
        }
    }
}

inline double phredToErrorProb(char phred)
{
    if (phred < kPhredBase) {
        return kMaxErrorProb;
    }

    int q = static_cast<int>(phred) - kPhredBase;
    if (q < 0) {
        q = 0;
    }

    double errorProb = std::pow(10.0, -q / 10.0);
    if (errorProb < kMinErrorProb) {
        errorProb = kMinErrorProb;
    } else if (errorProb > kMaxErrorProb) {
        errorProb = kMaxErrorProb;
    }
    return errorProb;
}

inline void accumulateCbQualityEvidence(const std::string &cbSeq,
                                        std::string cbQual,
                                        std::vector<double> &logLikMatch,
                                        std::vector<double> &logLikMismatch,
                                        uint32_t &evidenceReads)
{
    normalizeCbQual(cbQual, cbSeq);
    if (logLikMatch.size() != cbSeq.size()) {
        logLikMatch.assign(cbSeq.size(), 0.0);
    }
    if (logLikMismatch.size() != cbSeq.size()) {
        logLikMismatch.assign(cbSeq.size(), 0.0);
    }

    for (size_t i = 0; i < cbSeq.size(); ++i) {
        if (cbSeq[i] == 'N' || cbSeq[i] == 'n') {
            const double ll = std::log(0.25 + kMinErrorProb);
            logLikMatch[i] += ll;
            logLikMismatch[i] += ll;
            continue;
        }

        const double errorProb = phredToErrorProb(cbQual[i]);
        const double logMatch = std::log((1.0 - errorProb) + kMinErrorProb);
        const double logMismatch = std::log((errorProb / 3.0) + kMinErrorProb);
        logLikMatch[i] += logMatch;
        logLikMismatch[i] += logMismatch;
    }

    ++evidenceReads;
}

inline void mergeCbQualityEvidence(const std::vector<double> &srcMatch,
                                   const std::vector<double> &srcMismatch,
                                   uint32_t srcEvidenceReads,
                                   std::vector<double> &dstMatch,
                                   std::vector<double> &dstMismatch,
                                   uint32_t &dstEvidenceReads)
{
    if (srcMatch.empty() || srcMismatch.empty() || srcEvidenceReads == 0) {
        return;
    }

    if (dstMatch.empty()) {
        dstMatch = srcMatch;
        dstMismatch = srcMismatch;
        dstEvidenceReads = srcEvidenceReads;
        return;
    }

    if (dstMatch.size() != srcMatch.size() || dstMismatch.size() != srcMismatch.size()) {
        return;
    }

    for (size_t i = 0; i < dstMatch.size(); ++i) {
        dstMatch[i] += srcMatch[i];
        dstMismatch[i] += srcMismatch[i];
    }
    dstEvidenceReads += srcEvidenceReads;
}
} // namespace cb_bayesian

// CB Context: sequence and quality scores
struct CBContext {
    std::string cbSeq;      // Observed CB sequence (16 bases)
    std::string cbQual;     // Phred quality scores (same length as cbSeq)
    const std::vector<double> *aggLogLikMatch = nullptr;
    const std::vector<double> *aggLogLikMismatch = nullptr;
    uint32_t evidenceReads = 0;
    
    CBContext() {}
    CBContext(const std::string &seq, const std::string &qual)
        : cbSeq(seq), cbQual(qual) {}
    CBContext(const std::string &seq, const std::string &qual,
              const std::vector<double> &match,
              const std::vector<double> &mismatch,
              uint32_t reads)
        : cbSeq(seq), cbQual(qual), aggLogLikMatch(&match), aggLogLikMismatch(&mismatch), evidenceReads(reads) {}

    bool hasAggregatedEvidence() const
    {
        return aggLogLikMatch != nullptr
            && aggLogLikMismatch != nullptr
            && evidenceReads > 0
            && aggLogLikMatch->size() == cbSeq.size()
            && aggLogLikMismatch->size() == cbSeq.size();
    }
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

    // Resolve the legacy STARsolo multi-candidate payload exactly: each
    // candidate is weighted by its exact-barcode abundance and by the Phred
    // probability that its mismatching base is an error. candidateIdx is
    // 1-based; candidateQual contains the corresponding raw quality byte.
    BayesianResult resolveStarSolo(const std::vector<uint32_t> &candidateIdx,
                                    const std::vector<uint8_t> &candidateQual,
                                    const std::vector<uint32_t> &exactCbReadCount,
                                    int qsBase,
                                    uint32_t qsMax,
                                    double minPosterior) const;
    
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
    static constexpr char PHRED_BASE = cb_bayesian::kPhredBase;            // Phred base offset (ASCII '!')
    static constexpr double MIN_ERROR_PROB = cb_bayesian::kMinErrorProb;   // Minimum error probability (avoid log(0))
    static constexpr double MAX_ERROR_PROB = cb_bayesian::kMaxErrorProb;    // Maximum error probability (cap at 75%)
};

#endif
