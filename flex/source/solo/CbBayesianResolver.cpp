#include "solo/CbBayesianResolver.h"
#include <algorithm>
#include <cmath>
#include <limits>

CbBayesianResolver::CbBayesianResolver(size_t whitelistSize, 
                                       const std::vector<std::string> *whitelistSeqs)
    : whitelistSize_(whitelistSize), whitelistSeqs_(whitelistSeqs)
{
}

double CbBayesianResolver::phredToErrorProb(char phred) const {
    // Convert Phred quality score to error probability
    // Phred Q = -10 * log10(P_error)
    // P_error = 10^(-Q/10)
    // ASCII offset: '!' (33) = Q0, 'A' (65) = Q33, etc.
    
    if (phred < PHRED_BASE) {
        // Invalid quality: use maximum error probability
        return MAX_ERROR_PROB;
    }
    
    int q = static_cast<int>(phred) - PHRED_BASE;
    if (q < 0) {
        q = 0;
    }
    
    // Convert to error probability: 10^(-Q/10)
    double errorProb = std::pow(10.0, -q / 10.0);
    
    // Clamp to reasonable range
    if (errorProb < MIN_ERROR_PROB) {
        errorProb = MIN_ERROR_PROB;
    } else if (errorProb > MAX_ERROR_PROB) {
        errorProb = MAX_ERROR_PROB;
    }
    
    return errorProb;
}

double CbBayesianResolver::computeLogLikelihood(const CBContext &context, 
                                                 const std::string &candidateSeq) const {
    // Compute log-likelihood: log P(observed_seq | candidate_seq)
    // This matches Cell Ranger's compute_log_likelihood() logic
    
    if (context.cbSeq.length() != candidateSeq.length() || 
        context.cbSeq.length() != context.cbQual.length()) {
        // Length mismatch: return very low likelihood
        return -1e10;
    }
    
    double logLikelihood = 0.0;
    size_t len = context.cbSeq.length();
    
    // Accumulate log-likelihood over each base position
    // This matches process_features: compare observed sequence (may contain Ns)
    // against whitelist candidate sequence, using quality scores for per-base error probabilities
    for (size_t i = 0; i < len; ++i) {
        char observedBase = context.cbSeq[i];
        char candidateBase = candidateSeq[i];
        char phred = context.cbQual[i];
        
        // Handle N bases: treat as unknown with high error probability
        // N in observed sequence means we don't know the true base, so use maximum error prob
        if (observedBase == 'N' || observedBase == 'n') {
            // N base: use maximum error probability (uniform over all bases)
            // P(observed N | candidate) = 0.25 (uniform probability over A/C/G/T)
            logLikelihood += std::log(0.25 + MIN_ERROR_PROB);
            continue;
        }
        
        // Convert Phred to error probability
        double errorProb = phredToErrorProb(phred);
        double correctProb = 1.0 - errorProb;
        
        // Probability of observing this base given the candidate
        // If bases match: P(observed | candidate) = 1 - error_prob
        // If bases mismatch: P(observed | candidate) = error_prob / 3 (assuming uniform substitution)
        double baseProb;
        if (observedBase == candidateBase) {
            baseProb = correctProb;
        } else {
            // Mismatch: error probability divided by 3 (uniform substitution model)
            // This assumes equal probability of substituting to any of the 3 other bases
            baseProb = errorProb / 3.0;
        }
        
        // Add to log-likelihood (use log-sum to avoid underflow)
        logLikelihood += std::log(baseProb + MIN_ERROR_PROB);
    }
    
    return logLikelihood;
}

double CbBayesianResolver::computePrior(const Candidate &candidate) const {
    // Compute prior probability for a candidate CB
    // If frequency is provided (> 0), use it; otherwise use uniform prior
    
    if (candidate.frequency > 0.0) {
        // Use provided frequency as prior
        return candidate.frequency;
    } else {
        // Uniform prior: 1 / whitelist_size
        return 1.0 / whitelistSize_;
    }
}

double CbBayesianResolver::umiWeight(uint32_t /*umi24*/, uint32_t count) const {
    // UMI weighting: matching process_features' weight_table if applicable
    // For now: uniform weighting (return count as-is)
    // In Cell Ranger, this can be used to down-weight low-quality UMIs
    return static_cast<double>(count);
}

BayesianResult CbBayesianResolver::resolve(const CBContext &context,
                                          const std::vector<Candidate> &candidates,
                                          const std::unordered_map<uint32_t, uint32_t> &umiCounts) const {
    BayesianResult result;
    
    if (candidates.empty()) {
        result.status = BayesianResult::Unresolved;
        return result;
    }
    
    if (candidates.size() == 1) {
        // Single candidate: no ambiguity to resolve
        result.status = BayesianResult::Unresolved;
        result.bestIdx = candidates[0].whitelistIdx;
        result.posteriorBest = 1.0;
        result.posteriorRunner = 0.0;
        return result;
    }
    
    // Step 1: Compute log-posterior for each candidate
    // log P(candidate | observed) = log P(observed | candidate) + log P(candidate) + log P(UMIs | candidate)
    // We'll normalize later, so we can drop the P(observed) normalization constant
    
    std::vector<double> logPost(candidates.size(), 0.0);
    
    // Compute total UMI count for normalization
    double totalUmiWeight = 0.0;
    for (const auto &umi : umiCounts) {
        totalUmiWeight += umiWeight(umi.first, umi.second);
    }
    
    if (totalUmiWeight == 0.0) {
        // No UMI counts: cannot resolve
        result.status = BayesianResult::Unresolved;
        return result;
    }
    
    // For each candidate, compute log-posterior
    for (size_t i = 0; i < candidates.size(); ++i) {
        const Candidate &cand = candidates[i];
        
        // Get candidate sequence (from candidate struct or whitelist)
        std::string candidateSeq = cand.whitelistSeq;
        if (candidateSeq.empty() && whitelistSeqs_ && 
            cand.whitelistIdx > 0 && cand.whitelistIdx <= whitelistSeqs_->size()) {
            // Look up from whitelist if not provided
            candidateSeq = (*whitelistSeqs_)[cand.whitelistIdx - 1];
        }
        
        if (candidateSeq.empty()) {
            // Cannot compute likelihood without candidate sequence
            logPost[i] = -1e10;
            continue;
        }
        
        // Compute log-likelihood for CB sequence (using quality scores)
        double logLikelihoodCB = computeLogLikelihood(context, candidateSeq);
        
        // Compute log-likelihood for UMIs
        // For uniform UMI model: log P(UMIs | candidate) = sum over UMIs: weight * log(1/N)
        // where N is the number of possible UMIs (we normalize by total weight)
        // In practice, we use: sum over UMIs: weight * log(weight / total_weight)
        double logLikelihoodUMI = 0.0;
        for (const auto &umi : umiCounts) {
            double weight = umiWeight(umi.first, umi.second);
            double prob = weight / totalUmiWeight;
            logLikelihoodUMI += weight * std::log(prob + MIN_ERROR_PROB);
        }
        
        // Compute prior
        double prior = computePrior(cand);
        
        // Log-posterior = log(prior) + log_likelihood_CB + log_likelihood_UMI
        logPost[i] = std::log(prior + MIN_ERROR_PROB) + logLikelihoodCB + logLikelihoodUMI;
    }
    
    // Step 2: Normalize to posteriors (log-sum-exp trick for numerical stability)
    double maxLog = *std::max_element(logPost.begin(), logPost.end());
    double sum = 0.0;
    std::vector<double> posterior(logPost.size());
    
    for (size_t i = 0; i < logPost.size(); ++i) {
        posterior[i] = std::exp(logPost[i] - maxLog);
        sum += posterior[i];
    }
    
    // Normalize
    if (sum > 0.0) {
        for (double &p : posterior) {
            p /= sum;
        }
    } else {
        // All posteriors are zero (shouldn't happen, but handle gracefully)
        result.status = BayesianResult::Unresolved;
        return result;
    }
    
    // Step 3: Find best and runner-up
    size_t bestIdxPos = 0;
    size_t runnerPos = (candidates.size() > 1) ? 1 : 0;
    
    for (size_t i = 0; i < posterior.size(); ++i) {
        if (posterior[i] > posterior[bestIdxPos]) {
            runnerPos = bestIdxPos;
            bestIdxPos = i;
        } else if (i != bestIdxPos && posterior[i] > posterior[runnerPos]) {
            runnerPos = i;
        }
    }
    
    double bestP = posterior[bestIdxPos];
    double runnerP = (candidates.size() > 1) ? posterior[runnerPos] : 0.0;
    
    // Step 4: Apply confidence rule (>= 0.9 and >= 2x runner-up)
    if (bestP >= MIN_POSTERIOR && bestP >= MIN_RATIO * runnerP) {
        result.status = BayesianResult::Resolved;
        result.bestIdx = candidates[bestIdxPos].whitelistIdx;
        result.posteriorBest = bestP;
        result.posteriorRunner = runnerP;
    } else if (candidates.size() > 1) {
        result.status = BayesianResult::Ambiguous;
        result.bestIdx = candidates[bestIdxPos].whitelistIdx; // Still report best, but mark as ambiguous
        result.posteriorBest = bestP;
        result.posteriorRunner = runnerP;
    } else {
        result.status = BayesianResult::Unresolved;
    }
    
    return result;
}
