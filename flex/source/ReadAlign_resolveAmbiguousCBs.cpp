#include "ReadAlign.h"
#include "solo/CbBayesianResolver.h"
#include "solo/CbCorrector.h"
#include "ParametersSolo.h"

// FNV-1a hash for CB sequence -> AmbigKey
ReadAlign::AmbigKey ReadAlign::hashCbSeq(const std::string &cbSeq) {
    // FNV-1a 64-bit hash
    constexpr uint64_t FNV_OFFSET = 14695981039346656037ULL;
    constexpr uint64_t FNV_PRIME = 1099511628211ULL;
    
    uint64_t hash = FNV_OFFSET;
    for (char c : cbSeq) {
        hash ^= static_cast<uint64_t>(c);
        hash *= FNV_PRIME;
    }
    return hash;
}

// Phase 2: Resolve accumulated ambiguous CBs using Bayesian inference
void ReadAlign::resolveAmbiguousCBs() {
    if (pendingAmbiguous_.empty()) {
        return;
    }
    
    // Need whitelist sequences for Bayesian comparison
    if (!P.pSolo.cbCorrector) {
        P.inOut->logMain << "resolveAmbiguousCBs: CbCorrector not available, skipping " 
                         << pendingAmbiguous_.size() << " ambiguous CBs" << endl;
        cbResolutionStats_.stillAmbiguous += pendingAmbiguous_.size();
        return;
    }
    
    // Get whitelist sequences from CbCorrector
    const std::vector<std::string> &whitelistSeqs = P.pSolo.cbCorrector->whitelist();
    
    // Create Bayesian resolver
    CbBayesianResolver resolver(whitelistSeqs.size(), &whitelistSeqs);
    
    uint64_t resolved = 0;
    uint64_t stillAmbiguous = 0;
    
    for (auto &kv : pendingAmbiguous_) {
        const AmbigKey key = kv.first;
        AmbiguousEntry &entry = kv.second;
        
        if (entry.candidateIdx.empty() || entry.umiCounts.empty()) {
            stillAmbiguous++;
            continue;
        }
        
        // Build context and candidates for resolver
        CBContext context(entry.cbSeq, entry.cbQual);
        
        std::vector<Candidate> candidates;
        candidates.reserve(entry.candidateIdx.size());
        for (uint32_t idx : entry.candidateIdx) {
            if (idx > 0 && idx <= whitelistSeqs.size()) {
                candidates.emplace_back(idx, whitelistSeqs[idx - 1], 0.0);
            }
        }
        
        if (candidates.empty()) {
            stillAmbiguous++;
            continue;
        }
        
        // Run Bayesian resolution
        BayesianResult result = resolver.resolve(context, candidates, entry.umiCounts);
        
        if (result.status == BayesianResult::Resolved && result.bestIdx > 0) {
            // Store resolved CB index for later lookup
            resolvedCbByKey_[key] = result.bestIdx;
            resolved++;
        } else {
            stillAmbiguous++;
        }
    }
    
    cbResolutionStats_.resolvedBayesian += resolved;
    cbResolutionStats_.stillAmbiguous += stillAmbiguous;
    
    if (resolved > 0 || stillAmbiguous > 0) {
        P.inOut->logMain << "resolveAmbiguousCBs: resolved " << resolved 
                         << ", still ambiguous " << stillAmbiguous << endl;
    }
}

