#include "ReadAlign.h"
#include "solo/CbBayesianResolver.h"
#include "solo/CbCorrector.h"
#include "ParametersSolo.h"
#include <atomic>
#include <cstdlib>

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
    static const bool g_disableAmbigResolve =
        (std::getenv("STAR_DISABLE_AMBIG_CB_RESOLVE") != nullptr);
    static const bool g_debugAmbigResolve =
        (std::getenv("STAR_DEBUG_AMBIG_CB_RESOLVE") != nullptr);
    static std::atomic<uint64_t> g_debugAmbigResolveCount{0};
    const uint64_t kMaxDebugEntries = 50;
    if (g_disableAmbigResolve) {
        pendingAmbiguous_.clear();
        return;
    }

    if (pendingAmbiguous_.empty()) {
        return;
    }
    
    // Need whitelist sequences for Bayesian comparison
    if (!P.pSolo.cbCorrector) {
        if (P.inOut && P.inOut->logMain.good()) {
            P.inOut->logMain << "resolveAmbiguousCBs: CbCorrector not available, skipping " 
                             << pendingAmbiguous_.size() << " ambiguous CBs" << endl;
        }
        cbResolutionStats_.stillAmbiguous += pendingAmbiguous_.size();
        return;
    }

    if (g_debugAmbigResolve && P.inOut && P.inOut->logMain.good()) {
        P.inOut->logMain << "[AMBIG-CB-DEBUG] resolveAmbiguousCBs: pending="
                         << pendingAmbiguous_.size() << endl;
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

        bool badSeq = entry.cbSeq.empty() || entry.cbQual.empty() ||
                      (entry.cbSeq.size() != entry.cbQual.size());
        bool badCandidates = false;
        uint32_t minIdx = UINT32_MAX;
        uint32_t maxIdx = 0;
        for (uint32_t idx : entry.candidateIdx) {
            if (idx == 0 || idx > whitelistSeqs.size()) {
                badCandidates = true;
            }
            if (idx < minIdx) minIdx = idx;
            if (idx > maxIdx) maxIdx = idx;
        }
        bool badUmi = entry.umiCounts.empty();

        if (g_debugAmbigResolve) {
            uint64_t entryNo = g_debugAmbigResolveCount.fetch_add(1);
            if ((entryNo < kMaxDebugEntries) || badSeq || badCandidates || badUmi) {
                if (P.inOut && P.inOut->logMain.good()) {
                    P.inOut->logMain << "[AMBIG-CB-DEBUG] entry=" << entryNo
                                     << " key=0x" << std::hex << key << std::dec
                                     << " cbSeqLen=" << entry.cbSeq.size()
                                     << " cbQualLen=" << entry.cbQual.size()
                                     << " cbSeq=" << entry.cbSeq
                                     << " cbQual=" << entry.cbQual
                                     << " candidates=" << entry.candidateIdx.size()
                                     << " minIdx=" << (entry.candidateIdx.empty() ? 0 : minIdx)
                                     << " maxIdx=" << maxIdx
                                     << " umiCounts=" << entry.umiCounts.size()
                                     << " badSeq=" << badSeq
                                     << " badCandidates=" << badCandidates
                                     << " badUmi=" << badUmi
                                     << endl;
                }
            }
        }
        
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
        if (P.inOut && P.inOut->logMain.good()) {
            P.inOut->logMain << "resolveAmbiguousCBs: resolved " << resolved 
                             << ", still ambiguous " << stillAmbiguous << endl;
        }
    }
}
