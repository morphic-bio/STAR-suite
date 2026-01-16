#include "GeneResolver.h"
#include <unordered_set>
#include <unordered_map>

uint16_t resolveGeneFromCandidates(const std::vector<CandidateView>& candidates) {
    // Helper to compute a ranking score: prefer AS-like score, fall back to -NM.
    auto scoreOf = [](const CandidateView* cv) -> int {
        if (cv->asScore != 0) {
            return cv->asScore;
        } else if (cv->nm >= 0) {
            return -cv->nm;
        } else {
            return 0;
        }
    };

    // Partition by isGenomic
    std::vector<const CandidateView*> probeCandidates;
    std::vector<const CandidateView*> genomicCandidates;
    
    for (const CandidateView &cand : candidates) {
        if (!cand.isGenomic) {
            probeCandidates.push_back(&cand);
        } else {
            genomicCandidates.push_back(&cand);
        }
    }
    
    // 1. Probe path: apply probe CIGAR filter, then rank by AS/NM.
    std::vector<const CandidateView*> goodProbe;
    for (const CandidateView *pc : probeCandidates) {
        // Require at least one gene and a passing probe CIGAR
        if (!pc->probeCigarOk) continue;
        bool hasGene = (pc->geneIdx15 != 0);
        for (uint16_t g : pc->zgGeneIdx15) {
            if (g != 0) { hasGene = true; break; }
        }
        if (!hasGene) continue;
        goodProbe.push_back(pc);
    }
    
    if (!goodProbe.empty()) {
        // Find best score among good probes
        int bestScore = scoreOf(goodProbe.front());
        for (const CandidateView *pc : goodProbe) {
            int s = scoreOf(pc);
            if (s > bestScore) bestScore = s;
        }
        
        // Collect genes from top-scoring probes
        std::unordered_set<uint16_t> genes;
        for (const CandidateView *pc : goodProbe) {
            if (scoreOf(pc) != bestScore) continue;
            if (pc->geneIdx15 != 0) {
                genes.insert(pc->geneIdx15);
            }
            for (uint16_t g : pc->zgGeneIdx15) {
                if (g != 0) {
                    genes.insert(g);
                }
            }
        }
        
        if (genes.size() == 1) {
            return *genes.begin();
        } else {
            // Multiple distinct probe genes at best score → drop
            return 0;
        }
    }
    
    // 2. Genomic fallback: rank genomic hits by AS/NM, no CIGAR filter.
    std::vector<const CandidateView*> goodGenomic;
    for (const CandidateView *pc : genomicCandidates) {
        bool hasGene = false;
        for (uint16_t g : pc->zgGeneIdx15) {
            if (g != 0) { hasGene = true; break; }
        }
        if (!hasGene) continue;
        goodGenomic.push_back(pc);
    }
    
    if (!goodGenomic.empty()) {
        int bestScore = scoreOf(goodGenomic.front());
        for (const CandidateView *pc : goodGenomic) {
            int s = scoreOf(pc);
            if (s > bestScore) bestScore = s;
        }
        
        // Examine all top-scoring genomic candidates: all of their probe-mapped
        // genes (zgGeneIdx15) must agree on a single probe gene. If any
        // top-scoring candidate has multiple distinct probe genes, or top-tier
        // candidates disagree, we drop.
        bool firstSet = false;
        uint16_t chosen = 0;
        
        for (const CandidateView *pc : goodGenomic) {
            if (scoreOf(pc) != bestScore) continue;
            
            // Collect distinct genes for this candidate
            std::unordered_set<uint16_t> local;
            for (uint16_t g : pc->zgGeneIdx15) {
                if (g != 0) {
                    local.insert(g);
                }
            }
            
            if (local.empty()) {
                // No probe-mapped gene for this top-tier genomic candidate -> drop
                return 0;
            }
            if (local.size() > 1) {
                // Candidate maps to multiple distinct probe genes -> ambiguous
                return 0;
            }
            
            uint16_t gLocal = *local.begin();
            if (!firstSet) {
                chosen = gLocal;
                firstSet = true;
            } else if (gLocal != chosen) {
                // Different top-tier genomic candidates disagree on probe gene
                return 0;
            }
        }
        
        return firstSet ? chosen : 0;
    }
    
    // 3. No usable candidates
    return 0;
}
