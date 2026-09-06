#include "GeneResolver.h"
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

    auto hasGene = [](const CandidateView& cv) {
        if (cv.geneIdx15 != 0) return true;
        for (uint16_t gene : cv.zgGeneIdx15) {
            if (gene != 0) return true;
        }
        return false;
    };

    // Probe-contig and genomic alignments use the same STAR score scale. Find
    // the best eligible alignment tier across both kinds so an equal-scoring
    // cross-gene genomic hit cannot be hidden by a probe hit.
    bool haveBest = false;
    int bestScore = 0;
    for (const CandidateView& candidate : candidates) {
        if (!hasGene(candidate)) continue;
        const int score = scoreOf(&candidate);
        if (!haveBest || score > bestScore) {
            bestScore = score;
            haveBest = true;
        }
    }
    if (!haveBest) return 0;

    uint16_t chosen = 0;
    auto addGene = [&chosen](uint16_t gene) {
        if (gene == 0) return true;
        if (chosen == 0) {
            chosen = gene;
            return true;
        }
        return chosen == gene;
    };

    for (const CandidateView& candidate : candidates) {
        if (!hasGene(candidate) || scoreOf(&candidate) != bestScore) continue;
        if (!addGene(candidate.geneIdx15)) return 0;
        for (uint16_t gene : candidate.zgGeneIdx15) {
            if (!addGene(gene)) return 0;
        }
    }
    return chosen;
}
