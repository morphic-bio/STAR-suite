#include "alignment_filter.h"
#include <algorithm>
#include <cmath>

// Update per-transcript best scores (called per alignment)
// Matches Salmon's updateRefMappings logic from SalmonMappingUtils.hpp lines 213-269
void updateRefMappings(
    uint32_t tid, 
    int32_t score, 
    bool is_compat, 
    size_t idx,
    const std::vector<bool>& is_decoy,
    MappingScoreInfo& msi,
    std::vector<int32_t>& scores
) {
    constexpr int32_t invalid_score = std::numeric_limits<int32_t>::min();
    
    // Decoy handling
    if (tid < is_decoy.size() && is_decoy[tid]) {
        if (score > msi.best_decoy_score) {
            msi.best_decoy_score = score;
        }
        return;
    }
    
    // Decoy threshold check
    // If bestDecoyScore is invalid, set it to a value that won't filter everything
    if (msi.best_decoy_score == invalid_score) {
        msi.best_decoy_score = invalid_score + 1;
    }
    int32_t decoy_cutoff = static_cast<int32_t>(msi.decoy_thresh * msi.best_decoy_score);
    
    if (score < decoy_cutoff || score == invalid_score) {
        scores[idx] = invalid_score;
        return;
    }
    
    // Per-transcript best hit tracking (Salmon keeps only best hit per transcript)
    auto it = msi.best_score_per_transcript.find(tid);
    if (it == msi.best_score_per_transcript.end()) {
        // First alignment for this transcript
        msi.best_score_per_transcript[tid] = {score, idx};
    } else if (score > it->second.first || (score == it->second.first && is_compat)) {
        // Better score or equal score with compatibility preference
        scores[it->second.second] = invalid_score; // Invalidate old
        it->second = {score, idx};
    } else {
        // Already have a better alignment for this transcript
        scores[idx] = invalid_score;
        return;
    }
    
    // Update global best scores
    if (score > msi.best_score) {
        msi.second_best_score = msi.best_score;
        msi.best_score = score;
    }
}

// Filter and collect alignments matching Salmon's filterAndCollectAlignments
// Matches SalmonMappingUtils.hpp lines 271-390
std::vector<RawAlignment> filterAndCollectAlignments(
    std::vector<RawAlignment>& alignments,
    const FilterParams& params,
    MappingScoreInfo& msi
) {
    constexpr int32_t invalid_score = std::numeric_limits<int32_t>::min();
    std::vector<RawAlignment> filtered;
    
    // Handle invalid bestDecoyScore
    if (msi.best_decoy_score == invalid_score) {
        msi.best_decoy_score = invalid_score + 1;
    }
    
    int32_t decoy_threshold = static_cast<int32_t>(msi.decoy_thresh * msi.best_decoy_score);
    
    // Apply hard vs soft filter
    // Hard filter: keep only alignments with score == best_score
    // Soft filter: keep alignments above decoy_threshold
    int32_t filter_threshold = params.hard_filter 
        ? msi.best_score 
        : decoy_threshold;
    
    double best_score_d = static_cast<double>(msi.best_score);
    
    // Apply min_score_fraction threshold (Salmon's score fraction gating)
    // Only active when min_score_fraction > 0.0 (disabled by default for STAR inline mode)
    // Set to 0.65 for Salmon-parity mode to drop alignments < 65% of best score
    int32_t min_score_threshold = invalid_score;
    if (params.min_score_fraction > 0.0 && msi.best_score != invalid_score) {
        min_score_threshold = static_cast<int32_t>(best_score_d * params.min_score_fraction);
    }
    
    // Filter alignments
    for (auto& aln : alignments) {
        // Score threshold check
        if (aln.score < filter_threshold) {
            continue;
        }
        
        // Min score fraction check (only active when min_score_fraction > 0.0)
        // Default is 0.0 (disabled) to preserve STAR's prior alignment retention behavior
        if (min_score_threshold != invalid_score && aln.score < min_score_threshold) {
            continue;
        }
        
        // Compute estAlnProb = exp(-scoreExp * (bestScore - currScore))
        double v = best_score_d - static_cast<double>(aln.score);
        aln.est_aln_prob = params.hard_filter ? 1.0 : std::exp(-params.score_exp * v);
        
        // Skip if below minAlnProb (soft filter only)
        if (!params.hard_filter && aln.est_aln_prob < params.min_aln_prob) {
            continue;
        }
        
        filtered.push_back(aln);
    }
    
    // Discard if too many hits (Salmon's maxReadOccs behavior)
    if (filtered.size() > params.max_read_occs) {
        filtered.clear();
    }
    
    return filtered;
}
