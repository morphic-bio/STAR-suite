#include "extended_pruning.h"
#include <algorithm>
#include <cmath>
#include <unordered_set>

// Remove low-probability transcripts from individual ECs
// NOT part of Salmon - use only when enable_local_pruning = true
void applyLocalPruning(EC& ec, double min_local_prob) {
    if (ec.transcript_ids.empty() || ec.weights.empty()) return;
    
    // Find transcripts to keep
    std::vector<uint32_t> kept_ids;
    std::vector<double> kept_weights;
    
    for (size_t i = 0; i < ec.transcript_ids.size(); ++i) {
        if (ec.weights[i] >= min_local_prob) {
            kept_ids.push_back(ec.transcript_ids[i]);
            kept_weights.push_back(ec.weights[i]);
        }
    }
    
    // If nothing kept, discard the EC entirely
    if (kept_ids.empty()) {
        ec.transcript_ids.clear();
        ec.weights.clear();
        ec.count = 0.0;
        return;
    }
    
    // Re-normalize weights to sum to 1.0
    double weight_sum = 0.0;
    for (double w : kept_weights) {
        weight_sum += w;
    }
    
    if (weight_sum > 0.0) {
        for (double& w : kept_weights) {
            w /= weight_sum;
        }
    }
    
    ec.transcript_ids = kept_ids;
    ec.weights = kept_weights;
}

// Remove transcripts with low total weight across all ECs
// NOT part of Salmon - use only when enable_global_pruning = true
void applyGlobalPruning(ECTable& ecs, double min_global_abundance) {
    if (ecs.ecs.empty() || ecs.n_transcripts == 0) return;
    
    // Accumulate total weight per transcript
    std::vector<double> global_weights(ecs.n_transcripts, 0.0);
    
    for (const EC& ec : ecs.ecs) {
        for (size_t i = 0; i < ec.transcript_ids.size(); ++i) {
            uint32_t tid = ec.transcript_ids[i];
            if (tid < ecs.n_transcripts) {
                double weight = ec.weights.empty() ? 1.0 / ec.transcript_ids.size() : ec.weights[i];
                global_weights[tid] += ec.count * weight;
            }
        }
    }
    
    // Identify dead transcripts
    std::unordered_set<uint32_t> dead_transcripts;
    for (size_t i = 0; i < global_weights.size(); ++i) {
        if (global_weights[i] < min_global_abundance) {
            dead_transcripts.insert(static_cast<uint32_t>(i));
        }
    }
    
    if (dead_transcripts.empty()) return;
    
    // Remove dead transcripts from all ECs
    for (EC& ec : ecs.ecs) {
        std::vector<uint32_t> kept_ids;
        std::vector<double> kept_weights;
        
        for (size_t i = 0; i < ec.transcript_ids.size(); ++i) {
            if (dead_transcripts.find(ec.transcript_ids[i]) == dead_transcripts.end()) {
                kept_ids.push_back(ec.transcript_ids[i]);
                if (!ec.weights.empty()) {
                    kept_weights.push_back(ec.weights[i]);
                }
            }
        }
        
        // If EC becomes empty, mark for removal
        if (kept_ids.empty()) {
            ec.count = 0.0;
            continue;
        }
        
        // Re-normalize weights if needed
        if (!kept_weights.empty()) {
            double weight_sum = 0.0;
            for (double w : kept_weights) {
                weight_sum += w;
            }
            if (weight_sum > 0.0) {
                for (double& w : kept_weights) {
                    w /= weight_sum;
                }
            }
        }
        
        ec.transcript_ids = kept_ids;
        ec.weights = kept_weights;
    }
    
    // Remove ECs with zero count
    ecs.ecs.erase(
        std::remove_if(ecs.ecs.begin(), ecs.ecs.end(),
            [](const EC& ec) { return ec.count == 0.0 || ec.transcript_ids.empty(); }),
        ecs.ecs.end()
    );
    
    ecs.n_ecs = ecs.ecs.size();
}

// Apply both pruning steps if enabled
void applyExtendedPruning(ECTable& ecs, const ExtendedPruningParams& params) {
    if (!params.enable_local_pruning && !params.enable_global_pruning) {
        return;
    }
    
    // Apply local pruning first
    if (params.enable_local_pruning) {
        for (EC& ec : ecs.ecs) {
            applyLocalPruning(ec, params.min_local_prob);
        }
        
        // Remove ECs that became empty
        ecs.ecs.erase(
            std::remove_if(ecs.ecs.begin(), ecs.ecs.end(),
                [](const EC& ec) { return ec.transcript_ids.empty() || ec.count == 0.0; }),
            ecs.ecs.end()
        );
    }
    
    // Apply global pruning second
    if (params.enable_global_pruning) {
        applyGlobalPruning(ecs, params.min_global_abundance);
    }
    
    ecs.n_ecs = ecs.ecs.size();
}
