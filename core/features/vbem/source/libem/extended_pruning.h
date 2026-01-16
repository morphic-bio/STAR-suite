#ifndef EXTENDED_PRUNING_H
#define EXTENDED_PRUNING_H

#include "em_types.h"
#include <vector>

// Configuration for extended pruning (non-Salmon features)
struct ExtendedPruningParams {
    bool enable_local_pruning = false;    // Default OFF for Salmon parity
    double min_local_prob = 1e-3;
    
    bool enable_global_pruning = false;   // Default OFF for Salmon parity
    double min_global_abundance = 1e-7;
};

// Remove low-probability transcripts from individual ECs
// NOT part of Salmon - use only when enable_local_pruning = true
void applyLocalPruning(EC& ec, double min_local_prob);

// Remove transcripts with low total weight across all ECs
// NOT part of Salmon - use only when enable_global_pruning = true
void applyGlobalPruning(ECTable& ecs, double min_global_abundance);

// Apply both pruning steps if enabled
void applyExtendedPruning(ECTable& ecs, const ExtendedPruningParams& params);

#endif // EXTENDED_PRUNING_H
