#ifndef VB_ENGINE_H
#define VB_ENGINE_H

#include "em_types.h"

// Salmon-style VB initialization, exported so a setup/sharding stage can seed
// shards with exactly what run_vb would have used. Sharing one implementation
// is deliberate: computing alpha separately in a splitter is how the two drift,
// and an initialization mismatch changes results silently.
std::vector<double> compute_initial_alpha(const ECTable& ecs,
                                          const TranscriptState& state,
                                          const EMParams& params);

// Run Variational Bayes algorithm on equivalence classes
EMResult run_vb(const ECTable& ecs, TranscriptState& state, const EMParams& params);

// Compute ELBO (Evidence Lower BOund) for VB convergence (with effective-length weighting)
double compute_elbo(const ECTable& ecs, const double* abundances, const double* eff_lengths, double vb_prior);


// Compute unique EC counts per transcript (for VB initialization)
std::vector<double> compute_unique_counts(const ECTable& ecs);

// Check if transcript has unique evidence (appears in single-transcript ECs)
std::vector<bool> compute_unique_evidence(const ECTable& ecs);

#endif // VB_ENGINE_H
