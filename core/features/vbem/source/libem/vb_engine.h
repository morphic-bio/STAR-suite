#ifndef VB_ENGINE_H
#define VB_ENGINE_H

#include "em_types.h"

// Run Variational Bayes algorithm on equivalence classes
EMResult run_vb(const ECTable& ecs, TranscriptState& state, const EMParams& params);

// Compute ELBO (Evidence Lower BOund) for VB convergence (with effective-length weighting)
double compute_elbo(const ECTable& ecs, const double* abundances, const double* eff_lengths, double vb_prior);


// Compute unique EC counts per transcript (for VB initialization)
std::vector<double> compute_unique_counts(const ECTable& ecs);

// Check if transcript has unique evidence (appears in single-transcript ECs)
std::vector<bool> compute_unique_evidence(const ECTable& ecs);

#endif // VB_ENGINE_H
