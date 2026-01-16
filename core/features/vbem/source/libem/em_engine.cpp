#include "em_engine.h"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <cstring>

// Compute projected counts from ECs (Salmon-style initialization)
// This distributes multi-mapper counts using EC weights * (1/effLen)
std::vector<double> compute_projected_counts_em(const ECTable& ecs, const std::vector<double>& eff_lengths) {
    std::vector<double> projected_counts(ecs.n_transcripts, 0.0);
    
    for (const EC& ec : ecs.ecs) {
        if (ec.count <= 0) continue;
        
        size_t groupSize = ec.transcript_ids.size();
        
        if (groupSize == 1) {
            // Single-transcript EC: full count
            projected_counts[ec.transcript_ids[0]] += ec.count;
        } else {
            // Multi-transcript EC: distribute by weights * (1/effLen)
            double total_weight = 0.0;
            std::vector<double> weights(groupSize);
            
            for (size_t i = 0; i < groupSize; ++i) {
                uint32_t tid = ec.transcript_ids[i];
                double effLen_factor = (eff_lengths[tid] > 0) ? 1.0 / eff_lengths[tid] : 0.0;
                if (ec.has_weights() && i < ec.weights.size()) {
                    weights[i] = ec.weights[i] * effLen_factor;
                } else {
                    weights[i] = effLen_factor;
                }
                total_weight += weights[i];
            }
            
            if (total_weight > 0) {
                for (size_t i = 0; i < groupSize; ++i) {
                    uint32_t tid = ec.transcript_ids[i];
                    projected_counts[tid] += ec.count * (weights[i] / total_weight);
                }
            }
        }
    }
    
    return projected_counts;
}

// Initialize abundances (Salmon-style: projected counts + fracObserved mixing)
void initialize_abundances(TranscriptState& state, const EMParams& params) {
    // Salmon uses projected counts with fracObserved mixing, not uniform
    // This is now handled in run_em() after we have access to ECs
    // For backward compatibility, this function does basic init
    if (params.init_by_length) {
        double total_length = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            total_length += state.lengths[i];
        }
        
        if (total_length > 0) {
            for (size_t i = 0; i < state.n; ++i) {
                state.abundances[i] = state.lengths[i] / total_length;
            }
        } else {
            double uniform = 1.0 / state.n;
            for (size_t i = 0; i < state.n; ++i) {
                state.abundances[i] = uniform;
            }
        }
    } else {
        double uniform = 1.0 / state.n;
        for (size_t i = 0; i < state.n; ++i) {
            state.abundances[i] = uniform;
        }
    }
}

// Compute log-likelihood using effective-length-weighted probabilities
// P(read | EC) = sum_i (abundance[i] / eff_length[i]) / sum_j (abundance[j] / eff_length[j])
// This matches Salmon's approach where theta_i / effLen_i gives the probability
double compute_log_likelihood(const ECTable& ecs, const double* abundances, const double* eff_lengths) {
    double ll = 0.0;
    for (const EC& ec : ecs.ecs) {
        double prob = 0.0;
        for (uint32_t tid : ec.transcript_ids) {
            // Weight by abundance / effective_length (Salmon's approach)
            if (eff_lengths[tid] > 0) {
                prob += abundances[tid] / eff_lengths[tid];
            }
        }
        if (prob > 0) {
            ll += ec.count * std::log(prob);
        }
    }
    return ll;
}

void zero_low_abundance(std::vector<double>& counts, TranscriptState& state, double threshold) {
    // Zero out transcripts with counts below threshold (Salmon-like behavior)
    size_t n_zeros = 0;
    double total_counts = 0.0;
    
    for (size_t i = 0; i < state.n; ++i) {
        if (counts[i] < threshold) {
            counts[i] = 0.0;
            state.abundances[i] = 0.0;
            n_zeros++;
        } else {
            total_counts += counts[i];
        }
    }
    
    // Renormalize abundances for non-zero transcripts
    if (total_counts > 0 && n_zeros > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            if (counts[i] > 0) {
                state.abundances[i] = counts[i] / total_counts;
            }
        }
    }
}

EMResult run_em(const ECTable& ecs, TranscriptState& state, const EMParams& params) {
    EMResult result;
    result.counts.resize(state.n);
    result.tpm.resize(state.n);
    
    // Salmon-style initialization: projected counts + fracObserved mixing
    // This matches CollapsedEMOptimizer.cpp lines 829-895
    std::vector<double> projectedCounts = compute_projected_counts_em(ecs, state.eff_lengths);
    
    // Compute total weight
    double totalWeight = 0.0;
    size_t numActive = state.n;  // Use all transcripts like Salmon
    for (size_t i = 0; i < state.n; ++i) {
        totalWeight += projectedCounts[i];
    }
    
    // Salmon's fracObserved calculation: min(0.999, totalWeight / numRequiredFragments)
    const double numRequiredFragments = 5e6;  // Salmon default
    double maxFrac = 0.999;
    double fracObserved = std::min(maxFrac, totalWeight / numRequiredFragments);
    
    // Salmon's uniformPrior: totalWeight / numActive
    double uniformPrior = (numActive > 0) ? (totalWeight / static_cast<double>(numActive)) : 1.0;
    
    // Initialize alphas/abundances with Salmon's linear combination
    // alphas[i] = projectedCounts[i] * fracObserved + uniformPrior * (1 - fracObserved)
    double total_abundance = 0.0;
    for (size_t i = 0; i < state.n; ++i) {
        double alpha_init = projectedCounts[i] * fracObserved + uniformPrior * (1.0 - fracObserved);
        state.abundances[i] = alpha_init;
        total_abundance += alpha_init;
    }
    
    // Normalize to sum to 1
    if (total_abundance > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            state.abundances[i] /= total_abundance;
        }
    }
    
    // Set number of threads
    int num_threads = params.threads;
    if (num_threads == 0) {
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);
    
    // Allocate expected counts array
    std::vector<double> expected_counts(state.n, 0.0);
    
    // Thread-local storage for expected_counts: flat layout [num_threads * n_transcripts] for cache friendliness
    std::vector<double> expected_counts_tls(num_threads * state.n, 0.0);
    
    // Store previous expected_counts (alphas) for convergence checking
    std::vector<double> prev_expected_counts(state.n, 0.0);
    // Initialize prev_expected_counts from projected counts (like Salmon's alphas)
    for (size_t i = 0; i < state.n; ++i) {
        prev_expected_counts[i] = projectedCounts[i] * fracObserved + uniformPrior * (1.0 - fracObserved);
    }
    
    // Compute initial log-likelihood (for logging only, not used for convergence)
    double initial_ll = compute_log_likelihood(ecs, state.abundances.data(), state.eff_lengths.data());
    result.final_ll = initial_ll;
    
    // EM iterations
    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        // E-step: compute expected counts using alpha * aux (Salmon's exact approach)
        // Clear thread-local buffers (reuse per iteration to avoid reallocation)
        std::fill(expected_counts_tls.begin(), expected_counts_tls.end(), 0.0);
        std::memset(expected_counts.data(), 0, state.n * sizeof(double));
        
        // Parallel EC loop: each thread writes to its own TLS buffer (no atomics)
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t ec_idx = 0; ec_idx < ecs.ecs.size(); ++ec_idx) {
            int thread_id = omp_get_thread_num();
            double* thread_counts = expected_counts_tls.data() + thread_id * state.n;
            const EC& ec = ecs.ecs[ec_idx];
            size_t groupSize = ec.transcript_ids.size();
            
            // Single-transcript fast path (Salmon behavior: full count)
            if (groupSize == 1) {
                uint32_t tid = ec.transcript_ids[0];
                thread_counts[tid] += ec.count;  // No atomic - each thread has its own buffer
                continue;
            }
            
            // Multi-transcript: compute denominator using alpha * aux
            // Salmon-style: aux = weight * (1/effLen) where weight is alignment likelihood
            // Salmon's combinedWeights = alignmentWeight * probStartPos where probStartPos = 1/effLen
            double denom = 0.0;
            for (size_t i = 0; i < groupSize; ++i) {
                uint32_t tid = ec.transcript_ids[i];
                // Compute aux: weight * (1/effLen), always include 1/effLen
                double effLen_factor = (state.eff_lengths[tid] > 0 ? 1.0 / state.eff_lengths[tid] : 0.0);
                double aux = ec.has_weights()
                    ? ec.weights[i] * effLen_factor  // Multiply weight by 1/effLen
                    : effLen_factor;
                denom += state.abundances[tid] * aux;
            }
            
            if (denom > 0) {
                double invDenom = ec.count / denom;
                for (size_t i = 0; i < groupSize; ++i) {
                    uint32_t tid = ec.transcript_ids[i];
                    // Compute aux: same as in denom (weight * 1/effLen)
                    double effLen_factor = (state.eff_lengths[tid] > 0 ? 1.0 / state.eff_lengths[tid] : 0.0);
                    double aux = ec.has_weights()
                        ? ec.weights[i] * effLen_factor
                        : effLen_factor;
                    double contribution = state.abundances[tid] * aux * invDenom;
                    thread_counts[tid] += contribution;  // No atomic - each thread has its own buffer
                }
            }
        }
        
        // Deterministic reduction: sum across thread-local buffers into expected_counts
        // Use fixed thread order (0..num_threads-1) for determinism
        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (size_t i = 0; i < state.n; ++i) {
            double sum = 0.0;
            for (int t = 0; t < num_threads; ++t) {
                sum += expected_counts_tls[t * state.n + i];
            }
            expected_counts[i] = sum;
        }
        
        // M-step: update abundances (abundances are proportional to counts)
        double total_counts = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            total_counts += expected_counts[i];
        }
        
        if (total_counts > 0) {
            for (size_t i = 0; i < state.n; ++i) {
                state.abundances[i] = expected_counts[i] / total_counts;
            }
        }
        
        // Check convergence using expected_counts relative change
        // Only check transcripts with count above cutoff (like Salmon's alphaCheckCutoff)
        // Very low-count transcripts can oscillate with high relative diff but negligible absolute diff
        // Also enforce minIter before checking convergence (Salmon uses minIter=100)
        result.iterations = iter + 1;
        bool converged = false;
        double max_rel_diff = 0.0;
        const double countCheckCutoff = 1e-2;  // Same as Salmon's alphaCheckCutoff
        const uint32_t minIter = 100;  // Salmon's default minIter
        
        // Only check convergence after minIter iterations
        if (iter + 1 >= minIter) {
            converged = true;
            for (size_t i = 0; i < state.n; ++i) {
                if (expected_counts[i] > countCheckCutoff) {
                    // Compare current vs previous: |current - previous| / current
                    double rel_diff = std::abs(expected_counts[i] - prev_expected_counts[i]) / expected_counts[i];
                    max_rel_diff = std::max(max_rel_diff, rel_diff);
                    if (rel_diff > params.tolerance) {
                        converged = false;
                    }
                }
            }
        }
        
        // Update prev_expected_counts for next iteration
        for (size_t i = 0; i < state.n; ++i) {
            prev_expected_counts[i] = expected_counts[i];
        }
        
        if (converged) {
            result.converged = true;
            // Compute final log-likelihood for logging
            result.final_ll = compute_log_likelihood(ecs, state.abundances.data(), state.eff_lengths.data());
            break;
        }
    }
    
    // Compute final log-likelihood if didn't converge
    if (!result.converged) {
        result.final_ll = compute_log_likelihood(ecs, state.abundances.data(), state.eff_lengths.data());
    }
    
    // Copy final expected counts (NumReads in Salmon terminology)
    for (size_t i = 0; i < state.n; ++i) {
        result.counts[i] = expected_counts[i];
    }
    
    // Apply Salmon-style truncation (minAlpha = 1e-8)
    // Zero out transcripts with counts <= minAlpha and renormalize abundances
    const double minAlpha = params.zero_threshold;  // Default is 1e-8 (Salmon's minAlpha)
    if (minAlpha > 0) {
        size_t n_zeros = 0;
        double total_counts = 0.0;
        
        for (size_t i = 0; i < state.n; ++i) {
            if (result.counts[i] <= minAlpha) {
                result.counts[i] = 0.0;
                state.abundances[i] = 0.0;
                n_zeros++;
            } else {
                total_counts += result.counts[i];
            }
        }
        
        // Renormalize abundances for non-zero transcripts (Salmon's truncateCountVector behavior)
        if (total_counts > 0 && n_zeros > 0) {
            for (size_t i = 0; i < state.n; ++i) {
                if (result.counts[i] > 0) {
                    state.abundances[i] = result.counts[i] / total_counts;
                }
            }
        }
    }
    
    // Compute TPM the Salmon way: TPM_i = (count_i / eff_length_i) / sum_j(count_j / eff_length_j) * 1e6
    double total_normalized = 0.0;
    for (size_t i = 0; i < state.n; ++i) {
        if (state.eff_lengths[i] > 0) {
            total_normalized += result.counts[i] / state.eff_lengths[i];
        }
    }
    
    if (total_normalized > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            if (state.eff_lengths[i] > 0) {
                result.tpm[i] = (result.counts[i] / state.eff_lengths[i]) / total_normalized * 1e6;
            } else {
                result.tpm[i] = 0.0;
            }
        }
    }
    
    return result;
}
