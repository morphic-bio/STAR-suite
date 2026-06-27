#include "vb_engine.h"
#include "em_engine.h"
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <unordered_map>

// A standalone, high-precision implementation of Digamma (Psi)
// Logic: Uses the recurrence relation psi(x) = psi(x+1) - 1/x to shift 
// small arguments to a stable region (x >= 6), then uses asymptotic expansion.
double robust_digamma(double x) {
    double result = 0.0;

    // 1. Handle the "Dirichlet Cliff" (x near 0) and negative numbers
    // The recurrence relation shifts us to the stable region.
    while (x < 6.0) {
        // If x is extremely close to zero, 1/x will be massive.
        // We rely on double precision to handle this naturally.
        if (x == 0.0) return -std::numeric_limits<double>::infinity(); // Singularity
        result -= 1.0 / x;
        x += 1.0;
    }

    // 2. Stable Region (x >= 6)
    // Use the asymptotic expansion:
    // psi(x) ~ ln(x) - 1/2x - 1/12x^2 + 1/120x^4 - 1/252x^6 ...
    
    // We compute 1/x^2 once to reuse it
    double inv_x = 1.0 / x;
    double inv_x2 = inv_x * inv_x;

    result += std::log(x);
    result -= 0.5 * inv_x;

    // Bernoulli number terms
    double term = inv_x2;       // 1/x^2
    result -= (1.0/12.0) * term; 
    
    term *= inv_x2;             // 1/x^4
    result += (1.0/120.0) * term; 
    
    term *= inv_x2;             // 1/x^6
    result -= (1.0/252.0) * term; 
    
    // Usually sufficient precision for double (term is now ~1e-15 for x=6)
    // Add one more term if strict parity with Boost is required:
    // term *= inv_x2;             // 1/x^8
    // result += (1.0/240.0) * term; 

    return result;
}

// Constants matching Salmon's implementation
constexpr double digammaMin = 1e-10;
constexpr double minEQClassWeight = std::numeric_limits<double>::min();
// Note: minAlpha and alphaCheckCutoff are now configurable via EMParams
// Defaults: minAlpha = zero_threshold = 1e-8, alphaCheckCutoff = 1e-2

// Compute unique EC counts per transcript (for fallback initialization)
std::vector<double> compute_unique_counts(const ECTable& ecs) {
    std::vector<double> unique_counts(ecs.n_transcripts, 0.0);
    
    for (const EC& ec : ecs.ecs) {
        if (ec.transcript_ids.size() == 1 && ec.count > 0) {
            // Single-transcript EC = unique evidence
            unique_counts[ec.transcript_ids[0]] += ec.count;
        }
    }
    
    return unique_counts;
}

namespace {

struct UnionFind {
    explicit UnionFind(size_t n) : parent(n), rank(n, 0) {
        for (size_t i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    size_t find(size_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    }

    void unite(size_t a, size_t b) {
        size_t ra = find(a);
        size_t rb = find(b);
        if (ra == rb) {
            return;
        }
        if (rank[ra] < rank[rb]) {
            parent[ra] = rb;
        } else if (rank[ra] > rank[rb]) {
            parent[rb] = ra;
        } else {
            parent[rb] = ra;
            ++rank[ra];
        }
    }

    std::vector<size_t> parent;
    std::vector<uint8_t> rank;
};

std::vector<uint32_t> distinct_transcript_ids(const EC& ec) {
    std::vector<uint32_t> ids = ec.transcript_ids;
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    return ids;
}

void project_cluster_to_polytope(
    const std::vector<uint32_t>& members,
    double cluster_count,
    const std::vector<double>& unique_counts,
    const std::vector<double>& total_counts,
    std::vector<double>& projected_counts
) {
    if (members.empty() || cluster_count <= 0.0) {
        return;
    }

    std::vector<bool> bound(members.size(), false);
    for (size_t round = 0; round <= 5000; ++round) {
        double unbound_counts = 0.0;
        double bound_counts = 0.0;

        for (size_t i = 0; i < members.size(); ++i) {
            const uint32_t tid = members[i];
            if (projected_counts[tid] > total_counts[tid]) {
                projected_counts[tid] = total_counts[tid];
                bound[i] = true;
            } else if (projected_counts[tid] < unique_counts[tid]) {
                projected_counts[tid] = unique_counts[tid];
                bound[i] = true;
            }

            if (bound[i]) {
                bound_counts += projected_counts[tid];
            } else {
                unbound_counts += projected_counts[tid];
            }
        }

        const double total = unbound_counts + bound_counts;
        if (std::abs(total - cluster_count) <= 1e-8 * std::max(1.0, cluster_count)) {
            return;
        }

        if (unbound_counts == 0.0) {
            std::fill(bound.begin(), bound.end(), false);
            unbound_counts = bound_counts;
            bound_counts = 0.0;
        }

        if (unbound_counts == 0.0) {
            return;
        }

        const double normalizer = (cluster_count - bound_counts) / unbound_counts;
        for (size_t i = 0; i < members.size(); ++i) {
            if (!bound[i]) {
                projected_counts[members[i]] *= normalizer;
            }
        }
    }
}

} // namespace

// Compute projected counts from ALL ECs for Salmon-style VB initialization.
// Salmon derives these from online transcript masses, then projects each
// transcript cluster into [uniqueCount, totalCount] bounds.  Reconstruct the
// same cluster-level initialization from the finalized EC table: EC weights
// approximate the online mass split, while distinct transcript IDs reproduce
// Salmon's per-fragment total/unique-count side effects.
std::vector<double> compute_projected_counts(const ECTable& ecs, const std::vector<double>& eff_lengths) {
    (void)eff_lengths;
    std::vector<double> projected_counts(ecs.n_transcripts, 0.0);
    std::vector<double> unique_counts(ecs.n_transcripts, 0.0);
    std::vector<double> total_counts(ecs.n_transcripts, 0.0);
    std::vector<double> transcript_mass(ecs.n_transcripts, 0.0);
    UnionFind clusters(ecs.n_transcripts);

    std::vector<std::vector<uint32_t>> ec_distinct_ids;
    ec_distinct_ids.reserve(ecs.ecs.size());

    for (const EC& ec : ecs.ecs) {
        std::vector<uint32_t> ids = distinct_transcript_ids(ec);
        ec_distinct_ids.push_back(ids);
        if (ec.count <= 0.0 || ids.empty()) {
            continue;
        }

        for (uint32_t tid : ids) {
            if (tid < total_counts.size()) {
                total_counts[tid] += ec.count;
            }
        }

        if (ids.size() == 1) {
            unique_counts[ids[0]] += ec.count;
        } else {
            for (size_t i = 1; i < ids.size(); ++i) {
                clusters.unite(ids[0], ids[i]);
            }
        }

        if (ec.has_weights() && ec.weights.size() == ec.transcript_ids.size()) {
            for (size_t i = 0; i < ec.transcript_ids.size(); ++i) {
                const uint32_t tid = ec.transcript_ids[i];
                if (tid < transcript_mass.size() && std::isfinite(ec.weights[i]) && ec.weights[i] > 0.0) {
                    transcript_mass[tid] += ec.count * ec.weights[i];
                }
            }
        } else {
            const double share = ec.count / static_cast<double>(ids.size());
            for (uint32_t tid : ids) {
                if (tid < transcript_mass.size()) {
                    transcript_mass[tid] += share;
                }
            }
        }
    }

    std::unordered_map<size_t, std::vector<uint32_t>> cluster_members;
    cluster_members.reserve(ecs.n_transcripts);
    for (uint32_t tid = 0; tid < ecs.n_transcripts; ++tid) {
        cluster_members[clusters.find(tid)].push_back(tid);
    }

    std::unordered_map<size_t, double> cluster_counts;
    cluster_counts.reserve(cluster_members.size());
    for (size_t ec_idx = 0; ec_idx < ecs.ecs.size(); ++ec_idx) {
        const EC& ec = ecs.ecs[ec_idx];
        const std::vector<uint32_t>& ids = ec_distinct_ids[ec_idx];
        if (ec.count <= 0.0 || ids.empty()) {
            continue;
        }
        cluster_counts[clusters.find(ids[0])] += ec.count;
    }

    for (const auto& kv : cluster_members) {
        const size_t rep = kv.first;
        const std::vector<uint32_t>& members = kv.second;
        const double cluster_count = cluster_counts[rep];
        if (cluster_count <= 0.0) {
            continue;
        }

        double cluster_mass = 0.0;
        for (uint32_t tid : members) {
            cluster_mass += transcript_mass[tid];
        }

        if (cluster_mass > 0.0) {
            const double scale = cluster_count / cluster_mass;
            for (uint32_t tid : members) {
                projected_counts[tid] = transcript_mass[tid] * scale;
            }
        } else {
            const double share = cluster_count / static_cast<double>(members.size());
            for (uint32_t tid : members) {
                projected_counts[tid] = share;
            }
        }

        if (members.size() > 1) {
            project_cluster_to_polytope(members, cluster_count, unique_counts, total_counts, projected_counts);
        }
    }

    return projected_counts;
}

// Check if transcript has unique evidence (appears in single-transcript ECs)
std::vector<bool> compute_unique_evidence(const ECTable& ecs) {
    std::vector<bool> has_unique(ecs.n_transcripts, false);
    
    for (const EC& ec : ecs.ecs) {
        if (ec.transcript_ids.size() == 1 && ec.count > 0) {
            has_unique[ec.transcript_ids[0]] = true;
        }
    }
    
    return has_unique;
}

double compute_elbo(const ECTable& ecs, const double* abundances, const double* eff_lengths, double vb_prior) {
    double ll = compute_log_likelihood(ecs, abundances, eff_lengths);
    
    // Add Dirichlet prior term: sum((alpha-1) * log(theta))
    // For uniform Dirichlet: alpha_i = vb_prior for all i
    double prior_term = 0.0;
    for (size_t i = 0; i < ecs.n_transcripts; ++i) {
        if (abundances[i] > 0) {
            prior_term += (vb_prior - 1.0) * std::log(abundances[i]);
        }
    }
    
    return ll + prior_term;
}

EMResult run_vb(const ECTable& ecs, TranscriptState& state, const EMParams& params) {
    EMResult result;
    result.counts.resize(state.n);
    result.tpm.resize(state.n);
    
    // Set number of threads
    int num_threads = params.threads;
    if (num_threads == 0) {
        num_threads = omp_get_max_threads();
    }
    omp_set_num_threads(num_threads);
    
    // Initialize priorAlphas based on mode (per-transcript vs per-nucleotide)
    std::vector<double> priorAlphas(state.n, params.vb_prior);
    if (!params.per_transcript_prior) {
        for (size_t i = 0; i < state.n; ++i) {
            priorAlphas[i] = params.vb_prior * state.eff_lengths[i];
        }
    }
    
    // Salmon-style VB initialization using projected counts
    // This distributes multi-mapper counts using EC weights to break symmetry
    std::vector<double> projectedCounts = compute_projected_counts(ecs, state.eff_lengths);
    
    // Compute total weight (sum of projected counts)
    // Salmon uses ALL transcripts for numActive, not just those with projected counts > 0
    double totalWeight = 0.0;
    size_t numActive = state.n;  // Use all transcripts like Salmon
    for (size_t i = 0; i < state.n; ++i) {
        totalWeight += projectedCounts[i];
    }
    
    // Salmon's fracObserved calculation: min(0.999, totalWeight / numRequiredFragments)
    double numRequiredFragments = params.num_required_fragments;
    if (numRequiredFragments <= 0) {
        numRequiredFragments = 5e6;  // Salmon default
    }
    double maxFrac = 0.999;  // Salmon uses 0.999 as max
    double fracObserved = std::min(maxFrac, totalWeight / numRequiredFragments);
    
    // Salmon's uniformPrior: totalWeight / numActive
    double uniformPrior = (numActive > 0) ? totalWeight / numActive : 0.0;
    
    // Initialize alpha: Salmon's formula (NO prior added here!)
    // alpha[i] = projectedCounts[i] * fracObserved + uniformPrior * (1 - fracObserved)
    // The prior is added only during digamma calculation, not stored in alpha
    std::vector<double> alpha(state.n);
    for (size_t i = 0; i < state.n; ++i) {
        alpha[i] = projectedCounts[i] * fracObserved + uniformPrior * (1.0 - fracObserved);
    }
    
    // Initialize abundances from alpha (normalized for E-step)
    double total_alpha = 0.0;
    for (size_t i = 0; i < state.n; ++i) {
        total_alpha += alpha[i];
    }
    if (total_alpha > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            state.abundances[i] = alpha[i] / total_alpha;
        }
    } else {
        // Fallback to uniform if all zero
        double uniform = 1.0 / state.n;
        for (size_t i = 0; i < state.n; ++i) {
            state.abundances[i] = uniform;
            alpha[i] = params.vb_prior; // Reset alpha to prior
        }
    }
    
    // Allocate expected counts array
    std::vector<double> expected_counts(state.n, 0.0);
    
    // Debug instrumentation: track selected transcripts
    std::vector<size_t> debug_indices;
    std::ofstream debug_stream;
    if (params.debug_trace && !params.debug_file.empty()) {
        debug_stream.open(params.debug_file);
        if (debug_stream.is_open()) {
            debug_stream << "iter\ttranscript\talpha\tlogNorm\texpTheta\texpected_count\n";
            debug_stream << "# EC-level trace: EC\titer\tec_id\ttranscript\tdenom\texpTheta\taux\tcontribution\n";
            // Find indices for debug transcripts
            for (const std::string& txp_name : params.debug_transcripts) {
                for (size_t i = 0; i < state.n; ++i) {
                    if (state.names[i] == txp_name) {
                        debug_indices.push_back(i);
                        break;
                    }
                }
            }
            // Log initial values (iter -1) - BEFORE first iteration
            if (!debug_indices.empty()) {
                // Compute logNorm for debug output (add priors like Salmon does)
                double alphaSum = 0.0;
                for (size_t i = 0; i < state.n; ++i) {
                    alphaSum += alpha[i] + priorAlphas[i];  // Add prior for digamma
                }
                double logNorm = robust_digamma(alphaSum);
                
                for (size_t idx : debug_indices) {
                    double ap = alpha[idx] + priorAlphas[idx];  // Add prior for digamma
                    double expTheta_val = (ap > digammaMin) ? std::exp(robust_digamma(ap) - logNorm) : 0.0;
                    debug_stream << "-1\t" << state.names[idx] << "\t" 
                                << ap << "\t" << logNorm << "\t"
                                << expTheta_val << "\t0.0\n";
                    debug_stream.flush(); // Ensure it's written
                }
            }
        }
    }
    
    // Store previous alpha for convergence checking
    std::vector<double> prev_alpha(state.n);
    for (size_t i = 0; i < state.n; ++i) {
        prev_alpha[i] = alpha[i];
    }
    
    // NOTE: VB/EM runs after alignment is complete, so all threads are available.
    // This avoids nested OpenMP parallelism (alignment uses threads, VB uses threads separately).
    
    // Thread-local storage for expected_counts: flat layout [num_threads * n_transcripts] for cache friendliness
    // Each thread writes to its own buffer, then we reduce deterministically
    std::vector<double> expected_counts_tls(num_threads * state.n, 0.0);
    bool effective_lengths_updated = false;
    
    // VB iterations
    for (uint32_t iter = 0; iter < params.max_iters; ++iter) {
        if (!effective_lengths_updated && params.effective_length_update &&
            iter > params.effective_length_update_target_iter) {
            params.effective_length_update(iter, state, alpha);
            effective_lengths_updated = true;
            if (!params.per_transcript_prior) {
                for (size_t i = 0; i < state.n; ++i) {
                    priorAlphas[i] = params.vb_prior * state.eff_lengths[i];
                }
            }
        }

        // Compute logNorm = digamma(sum of all alphas + priors) ONCE before E-step (Salmon's approach)
        // Salmon adds priors only for digamma calculation, not stored in alpha
        double alphaSum = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            alphaSum += alpha[i] + priorAlphas[i];  // Add prior for digamma
        }
        double logNorm = robust_digamma(alphaSum);
        
        // Precompute expTheta vector with digammaMin guard (Salmon's approach)
        // ap = alpha + prior for digamma computation
        std::vector<double> expTheta(state.n, 0.0);
        for (size_t i = 0; i < state.n; ++i) {
            double ap = alpha[i] + priorAlphas[i];  // Add prior for digamma
            if (ap > digammaMin) {
                expTheta[i] = std::exp(robust_digamma(ap) - logNorm);
            } else {
                expTheta[i] = 0.0;  // Zero out very small alphas
            }
        }
        
        // E-step: compute expected counts using expTheta * aux (Salmon's exact approach)
        // Clear thread-local buffers (reuse per iteration to avoid reallocation)
        std::fill(expected_counts_tls.begin(), expected_counts_tls.end(), 0.0);
        std::memset(expected_counts.data(), 0, state.n * sizeof(double));
        
        // Check if we need to log EC details for debug transcripts
        bool log_ec_details = params.debug_trace && !params.debug_file.empty() && !debug_indices.empty();
        
        // Parallel EC loop: each thread writes to its own TLS buffer (no atomics)
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t ec_idx = 0; ec_idx < ecs.ecs.size(); ++ec_idx) {
            int thread_id = omp_get_thread_num();
            double* thread_counts = expected_counts_tls.data() + thread_id * state.n;
            const EC& ec = ecs.ecs[ec_idx];
            size_t groupSize = ec.transcript_ids.size();
            
            // Check if this EC contains any debug transcripts
            bool ec_contains_debug = false;
            if (log_ec_details) {
                for (size_t i = 0; i < groupSize; ++i) {
                    uint32_t tid = ec.transcript_ids[i];
                    for (size_t debug_idx : debug_indices) {
                        if (tid == debug_idx) {
                            ec_contains_debug = true;
                            break;
                        }
                    }
                    if (ec_contains_debug) break;
                }
            }
            
            // Single-transcript fast path (Salmon behavior: full count, no expTheta guard)
            if (groupSize == 1) {
                uint32_t tid = ec.transcript_ids[0];
                thread_counts[tid] += ec.count;  // No atomic - each thread has its own buffer
                
                // Log single-transcript EC if it's a debug transcript
                if (log_ec_details && ec_contains_debug) {
                    #pragma omp critical
                    {
                        debug_stream << "EC\t" << iter << "\t" << ec_idx << "\t" << state.names[tid] 
                                    << "\tSINGLE\t" << ec.count << "\t" << ec.count << "\t" << ec.count << "\n";
                        debug_stream.flush();
                    }
                }
                continue;
            }
            
            // Multi-transcript: compute denominator using expTheta * aux
            double denom = 0.0;
            for (size_t i = 0; i < groupSize; ++i) {
                uint32_t tid = ec.transcript_ids[i];
                if (expTheta[tid] > 0.0) {
                    // Use expTheta * aux (Salmon's VBEMOptimizer approach)
                    // aux = weight * (1/effLen), where weight is alignment likelihood
                    // Salmon's combinedWeights = weight * probStartPos where probStartPos = 1/effLen
                    double effLen_factor = (state.eff_lengths[tid] > 0 ? 1.0 / state.eff_lengths[tid] : 0.0);
                    double aux = ec.has_weights() 
                        ? ec.weights[i] * effLen_factor  // Multiply weight by 1/effLen
                        : effLen_factor;
                    denom += expTheta[tid] * aux;
                }
            }
            
            // Skip EC if denom too small (minEQClassWeight guard)
            if (denom <= minEQClassWeight) {
                continue;
            }
            
            // Distribute counts proportionally
            double invDenom = ec.count / denom;
            for (size_t i = 0; i < groupSize; ++i) {
                uint32_t tid = ec.transcript_ids[i];
                if (expTheta[tid] > 0.0) {
                    // aux = weight * (1/effLen), same as in denom calculation
                    double effLen_factor = (state.eff_lengths[tid] > 0 ? 1.0 / state.eff_lengths[tid] : 0.0);
                    double aux = ec.has_weights()
                        ? ec.weights[i] * effLen_factor
                        : effLen_factor;
                    double contribution = expTheta[tid] * aux * invDenom;
                    thread_counts[tid] += contribution;  // No atomic - each thread has its own buffer
                    
                    // Log per-transcript contribution if this is a debug transcript
                    if (log_ec_details && ec_contains_debug) {
                        bool is_debug_txp = false;
                        for (size_t debug_idx : debug_indices) {
                            if (tid == debug_idx) {
                                is_debug_txp = true;
                                break;
                            }
                        }
                        if (is_debug_txp) {
                            #pragma omp critical
                            {
                                debug_stream << "EC\t" << iter << "\t" << ec_idx << "\t" << state.names[tid]
                                            << "\t" << denom << "\t" << expTheta[tid] << "\t" << aux 
                                            << "\t" << contribution << "\n";
                                debug_stream.flush();
                            }
                        }
                    }
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
        
        // VB M-step: update alpha = new expected counts (NO prior added!)
        // Prior is only added during digamma calculation in E-step (Salmon's approach)
        for (size_t i = 0; i < state.n; ++i) {
            alpha[i] = expected_counts[i];
        }
        
        // Normalize abundances from alpha (for ELBO computation and next iteration)
        double total_alpha = 0.0;
        for (size_t i = 0; i < state.n; ++i) {
            total_alpha += alpha[i];
        }
        
        if (total_alpha > 0) {
            for (size_t i = 0; i < state.n; ++i) {
                state.abundances[i] = alpha[i] / total_alpha;
            }
        }
        
        // Debug instrumentation: log alpha/weights for selected transcripts
        // Log AFTER M-step (alpha updated) but BEFORE convergence check
        if (debug_stream.is_open()) {
            // Recompute logNorm and expTheta for debug output (add priors like Salmon)
            double debug_alphaSum = 0.0;
            for (size_t i = 0; i < state.n; ++i) {
                debug_alphaSum += alpha[i] + priorAlphas[i];  // Add prior
            }
            double debug_logNorm = robust_digamma(debug_alphaSum);
            
            for (size_t idx : debug_indices) {
                double ap = alpha[idx] + priorAlphas[idx];  // Add prior
                double expTheta_val = (ap > digammaMin) ? std::exp(robust_digamma(ap) - debug_logNorm) : 0.0;
                debug_stream << iter << "\t" << state.names[idx] << "\t"
                            << ap << "\t" << debug_logNorm << "\t"
                            << expTheta_val << "\t" << expected_counts[idx] << "\n";
                debug_stream.flush();
            }
        }
        
        // Check convergence using alpha relative difference
        // Only check after min_iters (Salmon behavior)
        bool converged = false;
        double maxRelDiff = 0.0;
        
        if (iter + 1 >= params.min_iters) {
            converged = true;
            double alphaCheckCutoff = params.alpha_check_cutoff;
            
            for (size_t i = 0; i < state.n; ++i) {
                // Salmon checks alpha + prior (called alphaPrime in their code)
                double alphaPrime = alpha[i];  // In Salmon, this is the expected count
                // Only check transcripts with alpha above cutoff (Salmon's alphaCheckCutoff)
                if (alphaPrime > alphaCheckCutoff) {
                    double relDiff = std::abs(alphaPrime - prev_alpha[i]) / alphaPrime;
                    maxRelDiff = std::max(maxRelDiff, relDiff);
                    if (relDiff > params.tolerance) {
                        converged = false;
                    }
                }
            }
        }
        
        // Update prev_alpha for next iteration
        for (size_t i = 0; i < state.n; ++i) {
            prev_alpha[i] = alpha[i];
        }
        
        // Compute ELBO only at end or for debug (not every iteration - expensive)
        // Store iteration count
        result.iterations = iter + 1;
        
        // Only compute ELBO if converged, at max iterations, or debug mode
        if (converged || iter + 1 >= params.max_iters || params.debug_trace) {
            double curr_elbo = compute_elbo(ecs, state.abundances.data(), state.eff_lengths.data(), params.vb_prior);
            result.final_ll = curr_elbo;
        }
        
        if (converged) {
            result.converged = true;
            break;
        }
    }
    
    if (debug_stream.is_open()) {
        debug_stream.close();
    }
    
    // Compute ELBO at end if not already computed (for final result)
    // ELBO is computed during iterations only if converged, at max_iters, or debug mode
    // Always compute at end to ensure final_ll is set
    if (result.final_ll == 0.0) {
        result.final_ll = compute_elbo(ecs, state.abundances.data(), state.eff_lengths.data(), params.vb_prior);
    }
    
    // Copy final expected counts: alpha is now just expected_counts (no prior subtraction needed)
    for (size_t i = 0; i < state.n; ++i) {
        result.counts[i] = alpha[i];
    }
    
    // Post-convergence truncation: zero out counts below zero_threshold (configurable via EMParams)
    double zero_threshold = params.zero_threshold;  // Use configurable threshold (default: 1e-8)
    for (size_t i = 0; i < state.n; ++i) {
        if (result.counts[i] <= zero_threshold) {
            result.counts[i] = 0.0;
        }
    }
    
    // Zero out transcripts where alpha_i <= epsilon AND no unique support
    // alpha now stores just expected_counts (no prior included)
    // This matches Salmon's implicit zeroing behavior
    std::vector<bool> has_unique = compute_unique_evidence(ecs);
    const double epsilon = 1e-8;
    double total_counts = 0.0;
    
    for (size_t i = 0; i < state.n; ++i) {
        if (alpha[i] <= epsilon && !has_unique[i]) {
            // Zero out: very small expected count AND no unique evidence
            result.counts[i] = 0.0;
            state.abundances[i] = 0.0;
        } else {
            total_counts += result.counts[i];
        }
    }
    
    // Renormalize abundances for non-zero transcripts
    if (total_counts > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            if (result.counts[i] > 0) {
                state.abundances[i] = result.counts[i] / total_counts;
            }
        }
    }
    
    // Compute TPM the Salmon way: TPM_i = (count_i / eff_length_i) / sum_j(count_j / eff_length_j) * 1e6
    double total_normalized = 0.0;
    for (size_t i = 0; i < state.n; ++i) {
        if (state.eff_lengths[i] > 0 && result.counts[i] > 0) {
            total_normalized += result.counts[i] / state.eff_lengths[i];
        }
    }
    
    if (total_normalized > 0) {
        for (size_t i = 0; i < state.n; ++i) {
            if (state.eff_lengths[i] > 0 && result.counts[i] > 0) {
                result.tpm[i] = (result.counts[i] / state.eff_lengths[i]) / total_normalized * 1e6;
            } else {
                result.tpm[i] = 0.0;
            }
        }
    }
    
    return result;
}
