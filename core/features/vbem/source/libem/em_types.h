#ifndef EM_TYPES_H
#define EM_TYPES_H

#include <vector>
#include <string>
#include <cstdint>
#include <cstddef>

// Structure of Arrays (SoA) for transcript state - cache/SIMD friendly
struct TranscriptState {
    std::vector<std::string> names;      // transcript IDs
    std::vector<double> lengths;         // raw transcript lengths
    std::vector<double> eff_lengths;     // effective lengths (= lengths for v1)
    std::vector<double> abundances;      // current abundance estimates
    std::vector<double> counts;          // expected counts (after EM)
    size_t n;                            // number of transcripts
    
    TranscriptState() : n(0) {}
    
    void resize(size_t size) {
        n = size;
        names.resize(size);
        lengths.resize(size);
        eff_lengths.resize(size);
        abundances.resize(size);
        counts.resize(size);
    }
};

// Array of Structures (AoS) for equivalence classes
struct EC {
    std::vector<uint32_t> transcript_ids;  // indices into TranscriptState (order preserved from file)
    std::vector<double> weights;           // per-transcript combinedWeights (auxs) from Salmon
                                           // Empty if file was generated without --dumpEqWeights
    double count;                          // observed fragment count
    
    EC() : count(0.0) {}
    
    // Check if this EC has explicit weights
    bool has_weights() const { return !weights.empty(); }
};

struct ECTable {
    std::vector<EC> ecs;
    size_t n_ecs;
    size_t n_transcripts;
    
    ECTable() : n_ecs(0), n_transcripts(0) {}
};

// EM algorithm parameters
struct EMParams {
    uint32_t max_iters = 10000;       // Maximum iterations
    double tolerance = 0.01;          // Relative change tolerance for convergence
    double vb_prior = 0.01;           // Dirichlet concentration
    bool use_vb = false;              // --vb flag
    bool init_by_length = false;      // weight initial abundances by length
    int threads = 0;                  // 0 = OMP default
    double zero_threshold = 1e-8;     // Truncation threshold
    bool debug_trace = false;         // enable debug tracing for selected transcripts
    std::vector<std::string> debug_transcripts; // transcript IDs to trace
    std::string debug_file;           // output file for debug traces
    bool per_transcript_prior = true; // default: per-transcript
                                      // If false: prior_i = vb_prior * eff_length_i (per-nucleotide)
    
    // VB-specific parameters (Salmon-style initialization)
    uint32_t min_iters = 100;         // Minimum iterations before checking convergence (VB only)
    double num_required_fragments = 5e6; // For VB init: fracObserved = min(1, totalWeight/numRequired)
    double alpha_check_cutoff = 1e-2; // Only check convergence for transcripts with alpha > cutoff
};

// EM algorithm results
struct EMResult {
    std::vector<double> counts;       // estimated counts per transcript
    std::vector<double> tpm;          // TPM values
    double final_ll;                  // final log-likelihood (or ELBO for VB)
    uint32_t iterations;              // iterations until convergence
    bool converged;
    
    EMResult() : final_ll(0.0), iterations(0), converged(false) {}
};

#endif // EM_TYPES_H
