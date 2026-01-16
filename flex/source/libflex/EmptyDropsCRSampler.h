#ifndef H_EmptyDropsCRSampler
#define H_EmptyDropsCRSampler

#include "IncludeDefine.h"
#include <vector>

// Forward declarations
class Parameters;

/**
 * DropletUtils-style nested multinomial sampler for EmptyDrops.
 * 
 * Ports the logic from DropletUtils montecarlo_pval.cpp:
 * - Uses run-length encoded totals/logprobs (sorted by total ascending)
 * - Per-UMI discrete draws, accumulating as totals increase
 * - PCG32 RNG with seed ^ iteration per simulation
 * - Returns counts of sims with logprob <= observed
 */
class EmptyDropsCRSampler {
public:
    /**
     * Run Monte Carlo simulations using DropletUtils-style nested multinomial sampling.
     * 
     * @param totalval Distinct totals (ascending order)
     * @param totallen Run lengths (number of candidates with each total)
     * @param prob Observed log probabilities (in order matching run-length encoding)
     * @param ambient Dense ambient probability vector (normalized, all genes)
     * @param iterations Number of Monte Carlo simulations
     * @param seed Base seed (will be combined with iteration: seed ^ iteration)
     * @param numThreads Number of threads (0 = single-threaded, >0 = parallel)
     * @return Vector of counts: for each candidate, number of sims with logprob <= observed
     */
    static std::vector<uint32_t> montecarloPval(
        const std::vector<uint32_t>& totalval,
        const std::vector<uint32_t>& totallen,
        const std::vector<double>& prob,
        const std::vector<double>& ambient,
        uint32_t iterations,
        uint64_t seed,
        uint32_t numThreads = 0
    );
};

#endif
