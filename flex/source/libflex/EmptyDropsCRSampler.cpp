#include "EmptyDropsCRSampler.h"
#include "pcg_random.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <random>
#include <thread>

using namespace std;

vector<uint32_t> EmptyDropsCRSampler::montecarloPval(
    const vector<uint32_t>& totalval,
    const vector<uint32_t>& totallen,
    const vector<double>& prob,
    const vector<double>& ambient,
    uint32_t iterations,
    uint64_t seed,
    uint32_t numThreads)
{
    // Validate inputs
    if (totalval.size() != totallen.size()) {
        throw runtime_error("length of run value/length vectors are not the same");
    }
    
    const size_t nvalues = accumulate(totallen.begin(), totallen.end(), 0ULL);
    if (nvalues != prob.size()) {
        throw runtime_error("sum of run lengths does not equal length of prob vector");
    }
    
    if (iterations == 0) {
        return vector<uint32_t>(nvalues, 0);
    }
    
    const size_t ngenes = ambient.size();
    if (ngenes == 0) {
        return vector<uint32_t>(nvalues, 0);
    }
    
    // Precompute log probabilities for ambient (multinomial path, alpha = Inf)
    vector<double> logprob;
    logprob.reserve(ngenes);
    for (double a : ambient) {
        logprob.push_back(log(a));
    }
    
    // Determine actual thread count
    uint32_t actualThreads = (numThreads == 0) ? 1 : numThreads;
    actualThreads = min(actualThreads, iterations);
    
    // Precompute cumulative probabilities for fast discrete sampling
    // This avoids creating discrete_distribution inside each iteration
    vector<double> cumProb(ngenes);
    double cumSum = 0.0;
    for (size_t i = 0; i < ngenes; i++) {
        cumSum += ambient[i];
        cumProb[i] = cumSum;
    }
    // Normalize
    if (cumSum > 0.0) {
        for (size_t i = 0; i < ngenes; i++) {
            cumProb[i] /= cumSum;
        }
    }
    
    // Single-threaded path (original algorithm)
    if (actualThreads <= 1) {
        vector<uint32_t> above(nvalues, 0);
        vector<uint32_t> tracker(ngenes, 0);
        
        for (uint32_t it = 0; it < iterations; ++it) {
            uint64_t iterSeed = seed ^ static_cast<uint64_t>(it);
            pcg32 generator(iterSeed, it);
            
            fill(tracker.begin(), tracker.end(), 0);
            uint32_t curtotal = 0;
            double curp = 0.0;
            
            auto abIt = above.begin();
            auto tvIt = totalval.begin();
            auto tlIt = totallen.begin();
            auto pIt = prob.begin();
            
            while (tvIt != totalval.end()) {
                const uint32_t curlen = *tlIt;
                const uint32_t curval = *tvIt;
                
                while (curtotal < curval) {
                    // Fast discrete sampling using binary search on cumulative probs
                    double u = ldexp(static_cast<double>(generator()), -32);
                    auto it = lower_bound(cumProb.begin(), cumProb.end(), u);
                    size_t chosen = (it != cumProb.end()) ? (it - cumProb.begin()) : (ngenes - 1);
                    
                    uint32_t& curnum = tracker[chosen];
                    curp += logprob[chosen];
                    curp -= log(static_cast<double>(++curnum));
                    ++curtotal;
                }
                
                size_t higher = lower_bound(pIt, pIt + curlen, curp) - pIt;
                if (higher < curlen) {
                    ++(*(abIt + higher));
                }
                
                ++tlIt;
                ++tvIt;
                pIt += curlen;
                abIt += curlen;
            }
        }
        
        // Cumulative counts
        auto abIt2 = above.begin();
        for (uint32_t curlen : totallen) {
            for (uint32_t i = 1; i < curlen; ++i) {
                uint32_t prev = *abIt2;
                ++abIt2;
                *abIt2 += prev;
            }
            ++abIt2;
        }
        return above;
    }
    
    // Multi-threaded path: each thread has its own local_above vector
    vector<vector<uint32_t>> threadAbove(actualThreads);
    for (auto& v : threadAbove) {
        v.resize(nvalues, 0);
    }
    
    // Worker function for each thread
    auto worker = [&](uint32_t threadId, uint32_t startIter, uint32_t endIter) {
        vector<uint32_t>& local_above = threadAbove[threadId];
        vector<uint32_t> tracker(ngenes, 0);
        
        for (uint32_t it = startIter; it < endIter; ++it) {
            uint64_t iterSeed = seed ^ static_cast<uint64_t>(it);
            pcg32 generator(iterSeed, it);
            
            fill(tracker.begin(), tracker.end(), 0);
            uint32_t curtotal = 0;
            double curp = 0.0;
            
            auto abIt = local_above.begin();
            auto tvIt = totalval.begin();
            auto tlIt = totallen.begin();
            auto pIt = prob.begin();
            
            while (tvIt != totalval.end()) {
                const uint32_t curlen = *tlIt;
                const uint32_t curval = *tvIt;
                
                while (curtotal < curval) {
                    // Fast discrete sampling using binary search on cumulative probs
                    double u = ldexp(static_cast<double>(generator()), -32);
                    auto it = lower_bound(cumProb.begin(), cumProb.end(), u);
                    size_t chosen = (it != cumProb.end()) ? (it - cumProb.begin()) : (ngenes - 1);
                    
                    uint32_t& curnum = tracker[chosen];
                    curp += logprob[chosen];
                    curp -= log(static_cast<double>(++curnum));
                    ++curtotal;
                }
                
                size_t higher = lower_bound(pIt, pIt + curlen, curp) - pIt;
                if (higher < curlen) {
                    ++(*(abIt + higher));
                }
                
                ++tlIt;
                ++tvIt;
                pIt += curlen;
                abIt += curlen;
            }
        }
    };
    
    // Launch threads with balanced iteration ranges
    vector<thread> threads;
    uint32_t iterPerThread = iterations / actualThreads;
    uint32_t remainder = iterations % actualThreads;
    uint32_t startIter = 0;
    
    for (uint32_t t = 0; t < actualThreads; ++t) {
        uint32_t count = iterPerThread + (t < remainder ? 1 : 0);
        uint32_t endIter = startIter + count;
        threads.emplace_back(worker, t, startIter, endIter);
        startIter = endIter;
    }
    
    // Wait for all threads
    for (auto& th : threads) {
        th.join();
    }
    
    // Merge thread-local results
    vector<uint32_t> above(nvalues, 0);
    for (const auto& local : threadAbove) {
        for (size_t i = 0; i < nvalues; ++i) {
            above[i] += local[i];
        }
    }
    
    // Cumulative counts
    auto abIt = above.begin();
    for (uint32_t curlen : totallen) {
        for (uint32_t i = 1; i < curlen; ++i) {
            uint32_t prev = *abIt;
            ++abIt;
            *abIt += prev;
        }
        ++abIt;
    }
    
    return above;
}
