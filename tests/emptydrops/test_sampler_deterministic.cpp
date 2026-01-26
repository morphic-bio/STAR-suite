/**
 * @file test_sampler_deterministic.cpp
 * @brief Unit test: EmptyDropsCRSampler determinism
 * 
 * Target: EmptyDropsCRSampler::montecarloPval
 * Test: Same seed produces identical results (deterministic behavior)
 */

#include <iostream>
#include <vector>
#include <cstdint>

#include "EmptyDropsCRSampler.h"

using namespace std;

int main() {
    cout << "Test: Sampler determinism (EmptyDropsCRSampler::montecarloPval)" << endl;
    
    // Tiny test data
    vector<uint32_t> totalval = {100, 200, 300};
    vector<uint32_t> totallen = {2, 3, 2};  // 7 total candidates
    vector<double> prob = {-50.0, -55.0, -60.0, -65.0, -70.0, -75.0, -80.0};
    
    // Uniform ambient profile (100 genes)
    vector<double> ambient(100, 0.01);
    
    uint32_t iterations = 1000;
    uint64_t seed = 1;  // Default seed per emptydrops_refactor_plan.md
    
    // Run twice with same seed
    vector<uint32_t> result1 = EmptyDropsCRSampler::montecarloPval(
        totalval, totallen, prob, ambient, iterations, seed, 0);
    
    vector<uint32_t> result2 = EmptyDropsCRSampler::montecarloPval(
        totalval, totallen, prob, ambient, iterations, seed, 0);
    
    // Verify identical results
    if (result1.size() != result2.size()) {
        cerr << "FAIL: Result sizes differ (" << result1.size() << " vs " << result2.size() << ")" << endl;
        return 1;
    }
    
    // Verify expected size (sum of totallen)
    size_t expectedSize = 0;
    for (auto len : totallen) expectedSize += len;
    if (result1.size() != expectedSize) {
        cerr << "FAIL: Result size " << result1.size() << " != expected " << expectedSize << endl;
        return 1;
    }
    
    // Verify identical values
    for (size_t i = 0; i < result1.size(); i++) {
        if (result1[i] != result2[i]) {
            cerr << "FAIL: Results differ at index " << i << " (" << result1[i] << " vs " << result2[i] << ")" << endl;
            return 1;
        }
    }
    
    cout << "PASS: Sampler produces identical results with seed=" << seed << " (size=" << result1.size() << ")" << endl;
    return 0;
}
