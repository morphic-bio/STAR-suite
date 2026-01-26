/**
 * @file test_ambient_profile.cpp
 * @brief Unit test: EmptyDropsMultinomial ambient profile computation
 * 
 * Target: EmptyDropsMultinomial::computeAmbientProfile
 * Test: Verify sum≈1, nonzero genes tracked, stable output on tiny counts vector
 */

#include <iostream>
#include <vector>
#include <cstdint>
#include <cmath>

#include "EmptyDropsMultinomial.h"

using namespace std;

int main() {
    cout << "Test: Ambient profile sanity (EmptyDropsMultinomial::computeAmbientProfile)" << endl;
    
    // Create ambient counts (simulating gene expression)
    uint32_t nFeatures = 1000;
    vector<uint32_t> ambCount(nFeatures, 0);
    
    // Set up realistic distribution
    // High expressers (50 genes)
    for (uint32_t i = 0; i < 50; i++) {
        ambCount[i] = 100 + (i * 10);
    }
    // Medium expressers (150 genes)
    for (uint32_t i = 50; i < 200; i++) {
        ambCount[i] = 10 + (i % 20);
    }
    // Low expressers (300 genes)
    for (uint32_t i = 200; i < 500; i++) {
        ambCount[i] = 1 + (i % 5);
    }
    // Rest remain at 0 (undetected)
    
    // Build detected features list
    vector<uint32_t> featDetVec;
    for (uint32_t i = 0; i < nFeatures; i++) {
        if (ambCount[i] > 0) featDetVec.push_back(i);
    }
    
    // Compute profile
    AmbientProfile profile = EmptyDropsMultinomial::computeAmbientProfile(
        ambCount, nFeatures, featDetVec, featDetVec.size());
    
    // Verify profile size
    if (profile.ambProfileLogP.size() != nFeatures) {
        cerr << "FAIL: Profile size " << profile.ambProfileLogP.size() << " != nFeatures " << nFeatures << endl;
        return 1;
    }
    
    // Verify features detected
    if (profile.featuresDetected == 0) {
        cerr << "FAIL: No features detected" << endl;
        return 1;
    }
    
    // Expected: 50 + 150 + 300 = 500 detected features
    if (profile.featuresDetected != 500) {
        cerr << "FAIL: Expected 500 detected features, got " << profile.featuresDetected << endl;
        return 1;
    }
    
    // Verify sum of exp(logP) ≈ 1 (probability distribution)
    double probSum = 0.0;
    for (size_t i = 0; i < profile.ambProfileLogP.size(); i++) {
        if (profile.ambProfileLogP[i] > -500) {  // Avoid underflow
            probSum += exp(profile.ambProfileLogP[i]);
        }
    }
    
    if (fabs(probSum - 1.0) > 0.01) {
        cerr << "FAIL: Probability sum = " << probSum << " (expected ~1.0)" << endl;
        return 1;
    }
    
    // Verify nonzero probabilities are tracked
    if (profile.ambProfilePnon0.empty()) {
        cerr << "FAIL: ambProfilePnon0 is empty" << endl;
        return 1;
    }
    
    cout << "PASS: Ambient profile computed correctly (sum=" << probSum 
         << ", detected=" << profile.featuresDetected 
         << ", nonzero=" << profile.ambProfilePnon0.size() << ")" << endl;
    return 0;
}
