#include "OrdMagStage.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

// Frozen pre-optimization STAR oracle. It deliberately retains the repeated
// sort; do not share the optimized helper or alter its float trial arithmetic.
static uint32 legacyFind(vector<uint32> values, uint32 index) {
    if (values.empty()) return 0;
    std::sort(values.begin(), values.end());
    const uint32 n = values.size();
    if (index >= n) index = n - 1;
    const uint32 cutoff = std::max(uint32(1), uint32(std::round(.1 * values[n-index-1])));
    return n - (std::lower_bound(values.begin(), values.end(), cutoff) - values.begin());
}
static std::pair<uint32,double> legacyEstimate(const vector<uint32>& values, uint32 maximum, double q) {
    if (values.empty()) return {0,0};
    vector<uint32> options;
    const float logMax = log2f(float(maximum)), quantile = float(q);
    for (int i=0;i<2000;++i) {
        const float x = 1.f + (logMax-1.f)*i/1999.f;
        const uint32 n = uint32(roundf(powf(2.f,x)));
        if (options.empty() || n != options.back()) options.push_back(n);
    }
    uint32 best=1; float loss=1e30f;
    for (uint32 n : options) {
        uint32 index=uint32(roundf(n*(1.f-quantile)));
        if (index>=values.size()) index=values.size()-1;
        const float diff=float(legacyFind(values,index))-float(n);
        const float candidate=(diff*diff)/float(n);
        if (candidate<loss) { loss=candidate;best=n; }
    }
    return {best,double(loss)};
}
int main() {
    vector<vector<uint32>> cases={{},{0},{1},{4,5,6,14,15,16,24,25,26},
        {UINT32_MAX,UINT32_MAX-1,429496729,500,499,1,0},vector<uint32>(311,1),vector<uint32>(311,500)};
    std::mt19937 rng(3907);
    vector<uint32> random(997);
    for (auto& n:random) n=rng()%17000;
    cases.push_back(random);
    std::sort(random.begin(),random.end()); cases.push_back(random);
    std::reverse(random.begin(),random.end()); cases.push_back(random);
    size_t compared=0;
    for (const auto& values:cases) {
        const auto original=values;
        for (uint32 maximum:{1u,50u,45000u}) for (double q:{.99,.985,1.0}) {
            if (legacyEstimate(values,maximum,q)!=OrdMagStage::estimateRecoveredCellsOrdmag(values,maximum,q))
                throw std::runtime_error("estimator differs from frozen STAR oracle");
            ++compared;
        }
        for (uint32 index:{0u,1u,20u,100000u})
            if (legacyFind(values,index)!=OrdMagStage::findWithinOrdmag(values,index))
                throw std::runtime_error("inclusive threshold or index clamping changed");
        if (values!=original) throw std::runtime_error("input counts modified");
    }
    std::cout << "PASS: " << compared << " exact frozen-estimator comparisons, boundary cutoffs and immutable inputs\n";
}
