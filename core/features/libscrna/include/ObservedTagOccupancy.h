#ifndef SCRNA_OBSERVED_TAG_OCCUPANCY_H
#define SCRNA_OBSERVED_TAG_OCCUPANCY_H

#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

// Occupancy is the number of distinct observed TAG8 calls in a CB16 GEM.
// Fit the positive observations with a zero-truncated Poisson; absent tags
// and absent GEMs must not be manufactured in the fit.
struct ObservedTagOccupancy {
    double lambda = 0;
    double occupiedMean = 0;
    unsigned cutoff = 1;
    size_t occupiedGems = 0;
    std::set<std::string> rejectedGems;
};

inline ObservedTagOccupancy fitObservedTagOccupancy(
    const std::vector<std::string>& calls, double percentile = .999)
{
    if (!(percentile > 0 && percentile < 1))
        throw std::runtime_error("Occupancy percentile must be between zero and one");
    std::map<std::string, std::set<std::string>> gems;
    for (auto bc : calls) {
        if (bc.size() == 26 && bc.substr(24) == "-1") bc.resize(24);
        if (bc.size() != 24 || bc.find_first_not_of("ACGT") != std::string::npos)
            throw std::runtime_error("Occupancy requires a full CB16+TAG8 barcode");
        gems[bc.substr(0, 16)].insert(bc.substr(16));
    }
    ObservedTagOccupancy result;
    result.occupiedGems = gems.size();
    if (gems.empty()) return result;
    size_t n = 0;
    for (const auto& gem : gems) n += gem.second.size();
    result.occupiedMean = double(n) / gems.size();
    // At mean=1 the ZTP MLE is the boundary lambda=0, with all positive
    // observations equal to one. Keep those singleton observations.
    if (result.occupiedMean > 1) {
        double lo = 0, hi = result.occupiedMean;
        for (unsigned i = 0; i < 100; ++i) {
            const double mid = (lo + hi) / 2;
            if (mid / -std::expm1(-mid) < result.occupiedMean) lo = mid;
            else hi = mid;
        }
        result.lambda = (lo + hi) / 2;
        double mass = std::exp(-result.lambda), cdf = mass;
        unsigned k = 0;
        while (cdf < percentile) {
            mass *= result.lambda / ++k;
            cdf += mass;
            if (k > 10000) throw std::runtime_error("Occupancy CDF did not converge");
        }
        result.cutoff = k < 1 ? 1 : k;
    }
    for (const auto& gem : gems)
        if (gem.second.size() > result.cutoff) result.rejectedGems.insert(gem.first);
    return result;
}
#endif
