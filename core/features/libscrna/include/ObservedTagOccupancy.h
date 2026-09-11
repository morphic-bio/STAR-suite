#ifndef SCRNA_OBSERVED_TAG_OCCUPANCY_H
#define SCRNA_OBSERVED_TAG_OCCUPANCY_H

#include <cmath>
#include <algorithm>
#include <cstdint>
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
    // One 48-bit CB16+TAG8 key per call, then one contiguous unique-key array.
    // The tag space is the full TAG8 space, not a fixed sixteen-tag mask.
    std::vector<uint64_t> keys;
    keys.reserve(calls.size());
    for (const auto& bc : calls) {
        if (bc.size() != 24 && !(bc.size() == 26 && bc.compare(24, 2, "-1") == 0))
            throw std::runtime_error("Occupancy requires a full CB16+TAG8 barcode");
        uint64_t key = 0;
        for (size_t i = 0; i < 24; ++i) {
            unsigned base;
            switch (bc[i]) {
                case 'A': base = 0; break;
                case 'C': base = 1; break;
                case 'G': base = 2; break;
                case 'T': base = 3; break;
                default: throw std::runtime_error("Occupancy requires a full CB16+TAG8 barcode");
            }
            key = (key << 2) | base;
        }
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    ObservedTagOccupancy result;
    if (keys.empty()) return result;
    for (size_t i = 0; i < keys.size(); ++i)
        if (i == 0 || (keys[i] >> 16) != (keys[i - 1] >> 16)) ++result.occupiedGems;
    result.occupiedMean = double(keys.size()) / result.occupiedGems;
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
    for (size_t first = 0; first < keys.size();) {
        const uint64_t gem = keys[first] >> 16;
        size_t end = first + 1;
        while (end < keys.size() && (keys[end] >> 16) == gem) ++end;
        if (end - first > result.cutoff) {
            std::string barcode(16, 'A');
            uint64_t remaining = gem;
            for (size_t p = 16; p > 0; --p) {
                barcode[p - 1] = "ACGT"[remaining & 3];
                remaining >>= 2;
            }
            result.rejectedGems.insert(std::move(barcode));
        }
        first = end;
    }
    return result;
}
#endif
