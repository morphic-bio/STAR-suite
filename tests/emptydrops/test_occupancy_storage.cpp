#include "ObservedTagOccupancy.h"
#include "OccupancyGuard.h"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <random>

static std::string dna(uint64_t n, unsigned length) {
    std::string s(length, 'A');
    for (auto& c : s) { c = "ACGT"[n % 4]; n /= 4; }
    return s;
}

int main() {
    std::mt19937 rng(32143);
    for (unsigned trial = 0; trial < 30; ++trial) {
        std::vector<std::string> calls;
        for (unsigned gem = 0; gem < trial * 13; ++gem) {
            unsigned tags = 1 + rng() % (trial % 3 ? 4 : 80);
            for (unsigned t = 0; t < tags; ++t) {
                calls.push_back(dna(gem, 16) + dna(t, 8));
                if (rng() % 4 == 0) calls.push_back(calls.back() + "-1");
            }
        }
        std::shuffle(calls.begin(), calls.end(), rng);
        for (double p : {.95, .99, .999}) {
            const auto r = fitObservedTagOccupancy(calls, p);
            std::cout << trial << ' ' << std::hexfloat << r.lambda << ' ' << r.occupiedMean
                      << std::defaultfloat << ' ' << r.cutoff << ' ' << r.occupiedGems;
            for (const auto& gem : r.rejectedGems) std::cout << ' ' << gem;
            std::cout << '\n';
        }
    }
    for (const std::string& invalid : {std::string(), dna(0, 23), dna(0, 25),
         dna(0, 24) + "-2", dna(0, 23) + "N", dna(0, 23) + "a"}) {
        bool failed = false;
        try { fitObservedTagOccupancy({invalid}); } catch (const std::runtime_error&) { failed = true; }
        assert(failed);
    }
    for (double p : {0., 1., -1., std::nan("")}) {
        bool failed = false;
        try { fitObservedTagOccupancy({}, p); } catch (const std::runtime_error&) { failed = true; }
        assert(failed);
    }
    // Exercise general (non-DNA and variable-length) GEM/tag grouping as well.
    for (unsigned tagLength : {0, 1, 8, 20}) {
        std::vector<std::string> barcodes;
        for (unsigned gem = 0; gem < 80; ++gem) {
            const std::string prefix = "gem-prefix-" + dna(gem, 6);
            for (unsigned t = 0; t < 1 + gem % 23; ++t) {
                const std::string bc = prefix + dna(t, tagLength);
                barcodes.push_back(bc);
                if (t % 3 == 0) barcodes.push_back(bc);
            }
        }
        barcodes.push_back("x"); // skipped by GEM grouping; counted by old tag distribution
        barcodes.push_back(std::string(60, 'x')); // preserves last-tag semantics
        const std::vector<uint32_t> counts(barcodes.size(), 900);
        OccupancyGuard::Config config;
        for (double lambda : {0., .1, 2., 5.}) {
            auto removed = OccupancyGuard::filterHighOccupancy(counts, barcodes, config,
                OccupancyMode::MonteCarlo, lambda, tagLength, 10000, 8831);
            std::sort(removed.begin(), removed.end());
            std::cout << "MC " << tagLength << ' ' << lambda;
            for (auto i : removed) std::cout << ' ' << i;
            std::cout << '\n';
        }
    }
}
