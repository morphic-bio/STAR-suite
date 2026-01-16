#include "OccupancyGuard.h"
#include "pcg_random.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <random>
#include <set>
#include <fstream>
#include <thread>
#include <vector>

using namespace std;

uint32_t OccupancyGuard::computePartitions(const string& barcode, uint32_t totalPartitions) {
    // Compute partitions from barcode sequence
    // For CB16 barcodes, we use a hash-based approach to map to partition space
    // This matches the Python implementation's partition calculation
    
    if (barcode.empty()) return 0;
    
    // Simple hash of barcode to partition space
    // Use first 16 characters (CB16) if available
    string cb16 = barcode.substr(0, min(static_cast<size_t>(16), barcode.size()));
    
    // Hash function (simple polynomial hash)
    uint64_t hash = 0;
    for (char c : cb16) {
        hash = hash * 31 + static_cast<unsigned char>(c);
    }
    
    // Map to partition space [0, totalPartitions)
    return static_cast<uint32_t>(hash % totalPartitions);
}

double OccupancyGuard::computePercentileThreshold(const vector<double>& occupancies, double percentile) {
    if (occupancies.empty()) return 0.0;
    
    vector<double> sorted = occupancies;
    sort(sorted.begin(), sorted.end());
    
    // Compute percentile index
    double index = percentile * (sorted.size() - 1);
    size_t lower = static_cast<size_t>(floor(index));
    size_t upper = static_cast<size_t>(ceil(index));
    
    if (lower == upper) {
        return sorted[lower];
    }
    
    // Linear interpolation
    double weight = index - lower;
    return sorted[lower] * (1.0 - weight) + sorted[upper] * weight;
}

vector<uint32_t> OccupancyGuard::filterHighOccupancy(
    const vector<uint32_t>& umiCounts,
    const vector<string>& barcodes,
    const Config& config,
    OccupancyMode mode,
    double lambdaEstimate,
    uint32_t tagLength,
    uint32_t simulatedGems,
    uint64_t seed) {
    
    vector<uint32_t> toRemove;
    
    if (umiCounts.size() != barcodes.size()) {
        return toRemove; // Size mismatch
    }
    
    if (umiCounts.empty()) {
        return toRemove;
    }
    
    if (mode == OccupancyMode::MonteCarlo) {
        // Monte Carlo mode: Port Python's occupancy_threshold() logic (lines 585-615)
        // Estimate tag distribution from observed barcodes
        unordered_map<string, uint32_t> tagCounts;
        for (const auto& barcode : barcodes) {
            if (barcode.size() >= tagLength) {
                string tag = barcode.substr(barcode.size() - tagLength);
                tagCounts[tag]++;
            }
        }
        
        if (tagCounts.empty() || lambdaEstimate <= 0.0) {
            // Fallback to hash-based if no tags or invalid lambda
            return filterHighOccupancy(umiCounts, barcodes, config, OccupancyMode::HashBased, 0.0, tagLength, simulatedGems, seed);
        }
        
        // Extract tags and probabilities
        vector<string> tags;
        vector<double> probs;
        uint32_t totalCount = 0;
        for (const auto& pair : tagCounts) {
            tags.push_back(pair.first);
            totalCount += pair.second;
        }
        for (const auto& tag : tags) {
            probs.push_back(static_cast<double>(tagCounts[tag]) / static_cast<double>(totalCount));
        }
        
        // Simulate GEMs: draw Poisson samples and count unique tags per GEM (parallel with std::thread)
        vector<uint32_t> uniqueCounts(simulatedGems, 0);
        uint64_t baseSeed = (seed == 0) ? 19760110LLU : seed;

        // Use std::threads - default to hardware_concurrency, but respect OMP_NUM_THREADS if set
        uint32_t nThreads = thread::hardware_concurrency();
        if (nThreads == 0) nThreads = 4;  // Fallback
        const char* ompEnv = getenv("OMP_NUM_THREADS");
        if (ompEnv != nullptr) {
            int envThreads = atoi(ompEnv);
            if (envThreads > 0) nThreads = static_cast<uint32_t>(envThreads);
        }
        nThreads = min(nThreads, simulatedGems);
        vector<thread> threads;
        
        // Precompute cumulative probabilities for fast discrete sampling (float precision is sufficient)
        vector<float> cumProbs(probs.size());
        float cumSum = 0.0f;
        for (size_t i = 0; i < probs.size(); i++) {
            cumSum += static_cast<float>(probs[i]);
            cumProbs[i] = cumSum;
        }
        // Normalize to ensure last entry is exactly 1.0
        if (cumSum > 0.0f) {
            for (size_t i = 0; i < cumProbs.size(); i++) {
                cumProbs[i] /= cumSum;
            }
        }
        
        auto simulateWorker = [&](uint32_t startGem, uint32_t endGem) {
            // PCG32 is ~3x faster than mt19937
            pcg32 rng(baseSeed + startGem, startGem);
            // Use float for Poisson parameter (sufficient precision for lambda < 100)
            float lambdaF = static_cast<float>(lambdaEstimate);
            
            for (uint32_t gem = startGem; gem < endGem; gem++) {
                // Fast Poisson sampling using Knuth's algorithm (sufficient for small lambda)
                // For lambda < 30, this is faster than std::poisson_distribution
                uint32_t poissonDraw = 0;
                if (lambdaF > 0.0f) {
                    float L = expf(-lambdaF);
                    float p = 1.0f;
                    do {
                        poissonDraw++;
                        // PCG's ldexp gives uniform float in [0,1)
                        p *= ldexpf(static_cast<float>(rng()), -32);
                    } while (p > L);
                    poissonDraw--;
                }
                
                if (poissonDraw == 0) {
                    uniqueCounts[gem] = 0;
                    continue;
                }

                // Use bitset or small set for unique tag tracking (typically < 10 tags)
                unordered_set<uint32_t> uniqueTagIdx;
                uniqueTagIdx.reserve(min(poissonDraw, static_cast<uint32_t>(tags.size())));
                
                for (uint32_t cell = 0; cell < poissonDraw; cell++) {
                    // Fast discrete sampling using binary search on cumulative probs
                    float u = ldexpf(static_cast<float>(rng()), -32);
                    auto it = lower_bound(cumProbs.begin(), cumProbs.end(), u);
                    uint32_t tagIdx = static_cast<uint32_t>(it - cumProbs.begin());
                    if (tagIdx >= tags.size()) tagIdx = static_cast<uint32_t>(tags.size() - 1);
                    uniqueTagIdx.insert(tagIdx);
                }
                uniqueCounts[gem] = static_cast<uint32_t>(uniqueTagIdx.size());
            }
        };
        
        uint32_t chunkSize = (simulatedGems + nThreads - 1) / nThreads;
        for (uint32_t t = 0; t < nThreads; t++) {
            uint32_t startGem = t * chunkSize;
            uint32_t endGem = min(startGem + chunkSize, simulatedGems);
            if (startGem < endGem) {
                threads.emplace_back(simulateWorker, startGem, endGem);
            }
        }
        for (auto& th : threads) th.join();

        // Compact non-zero entries (matches previous behavior of skipping zero-draw GEMs)
        vector<uint32_t> filteredUniqueCounts;
        filteredUniqueCounts.reserve(simulatedGems);
        for (uint32_t v : uniqueCounts) {
            if (v > 0) {
                filteredUniqueCounts.push_back(v);
            }
        }

        if (filteredUniqueCounts.empty()) {
            // Fallback to hash-based if no valid simulations
            return filterHighOccupancy(umiCounts, barcodes, config, OccupancyMode::HashBased, 0.0, tagLength, simulatedGems, seed);
        }
        
        // Compute percentile threshold from simulated unique counts
        // This threshold represents the maximum number of unique tags per GEM
        sort(filteredUniqueCounts.begin(), filteredUniqueCounts.end());
        double index = config.percentile * (filteredUniqueCounts.size() - 1);
        size_t lower = static_cast<size_t>(floor(index));
        size_t upper = static_cast<size_t>(ceil(index));
        
        uint32_t occupancyThreshold;
        if (lower == upper) {
            occupancyThreshold = filteredUniqueCounts[lower];
        } else {
            double weight = index - lower;
            occupancyThreshold = static_cast<uint32_t>(ceil(filteredUniqueCounts[lower] * (1.0 - weight) + filteredUniqueCounts[upper] * weight));
        }
        
        // Per-tag filtering: Group cells by GEM (CB16) and count UNIQUE TAGS per GEM
        // Remove ALL cells from GEMs that have more unique tags than the MC threshold
        // (This matches what the MC simulation computes: unique tags per GEM)
        
        // Determine CB16 length (barcode length - tag length)
        size_t cb16Length = 16;  // Default CB16 length
        if (!barcodes.empty() && barcodes[0].size() > tagLength) {
            cb16Length = barcodes[0].size() - tagLength;
        }
        
        // Group cells by GEM (CB16) - track unique tags AND cell indices per GEM
        unordered_map<string, unordered_set<string>> gemToUniqueTags;  // CB16 -> unique TAG8s
        unordered_map<string, vector<uint32_t>> gemToCellIndices;      // CB16 -> cell indices
        
        for (size_t i = 0; i < barcodes.size(); i++) {
            if (barcodes[i].size() < cb16Length + tagLength) {
                continue;  // Skip malformed barcodes
            }
            string gem = barcodes[i].substr(0, cb16Length);  // CB16
            string tag = barcodes[i].substr(barcodes[i].size() - tagLength);  // TAG8
            gemToUniqueTags[gem].insert(tag);
            gemToCellIndices[gem].push_back(static_cast<uint32_t>(i));
        }
        
        // Find GEMs where UNIQUE TAG COUNT exceeds threshold (matches MC simulation)
        uint32_t highOccupancyGems = 0;
        for (const auto& gemEntry : gemToUniqueTags) {
            const string& gem = gemEntry.first;
            const unordered_set<string>& uniqueTags = gemEntry.second;
            if (uniqueTags.size() > occupancyThreshold) {
                highOccupancyGems++;
                // Remove ALL cells in this high-occupancy GEM
                const vector<uint32_t>& cellIndices = gemToCellIndices[gem];
                for (uint32_t idx : cellIndices) {
                    toRemove.push_back(idx);
                }
            }
        }
        
        // Debug output
        {
            ofstream debugFile("/tmp/occupancy_debug.txt", ios::app);
            debugFile << "Per-tag occupancy filter (unique tags per GEM):\n";
            debugFile << "  Total cells: " << barcodes.size() << "\n";
            debugFile << "  Total GEMs (unique CB16): " << gemToUniqueTags.size() << "\n";
            debugFile << "  Simulated threshold (99.9% quantile): " << occupancyThreshold << " unique tags/GEM\n";
            debugFile << "  High-occupancy GEMs (>" << occupancyThreshold << " unique tags): " << highOccupancyGems << "\n";
            debugFile << "  Cells removed from high-occupancy GEMs: " << toRemove.size() << "\n";
            debugFile << "\n";
            debugFile.close();
        }
    } else {
        // Hash-based mode: Compute occupancy directly (existing logic)
        vector<double> occupancies;
        occupancies.reserve(umiCounts.size());
        
        for (size_t i = 0; i < umiCounts.size(); i++) {
            uint32_t partitions = computePartitions(barcodes[i], config.totalPartitions);
            if (partitions == 0) {
                // Skip cells with zero partitions (edge case)
                occupancies.push_back(0.0);
                continue;
            }
            
            double occupancy = (static_cast<double>(umiCounts[i]) / partitions) * config.recoveryFactor;
            occupancies.push_back(occupancy);
        }
        
        // Compute percentile threshold
        double threshold = computePercentileThreshold(occupancies, config.percentile);
        
        // Find cells exceeding threshold
        for (size_t i = 0; i < occupancies.size(); i++) {
            if (occupancies[i] > threshold) {
                toRemove.push_back(static_cast<uint32_t>(i));
            }
        }
    }
    
    return toRemove;
}

vector<uint32_t> OccupancyGuard::filterLowUMI(
    const vector<uint32_t>& umiCounts,
    uint32_t threshold,
    double fracMedian) {
    
    vector<uint32_t> toRemove;
    
    if (umiCounts.empty()) {
        return toRemove;
    }
    
    uint32_t actualThreshold = threshold;
    
    // If threshold is 0, compute from median
    if (threshold == 0) {
        vector<uint32_t> sorted = umiCounts;
        sort(sorted.begin(), sorted.end());
        
        uint32_t median = sorted[sorted.size() / 2];
        actualThreshold = static_cast<uint32_t>(round(median * fracMedian));
    }
    
    // Find cells at or below threshold (FRP guard: <= threshold)
    // Note: For FRP guard with threshold=12, we remove cells with UMI <= 12
    for (size_t i = 0; i < umiCounts.size(); i++) {
        if (umiCounts[i] <= actualThreshold) {
            toRemove.push_back(static_cast<uint32_t>(i));
        }
    }
    
    return toRemove;
}
