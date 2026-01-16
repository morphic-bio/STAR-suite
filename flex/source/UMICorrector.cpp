#include "UMICorrector.h"
#include "UmiCodec.h"
#include <algorithm>
#include <queue>

UMICorrectionResult UMICorrector::correctClique(const std::vector<UMICount>& counts, 
                                                 const UMIParams& params) {
    UMICorrectionResult result;
    
    if (counts.empty()) {
        return result;
    }
    
    // Build UMI histogram (packed UMI -> count)
    std::unordered_map<uint32_t, uint32_t> hist;
    
    // First pass: aggregate counts for each UMI
    for (const auto& count : counts) {
        uint32_t packed = encodeUMI12(count.ur);
        if (packed == UINT32_MAX) continue; // Invalid UMI
        
        hist[packed] += count.readCount;
    }
    
    // Second pass: filter by min_count (after aggregation)
    for (auto it = hist.begin(); it != hist.end();) {
        if (it->second < static_cast<uint32_t>(params.minCount)) {
            it = hist.erase(it);
        } else {
            ++it;
        }
    }
    
    // Process each unvisited UMI to find its connected component
    std::unordered_set<uint32_t> visited;
    std::vector<uint32_t> componentUmis;
    std::vector<uint32_t> componentCounts;
    
    for (const auto& histEntry : hist) {
        uint32_t packedUmi = histEntry.first;
        
        // Check if already visited
        if (visited.find(packedUmi) != visited.end()) {
            continue;
        }
        
        // Find connected component
        int componentSize = 0;
        componentUmis.clear();
        componentCounts.clear();
        findConnectedComponent(packedUmi, componentUmis, componentCounts, 
                             componentSize, hist, visited, params.maxComponentSize);
        
        if (componentSize == 0) continue;
        
        result.components++;
        
        // Track component size
        result.componentSizes.push_back(componentSize);
        
        // Check component size cap
        if (componentSize > params.maxComponentSize) {
            result.componentsCapped++;
            // Skip this component (don't merge)
            continue;
        }
        
        // Find winner (highest count UMI) and second-best
        uint32_t winner = findWinnerUmi(componentUmis, componentCounts, componentSize);
        if (winner == UINT32_MAX) continue;
        
        // Get winner count and second-best count
        uint32_t winnerCount = 0;
        uint32_t secondBestCount = 0;
        for (int i = 0; i < componentSize; i++) {
            if (componentUmis[i] == winner) {
                winnerCount = componentCounts[i];
            } else if (componentCounts[i] > secondBestCount) {
                secondBestCount = componentCounts[i];
            }
        }
        
        // Enforce ratio_thresh: winner must be >= ratio_thresh * second_best
        // If second_best is 0, always accept (only one UMI in component)
        if (secondBestCount > 0) {
            double ratio = static_cast<double>(winnerCount) / static_cast<double>(secondBestCount);
            if (ratio < params.ratioThresh) {
                // Ratio threshold not met - skip merging this component
                result.componentsBelowThreshold++;
                continue;
            }
        }
        
        // Decode winner to UB string
        std::string winnerUb = decodeUMI12(winner);
        
        // Create corrections for all UMIs in component
        for (int i = 0; i < componentSize; i++) {
            if (componentUmis[i] == winner) continue; // Winner maps to itself
            
            // Decode UR
            std::string ur = decodeUMI12(componentUmis[i]);
            
            // Add correction mapping
            result.urToUb[ur] = winnerUb;
            result.merges++;
        }
    }
    
    return result;
}

void UMICorrector::findNeighbors(uint32_t packedUmi, 
                                std::vector<uint32_t>& neighbors,
                                const std::unordered_map<uint32_t, uint32_t>& hist,
                                const std::unordered_set<uint32_t>& visited) {
    neighbors.clear();
    
    // For each position in the UMI, toggle each possible base (A,C,G,T)
    for (int pos = 0; pos < 12; pos++) {
        uint32_t posMask = 3u << (pos * 2);  // Mask for this position
        uint32_t currentBase = (packedUmi >> (pos * 2)) & 3;
        
        // Try each of the 3 other bases
        for (uint32_t newBase = 0; newBase < 4; newBase++) {
            if (newBase == currentBase) continue;
            
            uint32_t neighbor = (packedUmi & ~posMask) | (newBase << (pos * 2));
            
            // Check if neighbor exists in histogram and not visited
            if (hist.find(neighbor) != hist.end() && 
                visited.find(neighbor) == visited.end()) {
                neighbors.push_back(neighbor);
            }
        }
    }
}

void UMICorrector::findConnectedComponent(uint32_t startUmi,
                                         std::vector<uint32_t>& componentUmis,
                                         std::vector<uint32_t>& componentCounts,
                                         int& componentSize,
                                         const std::unordered_map<uint32_t, uint32_t>& hist,
                                         std::unordered_set<uint32_t>& visited,
                                         int maxSize) {
    componentSize = 0;
    componentUmis.clear();
    componentCounts.clear();
    
    std::queue<uint32_t> queue;
    queue.push(startUmi);
    
    std::vector<uint32_t> neighbors;
    
    while (!queue.empty() && componentSize < maxSize) {
        uint32_t currentUmi = queue.front();
        queue.pop();
        
        // Check if already visited
        if (visited.find(currentUmi) != visited.end()) {
            continue;
        }
        
        // Check if exists in histogram
        auto histIt = hist.find(currentUmi);
        if (histIt == hist.end()) {
            continue;
        }
        
        // Mark as visited
        visited.insert(currentUmi);
        
        // Add to component
        componentUmis.push_back(currentUmi);
        componentCounts.push_back(histIt->second);
        componentSize++;
        
        // Find neighbors and enqueue them
        findNeighbors(currentUmi, neighbors, hist, visited);
        for (uint32_t neighbor : neighbors) {
            if (componentSize < maxSize) {
                queue.push(neighbor);
            }
        }
    }
}

uint32_t UMICorrector::findWinnerUmi(const std::vector<uint32_t>& componentUmis,
                                     const std::vector<uint32_t>& componentCounts,
                                     int componentSize) {
    if (componentSize == 0) return UINT32_MAX;
    
    uint32_t winner = componentUmis[0];
    uint32_t maxCount = componentCounts[0];
    
    for (int i = 1; i < componentSize; i++) {
        if (componentCounts[i] > maxCount) {
            maxCount = componentCounts[i];
            winner = componentUmis[i];
        }
    }
    
    return winner;
}

