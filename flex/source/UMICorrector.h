#ifndef CODE_UMICorrector
#define CODE_UMICorrector

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <string>
#include <cstdint>
#include <queue>

// UMI correction result
struct UMICorrectionResult {
    std::unordered_map<std::string, std::string> urToUb;  // UR -> corrected UB mapping
    uint32_t merges;
    uint32_t components;
    uint32_t componentsCapped;
    uint32_t componentsBelowThreshold;
    std::vector<uint32_t> componentSizes;  // Size of each component processed
    
    UMICorrectionResult() : merges(0), components(0), componentsCapped(0), componentsBelowThreshold(0) {}
};

// UMI count input
struct UMICount {
    std::string ur;  // UR tag (12bp)
    uint32_t readCount;
    
    UMICount(const std::string& u, uint32_t c) : ur(u), readCount(c) {}
};

// UMI correction parameters
struct UMIParams {
    int minCount;
    double ratioThresh;
    int maxComponentSize;
    
    UMIParams(int min_c, double ratio, int max_size) 
        : minCount(min_c), ratioThresh(ratio), maxComponentSize(max_size) {}
};

class UMICorrector {
public:
    // Main correction function (clique/connected-component mode)
    static UMICorrectionResult correctClique(const std::vector<UMICount>& counts, 
                                             const UMIParams& params);
    
private:
    // Find Hamming-1 neighbors of a packed UMI
    static void findNeighbors(uint32_t packedUmi, 
                             std::vector<uint32_t>& neighbors,
                             const std::unordered_map<uint32_t, uint32_t>& hist,
                             const std::unordered_set<uint32_t>& visited);
    
    // Find connected component using BFS
    static void findConnectedComponent(uint32_t startUmi,
                                      std::vector<uint32_t>& componentUmis,
                                      std::vector<uint32_t>& componentCounts,
                                      int& componentSize,
                                      const std::unordered_map<uint32_t, uint32_t>& hist,
                                      std::unordered_set<uint32_t>& visited,
                                      int maxSize);
    
    // Find UMI with highest count in component
    static uint32_t findWinnerUmi(const std::vector<uint32_t>& componentUmis,
                                  const std::vector<uint32_t>& componentCounts,
                                  int componentSize);
};

#endif

