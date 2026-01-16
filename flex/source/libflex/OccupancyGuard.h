#ifndef H_OccupancyGuard
#define H_OccupancyGuard

#include "IncludeDefine.h"
#include <vector>
#include <string>

// Occupancy mode enum
enum class OccupancyMode {
    HashBased,   // Deterministic hash (parity mode - matches current C++)
    MonteCarlo   // Full simulation (production mode - more accurate)
};

// Occupancy/FRP guard helper class
// Matches Python defaults exactly (immutable CR defaults)
class OccupancyGuard {
public:
    // Configuration (matches CR defaults exactly)
    struct Config {
        uint32_t totalPartitions = 115000;  // CR default (immutable)
        double recoveryFactor = 1.0 / 1.65;  // CR default (immutable)
        double percentile = 0.999;  // CR default (immutable)
    };
    
    // Filter cells with high occupancy (exceeding percentile threshold)
    // Input:
    //   - umiCounts: UMI counts per cell (aligned with barcodes)
    //   - barcodes: Cell barcode strings (for computing partitions)
    //   - config: Occupancy configuration (uses CR defaults)
    //   - mode: OccupancyMode::HashBased (deterministic) or OccupancyMode::MonteCarlo (simulation)
    //   - lambdaEstimate: Poisson lambda estimate (for Monte Carlo mode)
    //   - tagLength: Length of tag suffix in barcodes (for Monte Carlo mode)
    //   - simulatedGems: Number of GEMs to simulate (for Monte Carlo mode, default 1M)
    //   - seed: RNG seed (for Monte Carlo mode, optional)
    // Returns: Vector of cell indices to remove (0-based indices into umiCounts/barcodes)
    static std::vector<uint32_t> filterHighOccupancy(
        const std::vector<uint32_t>& umiCounts,
        const std::vector<std::string>& barcodes,
        const Config& config,
        OccupancyMode mode = OccupancyMode::HashBased,
        double lambdaEstimate = 0.0,
        uint32_t tagLength = 8,
        uint32_t simulatedGems = 1000000,
        uint64_t seed = 0);
    
    // Filter cells with low UMI counts
    // Input:
    //   - umiCounts: UMI counts per cell
    //   - threshold: Absolute UMI threshold (if > 0, use this; if 0, use fracMedian logic)
    //   - fracMedian: Fraction of median UMI (used if threshold == 0)
    // Returns: Vector of cell indices to remove
    static std::vector<uint32_t> filterLowUMI(
        const std::vector<uint32_t>& umiCounts,
        uint32_t threshold,
        double fracMedian);
    
private:
    // Helper: Compute partitions per barcode from barcode string
    // Partitions are computed from barcode sequence (typically CB16)
    static uint32_t computePartitions(const std::string& barcode, uint32_t totalPartitions);
    
    // Helper: Compute percentile threshold from occupancy values
    static double computePercentileThreshold(const std::vector<double>& occupancies, double percentile);
};

#endif

