#ifndef H_OrdMagStage
#define H_OrdMagStage

#include "IncludeDefine.h"
#include <string>

// Forward declarations
class Parameters;
class ParametersSolo;
class SoloFeature;

// ============================================================================
// Simple EmptyDrops (formerly OrdMag)
// A fast knee-based filter that identifies high-confidence cells without
// Monte Carlo simulation. Used as a fallback when EmptyDrops has insufficient
// data (too few ambient cells, candidates, or rescues).
// ============================================================================

// Structure to hold Simple EmptyDrops filtering results
struct SimpleEmptyDropsResult {
    vector<uint32> passingIndices;      // Cell indices that pass simple filter
    vector<uint32> candidateIndices;    // Candidate cell indices for EmptyDrops
    vector<uint32> ambientIndices;      // Ambient cell indices (for ambient profile)
    uint32 retainThreshold;             // UMI threshold for simple filtering
    uint32 nCellsSimple;                // Number of cells passing simple filter
    uint32 minUMI;                      // Minimum UMI for candidates
    uint32 medianVal;                   // Median UMI value
    uint32 candidateLastRank;           // Last rank considered as candidate
    pair<uint32, uint32> ambientRange;  // Range of ambient indices (1-based)
};

// Type alias for backwards compatibility
using OrdMagResult = SimpleEmptyDropsResult;

// Structure for Simple EmptyDrops parameters
// Zero-initialization allows populateConfigWithDefaults() to detect unset values
struct SimpleEmptyDropsParams {
    uint32 nExpectedCells = 0;      // Expected number of cells (default: 3000)
    double maxPercentile = 0.0;      // Max percentile for robust max (default: 0.99)
    double maxMinRatio = 0.0;         // Max/min ratio (default: 10.0)
    uint32 umiMin = 0;              // Minimum UMI threshold (default: 500)
    double umiMinFracMedian = 0.0;   // Min UMI as fraction of median (default: 0.01)
    uint32 candMaxN = 0;            // Maximum candidates (default: 20000)
    uint32 indMin = 0;              // Min index for ambient cells (default: 45000)
    uint32 indMax = 0;              // Max index for ambient cells (default: 90000)
    
    // Bootstrap parameters (Cell Ranger style)
    uint32 nBootstrapSamples = 100;  // Number of bootstrap samples (default: 100)
    double recoveredCellsQuantile = 0.99;  // Quantile for baseline (default: 0.99, i.e., top 1%)
    bool useBootstrap = true;        // Whether to use bootstrap (default: true)
    uint32 maxExpectedCells = 0;     // Max expected cells for search (0 = auto from chemistry)
    uint32 minRecoveredCells = 50;   // Minimum recovered cells (default: 50)
    
    // Simple EmptyDrops is DISABLED by default (runs as fallback only)
    // Set to false via --use-simple-empty-drops to force enable
    bool disabled = true;            // If true, only runs as fallback when ED fails
    
    uint32 bootstrapSeed = 0;        // Seed for bootstrap RNG (0 = use default seeds 1,2,3...)
    uint32 maxThreads = 0;           // Max threads for bootstrap (0 = auto: hardware_concurrency or OMP_NUM_THREADS)
};

// Type alias for backwards compatibility
using OrdMagParams = SimpleEmptyDropsParams;

// Main class for Simple EmptyDrops (knee/rank) filtering
// Formerly known as OrdMag - renamed for clarity
class SimpleEmptyDropsStage {
public:
    // Run Cell Ranger-style simple filtering (replicates cr_simple_filter from Python)
    // Input:
    //   - nUMIperCB: UMI counts per cell barcode
    //   - nCB: number of cell barcodes
    //   - params: Simple EmptyDrops parameters
    // Returns: SimpleEmptyDropsResult with filtered indices and thresholds
    static SimpleEmptyDropsResult runCRSimpleFilter(
        const vector<uint32>& nUMIperCB,
        uint32 nCB,
        const SimpleEmptyDropsParams& params
    );
    
    // Run Cell Ranger-style filtering with bootstrap (matches Python filter_cellular_barcodes_ordmag)
    // This is the full CR algorithm with bootstrap sampling for robustness
    static SimpleEmptyDropsResult runCRSimpleFilterBootstrap(
        const vector<uint32>& nUMIperCB,
        uint32 nCB,
        SimpleEmptyDropsParams& params  // Non-const: may update nExpectedCells if estimated
    );
    
    // Find number of cells within order of magnitude of baseline (matches Python find_within_ordmag)
    // Input:
    //   - counts: UMI counts (can be bootstrap sample)
    //   - baselineIdx: index of baseline cell (from top, 0-based)
    // Returns: number of cells with UMI >= cutoff (where cutoff = baseline/10)
    static uint32 findWithinOrdmag(
        const vector<uint32>& counts,
        uint32 baselineIdx
    );
    
    // Estimate recovered cells by minimizing loss (matches Python estimate_recovered_cells_ordmag)
    // Input:
    //   - nonzeroCounts: non-zero UMI counts
    //   - maxExpectedCells: maximum cells to search
    //   - recoveredCellsQuantile: quantile for baseline (default 0.99)
    // Returns: (estimated recovered cells, loss)
    static pair<uint32, double> estimateRecoveredCellsOrdmag(
        const vector<uint32>& nonzeroCounts,
        uint32 maxExpectedCells,
        double recoveredCellsQuantile = 0.99
    );
    
    // Write simple emptydrops outputs to directory
    // Outputs:
    //   - passing_barcodes.txt: barcodes that pass simple filter
    //   - filter_summary.json: JSON with counts, thresholds, parameters
    //   - filtered matrix (if requested)
    static void writeOutputs(
        const SimpleEmptyDropsResult& result,
        const vector<string>& barcodes,
        const string& outputDir,
        const SimpleEmptyDropsParams& params,
        bool writeFilteredMatrix = false
    );
    
    // Load simple emptydrops parameters from ParametersSolo (for EmptyDrops_CR compatibility)
    static SimpleEmptyDropsParams loadParamsFromSolo(const ParametersSolo& pSolo);
};

// Type alias for backwards compatibility
using OrdMagStage = SimpleEmptyDropsStage;

#endif

