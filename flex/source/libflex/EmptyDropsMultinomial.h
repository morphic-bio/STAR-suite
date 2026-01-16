#ifndef H_EmptyDropsMultinomial
#define H_EmptyDropsMultinomial

#include "IncludeDefine.h"
#include <vector>
#include <random>

// Forward declarations
class Parameters;
class ParametersSolo;

// Structure to hold p-value and FDR results for each candidate cell
struct EmptyDropsResult {
    uint32 cellIndex;
    double pValue;
    double pAdjusted;
    bool passesRawP;    // True if pValue <= rawPvalueThreshold
    bool passesFDR;     // True if pAdjusted <= FDR threshold
    double obsLogProb;  // Observed log probability (for diagnostics)
    uint32 monteCarloRank;  // Rank in Monte Carlo distribution (for diagnostics)
};

// Structure to hold ambient profile data
struct AmbientProfile {
    vector<double> ambProfileLogP;      // Log probabilities for all features
    vector<double> ambProfilePnon0;      // Non-zero probabilities (for sampling)
    vector<double> ambProfileLogPnon0;   // Log of non-zero probabilities
    uint32 featuresDetected;                // Number of detected features
};

// Structure for multinomial computation parameters
// Zero-initialization allows populateConfigWithDefaults() to detect unset values
struct EmptyDropsParams {
    uint32 indMin = 0;           // Min index for empty cells
    uint32 indMax = 0;           // Max index for empty cells
    uint32 umiMin = 0;           // Minimum UMI threshold
    double umiMinFracMedian = 0.0;   // Minimum UMI as fraction of median
    uint32 candMaxN = 0;         // Maximum number of candidates
    double FDR = 0.0;                // FDR threshold (for FDR-based filtering)
    double rawPvalueThreshold = 0.0; // Raw p-value threshold (default: 0.999, permissive)
    bool useFDRGate = true;          // If true, use FDR threshold; if false, use raw p-value threshold
    uint32 simN = 0;             // Number of Monte Carlo simulations
    uint64 seed = 0;             // Random seed (default: 19760110LLU)
    uint32 lowerTestingBound = 0;  // Lower UMI bound for testing: cells with UMI <= this are excluded from testing (default: 500, R's umi.min)
    uint32 ambientUmiMax = 0;      // Max UMI for ambient cells: cells with UMI <= this used for ambient profile (default: 100, R's lower)
    uint32 maxTotalBuckets = 0;    // Max buckets for total binning (0 = disabled, use exact totals)
    uint32 mcThreads = 0;        // Threads for Monte Carlo simulation (0 = single-threaded)
};

// Main class for EmptyDrops multinomial computations
class EmptyDropsMultinomial {
public:
    // Compute ambient profile from empty cells
    // Input:
    //   - ambCount: gene counts summed over empty cells (featuresNumber elements)
    //   - featuresNumber: total number of features
    //   - featDet: set of detected feature indices
    //   - featDetN: number of detected features
    // Returns: AmbientProfile structure
    static AmbientProfile computeAmbientProfile(
        const vector<uint32>& ambCount,
        uint32 featuresNumber,
        const vector<uint32>& featDetVec,
        uint32 featDetN
    );
    
    // Compute log multinomial PDF for sparse representation
    // This matches the existing logMultinomialPDFsparse function
    static double logMultinomialPDFsparse(
        const vector<double>& ambProfileLogP,
        const vector<uint32>& countCellGeneUMI,
        uint32 stride,
        uint32 shift,
        int64 start,
        uint32 nGenes,
        const vector<double>& logFactorial
    );
    
    // Run Monte Carlo simulations and compute p-values/FDR
    // NOTE: This function expects a pre-built COMPACT matrix and ambient profile
    //       (built by FlexFilter per compact_copy_plan.md). Gene IDs in countCellGeneUMI
    //       are already compact indices; ambProfile is the compact ambient profile.
    // Input:
    //   - ambProfile: COMPACT ambient profile from computeAmbientProfile (compact gene indices)
    //   - candidateIndices: cell indices to test (sorted by UMI count descending)
    //   - candidateCounts: UMI counts for each candidate
    //   - countCellGeneUMI: COMPACT matrix data (gene IDs are already compact indices)
    //   - countCellGeneUMIindex: index array for compact matrix
    //   - nGenePerCB: number of genes per cell (in compact matrix)
    //   - countMatStride: stride for count matrix (typically 2 for compact)
    //   - umiDedupCountIndMain: offset for main count (typically 1)
    //   - params: EmptyDrops parameters
    //   - nSimpleCells: number of simple cells (first nSimpleCells candidates) to skip simulation (default 0)
    //   - featureNames, nTotalCells, debugOutputDir: kept for backward compatibility, now unused
    //   - tagName: tag name for diagnostics (default empty)
    //   - enableInvariantChecks: if true, run invariant checks and throw on failure (default false)
    // Returns: vector of EmptyDropsResult with p-values and FDR
    static vector<EmptyDropsResult> computePValues(
        const AmbientProfile& ambProfile,
        const vector<uint32>& candidateIndices,
        const vector<uint32>& candidateCounts,
        const vector<uint32>& countCellGeneUMI,
        const vector<uint32>& countCellGeneUMIindex,
        const vector<uint32>& nGenePerCB,
        uint32 countMatStride,
        uint32 umiDedupCountIndMain,
        const EmptyDropsParams& params,
        uint32 nSimpleCells = 0,
        const vector<string>& featureNames = vector<string>(),
        uint32 nTotalCells = 0,
        const string& debugOutputDir = string(),
        const string& tagName = string(),
        bool enableInvariantChecks = false
    );
    
    // Helper: compute log factorial table
    static vector<double> computeLogFactorial(uint32 maxCount);
};

#endif
