#ifndef H_FlexFilter
#define H_FlexFilter

#include "IncludeDefine.h"
#include "SampleMatrixData.h"
#include "OrdMagStage.h"
#include "EmptyDropsMultinomial.h"
#include "OccupancyGuard.h"  // Includes OccupancyMode enum
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>

// Forward declarations
struct FilterSummary;  // Defined in FixtureLoader.h (for tests) or inline here

// FlexFilter: Shared core library for OrdMag/EmptyDrops filtering
// Used by both STAR and standalone tools
class FlexFilter {
public:
    // File-based inputs (for CLI/standalone tool)
    struct FileInputs {
        std::string matrixDir;              // Path to MEX directory
        std::string allowedTagsPath;        // Path to allowed tags TSV (optional)
        std::vector<std::string> allowedTags;  // In-memory allowed tags (strings)
        std::vector<std::string> sampleLabels;  // Sample labels (from whitelist, optional)
        std::vector<std::string> sampleTags;    // TAG8 sequences (from whitelist, optional)
        std::string filterSummaryPath;     // Path to filter_summary.json (optional)
        // If filterSummaryPath is empty, FlexFilter must recompute per-tag expectations from MEX
        std::map<std::string, uint32_t> tagExpectedCells;  // Override: tag -> expected cells (optional)
    };
    
    // In-memory inputs (for STAR integration)
    struct MemoryInputs {
        SampleMatrixData matrixData;
        std::vector<std::string> observedBarcodes;  // Full observed barcodes list (denseId = index + 1)
        std::map<std::string, uint32_t> tagExpectedCells;  // Optional override
        std::vector<std::string> sampleLabels;  // Sample labels (for output paths)
        std::vector<std::string> sampleTags;    // TAG8 sequences (for filtering)
        // sampleTags[i] maps to the tag that should be used to filter matrix slices for sampleLabels[i]
    };
    
    // Outputs (structured results)
    struct Outputs {
        struct TagResults {
            std::string tag;
            std::string sampleLabel;
            uint32_t expectedCells;
            std::vector<std::string> tagBarcodes;       // Full barcode list for this tag (for index mapping)
            std::vector<std::string> retainBarcodes;    // Barcodes in retain window (ED test candidates) - cellIndex maps here
            std::vector<std::string> filteredBarcodes;  // Final filtered barcodes (same as passingBarcodes)
            std::vector<std::string> passingBarcodes;   // EmptyDrops passing + occupancy post-filter
            std::vector<std::string> edPasserBarcodes;  // ED passers before occupancy filter (for combined filter)
            std::vector<EmptyDropsResult> emptydropsResults;
            OrdMagResult ordmagResult;                  // OrdMag result for this tag (for compatibility)
            // Summary stats
            uint32_t nRetainWindow = 0;       // Cells in retain window
            uint32_t nSimpleCells = 0;        // Simple cells (auto-pass, UMI > retain threshold)
            uint32_t nTailTested = 0;         // Tail cells tested by ED
            uint32_t nSimplePassers = 0;      // Simple cells passing
            uint32_t nTailPassers = 0;        // Tail cells passing
            uint32_t occupancyRemoved = 0;    // Removed by occupancy post-filter
        };
        
        std::vector<TagResults> tagResults;
        // FilterSummary summary;  // Will be populated from tagResults
    };
    
    struct Config {
        SimpleEmptyDropsParams simpleEmptyDropsParams;  // Simple EmptyDrops (formerly OrdMag) - disabled by default
        EmptyDropsParams emptydropsParams;
        // Occupancy post-filter parameters
        uint32_t totalPartitions = 115000;
        double recoveryFactor = 1.0 / 1.65;
        double occupancyPercentile = 0.999;
        uint32_t occupancySimulatedGems = 1000000; // Number of GEMs to simulate in MC occupancy filter
        uint32_t totalExpectedCells = 0;  // Total expected cells across all tags (0 = auto: 3000 * numTags)
        OccupancyMode occupancyMode = OccupancyMode::MonteCarlo;
    // Tag allocation: use Winsorized UMI-weighted expected cells (default: true)
    // If false, allocate expected cells equally across tags
    bool useUMIWeights = true;
    // Testing flag: disable occupancy filter for testing
        bool disableOccupancyFilter = false;
        // Optional cap on ED candidates (default 120000 matches R's n.expected.cells)
        uint32_t edRetainCount = 120000;
        // Low UMI threshold for filtering (legacy, not currently used)
        uint32_t lowUMIThreshold = 0;
        
        // Simple EmptyDrops fallback thresholds
        // If ANY of these conditions are met, Simple EmptyDrops runs as fallback
        uint32_t simpleEDMinRescues = 50;      // Min ED rescues before fallback triggers
        uint32_t simpleEDMinAmbient = 100;     // Min ambient cells (UMI <= ambientUmiMax)
        uint32_t simpleEDMinCandidates = 100;  // Min candidates above lowerTestingBound
        bool useSimpleEmptyDrops = false;      // Force enable Simple EmptyDrops (--use-simple-empty-drops)
        
        // Debug flags
        bool debugTagLog = false;
        string debugOutputDir;
        
        // Invariant check flag (for testing/debugging)
        bool enableInvariantChecks = false;  // Enable EmptyDrops invariant checks (default: off for production)
        
        // Output options
        bool keepCBTag = false;  // If true, keep full CB+TAG barcodes in per-sample MEX (default: strip to 16bp)
    };
    
    // Type alias for backwards compatibility
    using OrdMagParams = SimpleEmptyDropsParams;
    
    // Helper: Get external defaults for Simple EmptyDrops parameters
    // These match the defaults from /mnt/pikachu/flex/scripts/run_ordmag_nonambient.py
    static SimpleEmptyDropsParams getExternalSimpleEmptyDropsDefaults();
    
    // Helper: Get external defaults for EmptyDrops parameters
    // These match the defaults from /mnt/pikachu/flex/scripts/run_ordmag_nonambient.py and run_emptydrops_cr_canonical.R
    static EmptyDropsParams getExternalEmptyDropsDefaults();
    
    // Helper: Populate config with external defaults for any zero/unset fields
    // Only fills fields that are zero or unset, preserving explicit overrides
    static void populateConfigWithDefaults(Config& config);
    
    // File-based entry point
    // Returns 0 on success, non-zero on error
    int runFromFiles(const FileInputs& inputs, Outputs* outputs, const Config& config);
    
    // In-memory entry point
    // Returns 0 on success, non-zero on error
    int runFromMemory(const MemoryInputs& inputs, Outputs* outputs, const Config& config);

    // Helper: Auto-derive sample whitelist from composite barcodes
    static void deriveSampleWhitelist(
        const std::vector<std::string>& compositeBarcodes,
        std::vector<std::string>& sampleLabels,
        std::vector<std::string>& sampleTags);
    
private:
    // Internal implementation (shared by both entry points)
    int runInternal(
        const SampleMatrixData& matrixData,
        const std::vector<std::string>& sampleLabels,
        const std::vector<std::string>& sampleTags,
        const std::map<std::string, uint32_t>& tagExpectedCells,
        Outputs* outputs,
        const Config& config);
    
    // Helper: Compute per-tag expected cells allocation
    // If tagExpectedCells override provided, use that; otherwise compute from matrix
    std::map<std::string, uint32_t> computeTagAllocations(
        const SampleMatrixData& matrixData,
        const std::vector<std::string>& sampleTags,
        const std::map<std::string, uint32_t>& override,
        const SimpleEmptyDropsParams& simpleEDParams,
        bool useNcellsSimpleWeights = true,
        uint32_t totalExpectedCells = 0);
    
    // Helper: Filter matrix to tag-specific barcodes
    SampleMatrixData filterMatrixToTag(
        const SampleMatrixData& fullMatrix,
        const std::string& tag);
    
    // Helper: Extract tag from composite barcode
    std::string extractTag(const std::string& compositeBarcode);
};

#endif
