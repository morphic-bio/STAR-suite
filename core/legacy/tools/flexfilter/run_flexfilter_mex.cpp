/**
 * run_flexfilter_mex: Standalone wrapper for FlexFilter from composite MEX
 * 
 * Reads a composite MEX directory (matrix.mtx, barcodes.tsv with CB+TAG, features.tsv)
 * and invokes the libflexfilter pipeline to produce per-sample split MEX outputs.
 * 
 * Usage:
 *   run_flexfilter_mex \
 *     --mex-dir /path/to/composite_mex \
 *     --allowed-tags allowed_tags.tsv \
 *     --total-expected 12000 \
 *     --output-prefix /path/to/output \
 *     [--cells-per-tag-json cells_per_tag.json] \
 *     [--emptydrops-params ...]
 */

#include "libflex/FlexFilter.h"
#include "libflex/FlexFilterAdapters.h"
#include "MexWriter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>  // For symlink/unlink

// Forward declarations for POSIX functions
extern "C" {
    int symlink(const char *target, const char *linkpath);
    int unlink(const char *pathname);
}

// Command-line arguments
struct Arguments {
    std::string mexDir;
    std::string allowedTagsPath;
    std::string sampleWhitelistPath;
    std::string cellsPerTagJson;
    std::string outputPrefix;
    uint32_t totalExpected = 0;
    
    // FlexFilter config
    FlexFilter::Config config;
    
    // Output barcode handling
    bool keepCBTag = false;  // If true, keep full CB+TAG barcodes; if false (default), strip to 16bp
    
    bool valid = false;
};

// Helper: Print usage
void printUsage(const char* progName) {
    std::cerr << "Usage: " << progName << " [OPTIONS]\n\n";
    std::cerr << "Required:\n";
    std::cerr << "  --mex-dir <path>           Path to composite MEX directory\n";
    std::cerr << "  --output-prefix <path>     Output directory prefix\n";
    std::cerr << "  --total-expected <N>       Total expected cells across all tags\n";
    std::cerr << "\nOptional:\n";
    std::cerr << "  --allowed-tags <path>      TSV file with allowed tags (tag sequences)\n";
    std::cerr << "  --sample-whitelist <path>  TSV file: sample_name<TAB>tag_sequence\n";
    std::cerr << "  --cells-per-tag-json <path> JSON file with per-tag cell counts\n";
    std::cerr << "\nTag allocation:\n";
    std::cerr << "  --simple-ed-expected-cells <N> Expected cells for allocation (default: 3000)\n";
    std::cerr << "  --no-tag-weights-umi        Disable UMI weighting, use equal allocation (default: UMI-weighted)\n";
    std::cerr << "\nEmptyDrops parameters:\n";
    std::cerr << "  --ed-sim-n <N>              Monte Carlo sims (default: 10000)\n";
    std::cerr << "  --ed-fdr <F>                FDR threshold (default: 0.01)\n";
    std::cerr << "  --ed-use-fdr                Use FDR gate\n";
    std::cerr << "  --ed-use-rawp               Use raw p-value gate (default)\n";
    std::cerr << "  --ed-rawp-threshold <F>     Raw p-value threshold (default: 0.01)\n";
    std::cerr << "  --ed-lower-bound <N>        Lower UMI bound (cells <= N excluded; default: 500)\n";
    std::cerr << "  --ed-ambient-umi-max <N>    Max UMI for ambient (default: 100)\n";
    std::cerr << "  --ed-retain-count <N>       Retain window size (default: 120000)\n";
    std::cerr << "\nSimple EmptyDrops (fallback filter):\n";
    std::cerr << "  --use-simple-empty-drops   Force enable Simple EmptyDrops (by default only runs as fallback)\n";
    std::cerr << "  --simple-ed-min-rescues <N> Min ED rescues before fallback triggers (default: 50)\n";
    std::cerr << "  --simple-ed-min-ambient <N> Min ambient cells before fallback triggers (default: 100)\n";
    std::cerr << "  --simple-ed-min-candidates <N> Min candidates before fallback triggers (default: 100)\n";
    std::cerr << "  --bootstrap-seed <N>       Seed for bootstrap RNG (default: 0 = deterministic)\n";
    std::cerr << "\nOccupancy post-filter parameters:\n";
    std::cerr << "  --total-partitions <N>     Total partitions (default: 115000)\n";
    std::cerr << "  --recovery-factor <F>      Recovery factor (default: 0.606)\n";
    std::cerr << "  --occupancy-percentile <F> Percentile (default: 0.999)\n";
    std::cerr << "  --occupancy-sim-gems <N>   Monte Carlo GEMs (default: 1000000)\n";
    std::cerr << "\nOutput options:\n";
    std::cerr << "  --keep-cb-tag              Keep full CB+TAG barcodes in output (default: strip to 16bp)\n";
    std::cerr << "\nFlags:\n";
    std::cerr << "  --disable-occupancy        Disable occupancy filter (for testing)\n";
    std::cerr << "  --debug-tag-log            Enable detailed logging\n";
    std::cerr << "  --debug-output-dir <path>  Debug output directory\n";
    std::cerr << "  --help                     Show this help\n";
}

// Helper: Parse command-line arguments
Arguments parseArguments(int argc, char** argv) {
    Arguments args;
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            exit(0);
        } else if (arg == "--mex-dir" && i + 1 < argc) {
            args.mexDir = argv[++i];
        } else if (arg == "--allowed-tags" && i + 1 < argc) {
            args.allowedTagsPath = argv[++i];
        } else if (arg == "--sample-whitelist" && i + 1 < argc) {
            args.sampleWhitelistPath = argv[++i];
        } else if (arg == "--cells-per-tag-json" && i + 1 < argc) {
            args.cellsPerTagJson = argv[++i];
        } else if (arg == "--output-prefix" && i + 1 < argc) {
            args.outputPrefix = argv[++i];
        } else if (arg == "--total-expected" && i + 1 < argc) {
            args.totalExpected = std::atoi(argv[++i]);
        } else if ((arg == "--simple-ed-expected-cells" || arg == "--ordmag-expected-cells") && i + 1 < argc) {
            args.config.simpleEmptyDropsParams.nExpectedCells = std::atoi(argv[++i]);
        } else if (arg == "--no-tag-weights-umi") {
            args.config.useUMIWeights = false;
        } else if (arg == "--ed-sim-n" && i + 1 < argc) {
            args.config.emptydropsParams.simN = std::atoi(argv[++i]);
        } else if (arg == "--ed-fdr" && i + 1 < argc) {
            args.config.emptydropsParams.FDR = std::atof(argv[++i]);
        } else if (arg == "--ed-use-fdr") {
            args.config.emptydropsParams.useFDRGate = true;
        } else if (arg == "--ed-use-rawp") {
            args.config.emptydropsParams.useFDRGate = false;
        } else if (arg == "--ed-rawp-threshold" && i + 1 < argc) {
            args.config.emptydropsParams.rawPvalueThreshold = std::atof(argv[++i]);
        } else if (arg == "--ed-lower-bound" && i + 1 < argc) {
            args.config.emptydropsParams.lowerTestingBound = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (arg == "--ed-ambient-umi-max" && i + 1 < argc) {
            args.config.emptydropsParams.ambientUmiMax = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (arg == "--ed-retain-count" && i + 1 < argc) {
            args.config.edRetainCount = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (arg == "--total-partitions" && i + 1 < argc) {
            args.config.totalPartitions = std::atoi(argv[++i]);
        } else if (arg == "--recovery-factor" && i + 1 < argc) {
            args.config.recoveryFactor = std::atof(argv[++i]);
        } else if (arg == "--occupancy-percentile" && i + 1 < argc) {
            args.config.occupancyPercentile = std::atof(argv[++i]);
        } else if (arg == "--occupancy-sim-gems" && i + 1 < argc) {
            args.config.occupancySimulatedGems = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (arg == "--disable-occupancy") {
            args.config.disableOccupancyFilter = true;
        } else if (arg == "--keep-cb-tag") {
            args.keepCBTag = true;
        } else if (arg == "--use-simple-empty-drops") {
            args.config.useSimpleEmptyDrops = true;
            args.config.simpleEmptyDropsParams.disabled = false;
        } else if (arg == "--simple-ed-min-rescues" && i + 1 < argc) {
            args.config.simpleEDMinRescues = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (arg == "--simple-ed-min-ambient" && i + 1 < argc) {
            args.config.simpleEDMinAmbient = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (arg == "--simple-ed-min-candidates" && i + 1 < argc) {
            args.config.simpleEDMinCandidates = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (arg == "--bootstrap-seed" && i + 1 < argc) {
            args.config.simpleEmptyDropsParams.bootstrapSeed = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
        } else if (arg == "--debug-tag-log") {
            args.config.debugTagLog = true;
        } else if (arg == "--debug-output-dir" && i + 1 < argc) {
            args.config.debugOutputDir = argv[++i];
        } else {
            std::cerr << "ERROR: Unknown argument: " << arg << "\n";
            printUsage(argv[0]);
            return args;
        }
    }
    
    // Validate required arguments
    if (args.mexDir.empty()) {
        std::cerr << "ERROR: --mex-dir is required\n";
        printUsage(argv[0]);
        return args;
    }
    
    if (args.outputPrefix.empty()) {
        std::cerr << "ERROR: --output-prefix is required\n";
        printUsage(argv[0]);
        return args;
    }
    
    if (args.totalExpected == 0) {
        std::cerr << "ERROR: --total-expected is required\n";
        printUsage(argv[0]);
        return args;
    }
    
    // Populate defaults for any unset parameters
    FlexFilter::populateConfigWithDefaults(args.config);
    
    args.valid = true;
    return args;
}

// Helper: Load allowed tags from TSV
std::vector<std::string> loadAllowedTags(const std::string& path) {
    std::vector<std::string> tags;
    
    if (path.empty()) {
        return tags;  // No filter
    }
    
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "WARNING: Could not open allowed tags file: " << path << "\n";
        return tags;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (!line.empty() && line.size() >= 8) {
            // Extract first 8 characters as tag sequence
            tags.push_back(line.substr(0, 8));
        }
    }
    
    file.close();
    return tags;
}

// Helper: Load sample whitelist from TSV (sample_name<TAB>tag_sequence)
bool loadSampleWhitelist(
    const std::string& path,
    std::vector<std::string>& sampleLabels,
    std::vector<std::string>& sampleTags)
{
    if (path.empty()) return false;
    
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open sample whitelist file: " << path << "\n";
        return false;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Parse tab-separated fields
        size_t tabPos = line.find('\t');
        if (tabPos == std::string::npos) {
            std::cerr << "WARNING: Invalid line in sample whitelist (missing tab): " << line << "\n";
            continue;
        }
        
        std::string label = line.substr(0, tabPos);
        std::string tag = line.substr(tabPos + 1);
        
        // Trim fields
        label.erase(0, label.find_first_not_of(" \t"));
        label.erase(label.find_last_not_of(" \t") + 1);
        tag.erase(0, tag.find_first_not_of(" \t"));
        tag.erase(tag.find_last_not_of(" \t") + 1);
        
        if (!label.empty() && !tag.empty()) {
            sampleLabels.push_back(label);
            sampleTags.push_back(tag);
        }
    }
    
    file.close();
    return !sampleLabels.empty();
}

// Helper: Load per-tag cell counts from JSON
std::map<std::string, uint32_t> loadCellsPerTagJson(const std::string& path) {
    std::map<std::string, uint32_t> tagCounts;
    
    if (path.empty()) {
        return tagCounts;  // No override
    }
    
    // TODO: Implement JSON parsing
    // For now, just warn
    std::cerr << "WARNING: --cells-per-tag-json not yet implemented\n";
    
    return tagCounts;
}

int main(int argc, char** argv) {
    // Parse arguments
    Arguments args = parseArguments(argc, argv);
    if (!args.valid) {
        return 1;
    }

    // Emit a run timestamp to confirm binary/params used
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::cout << "run_flexfilter_mex run timestamp: "
              << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << "\n";
    
    std::cout << "run_flexfilter_mex: Starting\n";
    std::cout << "  MEX directory: " << args.mexDir << "\n";
    std::cout << "  Output prefix: " << args.outputPrefix << "\n";
    std::cout << "  Total expected cells: " << args.totalExpected << "\n";
    
    // Load sample whitelist or allowed tags (optional)
    std::vector<std::string> sampleLabels;
    std::vector<std::string> sampleTags;
    std::vector<std::string> allowedTags;
    
    bool hasExplicitSamples = false;
    if (!args.sampleWhitelistPath.empty()) {
        hasExplicitSamples = loadSampleWhitelist(args.sampleWhitelistPath, sampleLabels, sampleTags);
        if (hasExplicitSamples) {
            std::cout << "  Loaded " << sampleLabels.size() << " samples from whitelist:\n";
            for (size_t i = 0; i < sampleLabels.size(); i++) {
                std::cout << "    " << sampleLabels[i] << ": " << sampleTags[i] << "\n";
            }
        }
    } else {
        allowedTags = loadAllowedTags(args.allowedTagsPath);
        if (!args.allowedTagsPath.empty()) {
            std::cout << "  Loaded " << allowedTags.size() << " allowed tags\n";
        }
    }
    
    // Load per-tag cell counts (optional)
    std::map<std::string, uint32_t> tagExpectedCells = loadCellsPerTagJson(args.cellsPerTagJson);
    
    // Load MEX files first (needed for output writing)
    std::vector<std::string> barcodes;
    std::vector<std::string> features;
    SampleMatrixData matrixData;
    
    // Build FlexFilter inputs
    FlexFilter::FileInputs inputs;
    inputs.matrixDir = args.mexDir;
    
    // Use explicit sample whitelist if provided, otherwise use allowed tags
    if (hasExplicitSamples) {
        inputs.sampleLabels = sampleLabels;
        inputs.sampleTags = sampleTags;
    } else {
        inputs.allowedTags = allowedTags;
    }
    inputs.tagExpectedCells = tagExpectedCells;
    
    // Run FlexFilter
    FlexFilter filter;
    FlexFilter::Outputs outputs;

    // Set totalExpectedCells in config
    args.config.totalExpectedCells = args.totalExpected;
    // If debug tag log is requested but no debug dir provided, default under output prefix
    if (args.config.debugTagLog && args.config.debugOutputDir.empty()) {
        args.config.debugOutputDir = args.outputPrefix + "/debug";
    }
    
    // Verify MEX files exist
    std::string matrixPath = args.mexDir + "/matrix.mtx";
    std::string barcodesPath = args.mexDir + "/barcodes.tsv";
    std::string featuresPath = args.mexDir + "/features.tsv";
    
    // Also check for InlineHashDedup_ prefix
    bool hasPrefix = false;
    std::ifstream testMatrix(matrixPath);
    if (!testMatrix.is_open()) {
        matrixPath = args.mexDir + "/InlineHashDedup_matrix.mtx";
        barcodesPath = args.mexDir + "/InlineHashDedup_barcodes.tsv";
        featuresPath = args.mexDir + "/InlineHashDedup_features.tsv";
        hasPrefix = true;
    }
    testMatrix.close();
    
    // Verify all files exist
    std::ifstream checkMatrix(matrixPath);
    std::ifstream checkBarcodes(barcodesPath);
    std::ifstream checkFeatures(featuresPath);
    
    if (!checkMatrix.is_open()) {
        std::cerr << "ERROR: Cannot open matrix file: " << matrixPath << "\n";
        std::cerr << "Expected files in " << args.mexDir << ":\n";
        std::cerr << "  - matrix.mtx (or InlineHashDedup_matrix.mtx)\n";
        std::cerr << "  - barcodes.tsv (or InlineHashDedup_barcodes.tsv)\n";
        std::cerr << "  - features.tsv (or InlineHashDedup_features.tsv)\n";
        return 1;
    }
    if (!checkBarcodes.is_open()) {
        std::cerr << "ERROR: Cannot open barcodes file: " << barcodesPath << "\n";
        return 1;
    }
    if (!checkFeatures.is_open()) {
        std::cerr << "ERROR: Cannot open features file: " << featuresPath << "\n";
        return 1;
    }
    
    checkMatrix.close();
    checkBarcodes.close();
    checkFeatures.close();
    
    std::cout << "Found MEX files:\n";
    std::cout << "  Matrix: " << matrixPath << "\n";
    std::cout << "  Barcodes: " << barcodesPath << "\n";
    std::cout << "  Features: " << featuresPath << "\n";
    
    // If files have prefix, create symlinks without prefix for FlexFilter
    if (hasPrefix) {
        std::cout << "Creating temporary symlinks for FlexFilter...\n";
        std::string symlinkMatrix = args.mexDir + "/matrix.mtx";
        std::string symlinkBarcodes = args.mexDir + "/barcodes.tsv";
        std::string symlinkFeatures = args.mexDir + "/features.tsv";
        
        // Remove old symlinks if they exist
        unlink(symlinkMatrix.c_str());
        unlink(symlinkBarcodes.c_str());
        unlink(symlinkFeatures.c_str());
        
        // Create new symlinks
        if (symlink("InlineHashDedup_matrix.mtx", symlinkMatrix.c_str()) != 0) {
            std::cerr << "WARNING: Could not create symlink for matrix.mtx\n";
        }
        if (symlink("InlineHashDedup_barcodes.tsv", symlinkBarcodes.c_str()) != 0) {
            std::cerr << "WARNING: Could not create symlink for barcodes.tsv\n";
        }
        if (symlink("InlineHashDedup_features.tsv", symlinkFeatures.c_str()) != 0) {
            std::cerr << "WARNING: Could not create symlink for features.tsv\n";
        }
    }
    
    std::cout << "Running FlexFilter pipeline...\n";
    
    int result = filter.runFromFiles(inputs, &outputs, args.config);
    
    // Clean up symlinks
    if (hasPrefix) {
        std::string symlinkMatrix = args.mexDir + "/matrix.mtx";
        std::string symlinkBarcodes = args.mexDir + "/barcodes.tsv";
        std::string symlinkFeatures = args.mexDir + "/features.tsv";
        unlink(symlinkMatrix.c_str());
        unlink(symlinkBarcodes.c_str());
        unlink(symlinkFeatures.c_str());
    }
    
    if (result != 0) {
        std::cerr << "ERROR: FlexFilter pipeline failed with code: " << result << "\n";
        std::cerr << "Check that MEX files are valid and composite barcodes are 24 characters.\n";
        return result;
    }
    
    std::cout << "FlexFilter completed successfully\n";
    std::cout << "  Processed " << outputs.tagResults.size() << " tags\n";
    
    // Note: runFromFiles internally loads MEX; load again for output writing after freeing its memory
    std::cout << "Loading MEX for output generation...\n";
    
    // Load full MEX including matrix (needed for output filtering)
    {
        std::string barcodesPath = args.mexDir + "/barcodes.tsv";
        std::string featuresPath = args.mexDir + "/features.tsv";
        std::string matrixPath = args.mexDir + "/matrix.mtx";
        
        // Check for InlineHashDedup prefix
        std::ifstream testBc(barcodesPath);
        if (!testBc.is_open()) {
            barcodesPath = args.mexDir + "/InlineHashDedup_barcodes.tsv";
            featuresPath = args.mexDir + "/InlineHashDedup_features.tsv";
            matrixPath = args.mexDir + "/InlineHashDedup_matrix.mtx";
        }
        testBc.close();
        
        // Load barcodes
        std::ifstream bcFile(barcodesPath);
        std::string line;
        while (std::getline(bcFile, line)) {
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            if (!line.empty()) barcodes.push_back(line);
        }
        bcFile.close();
        
        // Load features (gene IDs only)
        std::ifstream ftFile(featuresPath);
        while (std::getline(ftFile, line)) {
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);
            if (!line.empty()) {
                size_t tabPos = line.find('\t');
                std::string geneId = (tabPos != std::string::npos) ? line.substr(0, tabPos) : line;
                features.push_back(geneId);
            }
        }
        ftFile.close();
        
        // Load matrix
        std::ifstream mtxFile(matrixPath);
        if (!mtxFile.is_open()) {
            std::cerr << "ERROR: Cannot open matrix file: " << matrixPath << "\n";
            return 1;
        }
        
        // Skip header
        while (std::getline(mtxFile, line)) {
            if (line.empty() || line[0] != '%') break;
        }
        
        // Parse dimensions
        uint32_t numGenes, numCells, numNonzero;
        std::stringstream dimStream(line);
        dimStream >> numGenes >> numCells >> numNonzero;
        
        // Initialize matrix
        matrixData.nCells = numCells;
        matrixData.nGenes = numGenes;
        matrixData.countMatStride = 3;
        matrixData.nUMIperCB.resize(numCells, 0);
        matrixData.nGenePerCB.resize(numCells, 0);
        matrixData.countCellGeneUMIindex.resize(numCells + 1, 0);
        
        // Load all triplets into temp storage
        std::vector<std::vector<std::pair<uint32_t, uint32_t>>> cellGeneMap(numCells);
        
        uint32_t geneRow, cellCol, count;
        while (mtxFile >> geneRow >> cellCol >> count) {
            uint32_t geneIdx = geneRow - 1;  // Convert to 0-based
            uint32_t cellIdx = cellCol - 1;
            
            if (cellIdx < numCells && geneIdx < numGenes) {
                cellGeneMap[cellIdx].push_back({geneIdx, count});
                matrixData.nUMIperCB[cellIdx] += count;
            }
        }
        mtxFile.close();
        
        // Build sparse matrix
        uint32_t offset = 0;
        for (uint32_t cellIdx = 0; cellIdx < numCells; cellIdx++) {
            matrixData.countCellGeneUMIindex[cellIdx] = offset;
            matrixData.nGenePerCB[cellIdx] = cellGeneMap[cellIdx].size();
            
            for (const auto& gc : cellGeneMap[cellIdx]) {
                matrixData.countCellGeneUMI.push_back(gc.first);   // geneIdx
                matrixData.countCellGeneUMI.push_back(gc.second);  // count
                matrixData.countCellGeneUMI.push_back(0);          // reserved
                offset += 3;
            }
        }
        matrixData.countCellGeneUMIindex[numCells] = offset;
        
        std::cout << "  Loaded " << barcodes.size() << " barcodes, " 
                  << features.size() << " features, " 
                  << numNonzero << " entries\n";
    }
    
    // Write outputs using our own MEX writer
    std::cout << "Writing per-sample MEX outputs...\n";
    
    // Create output directories and write MEX for each tag
    for (const auto& tagResult : outputs.tagResults) {
        // Create sample directory
        std::string sampleDir = args.outputPrefix + "/" + tagResult.sampleLabel + "/Gene/filtered";
        
        // Create directory structure (mkdir -p equivalent)
        std::string mkdirCmd = "mkdir -p \"" + sampleDir + "\"";
        if (system(mkdirCmd.c_str()) != 0) {
            std::cerr << "WARNING: Failed to create directory: " << sampleDir << "\n";
            continue;
        }
        
        // Filter matrix to passing barcodes
        // Build set of passing barcode indices (0-based)
        std::unordered_set<std::string> passingBarcodeSet(
            tagResult.passingBarcodes.begin(),
            tagResult.passingBarcodes.end()
        );
        
        // Map barcode strings to indices in original matrix
        std::map<std::string, uint32_t> barcodeToIdx;
        for (uint32_t i = 0; i < barcodes.size(); i++) {
            barcodeToIdx[barcodes[i]] = i;
        }
        
        // Build filtered barcode list and index mapping
        std::vector<std::string> filteredBarcodes;
        std::map<uint32_t, uint32_t> oldToNewIdx;  // old idx -> new idx
        
        for (const auto& bc : tagResult.passingBarcodes) {
            auto it = barcodeToIdx.find(bc);
            if (it != barcodeToIdx.end()) {
                uint32_t oldIdx = it->second;
                uint32_t newIdx = filteredBarcodes.size();
                oldToNewIdx[oldIdx] = newIdx;
                filteredBarcodes.push_back(bc);
            }
        }
        
        // Build filtered triplets
        std::vector<MexWriter::Triplet> filteredTriplets;
        
        // Iterate through original matrix and filter to passing barcodes
        for (uint32_t cellIdx = 0; cellIdx < matrixData.nCells; cellIdx++) {
            // Check if this cell is in passing set
            if (oldToNewIdx.find(cellIdx) == oldToNewIdx.end()) {
                continue;  // Skip non-passing cells
            }
            
            uint32_t newCellIdx = oldToNewIdx[cellIdx];
            uint32_t start = matrixData.countCellGeneUMIindex[cellIdx];
            uint32_t end = matrixData.countCellGeneUMIindex[cellIdx + 1];
            
            for (uint32_t i = start; i < end; i += matrixData.countMatStride) {
                uint32_t geneIdx = matrixData.countCellGeneUMI[i];
                uint32_t count = matrixData.countCellGeneUMI[i + 1];
                
                if (count > 0 && geneIdx < matrixData.nGenes) {
                    filteredTriplets.push_back({newCellIdx, geneIdx, count});
                }
            }
        }
        
        // Convert features to Feature format
        std::vector<MexWriter::Feature> featureList;
        featureList.reserve(features.size());
        for (const auto& geneId : features) {
            featureList.emplace_back(geneId, geneId, "Gene Expression");
        }
        
        // Write MEX files directly into the filtered directory using dir-mode prefix
        std::string outputPrefix = sampleDir;
        if (outputPrefix.back() != '/') {
            outputPrefix.push_back('/');
        }
        
        std::cout << "    Writing MEX to: " << outputPrefix << "\n";

        if (filteredBarcodes.empty()) {
            std::cout << "    Skipping " << tagResult.sampleLabel << " (no passing barcodes)\n";
            continue;
        }
        
        // Per-sample MEX: strip sample tag from barcodes (16bp output) unless --keep-cb-tag
        int cb_len = args.keepCBTag ? -1 : 16;
        int writeResult = MexWriter::writeMex(
            outputPrefix,
            filteredBarcodes,
            featureList,
            filteredTriplets,
            cb_len
        );
        
        if (writeResult != 0) {
            std::cerr << "  ERROR: MexWriter failed for " << tagResult.sampleLabel
                      << " (barcodes=" << filteredBarcodes.size()
                      << ", features=" << featureList.size()
                      << ", entries=" << filteredTriplets.size() << ")\n";
            continue;
        }
        
        if (writeResult == 0) {
            std::cout << "  " << tagResult.sampleLabel << " (" << tagResult.tag << "): "
                     << filteredBarcodes.size() << " cells, "
                     << filteredTriplets.size() << " entries\n";
        } else {
            std::cerr << "  ERROR: Failed to write MEX for " << tagResult.sampleLabel << "\n";
        }
    }
    
    std::cout << "Done! Outputs written to: " << args.outputPrefix << "\n";
    
    // Summary: cells after EmptyDrops pipeline stages
    auto sumUMIFromMatrix = [](const std::string& matrixPath) -> uint64_t {
        std::ifstream mtx(matrixPath);
        if (!mtx.is_open()) {
            return 0;
        }
        std::string line;
        // Skip comments
        while (std::getline(mtx, line)) {
            if (!line.empty() && line[0] != '%') {
                break; // dimension line read into 'line'
            }
        }
        // Sum counts
        uint64_t total = 0;
        uint32_t r, c, v;
        while (mtx >> r >> c >> v) {
            total += v;
        }
        return total;
    };

    std::string summaryPath = args.outputPrefix + "/flexfilter_summary.tsv";
    std::ofstream summaryFile(summaryPath);
    if (summaryFile.is_open()) {
        summaryFile << "Sample\tExpected\tRetain\tSimple\tTail_Tested\tSimple_Pass\tTail_Pass\tOcc_Removed\tFinal\tTotal_UMIs\n";
    }

    std::cout << "\nSummary (saved to " << summaryPath << "):\n";
    printf("%-12s %8s %8s %8s %10s %10s %10s %10s %8s %14s\n",
           "Sample", "Expected", "Retain", "Simple", "Tail_Test",
           "Simp_Pass", "Tail_Pass", "Occ_Rem", "Final", "Total_UMIs");

    uint32_t totalExpected = 0;
    uint32_t totalRetain = 0;
    uint32_t totalSimple = 0;
    uint32_t totalTailTested = 0;
    uint32_t totalSimplePass = 0;
    uint32_t totalTailPass = 0;
    uint32_t totalOccRemoved = 0;
    uint32_t totalFinal = 0;
    uint64_t totalUMIs = 0;

    for (const auto& tagResult : outputs.tagResults) {
        // Write EmptyDrops per-candidate diagnostics
        {
            std::string sampleDir = args.outputPrefix + "/" + tagResult.sampleLabel + "/Gene/filtered";
            std::string edDir = sampleDir + "/EmptyDrops";
            mkdir(edDir.c_str(), 0775);
            std::string edPath = edDir + "/emptydrops_results.tsv";
            std::ofstream edFile(edPath);
            if (edFile.is_open()) {
                edFile << "barcode\traw_p\tadj_p\tpasses_rawp\tpasses_fdr\n";
                for (const auto& res : tagResult.emptydropsResults) {
                    uint32_t idx = res.cellIndex;
                    // Use retainBarcodes (retain window) for cellIndex mapping, not tagBarcodes (all)
                    if (idx < tagResult.retainBarcodes.size()) {
                        std::string bc = tagResult.retainBarcodes[idx];
                        edFile << bc << '\t'
                               << res.pValue << '\t'
                               << res.pAdjusted << '\t'
                               << (res.passesRawP ? 1 : 0) << '\t'
                               << (res.passesFDR ? 1 : 0) << '\n';
                    }
                }
            }
        }

        std::string mtxPath = args.outputPrefix + "/" + tagResult.sampleLabel + "/Gene/filtered/matrix.mtx";
        uint64_t totalUMI = sumUMIFromMatrix(mtxPath);
        uint32_t finalCells = static_cast<uint32_t>(tagResult.passingBarcodes.size());

        printf("%-12s %8u %8u %8u %10u %10u %10u %10u %8u %14lu\n",
               tagResult.sampleLabel.c_str(),
               tagResult.expectedCells,
               tagResult.nRetainWindow,
               tagResult.nSimpleCells,
               tagResult.nTailTested,
               tagResult.nSimplePassers,
               tagResult.nTailPassers,
               tagResult.occupancyRemoved,
               finalCells,
               totalUMI);

        if (summaryFile.is_open()) {
            summaryFile << tagResult.sampleLabel << '\t'
                        << tagResult.expectedCells << '\t'
                        << tagResult.nRetainWindow << '\t'
                        << tagResult.nSimpleCells << '\t'
                        << tagResult.nTailTested << '\t'
                        << tagResult.nSimplePassers << '\t'
                        << tagResult.nTailPassers << '\t'
                        << tagResult.occupancyRemoved << '\t'
                        << finalCells << '\t'
                        << totalUMI << '\n';
        }

        totalExpected += tagResult.expectedCells;
        totalRetain += tagResult.nRetainWindow;
        totalSimple += tagResult.nSimpleCells;
        totalTailTested += tagResult.nTailTested;
        totalSimplePass += tagResult.nSimplePassers;
        totalTailPass += tagResult.nTailPassers;
        totalOccRemoved += tagResult.occupancyRemoved;
        totalFinal += finalCells;
        totalUMIs += totalUMI;
    }

    // Totals row
    printf("%-12s %8u %8u %8u %10u %10u %10u %10u %8u %14lu\n",
           "TOTAL",
           totalExpected,
           totalRetain,
           totalSimple,
           totalTailTested,
           totalSimplePass,
           totalTailPass,
           totalOccRemoved,
           totalFinal,
           totalUMIs);

    if (summaryFile.is_open()) {
        summaryFile << "TOTAL\t"
                    << totalExpected << '\t'
                    << totalRetain << '\t'
                    << totalSimple << '\t'
                    << totalTailTested << '\t'
                    << totalSimplePass << '\t'
                    << totalTailPass << '\t'
                    << totalOccRemoved << '\t'
                    << totalFinal << '\t'
                    << totalUMIs << '\n';
        summaryFile.close();
    }
    
    return 0;
}
