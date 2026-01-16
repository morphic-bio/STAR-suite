#ifndef H_FlexFilterAdapters
#define H_FlexFilterAdapters

#include "IncludeDefine.h"
#include "FlexFilter.h"
#include <string>

// Output adapters for FlexFilter results
// Separates file I/O from core filtering logic

// STAR adapter: writes to Solo.out/sample_outs/<sample>/ structure
class STARFlexFilterAdapter {
public:
    // Write outputs to STAR directory structure
    // outputPrefix: base path (e.g., "Solo.out/sample_outs/<sample>/Gene/filtered/")
    // runDirPerm: directory permissions for output directories
    static void writeOutputs(
        const FlexFilter::Outputs& outputs,
        const std::string& outputPrefix,
        mode_t runDirPerm
    );
    
private:
    // Helper: Ensure directory exists
    static bool ensureDirectory(const std::string& path, mode_t mode);
    
    // Helper: Write OrdMag outputs
    static void writeOrdMagOutputs(
        const FlexFilter::Outputs::TagResults& tagResult,
        const std::string& outputPrefix
    );
    
    // Helper: Write EmptyDrops outputs
    static void writeEmptyDropsOutputs(
        const FlexFilter::Outputs::TagResults& tagResult,
        const std::string& outputPrefix
    );
    
    // Helper: Write MEX files
    static void writeMEXFiles(
        const FlexFilter::Outputs::TagResults& tagResult,
        const SampleMatrixData& matrixData,
        const std::string& outputPrefix
    );
    
    // Helper: Write aggregate filter_summary.json
    static void writeAggregateFilterSummary(
        const FlexFilter::Outputs& outputs,
        const std::string& outputPrefix
    );
};

// Standalone adapter: writes to legacy Flex layout
class StandaloneFlexFilterAdapter {
public:
    // Write outputs to standalone directory structure
    // outputPrefix: base path (e.g., "output/prefilter/")
    static void writeOutputs(
        const FlexFilter::Outputs& outputs,
        const std::string& outputPrefix
    );
    
private:
    // Helper: Ensure directory exists
    static bool ensureDirectory(const std::string& path);
    
    // Helper: Write per-tag filtered barcodes
    static void writePerTagFilteredBarcodes(
        const FlexFilter::Outputs::TagResults& tagResult,
        const std::string& outputPrefix
    );
    
    // Helper: Write EmptyDrops passing barcodes
    static void writeEmptyDropsPassing(
        const FlexFilter::Outputs::TagResults& tagResult,
        const std::string& outputPrefix
    );
    
    // Helper: Write per_tag_expected.tsv
    static void writePerTagExpected(
        const FlexFilter::Outputs& outputs,
        const std::string& outputPrefix
    );
    
    // Helper: Write aggregate filter_summary.json
    static void writeAggregateFilterSummary(
        const FlexFilter::Outputs& outputs,
        const std::string& outputPrefix
    );
};

#endif

