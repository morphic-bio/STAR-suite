#ifndef PF_MULTI_CONFIG_H
#define PF_MULTI_CONFIG_H

#include "IncludeDefine.h"
#include "pf_split_read.h"
#include <map>
#include <vector>

/**
 * @file PfMultiConfig.h
 * @brief Parser for pf-multi config CSV files (Cell Ranger-style CSV input)
 * 
 * Parses multi config files with [libraries] section and other sections
 * like [feature] and [reference].
 */

namespace PfMultiConfig {

/**
 * @struct LibraryEntry
 * @brief Single library entry from [libraries] section
 */
/**
 * @struct SampleEntry
 * @brief OCM multiplexing sample from [samples] section
 */
struct SampleEntry {
    string sample_id;
    string ocm_barcode_ids;     // OB1, OB2, or pipe-union OB1|OB2
    string description;

    vector<string> resolvedOcmIds() const;
};

struct LibraryEntry {
    string fastqs;              // FASTQ directory path
    string feature_types;       // Feature type (e.g., "Gene Expression", "CRISPR Guide Capture")
    string sample;              // Sample name (optional)
    string library_type;        // Library type (optional)
    string gem_well;            // GEM well suffix (e.g., "1", "2"), defaults to "1" if absent
    string starChemistry;       // Per-library chemistry override: "TRU", "NXT", "auto", or empty (use global)
    string starWhitelist;       // Per-library barcode whitelist override (defaults to global --crWhitelist)
    string starFeatureRef;      // Per-library feature reference CSV path (skip global filtering when set)
    string starLibraryId;       // Stable output/provenance key (auto-generated if absent)
    string starInputFormat;     // fastq (default) or table (resolved count table in fastqs column)
    int starMaxHamming = -1;    // Per-library max Hamming distance override (-1 = use global)

    // Hash / HTO / CMO demux (assignBarcodes adt_mex extension)
    string starHashDemux;              // yes | no | auto | empty (auto)
    string starHashFeatureSelector;    // feature_type:HTO, id_prefix:hashtag, ...
    string starHashDemuxMethod;        // ratio (default)
    int starHashMinTotal = -1;
    int starHashMinTop = -1;
    double starHashMinRatio = -1.0;

    bool isTableBacked() const;

    // Split-read guide layout (CAT-ATAC and similar assays)
    string starLayout;              // Named preset, e.g. catatac_guide
    string starBarcodeRead;         // R1, R2, or R3
    string starBarcodeFormat;       // Chromap-style, e.g. bc:8:23:-
    string starUmiRead;             // R1, R2, or R3
    int starUmiStart = -1;
    int starUmiLength = -1;
    string starFeatureRead;         // R1, R2, or R3
    string starCaptureRead;         // R1, R2, or R3
    string starCaptureSequences;    // Pipe-separated capture kmers
    int starCaptureMaxHamming = 0;
    string starBarcodeOutputMap;    // Two-column ATAC->GEX map for MEX output
    string starFeatureSearchMode;   // free | anchor

    bool hasSplitReadLayout() const;
    
    // Normalized feature type for matching
    string normalizedFeatureType() const;
};

/**
 * @struct Config
 * @brief Parsed multi config structure
 */
struct Config {
    vector<LibraryEntry> libraries;    // All library entries
    vector<SampleEntry> samples;         // OCM samples from [samples] section
    string featureRef;                  // Feature reference path from [feature] section
    string referencePath;               // Reference path from [reference] section
    
    bool hasOcmSamples() const { return !samples.empty(); }

    // Classify libraries by feature type
    vector<LibraryEntry> getGexLibraries() const;
    vector<LibraryEntry> getFeatureLibraries(const string& featureType) const;
};

/**
 * @brief Parse multi config file
 * @param configPath Path to multi config CSV file
 * @return Parsed config structure
 * @throws runtime_error on parse errors
 */
Config parseConfig(const string& configPath);

/**
 * @brief Resolve FASTQ directory path using fastq_map and fastq_root
 * @param configPath Path from config file
 * @param fastqRoot Fallback root directory
 * @param fastqMap Map of config paths to actual paths
 * @return Resolved path
 */
string resolveFastqDir(const string& configPath, const string& fastqRoot, 
                        const map<string, string>& fastqMap);

/**
 * @brief Parse fastq_map vector (key=value pairs) into map
 * @param fastqMapVec Vector of "key=value" strings
 * @return Map of config paths to actual paths
 */
map<string, string> parseFastqMap(const vector<string>& fastqMapVec);

/**
 * @brief Resolve read token (R1/R2/R3) to 0-based index.
 */
int parseReadIndexToken(const string& token);

/**
 * @brief Apply named layout presets and validate split-read fields.
 */
void finalizeSplitReadLayout(LibraryEntry& entry);

pf_split_read_layout buildSplitReadLayout(const LibraryEntry& entry);

} // namespace PfMultiConfig

#endif // PF_MULTI_CONFIG_H
