#ifndef CR_MULTI_CONFIG_H
#define CR_MULTI_CONFIG_H

#include "IncludeDefine.h"
#include <map>
#include <vector>

/**
 * @file CrMultiConfig.h
 * @brief Parser for Cell Ranger multi config CSV files
 * 
 * Parses multi config files with [libraries] section and other sections
 * like [feature] and [reference].
 */

namespace CrMultiConfig {

/**
 * @struct LibraryEntry
 * @brief Single library entry from [libraries] section
 */
struct LibraryEntry {
    string fastqs;              // FASTQ directory path
    string feature_types;       // Feature type (e.g., "Gene Expression", "CRISPR Guide Capture")
    string sample;              // Sample name (optional)
    string library_type;        // Library type (optional)
    string gem_well;            // GEM well suffix (e.g., "1", "2"), defaults to "1" if absent
    
    // Normalized feature type for matching
    string normalizedFeatureType() const;
};

/**
 * @struct Config
 * @brief Parsed multi config structure
 */
struct Config {
    vector<LibraryEntry> libraries;    // All library entries
    string featureRef;                  // Feature reference path from [feature] section
    string referencePath;               // Reference path from [reference] section
    
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

} // namespace CrMultiConfig

#endif // CR_MULTI_CONFIG_H
