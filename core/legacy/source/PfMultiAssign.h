#ifndef PF_MULTI_ASSIGN_H
#define PF_MULTI_ASSIGN_H

#include "IncludeDefine.h"
#include "PfMultiConfig.h"

/**
 * @file PfMultiAssign.h
 * @brief Execute feature-barcode assignment via in-process pf_api calls.
 */

namespace PfMultiAssign {

struct WhitelistNormalizationResult {
    string normalizedPath;
    bool hasTwoColumnSource = false;
    string assignmentNamespace = "UNKNOWN"; // NXT | TRU | UNKNOWN
    string sourcePath;
    uint64 normalizedRowCount = 0;
    bool namespaceConfidence = false;
};

struct AssignOptions {
    int maxHammingDistance = -1;
    int featureConstantOffset = -1;
    int limitSearch = -2; // -2 = unset; -1 is valid ("search entire read")
    int minCounts = -1;
    int maxBarcodeMismatches = -1;
    int featureN = -1;
    int barcodeN = -1;
    int consumerThreadsPerSet = -1;
    int searchThreads = -1;
    double minPosterior = -1.0;
    long long maxReads = -1; // <=0 means unlimited
    bool legacyCbRescue = false;
    bool enableStarDynamicPermitHooks = false;
    bool translateNxt = false;
    string filteredBarcodesPath;
    bool autodetectChemistry = false;
    int autodetectChemistryReads = 10000;
    int autodetectChemistryMinHits = 50;
    bool probeOnly = false;
    bool skipQcOutputs = false;
    bool allowUnionWhitelist = false;
    string sourceNamespace = "UNKNOWN";
    string targetNamespace = "UNKNOWN";
};

struct AssignResult {
    int returnCode = 0;
    string detectedMatchMode = "UNKNOWN";
    WhitelistNormalizationResult whitelistNormalization;
};

/**
 * @brief Normalize assignment whitelist and infer assignment namespace metadata.
 *
 * For 2-column translation lists, this writes column-1 barcodes to a normalized
 * whitelist file in assignOut and reports assignment namespace metadata.
 */
WhitelistNormalizationResult normalizeWhitelistForAssign(const string& whitelistPath,
                                                         const string& assignOut);

/**
 * @brief Run assignBarcodes for a feature library
 * @param whitelist Cell barcode whitelist file
 * @param featureRef Feature reference CSV file
 * @param fastqDir FASTQ directory for the feature library
 * @param assignOut Output directory for assignBarcodes
 * @param options Optional assignBarcodes CLI overrides
 * @return AssignResult with returnCode and detectedMatchMode
 */
AssignResult runAssignBarcodes(const string& whitelist,
                     const string& featureRef, const string& fastqDir,
                     const string& assignOut,
                     const AssignOptions& options = AssignOptions());

/**
 * @brief Process all feature libraries from config
 * @param config Parsed multi config
 * @param whitelist Cell barcode whitelist (from config or override)
 * @param featureRef Feature reference (from config or override)
 * @param fastqMap FASTQ path mapping
 * @param fastqRoot FASTQ root directory
 * @param outPrefix Output prefix (for creating subdirectories)
 * @param featureTypes List of feature types to process (e.g., {"CRISPR Guide Capture", "Antibody Capture"})
 * @param options Optional assignBarcodes CLI overrides
 * @return 0 on success, non-zero on error
 */
int processFeatureLibraries(const PfMultiConfig::Config& config,
                            const string& whitelist,
                            const string& featureRef, const map<string, string>& fastqMap,
                            const string& fastqRoot, const string& outPrefix,
                            const vector<string>& featureTypes,
                            const AssignOptions& options = AssignOptions());

} // namespace PfMultiAssign

#endif // PF_MULTI_ASSIGN_H
