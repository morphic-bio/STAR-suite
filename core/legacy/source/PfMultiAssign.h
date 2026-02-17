#ifndef PF_MULTI_ASSIGN_H
#define PF_MULTI_ASSIGN_H

#include "IncludeDefine.h"
#include "PfMultiConfig.h"

/**
 * @file PfMultiAssign.h
 * @brief Execute feature-barcode assignment via in-process pf_api calls.
 */

namespace PfMultiAssign {

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
    string filteredBarcodesPath;
};

/**
 * @brief Run assignBarcodes for a feature library
 * @param whitelist Cell barcode whitelist file
 * @param featureRef Feature reference CSV file
 * @param fastqDir FASTQ directory for the feature library
 * @param assignOut Output directory for assignBarcodes
 * @param options Optional assignBarcodes CLI overrides
 * @return 0 on success, non-zero on error
 */
int runAssignBarcodes(const string& whitelist,
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
