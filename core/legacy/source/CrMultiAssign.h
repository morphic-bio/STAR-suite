#ifndef CR_MULTI_ASSIGN_H
#define CR_MULTI_ASSIGN_H

#include "IncludeDefine.h"
#include "CrMultiConfig.h"

/**
 * @file CrMultiAssign.h
 * @brief Execute assignBarcodes for feature barcode processing
 * 
 * Invokes assignBarcodes binary as subprocess for each feature type library.
 */

namespace CrMultiAssign {

/**
 * @brief Run assignBarcodes for a feature library
 * @param assignBin Path to assignBarcodes binary
 * @param whitelist Cell barcode whitelist file
 * @param featureRef Feature reference CSV file
 * @param fastqDir FASTQ directory for the feature library
 * @param assignOut Output directory for assignBarcodes
 * @return 0 on success, non-zero on error
 */
int runAssignBarcodes(const string& assignBin, const string& whitelist,
                     const string& featureRef, const string& fastqDir,
                     const string& assignOut);

/**
 * @brief Process all feature libraries from config
 * @param config Parsed multi config
 * @param assignBin Path to assignBarcodes binary
 * @param whitelist Cell barcode whitelist (from config or override)
 * @param featureRef Feature reference (from config or override)
 * @param fastqMap FASTQ path mapping
 * @param fastqRoot FASTQ root directory
 * @param outPrefix Output prefix (for creating subdirectories)
 * @param featureTypes List of feature types to process (e.g., {"CRISPR Guide Capture", "Antibody Capture"})
 * @return 0 on success, non-zero on error
 */
int processFeatureLibraries(const CrMultiConfig::Config& config,
                          const string& assignBin, const string& whitelist,
                          const string& featureRef, const map<string, string>& fastqMap,
                          const string& fastqRoot, const string& outPrefix,
                          const vector<string>& featureTypes);

} // namespace CrMultiAssign

#endif // CR_MULTI_ASSIGN_H
