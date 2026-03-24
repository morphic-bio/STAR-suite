#ifndef CR_MULTI_MEX_STUB_H
#define CR_MULTI_MEX_STUB_H

#include "IncludeDefine.h"

/**
 * @file CrMultiMexStub.h
 * @brief Generate features.tsv and barcodes.tsv from assignBarcodes outputs
 * 
 * Ports logic from assignbarcodes_mex_stub.py to C++.
 * Reads feature_reference.csv and writes features.tsv/barcodes.tsv,
 * validating against features.txt when present.
 */

namespace CrMultiMexStub {

/**
 * @struct FeatureRow
 * @brief Single feature entry from feature reference CSV
 */
struct FeatureRow {
    string id;          // Feature ID
    string name;       // Feature name
    string featureType; // Feature type (e.g., "CRISPR Guide Capture")
};

/**
 * @brief Load feature reference CSV
 * @param csvPath Path to feature_reference.csv
 * @return Vector of feature rows
 * @throws runtime_error on parse errors
 */
vector<FeatureRow> loadFeatureCsv(const string& csvPath);

/**
 * @brief Read features.txt (one feature name per line)
 * @param txtPath Path to features.txt
 * @return Vector of feature names
 */
vector<string> readFeaturesTxt(const string& txtPath);

/**
 * @brief Compare feature names from CSV and features.txt
 * @param featureRows Features from CSV
 * @param featuresTxt Features from features.txt
 * @return Empty string if match, error message otherwise
 */
string compareFeatureNames(const vector<FeatureRow>& featureRows, 
                          const vector<string>& featuresTxt);

/**
 * @brief Write features.tsv file
 * @param outPath Output path for features.tsv
 * @param featureRows Features to write
 * @param defaultType Default feature type if missing
 * @param force Overwrite existing file
 * @return true if file was written, false if skipped
 */
bool writeFeaturesTsv(const string& outPath, const vector<FeatureRow>& featureRows,
                     const string& defaultType, bool force);

/**
 * @brief Copy barcodes.txt to barcodes.tsv
 * @param barcodesTxt Path to barcodes.txt
 * @param barcodesTsv Output path for barcodes.tsv
 * @param force Overwrite existing file
 * @return true if file was copied, false if skipped
 */
bool copyBarcodesTsv(const string& barcodesTxt, const string& barcodesTsv, bool force);

/**
 * @brief Process assignBarcodes output directory
 * @param assignOutDir assignBarcodes output directory (may contain filtered/ subdirectory)
 * @param featureCsvPath Path to feature_reference.csv
 * @param defaultFeatureType Default feature type (default: "Custom")
 * @param force Overwrite existing files
 * @return 0 on success, 1 on error
 */
int processAssignOutput(const string& assignOutDir, const string& featureCsvPath,
                       const string& defaultFeatureType = "Custom", bool force = false);

} // namespace CrMultiMexStub

#endif // CR_MULTI_MEX_STUB_H
