#ifndef PF_MULTI_TABLE_IMPORT_H
#define PF_MULTI_TABLE_IMPORT_H

#include "IncludeDefine.h"

/**
 * @file PfMultiTableImport.h
 * @brief Table-backed pf-multi feature library import (barcode,feature_id,count)
 */

namespace PfMultiTableImport {

struct TableImportStats {
    uint64_t rowsRead = 0;
    uint64_t rowsRetained = 0;
    uint64_t rowsRejectedBarcode = 0;
    uint64_t rowsRejectedFeature = 0;
    uint64_t rowsRejectedCount = 0;
    uint64_t duplicatePairsCollapsed = 0;
    uint64_t rowsZeroSkipped = 0;
    uint64_t rowsSuffixNormalized = 0;
    uint64_t permitChunksProcessed = 0;
    uint64_t featurePermitAcquires = 0;
    string delimiter;
    string inputPath;
    string barcodeNamespaceInput;
    string barcodeNamespaceOutput;
};

struct TableImportOptions {
    bool enableStarDynamicPermitHooks = false;
    string filteredBarcodesPath;
    string featureTypeLabel;
    string sampleName;
    string assignmentWhitelistNamespace;
};

struct TableImportResult {
    int returnCode = 0;
    TableImportStats stats;
};

/**
 * @brief Parse a headered CSV/TSV count table and emit per-library MEX artifacts.
 */
TableImportResult runTableFeatureImport(const string& whitelistNormalizedPath,
                                        const string& featureRefPath,
                                        const string& tablePath,
                                        const string& assignOut,
                                        const TableImportOptions& options);

} // namespace PfMultiTableImport

#endif // PF_MULTI_TABLE_IMPORT_H
