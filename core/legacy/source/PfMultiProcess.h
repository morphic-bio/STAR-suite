#ifndef PF_MULTI_PROCESS_H
#define PF_MULTI_PROCESS_H

#include "IncludeDefine.h"
#include "Parameters.h"

/**
 * @file PfMultiProcess.h
 * @brief Main entry point for pf-multi config processing
 * 
 * Orchestrates the entire pf-multi workflow:
 * 1. Parse multi config
 * 2. Run assignBarcodes for each feature type
 * 3. Generate MEX stub files
 * 4. Merge with GEX MEX
 * 5. Write combined output
 */

class Solo;

/**
 * @brief Process pf-multi config and generate combined MEX
 * @param P Parameters object (contains pfMulti flags)
 * @param solo Optional Solo instance for in-memory filtered barcode access
 * @return 0 on success, non-zero on error
 */
int processPfMultiConfig(Parameters& P, const Solo* solo = nullptr);

#endif // PF_MULTI_PROCESS_H
