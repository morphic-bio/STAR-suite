#ifndef CR_MULTI_PROCESS_H
#define CR_MULTI_PROCESS_H

#include "IncludeDefine.h"
#include "Parameters.h"

/**
 * @file CrMultiProcess.h
 * @brief Main entry point for Cell Ranger multi config processing
 * 
 * Orchestrates the entire CR multi workflow:
 * 1. Parse multi config
 * 2. Run assignBarcodes for each feature type
 * 3. Generate MEX stub files
 * 4. Merge with GEX MEX
 * 5. Write combined output
 */

/**
 * @brief Process Cell Ranger multi config and generate combined MEX
 * @param P Parameters object (contains crMulti flags)
 * @return 0 on success, non-zero on error
 */
int processCrMultiConfig(Parameters& P);

#endif // CR_MULTI_PROCESS_H
