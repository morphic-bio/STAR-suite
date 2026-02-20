#ifndef PF_MULTI_PROCESS_H
#define PF_MULTI_PROCESS_H

#include "IncludeDefine.h"
#include "Parameters.h"
#include <memory>

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
class PfMultiPreloadHandle;

/**
 * @brief Start asynchronous pf-multi preload of config/paths/estimates.
 * @param P Parameters snapshot source
 * @return shared preload handle (nullptr when preload is not started)
 */
std::shared_ptr<PfMultiPreloadHandle> startPfMultiConfigPreload(const Parameters& P);

/**
 * @brief Process pf-multi config and generate combined MEX
 * @param P Parameters object (contains pfMulti flags)
 * @param solo Optional Solo instance for in-memory filtered barcode access
 * @param preload Optional async preload handle from startPfMultiConfigPreload()
 * @return 0 on success, non-zero on error
 */
int processPfMultiConfig(Parameters& P,
                         const Solo* solo = nullptr,
                         const std::shared_ptr<PfMultiPreloadHandle>& preload = nullptr);

#endif // PF_MULTI_PROCESS_H
