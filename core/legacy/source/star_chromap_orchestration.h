#ifndef STAR_CHROMAP_ORCHESTRATION_H_
#define STAR_CHROMAP_ORCHESTRATION_H_

class Parameters;

// Validates Chromap ATAC configuration before expensive STAR work starts.
// Returns true if disabled or valid; false if enabled and invalid.
bool preflightStarChromapAtacIfEnabled(Parameters &P, bool batchModeActive);

// Runs Chromap ATAC via star::multiome::runChromapAtac when --chromapAtacEnable 1.
// Returns true if Chromap was skipped because it is disabled, or if it ran OK.
// Returns false if Chromap is enabled but cannot run or failed.
bool runStarChromapAtacIfEnabled(Parameters &P, bool batchModeActive);

#endif  // STAR_CHROMAP_ORCHESTRATION_H_
