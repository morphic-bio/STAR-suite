#ifndef STAR_CHROMAP_ORCHESTRATION_H_
#define STAR_CHROMAP_ORCHESTRATION_H_

class Parameters;

class StarChromapAtacAsyncRun {
 public:
  StarChromapAtacAsyncRun();
  ~StarChromapAtacAsyncRun();

  StarChromapAtacAsyncRun(const StarChromapAtacAsyncRun &) = delete;
  StarChromapAtacAsyncRun &operator=(const StarChromapAtacAsyncRun &) = delete;

 private:
  friend bool startStarChromapAtacIfEnabled(Parameters &P,
                                            bool batchModeActive,
                                            StarChromapAtacAsyncRun &run);
  friend bool runStarChromapAtacIfEnabled(Parameters &P,
                                          bool batchModeActive,
                                          StarChromapAtacAsyncRun &run);
  struct Impl;
  Impl *impl_;
};

// Validates Chromap ATAC configuration before expensive STAR work starts.
// Returns true if disabled or valid; false if enabled and invalid.
bool preflightStarChromapAtacIfEnabled(Parameters &P, bool batchModeActive);

// Starts Chromap ATAC in the background when chromapAtacStartMode=concurrent.
// Returns true if disabled, postMapping mode, or successfully started.
bool startStarChromapAtacIfEnabled(Parameters &P,
                                   bool batchModeActive,
                                   StarChromapAtacAsyncRun &run);

// Runs Chromap ATAC via star::multiome::runChromapAtac when --chromapAtacEnable 1.
// Returns true if Chromap was skipped because it is disabled, or if it ran OK.
// Returns false if Chromap is enabled but cannot run or failed.
bool runStarChromapAtacIfEnabled(Parameters &P,
                                 bool batchModeActive,
                                 StarChromapAtacAsyncRun &run);

#endif  // STAR_CHROMAP_ORCHESTRATION_H_
