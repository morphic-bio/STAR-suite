#include "star_chromap_orchestration.h"

#include "IncludeDefine.h"
#include "Parameters.h"

StarChromapAtacAsyncRun::StarChromapAtacAsyncRun() : impl_(nullptr) {}
StarChromapAtacAsyncRun::~StarChromapAtacAsyncRun() {}

bool preflightStarChromapAtacIfEnabled(Parameters &P, bool /*batchModeActive*/) {
  if (P.chromapAtac.enabled == 0) {
    return true;
  }
  P.inOut->logMain
      << "ERROR: --chromapAtacEnable 1 requires a STAR binary built with "
         "Chromap support. Rebuild with WITH_CHROMAP=1 (see docs/LIBCHROMAP_CONTRACT.md).\n";
  return false;
}

bool startStarChromapAtacIfEnabled(Parameters &P,
                                   bool batchModeActive,
                                   StarChromapAtacAsyncRun & /*run*/) {
  return preflightStarChromapAtacIfEnabled(P, batchModeActive);
}

bool runStarChromapAtacIfEnabled(Parameters &P,
                                 bool batchModeActive,
                                 StarChromapAtacAsyncRun & /*run*/) {
  return preflightStarChromapAtacIfEnabled(P, batchModeActive);
}
