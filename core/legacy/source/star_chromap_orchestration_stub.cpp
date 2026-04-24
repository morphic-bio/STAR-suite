#include "star_chromap_orchestration.h"

#include "IncludeDefine.h"
#include "Parameters.h"

bool preflightStarChromapAtacIfEnabled(Parameters &P, bool /*batchModeActive*/) {
  if (P.chromapAtac.enabled == 0) {
    return true;
  }
  P.inOut->logMain
      << "ERROR: --chromapAtacEnable 1 requires a STAR binary built with "
         "Chromap support. Rebuild with WITH_CHROMAP=1 (see docs/LIBCHROMAP_CONTRACT.md).\n";
  return false;
}

bool runStarChromapAtacIfEnabled(Parameters &P, bool batchModeActive) {
  return preflightStarChromapAtacIfEnabled(P, batchModeActive);
}
