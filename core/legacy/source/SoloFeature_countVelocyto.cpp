#include "SoloFeature.h"
#include "TimeFunctions.h"
#include "systemFunctions.h"
#include <cstdlib>
#include <cstring>

// Default: per-thread stream order (stream path lives in SoloFeature_countVelocytoBridge.cpp).
// STAR_VELOCYTO_DETERMINISTIC_REPLAY=1 uses sorted replay + same finalize.
// STAR_VELOCYTO_INTEGRATED_HASH=1 with deterministic replay uses Stage 2 merge (default: disk spill shards
// under Solo.out/Velocyto/; optional INMEMORY=1). See SoloFeature_countVelocytoBridge.cpp.
// Harness: stream@1t vs deterministic@1t only proves equivalence inside this refactored tree; for true
// legacy parity compare stream output to a frozen pre-refactor run (UCSF_VELOCYTO_BASELINE_OUTDIR).

void SoloFeature::countVelocyto()
{
    const char *det = std::getenv("STAR_VELOCYTO_DETERMINISTIC_REPLAY");
    if (det != nullptr && det[0] != '\0' && std::strcmp(det, "0") != 0) {
        time_t rawTime;
        time(&rawTime);
        const char *ih = std::getenv("STAR_VELOCYTO_INTEGRATED_HASH");
        const bool useCbBuckets =
            ih != nullptr && ih[0] != '\0' && std::strcmp(ih, "0") != 0;
        if (useCbBuckets) {
            P.inOut->logMain << timeMonthDayTime(rawTime)
                             << " ... Velocyto: deterministic replay + STAR_VELOCYTO_INTEGRATED_HASH "
                                "(CB-bucket / integrated hash path)"
                             << endl;
            countVelocytoSortedReplayCBuckets();
        } else {
            P.inOut->logMain << timeMonthDayTime(rawTime)
                            << " ... Velocyto: STAR_VELOCYTO_DETERMINISTIC_REPLAY enabled (sorted replay path)"
                            << endl;
            countVelocytoSortedReplay();
        }
    } else {
        countVelocytoStreamThreads();
    }
}
