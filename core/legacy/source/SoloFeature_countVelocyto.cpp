#include "SoloFeature.h"
#include "TimeFunctions.h"
#include "systemFunctions.h"
#include "SoloMemoryProfile.h"
#include <cstdlib>
#include <cstring>
#include <sstream>

namespace {

static bool envFlagTrue(const char *name)
{
    const char *v = std::getenv(name);
    return v != nullptr && v[0] != '\0' && std::strcmp(v, "0") != 0;
}

} // namespace

// Default: per-thread stream order (stream path lives in SoloFeature_countVelocytoBridge.cpp).
// STAR_VELOCYTO_LOW_MEM=1 selects the range-spill path directly.
// STAR_VELOCYTO_DETERMINISTIC_REPLAY=1 uses sorted replay + same finalize.
// STAR_VELOCYTO_INTEGRATED_HASH=1 with deterministic replay uses Stage 2 merge (default: disk spill shards
// under Solo.out/Velocyto/; optional INMEMORY=1). See SoloFeature_countVelocytoBridge.cpp.
// Harness: stream@1t vs deterministic@1t only proves equivalence inside this refactored tree; for true
// legacy parity compare stream output to a frozen pre-refactor run (UCSF_VELOCYTO_BASELINE_OUTDIR).

void SoloFeature::countVelocyto()
{
    {
        std::ostringstream extra;
        extra << "packedReadInfo_words=" << packedReadInfo.data.size()
              << " nCB=" << nCB;
        soloMemoryProfileCheckpoint(P.inOut->logMain, "countVelocyto_enter", extra.str());
    }
    if (envFlagTrue("STAR_VELOCYTO_LOW_MEM") || envFlagTrue("STAR_VELOCYTO_RANGE_SPILL")) {
        time_t rawTime;
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime)
                         << " ... Velocyto: STAR_VELOCYTO_LOW_MEM enabled "
                            "(range-spill / bounded-CB-map path)"
                         << endl;
        countVelocytoSortedReplayCBuckets();
        return;
    }

    const bool det = envFlagTrue("STAR_VELOCYTO_DETERMINISTIC_REPLAY");
    if (det) {
        time_t rawTime;
        time(&rawTime);
        const bool useCbBuckets = envFlagTrue("STAR_VELOCYTO_INTEGRATED_HASH");
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
