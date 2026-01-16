#include "ReadAlign.h"
#include <cstdlib>
#include <atomic>

// Static counter for debug logging (guarded by STAR_TRIM_DEBUG_N env var)
// Use atomic for thread-safety when debug logging is enabled
static atomic<uint64_t> g_mapDebugCount(0);
static atomic<int64_t> g_mapDebugMax(-1);  // -1 means not initialized

void ReadAlign::mappedFilter() {//filter mapped read, add to stats
    unmapType=-1;//mark as mapped
    if ( nW==0 ) {//no good windows
        statsRA.unmappedOther++;
        unmapType=0;
    } else if ( (trBest->maxScore < P.outFilterScoreMin) || (trBest->maxScore < (intScore) (P.outFilterScoreMinOverLread*(Lread-1))) \
              || (trBest->nMatch < P.outFilterMatchNmin)  || (trBest->nMatch < (uint) (P.outFilterMatchNminOverLread*(Lread-1))) ) {//too short
        statsRA.unmappedShort++;
        unmapType=1;
        
        // Debug logging for unmapped-too-short (guarded by STAR_TRIM_DEBUG_N env var)
        int64_t debugMax = g_mapDebugMax.load();
        if (debugMax == -1) {
            const char* debugEnv = getenv("STAR_TRIM_DEBUG_N");
            int64_t newMax = debugEnv ? atol(debugEnv) : 0;
            int64_t expected = -1;
            if (g_mapDebugMax.compare_exchange_strong(expected, newMax)) {
                debugMax = newMax;
            } else {
                debugMax = g_mapDebugMax.load();
            }
        }
        uint64_t debugCount = g_mapDebugCount.fetch_add(1);
        if (debugMax > 0 && debugCount < (uint64_t)debugMax) {
            P.inOut->logMain << "MAP_SHORT: " << readNameMates[0]
                             << " Lread=" << Lread
                             << " nMatch=" << trBest->nMatch
                             << " outFilterMatchNmin=" << P.outFilterMatchNmin
                             << " readLen[0]=" << readLength[0]
                             << " readLen[1]=" << readLength[1]
                             << " readLenPairOriginal=" << readLengthPairOriginal
                             << " maxScore=" << trBest->maxScore
                             << " outFilterScoreMin=" << P.outFilterScoreMin
                             << " outFilterMatchNminOverLread=" << P.outFilterMatchNminOverLread
                             << endl;
            g_mapDebugCount++;
        }
    } else if ( (trBest->nMM > outFilterMismatchNmaxTotal) || (double(trBest->nMM)/double(trBest->rLength)>P.outFilterMismatchNoverLmax) ) {//too many mismatches
        statsRA.unmappedMismatch++;
        unmapType=2;
    } else if (nTr > P.outFilterMultimapNmax){//too multi
        statsRA.unmappedMulti++;
        unmapType=3;
    };

    return;
};