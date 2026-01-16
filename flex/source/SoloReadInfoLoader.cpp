#include "SoloReadInfoLoader.h"
#include "SoloReadFeature.h"
#include "SoloFeature.h"
#include "soloInputFeatureUMI.h"
#include "SoloCommon.h"
#include "serviceFuns.cpp"
#include <unordered_set>
#include <sstream>
#include <string>

#ifdef DEBUG_CB_UB_PARITY
// Optional tracing of specific readIds via STAR_DEBUG_TRACE_READS=1,2,3
static std::unordered_set<uint32_t> buildTraceReadSetLoader() {
    std::unordered_set<uint32_t> out;
    const char* env = std::getenv("STAR_DEBUG_TRACE_READS");
    if (!env || env[0]=='\0') return out;
    std::string s(env);
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        try {
            out.insert(static_cast<uint32_t>(std::stoul(tok)));
        } catch (...) {
            // ignore malformed entries
        }
    }
    return out;
}
static const std::unordered_set<uint32_t> g_traceReadsLoader = buildTraceReadSetLoader();
#endif

void SoloReadInfoLoader::load(SoloReadFeature &rf,
                              SoloReadInfoMode mode,
                              const RecordSink &sink,
                              std::vector<uint32_t> &cbReadCountTotal,
                              SoloReadFlagClass &readFlagCounts,
                              std::vector<uint32_t> &nReadPerCBunique1,
                              std::vector<uint32_t> &nReadPerCBmulti1) {

    rf.streamReads->flush();
    rf.streamReads->seekg(0,std::ios::beg);

    uint32 feature;
    uint64 umi, iread, prevIread=(uint64)-1;
    int32 cbmatch;
    int64 cb;
    std::vector<uint32> trIdDist;

    while (soloInputFeatureUMI(rf.streamReads, rf.featureType, rf.readIndexYes, rf.P.sjAll, iread, cbmatch, feature, umi, trIdDist, readFlagCounts)) {
        if (feature == (uint32)(-1) && !rf.readIndexYes) {
            rf.streamReads->ignore((uint32)-1, '\n');
            continue;
        };

        bool readIsCounted = false;
        bool featGood = ( feature != (uint32)(-1) );
        bool noMMtoWLwithoutExact = false;
        bool noTooManyWLmatches = false;
        const char* reason = "good";
#ifndef DEBUG_CB_UB_PARITY
        (void)reason;
#endif

        if (cbmatch<=1) {
            *rf.streamReads >> cb;
            if ( rf.pSolo.CBmatchWL.oneExact && cbmatch==1 && cbReadCountTotal[cb]==0 ) {
                noMMtoWLwithoutExact = true;
                reason = "noMMtoWLwithoutExact";
                // Missing CB due to oneExact requirement without exact match
                if (sink) {
                    uint32_t umi32 = (uint32_t)-1;
                    uint8_t status = 0;
#ifdef DEBUG_CB_UB_PARITY
                    if (!g_traceReadsLoader.empty() && g_traceReadsLoader.count((uint32_t)iread)) {
                        fprintf(stderr, "[TRACE loader] read=%llu cb=%lld feature=%u umi=%u status=%u reason=%s\n",
                                (unsigned long long)iread, (long long)cb, (uint32_t)feature, umi32, status, reason);
                    }
#endif
                    sink({(uint64_t)iread, (uint32_t)-1, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
#ifdef DEBUG_CB_UB_PARITY
                    if (rf.owner) rf.owner->noteDebugStatus(status, reason);
#endif
                }
            } else {
                if (!rf.pSolo.cbWLyes) {
                    cb=binarySearchExact<uintCB>(cb, rf.pSolo.cbWL.data(), rf.pSolo.cbWLsize);
                    if (cb+1 == 0) continue;
                };
                if (featGood) {
                    readIsCounted = true;
                    if (sink) {
                        uint32_t umi32 = (uint32_t)umi;
                        uint8_t status = (umi32==(uint32_t)-1) ? 2 : 1;
#ifdef DEBUG_CB_UB_PARITY
                        if (!g_traceReadsLoader.empty() && g_traceReadsLoader.count((uint32_t)iread)) {
                            fprintf(stderr, "[TRACE loader] read=%llu cb=%lld feature=%u umi=%u status=%u reason=%s\n",
                                    (unsigned long long)iread, (long long)cb, (uint32_t)feature, umi32, status, reason);
                        }
#endif
                        sink({(uint64_t)iread, (uint32_t)cb, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
#ifdef DEBUG_CB_UB_PARITY
                        if (rf.owner) rf.owner->noteDebugStatus(status, reason);
#endif
                    }
                } else {
                    // no-feature, still emit minimal data; status indicates umi validity
                    if (sink) {
                        uint32_t umi32 = (uint32_t)umi;
                        uint8_t status = (umi32==(uint32_t)-1) ? 2 : 1;
#ifdef DEBUG_CB_UB_PARITY
                        if (!g_traceReadsLoader.empty() && g_traceReadsLoader.count((uint32_t)iread)) {
                            fprintf(stderr, "[TRACE loader] read=%llu cb=%lld feature=%u umi=%u status=%u reason=%s\n",
                                    (unsigned long long)iread, (long long)cb, (uint32_t)feature, umi32, status, reason);
                        }
#endif
                        sink({(uint64_t)iread, (uint32_t)cb, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
#ifdef DEBUG_CB_UB_PARITY
                        if (rf.owner) rf.owner->noteDebugStatus(status, reason);
#endif
                    }
                };
            };
        } else {
            #ifdef MATCH_CellRanger
            double ptot=0.0, pmax=0.0, pin;
            #else
            float ptot=0.0, pmax=0.0, pin;
            #endif
            for (uint32 ii=0; ii<(uint32)cbmatch; ii++) {
                uint32 cbin; char qin; *rf.streamReads >> cbin >> qin;
                if (cbReadCountTotal[cbin]>0) {
                    qin -= rf.pSolo.QSbase; qin = qin < rf.pSolo.QSmax ? qin : rf.pSolo.QSmax;
                    pin=cbReadCountTotal[cbin]*std::pow(10.0,-qin/10.0);
                    ptot+=pin; if (pin>pmax) { cb=cbin; pmax=pin; };
                };
            };
            if (ptot>0.0 && pmax>=rf.pSolo.cbMinP*ptot) {
                if (featGood) {
                    readIsCounted = true;
                    if (sink) {
                        uint32_t umi32 = (uint32_t)umi;
                        uint8_t status = (umi32==(uint32_t)-1) ? 2 : 1;
#ifdef DEBUG_CB_UB_PARITY
                        if (!g_traceReadsLoader.empty() && g_traceReadsLoader.count((uint32_t)iread)) {
                            fprintf(stderr, "[TRACE loader] read=%llu cb=%lld feature=%u umi=%u status=%u reason=%s\n",
                                    (unsigned long long)iread, (long long)cb, (uint32_t)feature, umi32, status, reason);
                        }
#endif
                        sink({(uint64_t)iread, (uint32_t)cb, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
#ifdef DEBUG_CB_UB_PARITY
                        if (rf.owner) rf.owner->noteDebugStatus(status, reason);
#endif
                    }
                } else {
                    if (sink) {
                        uint32_t umi32 = (uint32_t)umi;
                        uint8_t status = (umi32==(uint32_t)-1) ? 2 : 1;
#ifdef DEBUG_CB_UB_PARITY
                        if (!g_traceReadsLoader.empty() && g_traceReadsLoader.count((uint32_t)iread)) {
                            fprintf(stderr, "[TRACE loader] read=%llu cb=%lld feature=%u umi=%u status=%u reason=%s\n",
                                    (unsigned long long)iread, (long long)cb, (uint32_t)feature, umi32, status, reason);
                        }
#endif
                        sink({(uint64_t)iread, (uint32_t)cb, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
#ifdef DEBUG_CB_UB_PARITY
                        if (rf.owner) rf.owner->noteDebugStatus(status, reason);
#endif
                    }
                };
            } else {
                noTooManyWLmatches = true;
                reason = "noTooManyWLmatches";
                // Emit missing CB status for multi-match failure
                if (sink) {
                    uint32_t umi32 = (uint32_t)-1;
                    uint8_t status = 0;
                    sink({(uint64_t)iread, (uint32_t)-1, umi32, status, (uint32_t)feature, rf.readIndexYes ? (uint32_t)iread : (uint32_t)-1});
#ifdef DEBUG_CB_UB_PARITY
                    if (rf.owner) rf.owner->noteDebugStatus(status, reason);
#endif
                }
            };
        };

        if ( !rf.readIndexYes || iread != prevIread ) {
            prevIread = iread;
            if (mode==SoloReadInfoMode::Counting && featGood) {
                if (cbmatch==0) { rf.stats.V[rf.stats.yessubWLmatchExact]++; }
                else if (noMMtoWLwithoutExact) { rf.stats.V[rf.stats.noMMtoWLwithoutExact]++; }
                else if (noTooManyWLmatches) { rf.stats.V[rf.stats.noTooManyWLmatches]++; };
            };
            if (mode==SoloReadInfoMode::Counting && readIsCounted) {
                if (feature<geneMultMark) nReadPerCBunique1[cb]++; else nReadPerCBmulti1[cb]++;
            };
            if ( mode==SoloReadInfoMode::Counting && rf.pSolo.readStatsYes[rf.featureType] ) {
                if ( readIsCounted ) {
                    if ( readFlagCounts.checkBit(readFlagCounts.featureU) ) readFlagCounts.setBit(readFlagCounts.countedU);
                    if ( readFlagCounts.checkBit(readFlagCounts.featureM) ) readFlagCounts.setBit(readFlagCounts.countedM);
                };
                readFlagCounts.setBit(rf.readFlag.cbMatch);
                if (cbmatch==0) { readFlagCounts.setBit(readFlagCounts.cbPerfect); readFlagCounts.countsAdd(cb); }
                else if (cbmatch==1 && !noMMtoWLwithoutExact) { readFlagCounts.setBit(readFlagCounts.cbMMunique); readFlagCounts.countsAdd(cb); }
                else if (cbmatch>1 && !noTooManyWLmatches) { readFlagCounts.setBit(readFlagCounts.cbMMmultiple); readFlagCounts.countsAdd(cb); }
                else { readFlagCounts.countsAddNoCB(); };
            };
        };
    };
}

void SoloReadInfoLoader::loadMinimal(SoloReadFeature &rf, const RecordSink &sink, std::vector<uint32_t> &cbReadCountTotal) {
    SoloReadFlagClass dummyFlags;
    std::vector<uint32_t> nU(cbReadCountTotal.size(), 0), nM(cbReadCountTotal.size(), 0);
    load(rf, SoloReadInfoMode::Minimal, sink, cbReadCountTotal, dummyFlags, nU, nM);
}
