#include <cmath>
#include "SoloReadFeature.h"
#include "SoloCommon.h"
#include "SoloFeature.h"
#include "soloInputFeatureUMI.h"
#include "serviceFuns.cpp"

// Legacy overload removed

void SoloReadFeature::inputRecords(uint32 **cbP, uint32 cbPstride, vector<uint32> &cbReadCountTotal, SoloReadFlagClass &readFlagCounts,
                                   vector<uint32> &nReadPerCBunique1, vector<uint32> &nReadPerCBmulti1,
                                   const std::function<void(uint64, uint32, uint32, uint8)> &recordSink)
{   
    streamReads->flush();
    streamReads->seekg(0,std::ios::beg);

    uint32 feature;
    uint64 umi, iread, prevIread=(uint64)-1;
    int32 cbmatch;
    int64 cb;
    vector<uint32> trIdDist;

    while (soloInputFeatureUMI(streamReads, featureType, readIndexYes, P.sjAll, iread, cbmatch, feature, umi, trIdDist, readFlagCounts)) {
        if (feature == (uint32)(-1) && !readIndexYes) { streamReads->ignore((uint32)-1, '\n'); continue; };

        bool readIsCounted = false;
        bool featGood = ( feature != (uint32)(-1) );
        bool noMMtoWLwithoutExact = false;
        bool noTooManyWLmatches = false;

        if (cbmatch<=1) {
            *streamReads >> cb;
            if ( pSolo.CBmatchWL.oneExact && cbmatch==1 && cbReadCountTotal[cb]==0 ) {
                noMMtoWLwithoutExact = true;
            } else {
                if (!pSolo.cbWLyes) {
                    cb=binarySearchExact<uintCB>(cb, pSolo.cbWL.data(), pSolo.cbWLsize);
                    if (cb+1 == 0) continue;
                };
                if (featGood) {
                    readIsCounted = true;
                    cbP[cb][0]=feature;
                    cbP[cb][1]=umi;
                    if (readIndexYes) cbP[cb][2]=iread;
                    cbP[cb]+=cbPstride;
        } else {
                    recordSink(iread, (uint32)cb, (uint32)umi, 1);
                };
            };
        } else {
            #ifdef MATCH_CellRanger
            double ptot=0.0, pmax=0.0, pin;
            #else
            float ptot=0.0, pmax=0.0, pin;
            #endif
            for (uint32 ii=0; ii<(uint32)cbmatch; ii++) {
                uint32 cbin; char qin; *streamReads >> cbin >> qin;
                if (cbReadCountTotal[cbin]>0) {
                    qin -= pSolo.QSbase; qin = qin < pSolo.QSmax ? qin : pSolo.QSmax;
                    pin=cbReadCountTotal[cbin]*std::pow(10.0,-qin/10.0);
                    ptot+=pin; if (pin>pmax) { cb=cbin; pmax=pin; };
                };
            };
            if (ptot>0.0 && pmax>=pSolo.cbMinP*ptot) {
                if (featGood) {
                    readIsCounted = true;
                    cbP[cb][0]=feature; cbP[cb][1]=umi; if (readIndexYes) cbP[cb][2]=iread; cbP[cb]+=cbPstride;
                } else {
                    recordSink(iread, (uint32)cb, (uint32)umi, 1);
                };
            } else {
                noTooManyWLmatches = true;
            };
        };

        if ( !readIndexYes || iread != prevIread ) {
            prevIread = iread;
            if (featGood) {
                if (cbmatch==0) { stats.V[stats.yessubWLmatchExact]++; }
                else if (noMMtoWLwithoutExact) { stats.V[stats.noMMtoWLwithoutExact]++; }
                else if (noTooManyWLmatches) { stats.V[stats.noTooManyWLmatches]++; };
            };
            if (readIsCounted) {
                if (feature<geneMultMark) nReadPerCBunique1[cb]++; else nReadPerCBmulti1[cb]++;
            };
            if ( pSolo.readStatsYes[featureType] ) {
                if ( readIsCounted ) {
                    if ( readFlagCounts.checkBit(readFlagCounts.featureU) ) readFlagCounts.setBit(readFlagCounts.countedU);
                    if ( readFlagCounts.checkBit(readFlagCounts.featureM) ) readFlagCounts.setBit(readFlagCounts.countedM);
                };
                readFlagCounts.setBit(readFlag.cbMatch);
                if (cbmatch==0) { readFlagCounts.setBit(readFlagCounts.cbPerfect); readFlagCounts.countsAdd(cb); }
                else if (cbmatch==1 && !noMMtoWLwithoutExact) { readFlagCounts.setBit(readFlagCounts.cbMMunique); readFlagCounts.countsAdd(cb); }
                else if (cbmatch>1 && !noTooManyWLmatches) { readFlagCounts.setBit(readFlagCounts.cbMMmultiple); readFlagCounts.countsAdd(cb); }
                else { readFlagCounts.countsAddNoCB(); };
            };
        };
    };
};

