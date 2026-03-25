#include "SoloFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "Stats.h"
#include "GlobalVariables.h"
#include "SoloReadFeature.h"
#include "ErrorWarning.h"
#include <sstream>

void SoloFeature::sumThreads()
{   
    //stats
    nReadsInput=g_statsAll.readN+1; //reserve 1 extra

    auto releaseMergedThreadState = [](SoloReadFeature *rf, bool keepInlineHashAndBridgeMaps) {
        if (rf == nullptr) {
            return;
        }
        if (!keepInlineHashAndBridgeMaps && rf->inlineHash_) {
            kh_destroy(cg_agg, rf->inlineHash_);
            rf->inlineHash_ = nullptr;
        }
        if (rf->readIdTracker_) {
            kh_destroy(readid_cbumi, rf->readIdTracker_);
            rf->readIdTracker_ = nullptr;
        }
        decltype(rf->pendingAmbiguous_)().swap(rf->pendingAmbiguous_);
        decltype(rf->bridgeAmbigReadInfoOrphan_)().swap(rf->bridgeAmbigReadInfoOrphan_);
        decltype(rf->bridgeImmediateReadCounts_)().swap(rf->bridgeImmediateReadCounts_);
        rf->bridgePinNreadUnique_.clear();
        rf->bridgePinNreadMulti_.clear();
        decltype(rf->readFlag.flagCounts)().swap(rf->readFlag.flagCounts);
        rf->readFlag.flagCountsNoCB = {};
        std::vector<uint32_t>().swap(rf->cbReadCount);
        decltype(rf->cbReadCountMap)().swap(rf->cbReadCountMap);
    };

    const bool nonFlexDirectBridge = pSolo.inlineHashMode
        && !pSolo.flexMode
        && std::getenv("STAR_SOLO_NONFLEX_HASH_BRIDGE") != nullptr;
    
    ///////////////////////////// collect RAchunk->RA->soloRead->readFeat            
    for (int ii=0; ii<P.runThreadN; ii++) {//point to
        readFeatAll[ii]= RAchunk[ii]->RA->soloRead->readFeat[pSolo.featureInd[featureType]];
        readFeatAll[ii]->setOwner(this);
        if (readFeatAll[ii]->streamReads) {
            readFeatAll[ii]->streamReads->flush();
        }
        
        // Merge inline hash if enabled
        if (pSolo.inlineHashMode) {
            if (nonFlexDirectBridge) {
                // Keep the bulky non-ambiguous thread-local hashes for direct draining later;
                // merge only the smaller ambiguous/deferred sidecars into readFeatSum.
                readFeatSum->mergePendingAmbiguous(*readFeatAll[ii]);
            } else {
                readFeatSum->mergeInlineHash(*readFeatAll[ii]);
            }
            readFeatSum->addStats(*readFeatAll[ii]);
        }
        
        readFeatSum->addCounts(*readFeatAll[ii]);

        if (pSolo.inlineHashMode) {
            releaseMergedThreadState(readFeatAll[ii], nonFlexDirectBridge);
        }
    };       
    
    // if WL was not defined
    if (!pSolo.cbWLyes) {//now we can define WL and counts ??? we do not need to do it for every feature???
        pSolo.cbWLsize=readFeatSum->cbReadCountMap.size();
        pSolo.cbWL.resize(pSolo.cbWLsize);
        pSolo.cbWLstr.resize(pSolo.cbWLsize);
        readFeatSum->cbReadCount.resize(pSolo.cbWLsize);
        readBarSum->cbReadCountExact.resize(pSolo.cbWLsize);

        if (pSolo.CBtype.type==1) {//sequence cb
            uint64 icb=0;
            for (auto &cb : readFeatSum->cbReadCountMap) {
                pSolo.cbWL[icb] = cb.first;
                pSolo.cbWLstr[icb] = convertNuclInt64toString(pSolo.cbWL[icb],pSolo.cbL);
                readFeatSum->cbReadCount[icb]=cb.second;
                readBarSum->cbReadCountExact[icb]=cb.second;
                icb++;
            };
        } else if (pSolo.CBtype.type==2) {//string cb
            vector< std::unordered_map<string,uint32>::iterator > cbiter(pSolo.CBtype.strMap.size());
            for (auto cbi=pSolo.CBtype.strMap.begin(); cbi!=pSolo.CBtype.strMap.end(); cbi++)
                cbiter[cbi->second] = cbi;

            uint64 icb=0;
            for (auto &cb : readFeatSum->cbReadCountMap) {
                pSolo.cbWL[icb] = cb.first;
                pSolo.cbWLstr[icb] = cbiter[cb.first]->first;
                readFeatSum->cbReadCount[icb]=cb.second;
                readBarSum->cbReadCountExact[icb]=cb.second;
                icb++;
            };
        };
        pSolo.cbWLstrOut = pSolo.cbWLstr;

        //pseudocounts
        if (pSolo.CBmatchWL.mm1_multi_pc) {
            for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
                readBarSum->cbReadCountExact[ii]++;//add one to exact counts
            };
        };
    };

    // if restarting from _STARtmp/solo* file
    if (P.runRestart.type==1) {//this could happen if the run is restarted. Would be better to save/load cbReadCount, or recalculate it from
        for (int ii=0; ii<P.runThreadN; ii++) {
            if (readFeatAll[ii] != nullptr && readFeatAll[ii]->binarySpool) {
                ostringstream errOut;
                errOut << "EXITING because STAR_SOLO_BINARY_SPOOL does not support restart parsing from existing solo temp files.\n"
                       << "SOLUTION: rerun without restart files or disable STAR_SOLO_BINARY_SPOOL.\n";
                exitWithError(errOut.str(), std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
            }
        }
        if (pSolo.soloFlexMinimalMemory) {
            // Warn and skip restart even if inlineHashMode is false (flag will be ignored but restart still skipped)
            P.inOut->logMain << "WARNING: --soloFlexMinimalMemory is enabled; skipping restart logic that depends on stream files/packed read info" << endl;
        } else {
            // Existing restart logic (with null check for safety)
            for (int ii=0; ii<P.runThreadN; ii++) {
                if (readFeatAll[ii]->streamReads) {
                    readFeatAll[ii]->streamReads->clear(); //just in case EOF was reached in previous reading
                    readFeatAll[ii]->streamReads->seekg(0,ios::beg);
                    string line1;
                    while (std::getline(*readFeatAll[ii]->streamReads, line1)) {
                        istringstream line1stream(line1);
                        uint64 cb1;            
                        line1stream >> cb1 >> cb1 >> cb1;
                        if (featureType==SoloFeatureTypes::SJ)
                            line1stream >> cb1;
                        line1stream >> cb1;
                        //if (cb1>readFeatSum->cbReadCount.size())
                        //    continue;//this should not happen!
                        readFeatSum->cbReadCount[cb1]++;
                    };
                }
            };
        }
    };    
    
    //detected CBs
    nCB=0;nReadsMapped=0;
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            nCB++;
            nReadsMapped += readFeatSum->cbReadCount[ii];
        };
    };
    
    indCBwl.resize(pSolo.cbWLsize, (uint32) -1);
    indCB.resize(nCB);
    nCB=0;//will count it again below
    for (uint32 ii=0; ii<pSolo.cbWLsize; ii++) {
        if (readFeatSum->cbReadCount[ii]>0) {
            indCB[nCB]=ii;
            indCBwl[ii]=nCB;
            ++nCB;
        };
    };
    
};
    
