#include "SoloFeature.h"
#include "SoloReadFeature.h"
#include "streamFuns.h"
#include "TimeFunctions.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include "systemFunctions.h"

void SoloFeature::processRecords()
{
    if (pSolo.type==0)
        return;

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Starting Solo post-map for " <<SoloFeatureTypes::Names[featureType] <<endl;
    
    outputPrefix= P.outFileNamePrefix+pSolo.outFileNames[0];
    outputPrefix+= SoloFeatureTypes::Names[featureType] +'/';
    outputPrefixFiltered= outputPrefix + "filtered/";
    
    if (mkdir(outputPrefix.c_str(),P.runDirPerm)!=0 && errno!=EEXIST) {//create directory
        ostringstream errOut;
        errOut << "EXITING because of fatal OUTPUT FILE error: could not create Solo output directory"<<outputPrefix<<"\n";
        errOut << "SOLUTION: check the path and permisssions";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_PARAMETER, P);
    };    
     
    //prepare for feature-specific counting:
    if (featureType==SoloFeatureTypes::SJ && P.sjAll[0].empty()) {
        ifstream &sjIn = ifstrOpen(P.outFileTmp+"SJ.start_gap.tsv",  ERROR_OUT, "SOLUTION: re-run STAR", P);
        P.sjAll[0].reserve(10000000);
        P.sjAll[1].reserve(10000000);
        uint64 start1, gap1;
        while ( sjIn >> start1 >> gap1 ) {            
            P.sjAll[0].emplace_back(start1);
            P.sjAll[1].emplace_back(gap1);
        };
        sjIn.close();
        P.inOut->logMain <<"Read splice junctions for Solo SJ feature: "<< P.sjAll[0].size() <<endl;
    };

    SoloFeature::sumThreads();
    
    // Early exit for skipProcessing mode: populate readInfo but skip counting/matrices
    if (pSolo.skipProcessing) {
        // Call minimal processing to populate readInfo using collapseUMI in minimal mode
        prepareReadInfoOnly();
        
        // Write tag table if requested (readInfo is now available)
        writeTagTableIfRequested(false);

        // Sanity check: ensure readInfo is populated for all reads (packed)
        if (packedReadInfo.data.size() != nReadsInput) {
            P.inOut->logMain << "WARNING: packedReadInfo size (" << packedReadInfo.data.size() 
                             << ") does not equal nReadsInput (" << nReadsInput 
                             << ") in skipProcessing mode" << endl;
        }
        
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Solo: skipping counting and matrix output for " 
                         << SoloFeatureTypes::Names[featureType] << " (soloSkipProcessing=yes)" <<endl;
        return;
    }
    
    //call counting method
    if (featureType==SoloFeatureTypes::Velocyto) {
        countVelocyto();
    } else if (featureType==SoloFeatureTypes::Transcript3p) {
        quantTranscript();
        return;
    } else {//all others, standard processing
        if (pSolo.type==pSolo.SoloTypes::SmartSeq) {
            countSmartSeq();
        } else {
            countCBgeneUMI();
        };
    };
    
    // Inline hash path already wrote MEX directly; skip legacy matrix output to avoid
    // touching uninitialized Solo dense matrices.
    if (pSolo.inlineHashMode) {
        time(&rawTime);
        P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Solo: inline-hash mode completed (skipping legacy output)" <<endl;
        return;
    }
    
    // no need to include multi-gene in yesWLmatch:
    // for multi-gene options, it's already done in _collapseUMI_ function
    //if (!(pSolo.multiMap.yes.multi && (featureType==SoloFeatureTypes::Gene || featureType==SoloFeatureTypes::GeneFull)))
    //    readFeatSum->stats.V[readFeatSum->stats.yesWLmatch] += readFeatSum->stats.V[readFeatSum->stats.MultiFeature];   

    //output
    ofstream *statsStream = &ofstrOpen(outputPrefix+"Features.stats",ERROR_OUT, P);
    readFeatSum->statsOut(*statsStream);
    statsStream->close();
    
    // Write tag table if requested (after UMI collapse, before matrix output)
    writeTagTableIfRequested(false);
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Solo: writing raw matrix" <<endl;

    //output nU per gene per CB
    outputResults(false,  outputPrefix + "/raw/"); //unfiltered
    
    time(&rawTime);
    P.inOut->logMain << timeMonthDayTime(rawTime) << " ... Solo: cell filtering" <<endl;    
    cellFiltering();
    
    //summary stats output
    statsOutput();
    
    // Write UMI cleanup metrics and QC summary if UMI correction is enabled
    if (pSolo.umiCorrectionMode > 0 || !cellsAllowSet.empty() || !assignmentsMap.empty()) {
        writeUMICleanupMetrics();
        writeQCSummary();
    }
    
    //delete big arrays allocated in the previous functions
    clearLarge();
    //delete[] indCB;

    P.inOut->logMain << "RAM after completing solo:\n"
                     <<  linuxProcMemory() << flush;   
};
