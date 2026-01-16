#include "InOutStreams.h"
#include "GlobalVariables.h"

InOutStreams::InOutStreams() {
    logStdOut=NULL;
    outSAM=NULL;
    outBAMfileUnsorted=NULL;
    outQuantBAMfile=NULL;
    outBAMfileY=NULL;
    outBAMfileNoY=NULL;
};

InOutStreams::~InOutStreams() {

    if (logStdOut!=NULL) logStdOut->flush();
    if (outSAM!=NULL) outSAM->flush();

    logStdOutFile.flush();
    outSAMfile.flush();

    outChimSAM.flush();
    outChimJunction.flush();
    logProgress.flush();
    logMain.flush();
    logFinal.flush();
    outLocalChains.flush();
    
    // Solo tmp stream
    if (outBAMfileUnsortedSoloTmp.is_open()) {
        outBAMfileUnsortedSoloTmp.flush();
    };

    //logStdOutFile.close(); //do not want to close these log files, as some destructors (e.g. ~SharedMemory) might still write there
    //logMain.close();

    outSAMfile.close();
    outChimSAM.close();
    outChimJunction.close();
    logProgress.close();
    logFinal.close();
    outLocalChains.close();
    
    // Solo tmp stream
    if (outBAMfileUnsortedSoloTmp.is_open()) {
        outBAMfileUnsortedSoloTmp.close();
    };


    for (int ii=0;ii<2;ii++) {
        if (outUnmappedReadsStream[ii].is_open()) {
            outUnmappedReadsStream[ii].flush();
            outUnmappedReadsStream[ii].close();
        }
        if (outYFastqStream[ii].is_open()) {
            outYFastqStream[ii].flush();
            outYFastqStream[ii].close();
        }
        if (outNoYFastqStream[ii].is_open()) {
            outNoYFastqStream[ii].flush();
            outNoYFastqStream[ii].close();
        }
    };
};

