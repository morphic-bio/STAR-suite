#ifndef INOUTSTREAMS_DEF
#define INOUTSTREAMS_DEF

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H
#include <mutex>
#include <vector>

class InOutStreams {
    public:
    ostream *logStdOut, *outSAM;
    ofstream logStdOutFile, outSAMfile;
    BGZF *outBAMfileUnsorted, *outBAMfileCoord, *outQuantBAMfile;
    BGZF *outBAMfileY, *outBAMfileNoY;  // Y-chromosome split BAM handles

    ofstream outChimSAM, outChimJunction, logMain, logProgress, logFinal, outUnmappedReadsStream[MAX_N_MATES];
    ofstream outYFastqStream[MAX_N_MATES], outNoYFastqStream[MAX_N_MATES];  // Y/noY FASTQ output streams (uncompressed)
    ifstream readIn[MAX_N_MATES];

    //compilation-optional streams
    ofstream outLocalChains;

    // Track heap-allocated streams created via streamFuns (ifstrOpen/ofstrOpen/fstrOpen)
    // so they can be deleted at shutdown and not show up as LSAN leaks.
    std::mutex ownedStreamsMutex;
    std::vector<std::ifstream*> ownedIfstreams;
    std::vector<std::ofstream*> ownedOfstreams;
    std::vector<std::fstream*> ownedFstreams;

    InOutStreams();
    ~InOutStreams();
};

#endif
