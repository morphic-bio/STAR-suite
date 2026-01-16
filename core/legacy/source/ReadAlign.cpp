#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "ReadAlign.h"
#include "SampleDetector.h"
#include "GlobalVariables.h"
#include "TranscriptQuantEC.h"

ReadAlign::ReadAlign (Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk,
                      const libem::Transcriptome* libemTr)
                    : mapGen(genomeIn), genOut(*genomeIn.genomeOut.g), P(Pin), chunkTr(TrIn), quantEC(nullptr)
{
    readNmates=P.readNmates; //not readNends
    bamRecordIndexPtr = &g_bamRecordIndex;
    if (P.trimQcEnabled) {
        trimQc.init(static_cast<uint32_t>(P.readNmates), P.trimQcMaxReads, P.readQualityScoreBase);
    }
    //RNGs
    rngMultOrder.seed(P.runRNGseed*(iChunk+1));
    rngUniformReal0to1=std::uniform_real_distribution<double> (0.0, 1.0);
    //transcriptome
    if ( P.quant.trSAM.yes || P.quant.transcriptVB.yes ) {
        alignTrAll=new Transcript [P.alignTranscriptsPerReadNmax];
    };
    
    // Initialize EC table for transcript quantification
    if (P.quant.transcriptVB.yes && chunkTr != nullptr && chunkTr->nTr > 0) {
        // Create AlignmentModel per thread if error model is enabled
        std::unique_ptr<libem::AlignmentModel> alignment_model;
        if (P.quant.transcriptVB.errorModelMode != "off" && libemTr != nullptr) {
            alignment_model.reset(new libem::AlignmentModel(0.001, 4));  // alpha=0.001, readBins=4
        }
        
        quantEC = new TranscriptQuantEC(chunkTr->nTr, iChunk, 
                                         P.quant.transcriptVB.traceFile,
                                         P.quant.transcriptVB.traceLimit, P,
                                         alignment_model.release(),  // Transfer ownership
                                         libemTr);  // Shared read-only transcriptome
        // Set transcript lengths at construction time so startPosProb is enabled during EC building
        std::vector<int32_t> lens(chunkTr->nTr);
        for (uint i = 0; i < chunkTr->nTr; ++i) {
            lens[i] = static_cast<int32_t>(chunkTr->trLen[i]);
        }
        quantEC->setTranscriptLengths(lens);
    }

    if (P.pGe.gType==101) {//SuperTranscriptome
        splGraph = new SpliceGraph(*mapGen.superTr, P, this);
    } else {//standard map algorithm:
        winBin = new uintWinBin* [2];
        winBin[0] = new uintWinBin [P.winBinN];
        winBin[1] = new uintWinBin [P.winBinN];
        memset(winBin[0],255,sizeof(winBin[0][0])*P.winBinN);
        memset(winBin[1],255,sizeof(winBin[0][0])*P.winBinN);
        //split
        splitR=new uint*[3];
        splitR[0]=new uint[P.maxNsplit]; splitR[1]=new uint[P.maxNsplit]; splitR[2]=new uint[P.maxNsplit];
        //alignments
        PC=new uiPC[P.seedPerReadNmax];
        WC=new uiWC[P.alignWindowsPerReadNmax];
        nWA=new uint[P.alignWindowsPerReadNmax];
        nWAP=new uint[P.alignWindowsPerReadNmax];
        WALrec=new uint[P.alignWindowsPerReadNmax];
        WlastAnchor=new uint[P.alignWindowsPerReadNmax];
    
        WA=new uiWA*[P.alignWindowsPerReadNmax];
        for (uint ii=0;ii<P.alignWindowsPerReadNmax;ii++)
            WA[ii]=new uiWA[P.seedPerWindowNmax];
        WAincl = new bool [P.seedPerWindowNmax];        

        #ifdef COMPILE_FOR_LONG_READS
        swWinCov = new uint[P.alignWindowsPerReadNmax];
        scoreSeedToSeed = new intScore [P.seedPerWindowNmax*(P.seedPerWindowNmax+1)/2];
        scoreSeedBest = new intScore [P.seedPerWindowNmax];
        scoreSeedBestInd = new uint [P.seedPerWindowNmax];
        scoreSeedBestMM = new uint [P.seedPerWindowNmax];
        seedChain = new uint [P.seedPerWindowNmax];
        #endif
    };

    //aligns a.k.a. transcripts
    trAll = new Transcript**[P.alignWindowsPerReadNmax+1];
    nWinTr = new uint[P.alignWindowsPerReadNmax];
    trArray = new Transcript[P.alignTranscriptsPerReadNmax];
    trArrayPointer =  new Transcript*[P.alignTranscriptsPerReadNmax];
    for (uint ii=0;ii<P.alignTranscriptsPerReadNmax;ii++)
        trArrayPointer[ii]= &(trArray[ii]);
    trInit = new Transcript;
    
    if (mapGen.genomeOut.convYes) {//allocate output transcripts
        alignsGenOut.alMult = new Transcript*[P.outFilterMultimapNmax];
        for (uint32 ii=0; ii<P.outFilterMultimapNmax; ii++) 
            alignsGenOut.alMult[ii]=new Transcript;
    };
    
    //read
    Read0 = new char*[P.readNends];
    Qual0 = new char*[P.readNends];
    readNameMates=new char* [P.readNends];
    for (uint32 ii=0; ii<P.readNends; ii++) {
        readNameMates[ii]= new char [DEF_readNameLengthMax];
        Read0[ii]        = new char [DEF_readSeqLengthMax+1];
        Qual0[ii]        = new char [DEF_readSeqLengthMax+1];        
    };
    readNameExtra.resize(P.readNends);
    readName = readNameMates[0];

    Read1 = new char*[3];
    Read1[0]=new char[DEF_readSeqLengthMax+1]; 
    Read1[1]=new char[DEF_readSeqLengthMax+1]; 
    Read1[2]=new char[DEF_readSeqLengthMax+1];
    
    for (auto &q: qualHist)
        q.fill(0);
    
    //outBAM
    outBAMoneAlignNbytes = new uint [P.readNmates+2]; //extra piece for chimeric reads //not readNends: this is alignment
    outBAMoneAlign = new char* [P.readNmates+2]; //extra piece for chimeric reads //not readNends: this is alignment
    for (uint ii=0; ii<P.readNmates+2; ii++) {//not readNends: this is alignment
        outBAMoneAlign[ii]=new char [BAMoutput_oneAlignMaxBytes];
    };
    resetN();
    
    // Initialize Y/noY FASTQ streams
    for (uint32 imate = 0; imate < P.readNends; imate++) {
        chunkOutYFastqGz[imate] = nullptr;
        chunkOutNoYFastqGz[imate] = nullptr;
    }
    
    //chim
    chunkOutChimJunction = new fstream;
    chimDet = new ChimericDetection(P, trAll, nWinTr, Read1, mapGen, chunkOutChimJunction, this);
    
    //solo
    soloRead = new SoloRead (P, iChunk);
    
    // Sample detector initialization (mirror BAMoutputSoloTmp logic)
    sampleDet_ = nullptr;
    sampleDetReady_ = false;
    detectedSampleByte_ = 0xFFu; // Default: no sample detected
    extractedCbIdxPlus1_ = 0;    // Default: no CB extracted
    extractedUmi24_ = 0;          // Default: no UMI extracted
    if (!P.pSolo.sampleWhitelistPath.empty() && P.pSolo.sampleWhitelistPath != "-" &&
        !P.pSolo.sampleProbesPath.empty() && P.pSolo.sampleProbesPath != "-") {
        sampleDet_ = new SampleDetector(P.pSolo);
        if (sampleDet_->loadWhitelist(P.pSolo.sampleWhitelistPath) &&
            sampleDet_->loadProbes(P.pSolo.sampleProbesPath)) {
            sampleDetReady_ = true;
            P.inOut->logMain << "ReadAlign: SampleDetector initialized successfully (whitelist="
                             << P.pSolo.sampleWhitelistPath << ", probes="
                             << P.pSolo.sampleProbesPath << ")" << std::endl;
        } else {
            P.inOut->logMain << "WARNING: ReadAlign SampleDetector initialization failed (whitelist="
                             << P.pSolo.sampleWhitelistPath << ", probes="
                             << P.pSolo.sampleProbesPath << ")" << std::endl;
            delete sampleDet_;
            sampleDet_ = nullptr;
        }
    }
    
    //clipping
    P.pClip.initializeClipMates(clipMates);

    //debug
    {
    #ifdef DEBUG_OutputLastRead
        lastReadStream.open((P.outFileTmp+"/lastRead_"+to_string(iChunk)).c_str());
    #endif
    };
};

void ReadAlign::resetN () {//reset resets the counters to 0 for a new read
    mapMarker=0;
    nA=0; nP=0; nW=0;
    nTr=0;
    nUM[0]=0; nUM[1]=0;
    storedLmin=0; uniqLmax=0; uniqLmaxInd=0; multLmax=0; multLmaxN=0; multNminL=0; multNmin=0; multNmax=0; multNmaxL=0;
    chimN=0;
    detectedSampleByte_ = 0xFFu; // Reset detected sample for new read
    extractedCbIdxPlus1_ = 0;    // Reset extracted CB for new read
    extractedUmi24_ = 0;          // Reset extracted UMI for new read
    extractedCbSeq_.clear();      // Reset CB sequence for new read
    hasYAlignment_ = false;       // Reset Y-alignment flag for new read

    for (uint ii=0; ii<P.readNmates; ii++) {//not readNends: this is alignment
        maxScoreMate[ii]=0;
    };
};
