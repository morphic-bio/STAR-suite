#ifndef H_ReadAlign
#define H_ReadAlign

#include "IncludeDefine.h"
#include "Parameters.h"
#include "Transcript.h"
#include "Genome.h"
#include "Stats.h"
#include "OutSJ.h"
#include "Transcriptome.h"
#include "BAMoutput.h"
#include "Quantifications.h"
#include "ChimericDetection.h"
#include "SoloRead.h"
#include "ReadAnnotations.h"
#include "SpliceGraph.h"
#include "ClipMate.h"
#include "TrimQc.h"
#include <zlib.h>

// Forward declarations
class SampleDetector;
class TranscriptQuantEC;
class SlamQuant;
class SlamSnpMask;
class SlamCompat;
namespace libem {
class Transcriptome;
}

#include <time.h>
#include <random>
#include <atomic>
#include <unordered_map>
#include <string>
#include <vector>

class ReadAlign {
    public:
        ReadAlign (Parameters& Pin, Genome &genomeIn, Transcriptome *TrIn, int iChunk,
                   const libem::Transcriptome* libemTr = nullptr);//allocate arrays
        int oneRead();

        Genome &mapGen, &genOut; //mapped-to-genome structure

        uint64 iRead, iReadAll;
        char **Read1;
        char* readName; // Read name for current read

        Stats statsRA; //mapping statistics
        TrimQcCollector trimQc; // per-thread trim QC collector

        istream* readInStream[MAX_N_MATES];
        BAMoutput *outBAMcoord, *outBAMunsorted, *outBAMquant;//sorted by coordinate, unsorted, transcriptomic BAM structure
        std::atomic<uint64_t>* bamRecordIndexPtr; //pointer to global BAM record counter
        fstream chunkOutChimSAM, *chunkOutChimJunction, chunkOutUnmappedReadsStream[MAX_N_MATES], chunkOutFilterBySJoutFiles[MAX_N_MATES];
        fstream chunkOutYReadNames;
        fstream chunkOutYFastqStream[MAX_N_MATES], chunkOutNoYFastqStream[MAX_N_MATES];  // Y/noY FASTQ streams (uncompressed)
        gzFile chunkOutYFastqGz[MAX_N_MATES], chunkOutNoYFastqGz[MAX_N_MATES];  // Y/noY FASTQ streams (gzip-compressed)
        OutSJ *chunkOutSJ, *chunkOutSJ1;

        ostream* outSAMstream;
        uint outBAMbytes; //number of bytes output to SAM/BAM with oneRead
        char *outBAMarray;//pointer to the (last+1) position of the SAM/BAM output array

        uint outFilterMismatchNmaxTotal;
        uint Lread, readLength[MAX_N_MATES], readLengthOriginal[MAX_N_MATES], readLengthPair, readLengthPairOriginal;
        string readBarcodeSeq, readBarcodeQual; //Solo barcode sequence
        
        intScore maxScoreMate[MAX_N_MATES];

        uint32 readFilesIndex;
        
        //quality histogram
        array<array<uint64,256>,MAX_N_MATES> qualHist;
        
        //transcripts (aligns)
        uint nW;
        uint *nWinTr; //number of recorded transcripts per window
        Transcript trA, trA1, *trInit; //transcript, best tr, next best tr, initialized tr
        Transcript ***trAll; //all transcripts for all windows
        Transcript *trArray; //linear array of transcripts to store all of them from all windows
        Transcript **trArrayPointer; //linear array of transcripts to store all of them from all windows

        uint64 nTr; // number of transcripts called
        Transcript *trMult[MAX_N_MULTMAP];//selected alignments - to the mapGen        
        Transcript *trBest;//bset transcripts
        
        struct {
            bool yes;
            uint64 alN;
            Transcript **alMult;//multimapping transcripts - transformed to output genome
            Transcript *alBest;//best align
        } alignsGenOut;
        
        
        Transcript *alignTrAll;//alignments to transcriptome        

        ReadAlign *waspRA; //ReadAlign for alternative WASP alignment
        int waspType, waspType1; //alignment ASE-WASP type and

        ReadAlign *peMergeRA; //ReadAlign for merged PE mates

        ChimericDetection *chimDet;
        void peOverlapChimericSEtoPE(const Transcript *seTrIn1, const Transcript *seTrIn2, Transcript *peTrOut1, Transcript *peTrOut2);

        SoloRead *soloRead; //counts reads per CB per and outputs CB/UMI/gene into file, per thread
        
        SpliceGraph *splGraph;
        
        vector<vector<ClipMate>> clipMates;
        
        // Sample detection (per-thread detector for raw sequence detection)
        SampleDetector *sampleDet_; // Optional sample detector (owned here)
        bool sampleDetReady_;       // Whether detector resources loaded successfully
        uint8_t detectedSampleByte_; // Detected sample index (token) for current read, 0xFF if not detected
        uint32_t extractedCbIdxPlus1_; // Extracted CB whitelist index (1-based) from Solo structures, 0 if not available
        uint32_t extractedUmi24_;     // Extracted UMI (24-bit packed) from Solo structures, 0 if not available
        std::string extractedCbSeq_;   // CB sequence (for Phase 2: lookup resolved CB if cbIdxPlus1==0)
        
        // Phase 2: Ambiguous CB accumulation and resolution
        struct AmbiguousEntry {
            std::vector<uint32_t> candidateIdx;          // 1-based whitelist indices
            std::unordered_map<uint32_t, uint32_t> umiCounts; // UMI24 -> count
            std::string cbSeq;                           // Original CB sequence
            std::string cbQual;                          // CB quality scores (for Bayesian resolution)
            
            AmbiguousEntry() {}
        };
        
        using AmbigKey = uint64_t; // Hash of CB sequence
        std::unordered_map<AmbigKey, AmbiguousEntry> pendingAmbiguous_; // Per-thread accumulation
        std::unordered_map<AmbigKey, uint32_t> resolvedCbByKey_; // Resolved CB indices (AmbigKey -> 1-based whitelist idx)
        
        // Counters for CB resolution statistics
        struct CbResolutionStats {
            uint64_t resolvedImmediate;   // Resolved by Phase 1 (non-ambiguous)
            uint64_t resolvedBayesian;    // Resolved by Phase 2 (Bayesian)
            uint64_t stillAmbiguous;      // Still ambiguous after Bayesian
            uint64_t noMatch;             // No match found
            
            CbResolutionStats() : resolvedImmediate(0), resolvedBayesian(0), stillAmbiguous(0), noMatch(0) {}
        };
        CbResolutionStats cbResolutionStats_;
        
        // Helper: hash CB sequence to AmbigKey (FNV-1a hash)
        static AmbigKey hashCbSeq(const std::string &cbSeq);
        
        // Phase 2: Resolve accumulated ambiguous CBs using Bayesian inference
        void resolveAmbiguousCBs();

	//input,output
        char** outBAMoneAlign;
        uint* outBAMoneAlignNbytes;
        int alignBAM(Transcript const &trOut, uint nTrOut, uint iTrOut, uint trChrStart, uint mateChr, uint mateStart, char mateStrand, int unmapType, bool *mateMap, vector<int> outSAMattrOrder, char** outBAMarray, uint* outBAMarrayN);

    private:
        Parameters& P; //pointer to the parameters, will be initialized on construction

        //quantification
        Transcriptome *chunkTr;
    public:
        TranscriptQuantEC *quantEC;  // EC table for transcript quantification (TranscriptVB mode)
        // SLAM quantification (optional)
        SlamQuant* slamQuant = nullptr;
        const SlamSnpMask* slamSnpMask = nullptr;
        SlamCompat* slamCompat = nullptr;  // SLAM compatibility mode helper
    private:

        //mapping time
        time_t timeStart, timeFinish;

        //random number generators
        std::mt19937 rngMultOrder;//initialize in ReadAlign.cpp
        std::uniform_real_distribution<double> rngUniformReal0to1;//initialize in ReadAlign.cpp

        //input,output

        ostringstream samStreamCIGAR, samStreamSJmotif, samStreamSJintron;
        vector <string> matesCIGAR;

        intScore *scoreSeedToSeed, *scoreSeedBest;
        uint *scoreSeedBestInd, *seedChain, *scoreSeedBestMM;

        bool outFilterBySJoutPass; //true if alignment passed all filters and is output
        bool outSAMfilterPass; //true if alignment passed SAM filter

        //read
        uint64 iMate;
        char readFilter; //Illumina not passed Y/N
        bool revertStrand; //what to do with the strand, according to strandType and iMate
        int readFileType; //file type: 1=fasta; 2=fastq

        vector<string>readNameExtra;

        char dummyChar[4096];
        char** Read0;
        char** Qual0;
        char** readNameMates;

        uint readNmates;
        //split
        uint** splitR;
        uint Nsplit;

//         uint fragLength[MAX_N_FRAG], fragStart[MAX_N_FRAG]; //fragment Lengths and Starts in read space

        //binned alignments
        uintWinBin **winBin; //binned genome: window ID (number) per bin

        //alignments
        uiPC *PC; //pieces coordinates
        uiWC *WC; //windows coordinates
        uiWA **WA; //aligments per window

        int unmapType; //marker for why a read is unmapped
    public:
        bool hasYAlignment_;  // true if any alignment for this read touches Y-chromosome
    private:

        uint mapMarker; //alignment marker (typically, if there is something wrong)
        uint nA, nP, nWall, nUM[2]; //number of all alignments,  pieces, windows, U/M,
        uint *nWA, *nWAP, *WALrec, *WlastAnchor; //number of alignments per window, per window per piece, min recordable length per window
        bool *WAincl; //alginment inclusion mask

        uint *swWinCov, *swWinGleft, *swWinGright, *swWinRleft, *swWinRright; //read coverage per window
        char *swT;

        uint storedLmin, uniqLmax, uniqLmaxInd, multLmax, multLmaxN, multNmin, multNminL, multNmax, multNmaxL;
        intScore maxScore;//maximum alignment score
        bool mateMapped[2];
        
        //old chimeric detection
        uint chimN, chimRepeat, chimStr;
        int chimMotif;
        uint chimRepeat0, chimRepeat1, chimJ0, chimJ1;
        Transcript trChim[MAX_N_CHIMERAS];
        //new chimeric detection
        bool chimRecord; //true if chimeric aligment was detected

        struct {
            bool yes;
            uint nOv;//number of overlapping bases
            uint ovS;//first read base of the overlap
            uint mateStart[2];//mates starts in the merged read
        } peOv;//PE  mates overlap/merge/remap structure
        
        //read annotations
        ReadAnnotations readAnnot;

        #ifdef DEBUG_OutputLastRead
            ofstream lastReadStream;
        #endif

        /////////////////////////////////////////////////////////////////// METHODS
        void resetN();//resets the counters to 0
        void multMapSelect();
        int mapOneRead();
        void mapOneReadSpliceGraph();
        uint maxMappableLength2strands(uint pieceStart, uint pieceLength, uint iDir, uint iSA1, uint iSA2, uint& maxL, uint iFrag);
        void storeAligns (uint iDir, uint Shift, uint Nrep, uint L, uint indStartEnd[2], uint iFrag);

        bool outputTranscript(Transcript *trOut, uint nTrOut, ofstream *outBED);

        uint64 outputTranscriptSAM(Transcript const &trOut, uint nTrOut, uint iTrOut, uint mateChr, uint mateStart, char mateStrand, int unmapType, bool *mateMap, ostream *outStream);

        uint64 outputSpliceGraphSAM(Transcript const &trOut, uint nTrOut, uint iTrOut, ostream *outStream);

        void samAttrNM_MD (Transcript const &trOut, uint iEx1, uint iEx2, uint &tagNM, string &tagMD);

        void outputTranscriptSJ(Transcript const &trOut, uint nTrOut, OutSJ *outStream, uint sjReadStartN );
        string outputTranscriptCIGARp(Transcript const &trOut);
        int createExtendWindowsWithAlign(uint a1, uint aStr); //extends and windows with one alignment
        void assignAlignToWindow(uint a1, uint aLength, uint aStr, uint aNrep, uint aFrag, uint aRstart,bool aAnchor, uint sjA); //assigns one alignment to a window

        void mappedFilter();
        void chimericDetection();
        bool chimericDetectionOld();
        void chimericDetectionOldOutput();
        bool chimericDetectionMult();
        void chimericDetectionPEmerged(ReadAlign &seRa);

        void transformGenome();
        void outputAlignments();
        void calcCIGAR(Transcript const &trOut, uint nMates, uint iExMate, uint leftMate);

        void stitchWindowSeeds (uint iW, uint iWrec, bool *WAexcl, char *R);//stitches all seeds in one window: iW
        void stitchPieces(char **R, uint Lread);

        uint quantTranscriptome (Transcriptome *Tr, uint nAlignG, Transcript **alignG, Transcript *alignT);

        void copyRead(ReadAlign&);
        void waspMap();
        void peOverlapMergeMap();
        void peMergeMates();
        void peOverlapSEtoPE(ReadAlign &seRA);
        
        //output alignments functions
        void outFilterBySJout();
        void outReadsUnmapped();
        void writeFastxRecord(uint imate, bool isY);  // Write FASTQ/FASTA record to Y or noY stream
        void spliceGraphWriteSAM();
        void alignedAnnotation();
        bool slamCollect(const Transcript& trOut, const std::set<uint32_t>& geneIds, double weight, bool isIntronic);
        void writeSAM(uint64 nTrOutSAM, Transcript **trOutSAM, Transcript *trBestSAM);
        void recordSJ(uint64 nTrO, Transcript **trO, OutSJ *cSJ);

};

#endif
