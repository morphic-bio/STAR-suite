#ifndef PARAMETERS_DEF
#define PARAMETERS_DEF

#include "IncludeDefine.h"
#include "InOutStreams.h"
#include "ParameterInfo.h"
#include <map>
#include "TimeFunctions.h"
#include <unistd.h>
#include <signal.h>
#include "ParametersChimeric.h"
#include "ParametersSolo.h"
#include "ParametersClip.h"
#include "ParametersGenome.h"
#include <vector>
#include <array>
#include <unordered_set>
#include <atomic>

// Forward declaration for library format detector
class LibFormatDetector;
class SlamSnpMask;

class Parameters {

    public:
        vector <ParameterInfoBase*> parArray, parArrayInitial;
        vector <string> parameterInputName;

        string commandLine, commandLineFull;

        //version
        string versionGenome;

        //system parameters
        string sysShell; //shell for executing system commands

        // run parameters
        string runMode;
        vector<string> runModeIn;
        int runThreadN;
        mode_t runDirPerm;
        string runDirPermIn; //permission for directores created at run-time
        int runRNGseed; //random number generator seed

        struct {
            int32 type;//0 no restart, 1 no mapping - restart from _STARtmp files
        } runRestart; //restart options - in development
        
        //parameters
        vector <string> parametersFiles;

        //input
        string inputBAMfile;

        //genome
        char genomeNumToNT[6];
        ParametersGenome pGe, pGeOut;

        //binning,windows,anchors
        uint winBinChrNbits, winBinNbits, winAnchorDistNbins, winFlankNbins, winBinN;
        uint winAnchorMultimapNmax; //max number of alignments for anchors
        double winReadCoverageRelativeMin;
        uint winReadCoverageBasesMin;

        //read parameters
        vector <string> readFilesType;
        int readFilesTypeN;
        string readFilesPrefix, readFilesPrefixFinal;
        vector <string> readFilesIn, readFilesInTmp;
        uint32 readFilesN;
        vector <vector <string> > readFilesNames;
        vector <string> readFilesCommand;
        vector <string> readFilesManifest;
               
        string readFilesCommandString; //actual command string
        int readFilesIndex;
        pid_t readFilesCommandPID[MAX_N_MATES];

        uint readMapNumber;
        uint iReadAll;
        uint readNmates, readNends;
        string readMatesLengthsIn;
        uint32 readQualityScoreBase;

        vector <string> readNameSeparator;
        vector <char> readNameSeparatorChar;

        string outSAMreadID;
        bool outSAMreadIDnumber;
        
        //new: structure for readFiles parameters
        struct {
            vector<string> samAttrKeepIn; //input vector of SAM tags to keep, if readFilesType=SAMtag
            std::unordered_set<uint16_t> samAttrKeep;
            bool samAttrKeepAll, samAttrKeepNone;
        } readFiles;

        ParametersClip pClip;

        //cutadapt-style trimming
        string trimCutadapt;
        uint8 trimCutadaptQuality;
        uint32 trimCutadaptMinLength;
        vector <string> trimCutadaptAdapter;
        string trimCutadaptCompat;  // Compatibility mode: "-"/"Off" (default) or "Cutadapt3"

        // trim QC output (FastQC-like)
        string trimQcReport;
        string trimQcJson;
        string trimQcHtml;
        uint64 trimQcMaxReads = 0;
        bool trimQcEnabled = false;

        //align parameters
        uint alignSJoverhangMin,alignSJDBoverhangMin,alignSplicedMateMapLmin; //min SJ donor/acceptor length
        double alignSplicedMateMapLminOverLmate;
        uint alignWindowsPerReadNmax; //max number of alignment windows per read
        uint alignTranscriptsPerWindowNmax; //maximum number of transcripts recorded per window
        uint alignTranscriptsPerReadNmax;   //max number of alignments per read
        uint alignIntronMin;//min length to call a gap an intron
        uint alignIntronMax;//max length to call
        uint alignMatesGapMax;//max gap between the mates (if paired-end)
        vector <int32> alignSJstitchMismatchNmax;

        //         struct {
        //             string strandString;
        //             int32 strand;
        //         } pReads;

        struct {
            string in;
            bool yes;
        } alignSoftClipAtReferenceEnds;

        struct {
            string in;
            bool ext[2][2];
        } alignEndsType;

        struct {
            vector<string> in;
            int nBasesMax;
            bool concordantPair;
        } alignEndsProtrude;

        struct {
            string in;
            bool flushRight;
        } alignInsertionFlush;


        //seed parameters
        uint seedMultimapNmax; //max number of multiple alignments per piece
        uint seedSearchLmax; //max length of the seed
        uint seedPerReadNmax; //max number of pieces per Read
        uint seedPerWindowNmax; //max number of aligns per window
        uint seedNoneLociPerWindow; //max number of aligns from one piece per window
        uint seedSearchStartLmax;
        double seedSearchStartLmaxOverLread; //length of split start points
        uint64 seedSplitMin, seedMapMin;

        //chunk parameters
        uint chunkInSizeBytes,chunkInSizeBytesArray,chunkOutBAMsizeBytes;

        //output
        string outFileNamePrefix, outStd;
        string outTmpDir, outTmpKeep;
        string outLogFileName;

        //SAM output
        string outBAMfileCoordName, outBAMfileUnsortedName, outQuantBAMfileName;
        string outBAMfileUnsortedSoloTmpName;
        string samHeader, samHeaderHD, samHeaderSortedCoord, samHeaderExtra;
        string outSAMmode,  outSAMorder, outSAMprimaryFlag;
        vector<string> outSAMattributes, outSAMheaderHD, outSAMheaderPG;
        vector<string> outSAMattrRGline,outSAMattrRGlineSplit,outSAMattrRG;
        uint outSAMmultNmax,outSAMattrIHstart;
        string outSAMheaderCommentFile;
        int outSAMmapqUnique;

        struct {
            string in;
            uint32 type;
        } outSAMstrandField;

        int outSAMtlen;

        struct {bool NH,HI,AS,NM,MD,nM,jM,jI,RG,XS,rB,vG,vA,vW,ha,ch,MC,CR,CY,UR,UY,CB,UB,GX,GN,gx,gn,sM,sS,sQ,cN,sF,ZG,ZX,ZS;} outSAMattrPresent, outSAMattrPresentQuant;

        vector <int> outSAMattrOrder, outSAMattrOrderQuant;
        int outBAMcompression;
        vector <string> outSAMtype;
        bool outBAMunsorted, outBAMcoord, outSAMbool;
        bool outBAMunsortedUseSoloTmp; // whether to use solo tmp file for unsorted BAM output
        uint32 outBAMcoordNbins;
        uint32 outBAMsortingBinsN;//user-defined number of bins for sorting
        string outBAMsortTmpDir;
        string outBAMsortMethod;  // "star" (default) or "samtools"
        
        //Y-chromosome BAM split
        string emitNoYBAM;  // raw CLI: yes|no - split BAM output by Y-chromosome alignments
        bool emitNoYBAMyes;  // resolved: true if enabled
        string emitYReadNames;  // raw CLI: yes|no - emit Y-read name list for FASTQ filtering
        bool emitYReadNamesyes; // resolved: true if enabled
        string keepBAM;      // raw CLI: yes|no - when emitNoYBAM is enabled, also emit primary BAM
        bool keepBAMyes;    // resolved: true if enabled
        string outBAMfileYName, outBAMfileNoYName;  // derived output paths for Y/noY BAMs
        string noYOutput, YOutput;  // user-specified override paths
        string YReadNamesOutput; // user-specified override path for Y-read names list
        string outYReadNamesFile; // derived output path for Y-read names list
        
        //Y-chromosome FASTQ emission
        string emitYNoYFastq;  // raw CLI: yes|no - emit Y/noY FASTQ files
        bool emitYNoYFastqyes; // resolved: true if enabled
        string emitYNoYFastqCompression;  // raw CLI: gz|none - compression for FASTQ output
        string YFastqOutputPrefix, noYFastqOutputPrefix;  // user-specified output prefixes
        string outYFastqFile[MAX_N_MATES], outNoYFastqFile[MAX_N_MATES];  // derived output paths per mate

//         string bamRemoveDuplicatesType;
//         uint bamRemoveDuplicatesMate2basesN;
        struct {
            string mode;
            bool yes;
            bool markMulti;
            uint mate2basesN;
        } removeDuplicates;

        int outBAMsortingThreadN, outBAMsortingThreadNactual;
        uint64 *outBAMsortingBinStart; //genomic starts for bins for sorting BAM files
        uint16 outSAMflagOR, outSAMflagAND;

        struct {
            vector <string> mode;
            bool yes;
            bool within;//output unmapped reads within SAM/BAM files
            bool keepPairs;//keep mates together
        } outSAMunmapped;

        struct {
            vector <string> mode;
            bool yes;
            bool KeepOnlyAddedReferences;
            bool KeepAllAddedReferences;
        } outSAMfilter;

        struct {
            string mode;
            bool random;
        } outMultimapperOrder;

        struct {
            bool yes;
            uint NbasesMin;
            double MMp;
        } peOverlap;

        string outReadsUnmapped;
        int outQSconversionAdd;
        string outFileTmp;

        //output filtering
        uint outFilterMismatchNmax;
        double outFilterMismatchNoverLmax, outFilterMismatchNoverReadLmax; //max proportion of all MM within all bases

        uint outFilterMatchNmin,outFilterMultimapNmax;//min number of matches
        double outFilterScoreMinOverLread, outFilterMatchNminOverLread;//normalzied to read length
        intScore outFilterScoreMin,outFilterMultimapScoreRange;//min score to output
        string outFilterIntronMotifs,outFilterIntronStrands;
        string outFilterType; //type of filtering
        int outFilterBySJoutStage; //indicates the stage of filtering by SJout

        struct {
            vector<string> type;
            bool yes;
        } outSJ;
        
        //output filtering SJs
        string outSJfilterReads;
        vector <int32> outSJfilterCountUniqueMin, outSJfilterCountTotalMin;
        vector <int32> outSJfilterOverhangMin;
        vector <int32> outSJfilterDistToOtherSJmin; //min allowed distance to other SJ's donor/acceptor
        vector <int32> outSJfilterIntronMaxVsReadN;

        //wiggle output
        vector <string> outWigType, outWigStrand, outWigNorm;
        string outWigReferencesPrefix;
        struct {
            bool yes;
            bool strand;
            int type;
            int format;
            int norm;
        } outWigFlags;

        //2-pass
//         uint twoPass.pass1readsN, twoPass.sjLimit;
//         string twoPass.dir,twopassSJpass1file;
        struct {
            bool yes; //true in 2-pass mode
            bool pass2; //true if now running the 2nd pass
            uint pass1readsN;
            int pass1readsN_par;
            string dir;
            string pass1sjFile;
            string mode;
        } twoPass;

        //inserting junctions on the fly
        struct {
            bool yes; //insert?
            bool pass1;//insert on the 1st pass?
            bool pass2;//insert on the 2nd pass?
            string outDir;
        } sjdbInsert;

        //storage limits
        uint64 limitGenomeGenerateRAM;
        vector<uint64> limitIObufferSize; //max size of the in/out buffer, bytes
        uint64 limitOutSAMoneReadBytes;
        uint64 limitOutSJoneRead, limitOutSJcollapsed;
        uint64 limitBAMsortRAM;
        uint64 limitSjdbInsertNsj;
        uint64 limitNreadsSoft;

        // penalties
        intScore scoreGap, scoreGapNoncan, scoreGapGCAG, scoreGapATAC, scoreDelBase, scoreDelOpen, scoreInsBase, scoreInsOpen;
        intScore scoreStitchSJshift;//Max negative score when
        double scoreGenomicLengthLog2scale;

        //quantification parameters
        //input

        struct {
          bool yes=false; //if any quantification is done
          vector <string> mode; //quantification mode input string

            struct {
                bool yes=false;
                bool bamYes;
                bool indel;
                bool softClip;
                bool singleEnd;
                int bamCompression;
                string output;
            } trSAM;

            struct {
                bool yes=false;
                string outFile;
            } geCount;

            struct {
                bool yes=false;
                bool vb=true;           // Use VB (true) or EM (false)
                bool gcBias=false;      // Enable GC bias correction
                int gcBiasInt=0;        // Command-line flag for gcBias (0/1)
                int quantVBemInt=0;     // Command-line flag: if 1, use EM instead of VB
                double vbPrior=0.01;    // Dirichlet prior
                string outFile;         // Output file path
                // Gene-level output
                bool geneOutput=true;   // Output gene-level quant (default: yes)
                int geneOutputInt=1;    // CLI flag (0/1)
                string outFileGene;     // Gene output file path
                string outFileGeneTximport; // tximport-style gene output
                // Gene summarization mode
                string genesModeStr = "Legacy";  // "Legacy" or "Tximport"
                bool genesTximport = false;      // If true, use lengthScaledTPM counts
                
                // Library format detection and EC building parameters
                string libType = "A";      // A=auto, IU, ISF, ISR, U
                int autoDetectWindow = 1000;
                string traceFile;          // Empty = disabled
                int traceLimit = 0;        // 0 = unlimited (trace all reads); >0 = limit to N reads
                int preBurninFrags = 5000; // Pre-burn-in fragment count threshold (Salmon default: 5000)
                int miniBatchSize = 1000;  // Mini-batch size for processed reads counter (Salmon default: 1000)
                
                // Error model mode: auto|cigar|as|off (default: auto)
                string errorModelMode = "auto";  // auto=use CIGAR if available else AS, cigar=CIGAR only, as=AS only, off=disabled
                
                // Internal state (NOT CLI params, NOT LibraryFormat type)
                uint8_t detectedLibFormatId = 0;  // Stores LibraryFormat::typeId()
                bool detectionComplete = false;
                bool inDetectionMode = false;     // True during detection pass
                
                // Shared detector instance (created before detection, accessed by TranscriptQuantEC)
                LibFormatDetector* libFormatDetector = nullptr;  // Raw ptr, owned by STAR.cpp
            } transcriptVB;

            struct {
                bool yes=false;
                int modeInt=0;              // CLI flag (0/1)
                double errorRate=0.001;     // Background error rate
                double convRate=0.05;       // Conversion rate for labeled reads
                string snpBed;              // Optional SNP BED for prefilter
                int snpDetectInt=0;         // CLI flag (0/1) for internal SNP detection
                bool snpDetect=false;       // Enable internal SNP detection + masking
                double snpDetectFrac=-1.0;  // Mismatch fraction threshold: <=0=auto, >0=explicit
                string strandnessStr="Unspecific"; // Unspecific | Sense | Antisense
                uint8_t strandness=0;       // 0=Unspecific, 1=Sense, 2=Antisense
                string weightModeStr="Alignments"; // Alignments | Uniform
                uint8_t weightMode=0;       // 0=Alignments (1/nTr), 1=Uniform (1.0)
                string debugGeneList;       // File with gene IDs/names to instrument
                string debugReadList;       // File with read names to instrument
                string debugOutPrefix;      // Prefix for debug outputs
                int debugMaxReads=2000;     // Max debug read records (0 disables)
                string debugSnpLoc;         // Optional SNP site debug locus: <chrom>:<pos1> (1-based)
                int debugSnpWindow=2;       // Debug window around locus (pos +/- window)
                unordered_set<string> debugGenes; // Parsed gene list
                unordered_set<string> debugReads; // Parsed read list
                bool debugEnabled=false;    // Debug logging enabled
                string outFile;             // Output file path
                SlamSnpMask* snpMask=nullptr; // Loaded SNP mask (if any)
                
                // Compatibility mode
                string compatModeStr = "none";    // CLI: none|gedi
                int compatIntronicInt = -1;       // Parsing: -1=not set, 0=disable, 1=enable
                int compatLenientOverlapInt = -1; // Parsing: -1=not set, 0=disable, 1=enable
                int compatOverlapWeightInt = -1;  // Parsing: -1=not set, 0=disable, 1=enable
                int compatIgnoreOverlapInt = -1;  // Parsing: -1=not set, 0=disable, 1=enable
                bool compatIntronic = false;      // GEDI-style intronic predicate
                bool compatLenientOverlap = false;// 50% overlap + SJ concordance
                bool compatOverlapWeight = false; // Divide weight by gene count
                bool compatIgnoreOverlap = false; // Skip PE overlap positions
                int compatTrim5p = 0;             // 5' trim guard
                int compatTrim3p = 0;             // 3' trim guard
                
                // Auto-trim variance analysis (segmented regression on T→C stdev curve)
                string autoTrimMode = "";         // ""=disabled, "variance"=variance-based trimming
                string trimScope = "first";        // "first"=shared from first file, "per-file"=per-file
                string trimSource = "";           // Path to file for computing shared trims (overrides first file)
                int autoTrimMaxReads = 100000;    // Max reads to scan for variance (0=unlimited)
                int autoTrimMinReads = 1000;      // Minimum reads required for auto-trim
                int autoTrimSmoothWindow = 5;     // Median smoothing window for stdev curve
                int autoTrimSegMinLen = 3;        // Minimum segment length for segmented regression
                int autoTrimMaxTrim = 15;         // Maximum trim at either end
                string slamQcJson = "";           // Path for QC JSON output (empty=auto-generate)
                string slamQcHtml = "";            // Path for QC HTML output (empty=auto-generate)
                string slamQcReport = "";          // Prefix for comprehensive QC report (JSON + HTML, empty=disabled)
                int grandSlamOut = 1;              // --slamGrandSlamOut (0/1, default: 1)
                int autoTrim5p = 0;               // Auto-computed 5' trim (0=not computed)
                int autoTrim3p = 0;                // Auto-computed 3' trim (0=not computed)
                bool autoTrimComputed = false;    // Whether auto-trim has been computed
                uint32_t autoTrimFileIndex = 0;    // File index where auto-trim was computed
                std::vector<double> varianceStddevTcRate; // Variance curve from detection pass (optional)
                
                // Global SNP error rate estimation (from auto-trim detection pass)
                double snpErrEst = 0.0;           // Estimated T→C error rate (p_err)
                double snpErrUsed = 0.001;        // Error rate used for GEDI (with fallback)
                string snpErrFallbackReason = ""; // Reason for fallback (empty if no fallback)
                double snpErrMinThreshold = 0.001; // Minimum threshold for p_err (fallback if below)
                
                // Auto-trim detection mode (single-thread first pass with rewind)
                uint64_t autoTrimBufferReads = 1000000;  // Max reads for detection pass (default 1M) - legacy name for compat
                int autoTrimDetectionReads = 100000;     // Max reads for detection pass (used as readMapNumber limit)
                bool autoTrimReplayDone = false;         // Whether second pass has completed
                bool autoTrimDetectionPass = false;      // True during single-threaded detection pass
                bool perFileProcessing = false;          // True when processing files one at a time (for trimScope=per-file)
                int currentFileIndex = 0;                // Current file being processed in per-file mode
                int totalFileCount = 0;                  // Total number of input files
                int skipToFileIndex = -1;                // Skip reads until reaching this file index (-1=disabled)
                // Optional dump for external re-quant
                string dumpBinary;                       // --slamDumpBinary (path to dump)
                uint64_t dumpMaxReads = 1000000;         // --slamDumpMaxReads (max reads to dump)
                string dumpWeights;                      // --slamDumpWeights (path to weight sidecar)
                string dumpWeightsModeStr = "dump";      // --slamDumpWeightsMode (dump|vbGene)
                uint8_t dumpWeightsMode = 0;             // 0=dump, 1=vbGene
            } slam;

            struct {
                // Mask source flags (precedence: maskIn > vcfIn > buildFastqs)
                string maskIn;                // --slamSnpMaskIn (existing mask, highest priority)
                string vcfIn;                 // --slamSnpMaskVcfIn (VCF/gVCF SNP mask input)
                string vcfSample;             // --slamSnpMaskVcfSample (sample name for GT mode)
                string vcfMode = "gt";        // --slamSnpMaskVcfMode (gt|any)
                string vcfFilter = "pass";    // --slamSnpMaskVcfFilter (pass|all)
                string buildFastqsFofn;      // --slamSnpMaskBuildFastqs (FOFN for pre-pass)
                string buildBam;              // --slamSnpMaskBuildBam (BAM for pre-pass, alternative to FASTQ)
                int buildOnlyInt = 0;         // --slamSnpMaskOnly (int for parsing)
                bool buildOnly = false;       // derived from buildOnlyInt
                
                // Output paths
                string bedOut;                // --slamSnpMaskBedOut
                string summaryOut;           // --slamSnpMaskSummaryOut
                string bamOut;               // --slamSnpMaskBamOut (optional, default empty)
                
                // Compatibility mode: "gedi" mirrors GEDI thresholds, "" = default STAR
                string compat;                // --slamSnpMaskCompat (gedi|"")
                // What constitutes an "alt" observation for SNP masking:
                // - conv: only T->C (and A->G on opposite strand) conversions
                // - any: any mismatch (GEDI-style), regardless of base
                string kMode;                 // --slamSnpMaskKMode (conv|any)
                
                // Model selection
                string model = "binom";        // --slamSnpMaskModel (binom|em, default: binom)
                
                // Binomial model parameters
                // Initial values are sentinels (-2.0 for double, -1 for int); actual defaults are set
                // in finalization based on compat mode. This allows detecting user-set values.
                double pval = -2.0;           // --slamSnpMaskPval (sentinel: -2, STAR default: 0.001)
                double minTcRatio = -2.0;     // --slamSnpMaskMinTcRatio (sentinel: -2, STAR default: 0.3)
                double err = -1.0;            // --slamSnpMaskErr (-1 = use computed snp_err_used)
                
                // Coverage/count parameters (sentinel: -1, parsed from parametersDefault)
                int32_t minCov = -1;          // --slamSnpMaskMinCov (STAR default: 20, GEDI: 6)
                int32_t minAlt = -1;          // --slamSnpMaskMinAlt (STAR default: 3, GEDI: 1)
                
                // EM model parameters
                double posterior = 0.99;      // --slamSnpMaskPosterior
                uint32_t maxIter = 50;        // --slamSnpMaskMaxIter
                double convergeRelLL = 1e-7;  // --slamSnpMaskConvergeRelLL
                
                // Artifact filters (sentinel: -1, parsed from parametersDefault)
                int32_t junctionFlank = -1;   // --slamSnpMaskJunctionFlank (STAR: 6, GEDI: 0)
                int32_t indelFlank = -1;      // --slamSnpMaskIndelFlank (STAR: 3, GEDI: 0)
                int32_t minMapQ = -1;         // --slamSnpMaskMinMapQ (STAR: 20, GEDI: 0)
                int32_t minBaseQ = -1;        // --slamSnpMaskMinBaseQ (STAR: 20, GEDI: 0)
            } slamSnpMask;

            struct {
                bool yes=false;
            } geneFull;
          
            struct {
                bool yes=false;
            } geneFull_Ex50pAS;

            struct {
                bool yes=false;
            } geneFull_ExonOverIntron;

            struct {
                bool yes=false;
            } gene;
          
        } quant;
        
        // Global atomic counter for FLD observations (fragments that contribute to FLD)
        // Used for pre-burn-in gating: only increments when fragment length is valid
        // and alignment is compatible (matches Salmon's behavior)
        // Declared at class level because static members can't be in unnamed structs
        static std::atomic<uint64_t> global_fld_obs_count;

        //variation parameters
        struct {
            bool yes=false;
            bool heteroOnly=false;
            string vcfFile;
        } var;

        struct {
            bool yes=false;
            bool SAMtag;
            string outputMode;
        } wasp;

        //solo
        ParametersSolo pSolo;

        //chimeric
        ParametersChimeric pCh;

        //splitting
        uint maxNsplit;

        //not really parameters, but global variables:
        array<vector<uint64>,2> sjAll;
        uint64 sjNovelN, *sjNovelStart, *sjNovelEnd; //novel junctions collapased and filtered

    ////////////////////// CLEAN-UP needed
    InOutStreams *inOut; //main input output streams

    uint Lread;

    Parameters();
    int readParsFromFile(ifstream*, ofstream*, int); //read parameters from one file
    int readPars(); // read parameters from all files
    int scanOneLine (string &lineIn, int inputLevel, int inputLevelRequested);
    void scanAllLines (istream &streamIn, int inputLevel, int inputLevelRequested);
    void inputParameters (int argInN, char* argIn[]); //input parameters: default, from files, from command line
    void openReadsFiles();
    void readFilesInit();
    void closeReadsFiles();
    void readSAMheader(const string readFilesCommandString, const vector<string> readFilesNames);
    void samAttributes();
    void samAttrRequiresBAM(bool attrYes, string attrTag);
};
#endif  // Parameters.h
