#ifndef CODE_ParametersSolo
#define CODE_ParametersSolo

#include <array>
#include <unordered_map>
#include <mutex>
#include <memory>


#include "IncludeDefine.h"
#include "SoloBarcode.h"
#include "SoloFeatureTypes.h"

class Parameters;
class ParametersSolo;

class UMIdedup {
public:
    const static uint32 tN = 6;
    array<string,tN> typeNames { {"NoDedup", "Exact", "1MM_All", "1MM_Directional", "1MM_CR", "1MM_Directional_UMItools"} };
    enum typeI : int32 { NoDedup=0, Exact=1, All=2, Directional=3, CR=4, Directional_UMItools=5 };
    
    struct {
        uint32_t N;
        array<bool,tN> B;
        bool &NoDedup=B[0], &Exact=B[1], &All=B[2], &Directional=B[3], &CR=B[4], &Directional_UMItools=B[5]; 
    } yes;

    struct {
        //uint32_t N;
        array<uint32_t,tN> I;
        uint32_t &NoDedup=I[0], &Exact=I[1], &All=I[2], &Directional=I[3], &CR=I[4], &Directional_UMItools=I[5];
        uint32_t main; //index for SAM/stats/filtering output
    } countInd; //index in the countCellGennUMI
    
    vector<string> typesIn; //UMIdedup types from user options
    vector<int32> types; //the above converted to typeI numbers
    int32 typeMain; //the type to be used in SAM/stats/filtering output - for now just types[0]
    
    void initialize(ParametersSolo *pS);
    
//protected:
//    int it;
};

class MultiMappers {
public:
    const static uint32 tN = 5;
    array<string,tN> typeNames { {"Unique", "Uniform", "Rescue", "PropUnique", "EM"} };
    enum typeI : int32 { Unique=0, Uniform=1, Rescue=2, PropUnique=3, EM=4 };
    
    struct {
        bool multi; //if multimappers are requested
        uint32_t N;
        array<bool,tN> B;
        bool &Unique=B[0], &Uniform=B[1], &Rescue=B[2], &PropUnique=B[3], &EM=B[4] ;
    } yes;

    struct {
        //uint32_t N;
        array<uint32_t,tN> I;
        uint32_t &Unique=I[0], &Uniform=I[1], &Rescue=I[2], &PropUnique=I[3], &EM=I[4];
        uint32_t main; //index for SAM/stats/filtering output
    } countInd; //index in the countCellGennUMI
    
    vector<string> typesIn; //UMIdedup types from user options
    vector<int32> types; //the above converted to typeI numbers
    int32 typeMain; //the type to be used in SAM/stats/filtering output - for now just types[0]
    
    void initialize(ParametersSolo *pS);
};

class ParametersSolo {
public:
    Parameters *pP;
    bool yes;
    
    ~ParametersSolo();

    //chemistry, library etc
    string typeStr;
    enum SoloTypes : int32 {None=0, CB_UMI_Simple=1, CB_UMI_Complex=2, CB_samTagOut=3, SmartSeq=4};
    SoloTypes type;
    string strandStr;
    int32 strand;   
    
    uint32 barcodeRead, barcodeReadIn;//which read is the barcode read = 0,1,2?
    uint32 barcodeStart, barcodeEnd;//start/end of barcode sequence on barcodeRead
    bool barcodeReadSeparate;
    
    //simple barcodes
    uint32 cbS, cbL; //cell barcode start,length
    uint32 umiS, umiL; //umi start,length
    uint32 bL, cbumiL; //total barcode sequene length, CB+UMI length. Former does may not be equal to the latter

    vector<string> cbPositionStr;
    string umiPositionStr;
    
    //complex barcodes    
    vector<SoloBarcode> cbV;
    SoloBarcode umiV; //single UMI
    bool adapterYes; //anchor?  
    string adapterSeq; //anchor sequence
    uint32 adapterMismatchesNmax;//max number of mismatches in the anchor
    
    //input from SAM files
    vector<string> samAtrrBarcodeSeq, samAtrrBarcodeQual;
    
    //whitelist - general
    uint64 cbWLsize;
    bool cbWLyes;
    vector<string> soloCBwhitelist;
    vector <uint64> cbWL;    
    vector<string> cbWLstr;
    
    // CB Corrector instance (shared across threads)
    std::shared_ptr<class CbCorrector> cbCorrector;
    
    // Phase 2: Resolved ambiguous CB mappings (AmbigKey -> 1-based whitelist idx)
    // Populated by each thread after resolution, merged for global access
    std::unordered_map<uint64_t, uint32_t> resolvedCbByKey; // Thread-safe access via mutex if needed
    
    // Ambiguous CB mappings (hash of CB seq -> vector of 1-based whitelist neighbor indices)
    std::unordered_map<uint64_t, std::vector<uint32_t>> ambiguousCbByKey;
    
    MultiMappers multiMap;
    
    //features
    vector<string> featureIn;//string of requested features
    vector<uint32> features;
    uint32 nFeatures;//=features.size(), number of requested features
    
    int32 featureFirst; //which feature is the first on the list
    array<bool,SoloFeatureTypes::N> featureYes; //which features are requested
    array<bool,SoloFeatureTypes::N> readInfoYes;//which features will need readInfo (for now only Gene and GeneFull)
    array<bool,SoloFeatureTypes::N> readIndexYes;//which features will need recording of readIndex (for now only Gene and GeneFull, for multimappers)
    array<bool,SoloFeatureTypes::N> readStatsYes;//which features will need output of read statistics
    array<int32,SoloFeatureTypes::N> featureInd;//index of each feature - skips unrequested features
    
    //filtering
    char QSbase,QSmax;//quality score base and cutoff

    #ifdef MATCH_CellRanger
    double cbMinP;//for CBs with non-exact matching to WL, min posterior probability
    #else
    float cbMinP;//for CBs with non-exact matching to WL, min posterior probability
    #endif
    
    //cell filtering
    struct {
        vector<string> type;
        uint32 topCells;
        
        struct {
            double nExpectedCells;
            double maxPercentile;
            double maxMinRatio;
        } knee;
        
        struct {
            uint32 indMin, indMax; //min/max cell index, sorted by UMI counts,for empty cells
            uint32 umiMin;
            double umiMinFracMedian;
            uint32 candMaxN;
            double FDR;
            uint32 simN;
        } eDcr;//EmptyDrops-CellRanger
        
    } cellFilter;
      
    //CBtype
    struct {
        string typeString;
        int32 type;
        std::unordered_map<string,uint32> strMap;
        std::mutex *strMtx;
    } CBtype;

    //CB match
    struct {
        string type;
        bool mm1; //1 mismatch allowed
        bool mm1_multi; //1 mismatch, multiple matches to WL allowed
        bool oneExact; //CBs require at least one exact match
        bool mm1_multi_pc; //use psedocounts while calculating probabilities of multi-matches
        bool mm1_multi_Nbase; //allow multimatching to WL for CBs with N-bases
        bool EditDist_2; //allow EditDistance <=3
    } CBmatchWL;
    
    //UMIdedup
    UMIdedup umiDedup;
    
    //multi-gene umi
    struct {
        vector<string> type;
        bool MultiGeneUMI       = false;
        bool MultiGeneUMI_All   = false;
        bool yes                = false; //true for non-CR
        bool MultiGeneUMI_CR    = false;
    } umiFiltering;
    
    //clusters
    string clusterCBfile;
    
    //output
    vector<string> outFileNames;    
    struct {
    	string featuresGeneField3;
    } outFormat;

    bool samAttrYes;//post-processed SAM attributes: error-corrected CB and UMI
    int32 samAttrFeature;//which feature to use for error correction
    
    //two-pass unsorted CB/UB injection
    string addTagsToUnsortedStr;//string parameter input
    bool addTagsToUnsorted;//whether to add CB/UB tags to unsorted BAM via two-pass mode
    
    //tag table export (always enabled; legacy CLI removed)
    bool writeTagTableEnabled = false;//whether to export CB/UB tags to sidecar table
    string writeTagTablePath;//resolved absolute path
    
    // keys.bin (CR v1) emission
    string writeKeysBinStr = "no";   // raw CLI value: yes|no
    bool writeKeysBin = false;         // whether to emit keys.bin during pass-2
    string writeKeysBinPath;           // <prefix>Aligned.out.keys.bin
    string writeKeysBarcodesPath;      // <prefix>Aligned.out.keys.bin.barcodes

    // CR-compatible keys mode (handoff)
    string keysCompatStr = "raw"; // raw|cr
    enum KeysCompat : int32 { KeysRaw=0, KeysCR=1 } keysCompat = KeysRaw;

    // Probe list and sample detection resources
    string probeListPath;             // path to probe_list.txt
    string removeDeprecatedStr;      // --removeDeprecated Yes/No (remove deprecated entries from probe lists)
    bool removeDeprecated;           // Converted from removeDeprecatedStr
    string sampleWhitelistPath;       // path to sample whitelist TSV
    string sampleProbesPath;          // path to sample probes file
    uint32 sampleProbeOffset = 68;    // default offset
    string sampleSearchNearbyStr = "yes"; // yes|no
    bool sampleSearchNearby = true;
    string sampleStrictMatchStr = "no";   // yes|no
    bool sampleStrictMatch = false;
    bool sampleRequireMatch = false;       // drop reads with unmatched sample when true

    // Filters for CR-compatible keys
    int mapqThreshold = 255;          // MAPQ threshold (0-255)
    int nmMax = -1;                   // NM max; -1 disables
    double mmRateMax = 0.12;          // NM/read_len max
    // MAPQ filter mode: 0=off, 1=genomic-only, 2=all alignments
    string mapqModeStr = "off";
    enum MapqMode : int32 { MapqOff = 0, MapqGenomicOnly = 1, MapqAll = 2 } mapqMode = MapqOff;

    // Harmonization & barcodes mode
    string uniqueHarmonizeStr = "no";     // yes|no (default yes in cr)
    bool uniqueHarmonize = false;
    string barcodesObservedOnlyStr = "no"; // yes|no (default yes in cr)
    bool barcodesObservedOnly = false;

    // UMI correction (clique deduplication)
    string cellsAllowPath;                // path to cells_allow.tsv
    string sampleAssignmentsPath;         // path to assignments.tsv
    string umiCorrectionModeStr = "none"; // none|clique
    int umiCorrectionMode = 0;            // 0=none, 1=clique
    string umiCorrectionUseTagsStr = "no"; // yes|no (default: no, conservative grouping)
    bool umiCorrectionUseTags = false;     // resolved: include tag/sample in clique grouping
    int umiMinCount = 1;                  // minimum UMI count threshold
    double umiRatioThresh = 2.0;          // ratio threshold for clique merging
    int maxComponentSize = 4;             // maximum component size for clique
    string assignAmbiguousStr = "drop";   // drop|keep
    string assignAmbiguous = "drop";      // resolved value
    
    // Inline key replayer (embedded cr_key_replayer)
    string useInlineReplayerStr = "auto";    // raw CLI value: no|yes|auto (default: auto for "auto when skip + clique")
    bool useInlineReplayer = false;          // resolved value: true if enabled
    
    // Inline CB correction (process_features-style streaming correction)
    string inlineCBCorrectionStr = "no";     // raw CLI value: no|yes|auto (default: no)
    bool inlineCBCorrection = false;         // resolved value: true if enabled

    // Inline hash mode
    // Replaces temp solo* stream files with per-thread khash_t(cg_agg) hash tables
    // Default behavior: enabled when --flex yes, disabled when --flex no (default)
    // This ensures non-Flex runs use legacy STARsolo output (MEX format)
    string inlineHashModeStr = "";            // raw CLI value: yes|no|auto (default: auto - enabled for Flex, disabled otherwise)
    bool inlineHashMode = false;              // resolved value: true if enabled (gated behind Flex flag by default)

    //skip processing
    string skipProcessingStr = "no";//raw CLI value
    bool skipProcessing = false;//whether to skip Solo counting, matrix construction, and cell filtering

    // Flex omnibus flag - enables full Flex pipeline with production defaults
    string flexModeStr = "no";       // raw CLI: yes|no (default: no)
    bool flexMode = false;           // resolved: true if enabled
    
    // FlexFilter inline integration
    string runFlexFilterStr = "no";  // raw CLI: yes|no (default: no)
    bool runFlexFilter = false;      // resolved: true if enabled
    string flexFilterFatalOnErrorStr = "no";  // raw CLI: yes|no
    bool flexFilterFatalOnError = false;       // resolved: true if fail-fast

    // Expected cells configuration (one of these required when runFlexFilter=true)
    uint32_t flexFilterTotalExpected = 0;        // Total expected cells across all tags (deprecated, use below)
    uint32_t flexFilterExpectedCellsTotal = 0;   // --soloFlexExpectedCellsTotal: total cells across ALL tags
    uint32_t flexFilterExpectedCellsPerTag = 0;  // --soloFlexExpectedCellsPerTag: cells PER tag (multiplied by numTags)
    bool flexFilterExpectedPerTagMode = false;   // true if using per-tag mode

    // Optional paths
    string flexFilterAllowedTagsPath;      // Path to allowed tags TSV
    string flexFilterOutputPrefix;         // Output prefix (default: <outPrefix>/flexfilter/)

    // OrdMag parameters (map to FlexFilter::Config::ordmagParams)
    uint32_t flexFilterOrdmagNsamples = 0;        // Maps to nExpectedCells (use default if 0)
    double flexFilterOrdmagUmiMin = 0.0;          // Maps to umiMin (use default if 0)
    double flexFilterOrdmagTargetPct = 0.0;      // Maps to maxPercentile (use default if 0)

    // EmptyDrops parameters (map to FlexFilter::Config::emptydropsParams)
    uint32_t flexFilterEdLower = 0;               // Maps to indMin (use default if 0)
    uint32_t flexFilterEdNiters = 0;              // Maps to simN (use default if 0)
    double flexFilterEdFdrThreshold = 0.0;         // Maps to FDR (use default if 0)
    uint32_t flexFilterEdMaxTotalBuckets = 0;     // Maps to maxTotalBuckets (0 = disabled)

    // Occupancy parameters (map to FlexFilter::Config)
    uint32_t flexFilterTotalPartitions = 0;       // Use default if 0
    double flexFilterRecoveryFactor = 0.0;        // Use default if 0
    double flexFilterOccupancyPercentile = 0.0;   // Use default if 0
    uint32_t flexFilterLowUmiThreshold = 0;        // Use default if 0

    // Simple EmptyDrops (fallback filter) flags
    string flexFilterUseSimpleEDStr = "no";      // raw CLI: yes|no (default: no = fallback only)
    bool flexFilterUseSimpleED = false;          // resolved: if true, force enable Simple ED
    uint32_t flexFilterSimpleEDMinRescues = 50;  // Min ED rescues before fallback triggers
    uint32_t flexFilterSimpleEDMinAmbient = 100; // Min ambient cells before fallback triggers
    uint32_t flexFilterSimpleEDMinCandidates = 100; // Min candidates before fallback triggers

    // Debug flags
    string flexFilterDebugOutputDir;              // Debug output directory
    string flexFilterDebugTagLogStr = "no";      // raw CLI: yes|no
    bool flexFilterDebugTagLog = false;          // resolved
    string flexFilterDisableOccupancyStr = "no"; // raw CLI: yes|no
    bool flexFilterDisableOccupancy = false;     // resolved
    string flexFilterInvariantChecksStr = "no"; // raw CLI: yes|no
    bool flexFilterInvariantChecks = false;      // resolved

    // Output options
    string flexFilterKeepCBTagStr = "no";        // raw CLI: yes|no (default: no, strip to 16bp)
    bool flexFilterKeepCBTag = false;            // resolved: if true, keep full CB+TAG barcodes in per-sample MEX

    // Minimal memory mode (inline hash path only)
    string soloFlexMinimalMemoryStr = "no";      // raw CLI: yes|no (default: no)
    bool soloFlexMinimalMemory = false;          // resolved: true if enabled

    //processing
    uint32 redistrReadsNfiles; //numer of files to resditribute reads into
    
    //constants
    uint32 umiMaskLow, umiMaskHigh; //low/high half bit-mask or UMIs

    //readStats
    struct {
        string type; //input parameter
        bool yes=false;
    } readStats;

    void initialize(Parameters *pPin);
    void umiSwapHalves(uint32 &umi);
    void complexWLstrings();
    void cellFiltering();

    void init_CBmatchWL();
};
#endif
