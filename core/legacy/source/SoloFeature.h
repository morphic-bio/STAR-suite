#ifndef H_SoloFeature
#define H_SoloFeature

#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "IncludeDefine.h"
#include "ReadAlignChunk.h"
#include "Transcriptome.h"

#include "SoloCommon.h"
#include "SoloRead.h"
#include "ReadAlignChunk.h"
#include "SoloReadFeature.h"

#include "SoloFilteredCells.h"

#include "PackedReadInfo.h"
#include "UMICorrector.h"

#include "SampleMatrixData.h"
#include "MexWriter.h"

// Forward declaration
class ProbeListIndex;

class SoloFeature {
private:
    Parameters &P;
    ReadAlignChunk **RAchunk;    
    Transcriptome &Trans;
    
    SoloFeature **soloFeatAll;
    
    static const uint32 umiArrayStride=3;
    enum {rguG, rguU, rguR};
    uint32 rguStride;
    
public:
    // Expose stride accessor for sinks
    uint32 getRGUStride() const { return rguStride; }
    void setRGUStride(uint32 v) { rguStride = v; }
    // Expose Transcriptome for gene ID mapping (used by gene resolver)
    const Transcriptome& getTranscriptome() const { return Trans; }
    // Expose ProbeListIndex if available (stored in ReadAlignChunk)
    const ProbeListIndex* getProbeListIndex() const;
    // Global gene→probe index cache helpers
    static void initGeneProbeIdx(const Transcriptome& trans, const ProbeListIndex* probeIdx);
    static uint16_t getProbeIdxForGene(uint32_t geneIdx);
    // Global chromosome→genomic flag (true if chr name starts with "chr")
    static void initChrGenomicFlags(ReadAlignChunk **RAchunk);
    static bool isChrGenomic(uint32_t chrIdx);
    // Debug helper: size of gene→probe cache
    static size_t geneProbeCacheSize();
    ParametersSolo &pSolo;

    SoloReadFeature *readFeatSum, **readFeatAll;
    SoloReadBarcode *readBarSum;

    const int32 featureType;   
    
    uint64 nReadsMapped, nReadsInput; //total number of mapped reads
    uint32 nCB;
    uint32 featuresNumber; //number of features (i.e. genes, SJs, etc)

    uint32 *rGeneUMI = nullptr;//mapped reads sorted by CB
    uint32 *rCBn = nullptr;//number of reads for detected CBs in the whitelist
    uint32 **rCBp = nullptr;//array of pointers to each CB sub-array

    vector<uint32> indCB;//index of detected CBs in the whitelist
    vector<uint32> indCBwl; //reverse of indCB: index of WL CBs in detected CB list
    vector<uint32> nUMIperCB, nUMIperCBsorted;//number of UMIs per CB, and the same sorted (descendant)
    vector<uint32> nGenePerCB, nGenePerCBmulti;//number of genes (with >0 UMIs) per CB
    vector<uint32> nReadPerCB;//number of reads per CB. With multimappers: all aligns per CB
    vector<uint32> nReadPerCBunique, nReadPerCBtotal; //number of unique and multiple reads per CB
    uint32 nReadPerCBmax;

    vector<double> nUMIperCBmulti;

    vector<uint32> countCellGeneUMI;//sparsified matrix for the counts, each entry is: geneID count1 count2 ... countNcounts
    vector<uint32> countCellGeneUMIindex;//index of CBs in the count matrix
    uint32 countMatStride; //number of counts per entry in the count matrix
    
    struct {
        vector<double> m;
        vector<uint32> i;
        uint32 s;
    } countMatMult;
    
    vector<unordered_map<uint32, unordered_set<uint64>>> cbFeatureUMImap; //for SmartSeq counting
       
    string outputPrefix, outputPrefixFiltered;
    
    SoloFilteredCells filteredCells;
    
    array<vector<uint64>,2> sjAll;
    
    PackedReadInfo packedReadInfo;   // packed representation (8 bytes per read)
    SoloReadFlagClass readFlagCounts;

    
    vector<uint32> redistrFilesCBindex, redistrFilesCBfirst; //redistr file for each CB, CB boundaries in redistributed files
    vector<uint64> redistrFilesNreads; //number of reads in each file
    vector <fstream*> redistrFilesStreams;

    SoloFeature(Parameters &Pin, ReadAlignChunk **RAchunk, Transcriptome &inTrans, int32 feTy, SoloReadBarcode *readBarSumIn, SoloFeature **soloFeatAll);
    void clearLarge(); //clear large vectors
    void processRecords();
    void sumThreads();
    void countSmartSeq();
    void countCBgeneUMI();
    void countVelocyto();
    void quantTranscript();
    void prepareReadInfoOnly(); //minimal processing to populate readInfo without counting (for skipProcessing mode)
    
    void collapseUMI(uint32 iCB, uint32 *umiArray);
    void collapseUMI_CR(uint32 iCB, uint32 *umiArray);
    void collapseUMIall(bool minimalMode=false);
    void collapseUMIall_fromHash(); // Direct hash consumption (no materialization)
    void collapseUMIperCB(uint32 iCB, vector<uint32> &umiArray, vector<uint32> &gID,  vector<uint32> &gReadS, bool minimalMode);
    void materializeRGUFromHash(); // Materialize rGeneUMI/rCBp/rCBn from inlineHash_ (DEPRECATED)
    
    struct InlineMatrixBundle {
        SampleMatrixData matrixData;
        std::vector<MexWriter::Triplet> triplets;
    };

    // Build matrix bundle from inline hash dedup counts
    InlineMatrixBundle buildInlineMatrixFromHash(
        const std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>>& cbTagGeneCounts);

    // Write MEX directly from inline-hash dedup data (no Solo, no replayer)
    void writeMexFromInlineHashDedup(
        const std::string& outputPrefix,
        const InlineMatrixBundle& matrixBundle);

    // Run flexfilter inline after composite MEX is written
    void runFlexFilterInline(
        const InlineMatrixBundle& inlineMatrix,  // In-memory composite MEX data
        const std::string& outputPrefix          // FlexFilter output prefix (e.g., <outPrefix>/flexfilter/)
    );

    uint32 umiArrayCorrect_CR         (const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr);
    uint32 umiArrayCorrect_Directional(const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr, const int32 dirCountAdd);
    uint32 umiArrayCorrect_Graph      (const uint32 nU0, uintUMI *umiArr, const bool readInfoRec, const bool nUMIyes, unordered_map <uintUMI,uintUMI> &umiCorr);

    // Helpers to abstract legacy vs packed readInfo storage
    void resetPackedStorage(uint32_t nReads);
    void recordReadInfo(uint32_t readId, uint32_t cbIdx, uint32_t umiPacked, uint8_t status);
    uint32_t getPackedCB(uint32_t readId) const;
    uint32_t getPackedUMI(uint32_t readId) const;
    uint8_t  getPackedStatus(uint32_t readId) const;

    void outputResults(bool cellFilterYes, string outputPrefixMat);
    // Legacy: extracts iRead from trailing 8 bytes at bam0+size0 (includes Y-bit encoding)
    void addBAMtags(char *&bam0, uint32 &size0, char *bam1);
    // New: explicit readId for samtools sorter (uint32 readId, no shifting needed)
    void addBAMtags(char *&bam0, uint32 &size0, char *bam1, uint32_t readId);
    void statsOutput();
    void redistributeReadsByCB();
    
    void cellFiltering();
    void emptyDrops_CR();
    void loadRawMatrix();
    
    void writeTagTableIfRequested(bool) {} // stub - tag table not used in flex path
    void finalizeTagTableFromReadInfo();
    void writeReadIdTagTable();  // Export readId/CB/UB/status TSV table (env var gated)
    void initPackedReadInfo(uint32_t nReads) { packedReadInfo.init(nReads, pSolo.cbWLstr.size(), pSolo.umiL); }
    
    // UMI correction (clique deduplication) data structures
    struct AssignmentInfo {
        uint8_t sampleIdx;
        std::string sampleTag;
        std::string status;  // CONFIDENT, AMBIGUOUS, etc.
    };
    
    std::unordered_set<std::string> cellsAllowSet;  // CB16 (wildcard) or CB24 (exact match)
    std::unordered_map<std::string, AssignmentInfo> assignmentsMap;  // key: CB16, value: assignment info
    
    // UMI histogram collection: key is packed uint64 (cbIdx << 24 | sampleIdx << 16 | geneIdx)
    struct UMIHistogram {
        std::unordered_map<std::string, uint32_t> urCounts;  // UR string -> count
        uint32_t totalReads;
    };
    std::unordered_map<uint64_t, UMIHistogram> umiGroupHistograms;
    
    // UMI corrections: groupKey -> (UR -> corrected UB)
    std::unordered_map<uint64_t, std::unordered_map<std::string, std::string>> umiCorrections;
    
    // Load assignment files
    bool loadCellsAllowList(const std::string& path);
    bool loadSampleAssignments(const std::string& path);
    
    // UMI correction: collect UR histogram for a read
    void collectURHistogram(uint32_t readId, uint32_t cbIdx, uint32_t geneIdx);
    
    // UMI correction: run clique correction for all groups and apply corrections
    void runCliqueCorrection();
    
    // Ambiguous CB resolution: stub - not wired in flex path
    void resolveAmbiguousCBs() {}
    
    // Build UMI histograms from inline hash (for clique correction)
    void buildHistogramsFromHash();
    
    // Apply clique corrections back to inline hash (re-key entries with corrected UB)
    void applyCliqueCorrectionsToHash();
    
    // Metrics output
    void writeUMICleanupMetrics();
    void writeQCSummary();
    
    // UMI correction metrics
    uint64_t umisBeforeTotal, umisAfterTotal;
    uint64_t readsBeforeTotal, readsAfterTotal;
    uint32_t mergesTotal, componentsTotal, componentsCappedTotal, componentsBelowThresholdTotal;
    uint32_t maxComponentSeen;
    uint32_t componentSizeHist[5];  // [0]=size1, [1]=size2, [2]=size3, [3]=size4, [4]=size>4
    uint64_t readsURGrouped, readsURMissing;
    
    // Cell filtering metrics
    uint64_t cellsInput;           // Total unique CB16s seen
    uint64_t cellsKept;            // CB16s that passed allow-list filter
    uint64_t ambiguousCellsDropped; // Ambiguous cells dropped (per policy)
    
#ifdef DEBUG_CB_UB_PARITY
    // CB/UB Parity validation structures
    struct LegacyCBUB {
        uint32_t iCB;
        uint32_t iRead;
        std::string cb;
        std::string ub;
        uint32_t geneIdx;
        bool matched = false;
    };
    std::unordered_map<uint32_t, LegacyCBUB> legacyCBUBByRead; // Snapshot before clique correction, keyed by readId
    std::ofstream mismatchLog;           // File handle for Aligned.out.cbub_mismatches.txt
    uint64_t parityReadsLegacy, parityReadsInline, parityMatches, parityMismatches; // Counters
    bool parityEnabled;                  // Runtime flag from STAR_DEBUG_CB_UB_PARITY env var
    std::unordered_map<std::string,uint64_t> debugStatusCounters;
    // Stage instrumentation
    uint64_t dbgBufferedRecords = 0;     // records materialized into rGeneUMI
    uint64_t dbgBufferedCBs = 0;         // CBs with at least one buffered record
    uint64_t dbgWriteTotal = 0;          // calls to recordReadInfo
    uint64_t dbgWriteSiteUnique = 0;     // unique-gene path before clique correction
    uint64_t dbgWriteSiteCliqueGood = 0; // post-clique path kept UMIs
    uint64_t dbgWriteSiteCliqueDrop = 0; // post-clique path dropped/invalid UMIs
    uint64_t dbgWriteSiteMulti = 0;      // multigene tail path
    uint64_t dbgRAWPass = 0;             // read-after-write passes
    uint64_t dbgRAWFail = 0;             // read-after-write mismatches
    uint64_t dbgWriteStatus0 = 0;
    uint64_t dbgWriteStatus1 = 0;
    uint64_t dbgWriteStatus2 = 0;
    uint64_t dbgReadStatus0  = 0;        // inline status==0 observed at readback
    void resetDebugStatusCounters();
    void noteDebugStatus(uint8_t status, const char* reason);
    void logDebugStatusCounters();
    void logDebugStageCounters();
    
    // Parity validation methods
    void snapshotLegacyCBUB(uint32_t iCB, uint32_t *rGU, vector<uint32> &gReadS, uint32_t nGenes, vector<uint32> &gID, uint32_t rguStride);
    void compareCBUBParity(uint32_t iCB, uint32_t *rGU, vector<uint32> &gReadS, uint32_t nGenes, vector<uint32> &gID, uint32_t rguStride);
    void writeParityCounters();
#endif
};

#endif
