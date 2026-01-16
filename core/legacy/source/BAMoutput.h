#ifndef CODE_BAMoutput
#define CODE_BAMoutput

#include "IncludeDefine.h"
#include SAMTOOLS_BGZF_H
#include "Parameters.h"
#include <vector>
#include <cstdint>
#include <fstream>

// Forward declarations
class SampleDetector;

// Packed metadata for unsorted BAM output
struct PendingSoloMeta {
    uint32_t iReadIdx;        // Read index used for grouping/debug
    uint32_t cbIdxPlus1;      // 1-based whitelist index (0 => unmatched)
    uint32_t umi24;           // Packed 12bp UMI (LSB-first, 24 bits used)
    uint16_t geneIdx15;       // Harmonized gene index (0 if unset/multi)
    uint16_t sampleIdx;       // 1-based sample index, 0=no sample
    uint8_t dropFlags;        // MAPQ/NH/NM/gene/sample flags
    uint8_t newGroup;         // 1 => first record for a QNAME group
    uint8_t _pad;             // padding for alignment
    std::string cbSeq;        // CB sequence
    std::vector<uint16_t> zgGeneIdx15; // full ZG gene list for tie-break
    std::vector<uint32_t> neighborWL; // Ambiguous CB neighbor whitelist indices (1-based)
    bool urValid;             // UR validity flag (false if contains N)
    bool hasY;                // true if this read touches Y-chromosome
   
    PendingSoloMeta()
        : iReadIdx(0), cbIdxPlus1(0), umi24(0), geneIdx15(0),
          sampleIdx(0), dropFlags(0), newGroup(0), _pad(0), cbSeq(), zgGeneIdx15(), neighborWL(), urValid(true), hasY(false) {}
};

class BAMoutput {
public:
    //sorted output
    BAMoutput (int iChunk, string tmpDir, Parameters &Pin);
    void coordOneAlign (char *bamIn, uint bamSize, uint iRead, bool hasY = false);
    void coordBins ();
    void coordFlush ();
    //unsorted output
    BAMoutput (BGZF *bgzfBAMin, Parameters &Pin);
    ~BAMoutput();
    
    void unsortedOneAlign (char *bamIn, uint bamSize, uint bamSize2, uint64_t iReadAll, uint8_t sampleByte, 
                           uint32_t cbIdxPlus1 = 0, uint32_t umi24 = 0, const std::string &cbSeq = std::string(),
                           bool hasY = false);
    void unsortedFlush ();
    void coordUnmappedPrepareBySJout();

    uint32 nBins; //number of bins to split genome into
    uint* binTotalN; //total number of aligns in each bin
    uint* binTotalBytes;//total size of aligns in each bin
private:
    uint64 bamArraySize; //this size will be allocated
    char* bamArray; //large array to store the bam alignments, pre-sorted
    uint64 binSize, binSize1;//storage size of each bin
    uint64 binGlen;//bin genomic length
    char **binStart; //pointers to starts of the bins
    uint64 *binBytes, binBytes1;//number of bytes currently written to each bin
    ofstream **binStream;//output streams for each bin
    BGZF *bgzfBAM;
    BGZF *bgzfBAM_Y;         // _Y.bam handle (nullptr if !emitNoYBAM)
    BGZF *bgzfBAM_noY;       // _noY.bam handle (nullptr if !emitNoYBAM)
    bool suppressPrimary_;   // true if --keepBAM is not set
    Parameters &P;
    string bamDir;
    
    // Staging queues for flush-synchronized ledger (direct unsorted mode)
    std::vector<uint32_t> pendingReadIds_;
    std::vector<uint32_t> pendingAux_;
    
    // Staging for unsorted mode
    std::vector<PendingSoloMeta> pendingSoloMeta_;
    std::string lastQname_;  // Track QNAME changes for newGroup bit
    
    SampleDetector *sampleDet_;
    bool sampleDetReady_;
    
    // Ambiguous CB spillover file (optional)
    std::ofstream ambiguousCbSpilloverFile_;
    bool ambiguousCbSpilloverEnabled_;
    
    // Helper to flush buffer and append to ledger atomically
    void flushPendingToLedgerAndDisk();
};

#endif
