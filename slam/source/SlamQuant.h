#ifndef SLAM_QUANT_H
#define SLAM_QUANT_H

#include "SlamSolver.h"
#include "SlamVarianceAnalysis.h"
#include "SlamReadBuffer.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <vector>
#include <memory>

class Genome;
class Transcriptome;
class SlamCompat;

struct SlamSnpMaskVcfOptions {
    std::string sample;     // Sample name (optional; required in GT mode for multi-sample VCF)
    std::string mode;       // "gt" (use genotype) or "any" (ignore genotype)
    std::string filter;     // "pass" (PASS or .) or "all"
    std::string bedOut;     // Optional BED output path
    std::string summaryOut; // Optional summary TSV output path
};

struct SlamSnpMaskVcfStats {
    uint64_t recordsTotal = 0;
    uint64_t recordsFiltered = 0;
    uint64_t recordsNoSnpAlt = 0;
    uint64_t recordsNonSnpAlt = 0;
    uint64_t recordsMissingGt = 0;
    uint64_t recordsNoAltGt = 0;
    uint64_t recordsUnknownContig = 0;
    uint64_t sitesAdded = 0;
    uint64_t sitesDuplicate = 0;
};

class SlamSnpMask {
public:
    bool loadBed(const std::string& path, const Genome& genome, std::string* err);
    bool loadBedWithChrMap(const std::string& path,
                           const std::vector<std::string>& chrNames,
                           const std::vector<uint64_t>& chrStart,
                           std::string* err);
    bool loadVcf(const std::string& path,
                 const Genome& genome,
                 const SlamSnpMaskVcfOptions& opts,
                 SlamSnpMaskVcfStats* stats,
                 std::string* err);
    bool loadVcfWithChrMap(const std::string& path,
                           const std::vector<std::string>& chrNames,
                           const std::vector<uint64_t>& chrStart,
                           const SlamSnpMaskVcfOptions& opts,
                           SlamSnpMaskVcfStats* stats,
                           std::string* err);
    bool contains(uint64_t pos) const { return positions_.count(pos) > 0; }
    size_t size() const { return positions_.size(); }

private:
    std::unordered_set<uint64_t> positions_;
};

struct SlamGeneStats {
    MismatchHistogram histogram;
    double readCount = 0.0;
    double conversions = 0.0;
    double coverage = 0.0;
};

struct SlamDiagnostics {
    uint64_t readsDroppedSnpMask = 0;
    uint64_t readsDroppedStrandness = 0;
    uint64_t readsZeroGenes = 0;
    uint64_t readsProcessed = 0;
    uint64_t readsNAlignWithGeneZero = 0;  // reads where nAlignWithGene == 0
    uint64_t readsSumWeightLessThanOne = 0;  // reads where sumWeight < 1.0
    std::map<size_t, uint64_t> nTrDistribution;  // nTr -> count
    std::map<size_t, uint64_t> geneSetSizeDistribution;  // gene set size -> count
    std::map<size_t, uint64_t> nAlignWithGeneDistribution;  // nAlignWithGene -> count
    std::map<double, uint64_t> sumWeightDistribution;  // sumWeight bucket -> count (bucketed)
    
    // Compat mode counters (incremented by callers, not SlamCompat)
    // NOTE: Alignment-level counters (one per alignment, not per read)
    uint64_t compatAlignsReclassifiedIntronic = 0;  // Alignments reclassified as intronic
    uint64_t compatAlignsLenientAccepted = 0;       // Alignments accepted via lenient overlap
    uint64_t compatAlignsOverlapWeightApplied = 0;  // Alignments where weight was adjusted
    
    // Position-level counters (one per genomic position)
    uint64_t compatPositionsSkippedOverlap = 0;     // Positions skipped due to PE overlap
    uint64_t compatPositionsSkippedTrim = 0;        // Positions skipped due to trim guards
};

struct SlamSnpBufferStats {
    uint64_t bufferedReads = 0;
    uint64_t bufferedMismatches = 0;
    uint64_t bufferedMismatchesKept = 0;
    uint64_t bufferEntries = 0;
    uint64_t bufferBytes = 0;
    uint64_t maskEntries = 0;
    uint64_t blacklistEntries = 0;
    double avgMismatches = 0.0;
    double avgMismatchesKept = 0.0;
    // SNP threshold estimation stats
    double mismatchFracUsed = 0.0;      // Threshold actually used
    double mismatchFracAuto = 0.0;      // Auto-estimated threshold (0 if not computed)
    std::string mismatchFracMode;       // "auto" or "explicit"
    uint32_t kneeBin = 0;               // Bin index of detected knee (0 if fallback)
    uint64_t eligibleSites = 0;         // Sites with coverage >= threshold
};

struct SlamTransitionStats {
    double coverage[4] = {0.0, 0.0, 0.0, 0.0};
    double mismatches[4][4] = {{0.0}};
};

struct SlamPositionStats {
    double coverage[2][2][4] = {};
    double mismatches[2][2][4][4] = {};
};

enum class SlamDebugDropReason : uint8_t {
    None = 0,
    NoGenes = 1,
    SnpMask = 2,
    Strandness = 3,
    ZeroWeight = 4
};

struct SlamDebugGeneStats {
    double exonicWeight = 0.0;
    double intronicWeight = 0.0;
    double senseWeight = 0.0;
    double antisenseWeight = 0.0;
    double conversions = 0.0;
    double coverage = 0.0;
    uint64_t readsAssigned = 0;
    uint64_t readsAssignedIntronic = 0;
    uint64_t dropsSnpMask = 0;
    uint64_t dropsStrandness = 0;
};

struct SlamDebugReadRecord {
    std::string readName;
    std::string readLoc;
    uint32_t geneId = 0;
    bool intronic = false;
    bool oppositeStrand = false;
    double weight = 0.0;
    uint16_t nT = 0;
    uint16_t k = 0;
    uint32_t readLength = 0;
    SlamDebugDropReason status = SlamDebugDropReason::None;
    bool snpBuffered = false;
    std::string convReadPos;
    std::string convGenPos;
};

enum class SlamMismatchCategory : uint8_t {
    Exonic = 0,
    ExonicSense = 1,
    Intronic = 2,
    IntronicSense = 3,
    Count
};

constexpr size_t kSlamMismatchCategoryCount = static_cast<size_t>(SlamMismatchCategory::Count);
const char* slamMismatchCategoryName(SlamMismatchCategory cat);

class SlamQuant {
public:
    explicit SlamQuant(uint32_t nGenes, bool snpDetect = false, double snpMismatchFrac = -1.0, bool snpObsAnyMismatch = false);
    SlamQuant(uint32_t nGenes, std::vector<uint8_t> allowedGenes, bool snpDetect = false, double snpMismatchFrac = -1.0, bool snpObsAnyMismatch = false);

    void addRead(uint32_t geneId, uint16_t nT, uint8_t k, double weight);
    void addTransitionBase(SlamMismatchCategory category, uint32_t readPos, bool secondMate,
                           bool overlap, bool opposite, int genomicBase, int readBase, double weight);
    bool snpDetectEnabled() const { return snpDetectEnabled_; }
    // Record an observation for SNP masking. Depending on configuration, "alt" can mean:
    // - conversion mismatch (T->C/A->G) or
    // - any mismatch
    void recordSnpObservation(uint64_t pos, bool anyMismatch, bool convMismatch);
    void bufferSnpRead(uint32_t geneId, uint16_t nT,
                       const std::vector<uint32_t>& mismatchPositions, double weight);
    void finalizeSnpMask(SlamSnpBufferStats* outStats = nullptr);
    
    // Get raw SNP mask data (for mask build pre-pass)
    // Returns copy of internal snpMask_ map before finalization
    std::unordered_map<uint64_t, uint32_t> getSnpMaskData() const;
    
    // Variance analysis for auto-trim
    bool recordVarianceRead(); // Returns false if max reads reached
    void recordVariancePosition(uint32_t readPos, uint8_t qual, bool isT, bool isTc);
    SlamVarianceTrimResult computeVarianceTrim(uint32_t readLength);
    bool varianceAnalysisEnabled() const { return varianceAnalyzer_ != nullptr; }
    const SlamVarianceAnalyzer* varianceAnalyzer() const { return varianceAnalyzer_.get(); }
    void enableVarianceAnalysis(uint32_t maxReads, uint32_t minReads,
                                uint32_t smoothWindow = 5, uint32_t minSegLen = 3, uint32_t maxTrim = 15);
    void resetVarianceAnalysis();
    uint32_t getVarianceMaxReads() const;
    uint32_t getVarianceMinReads() const;
    
    // Read buffering for auto-trim replay
    void enableReadBuffer(uint64_t maxReads);
    bool readBufferEnabled() const { return readBuffer_ != nullptr; }
    bool readBufferFull() const { return readBuffer_ && readBuffer_->isFull(); }
    uint64_t readBufferSize() const { return readBuffer_ ? readBuffer_->size() : 0; }
    uint64_t readBufferCapacity() const { return readBuffer_ ? readBuffer_->capacity() : 0; }
    bool bufferRead(SlamBufferedRead&& read);
    const SlamReadBuffer* readBuffer() const { return readBuffer_.get(); }
    void clearReadBuffer();
    uint64_t readBufferMemoryBytes() const { return readBuffer_ ? readBuffer_->memoryBytes() : 0; }
    
    // Replay buffered reads with trim applied
    // Returns number of reads replayed
    uint64_t replayBufferedReads(SlamCompat* compat, const SlamSnpMask* snpMask, int strandness);

    // Dump buffer for external re-quantification
    void enableDumpBuffer(uint64_t maxReads);
    bool dumpEnabled() const { return dumpBuffer_ != nullptr; }
    bool dumpBufferFull() const { return dumpBuffer_ && dumpBuffer_->isFull(); }
    uint64_t dumpBufferSize() const { return dumpBuffer_ ? dumpBuffer_->size() : 0; }
    bool bufferDumpRead(SlamBufferedRead&& read);
    const SlamReadBuffer* dumpBuffer() const { return dumpBuffer_.get(); }
    void merge(const SlamQuant& other);
    void write(const Transcriptome& tr, const std::string& outFile,
               double errorRate, double convRate) const;
    void writeGrandSlam(const Transcriptome& tr, const std::string& outFile,
                        const std::string& outFileNamePrefix,
                        double errorRate, double convRate) const;
    void writeDiagnostics(const std::string& diagFile) const;
    void writeTransitions(const std::string& outFile) const;
    void writeMismatches(const std::string& outFile, const std::string& condition) const;
    void writeMismatchDetails(const std::string& outFile) const;
    void writeTopMismatches(const Transcriptome& tr, const std::string& refFile,
                           const std::string& mismatchFile, size_t topN) const;
    
    // Get position transition data for QC (ExonicSense category only)
    // Returns map: position -> {tc_cov, tc_mm, ta_cov, ta_mm}
    std::unordered_map<uint32_t, std::tuple<double, double, double, double>> getPositionTransitionData() const;
    void initDebug(const Transcriptome& tr,
                   const std::unordered_set<std::string>& debugGenes,
                   const std::unordered_set<std::string>& debugReads,
                   size_t maxReads,
                   const std::string& outPrefix);
    bool debugEnabled() const { return debugEnabled_; }
    bool debugGenesEnabled() const { return debugEnabled_ && !debugGeneMask_.empty(); }
    bool debugGeneEnabled(uint32_t geneId) const;
    bool debugReadMatch(const char* readName) const;
    void debugCountDrop(uint32_t geneId, SlamDebugDropReason reason);
    void debugAddAssignment(uint32_t geneId, double weight, bool intronic,
                            bool oppositeStrand, uint16_t nT, uint8_t k);
    void debugLogRead(const SlamDebugReadRecord& record);
    void writeDebug(const Transcriptome& tr, double errorRate, double convRate) const;
    // Enable/collect SNP-site debug for investigating counting parity.
    // absPos is STAR's 0-based genome-wide coordinate (same coordinate space as snpMask_ keys).
    void enableSnpSiteDebug(uint64_t absPos, int window, const std::string& locString);
    bool snpSiteDebugEnabled() const { return debugSnpEnabled_; }
    void debugSnpSiteObserve(uint64_t absPos, bool anyMismatch, bool convMismatch,
                             double weight, bool primaryFlag, int mapq);

    const std::vector<SlamGeneStats>& genes() const { return geneStats_; }
    SlamDiagnostics& diagnostics() { return diag_; }
    const SlamDiagnostics& diagnostics() const { return diag_; }

private:
    std::vector<SlamGeneStats> geneStats_;
    SlamDiagnostics diag_;
    std::array<SlamTransitionStats, kSlamMismatchCategoryCount> transitions_;
    std::array<SlamTransitionStats, kSlamMismatchCategoryCount> transitionsFirst_;
    std::array<SlamTransitionStats, kSlamMismatchCategoryCount> transitionsSecond_;
    std::array<std::vector<SlamPositionStats>, kSlamMismatchCategoryCount> positionTransitions_;
    bool snpDetectEnabled_ = false;
    bool snpFinalized_ = false;
    bool snpObsAnyMismatch_ = false; // if true, count any mismatch as alt; otherwise count conversions only
    double snpMismatchFrac_ = -1.0;  // <=0 means auto-estimate
    std::unordered_map<uint64_t, uint32_t> snpMask_;
    std::vector<uint32_t> snpReadBuffer_;
    std::vector<double> snpReadWeights_;
    std::vector<uint8_t> allowedGenes_;
    bool debugEnabled_ = false;
    size_t debugMaxReads_ = 0;
    std::string debugOutPrefix_;
    std::vector<uint8_t> debugGeneMask_;
    std::vector<SlamDebugGeneStats> debugGeneStats_;
    std::unordered_set<std::string> debugReadSet_;
    std::vector<SlamDebugReadRecord> debugReadRecords_;

    // SNP-site debug (separate from gene/read debug; shares debugOutPrefix_)
    bool debugSnpEnabled_ = false;
    uint64_t debugSnpAbsPos_ = 0;
    int debugSnpWindow_ = 0;
    std::string debugSnpLoc_;
    // Per-offset arrays (index = offset + debugSnpWindow_)
    std::vector<uint64_t> debugSnpCov_;
    std::vector<uint64_t> debugSnpAnyMis_;
    std::vector<uint64_t> debugSnpConvMis_;
    std::vector<uint64_t> debugSnpCovPrimary_;
    std::vector<uint64_t> debugSnpAnyMisPrimary_;
    std::vector<uint64_t> debugSnpConvMisPrimary_;
    std::vector<uint64_t> debugSnpCovNh1_;
    std::vector<uint64_t> debugSnpAnyMisNh1_;
    std::vector<uint64_t> debugSnpConvMisNh1_;
    std::vector<uint64_t> debugSnpCovNhGt1_;
    std::vector<uint64_t> debugSnpAnyMisNhGt1_;
    std::vector<uint64_t> debugSnpConvMisNhGt1_;
    std::vector<uint64_t> debugSnpCovMapq20_;
    std::vector<uint64_t> debugSnpAnyMisMapq20_;
    std::vector<uint64_t> debugSnpConvMisMapq20_;
    
    // Variance analysis
    std::unique_ptr<SlamVarianceAnalyzer> varianceAnalyzer_;
    uint32_t varianceMaxReads_ = 0;
    uint32_t varianceMinReads_ = 0;
    
    // Read buffer for auto-trim replay
    std::unique_ptr<SlamReadBuffer> readBuffer_;
    std::unique_ptr<SlamReadBuffer> dumpBuffer_;
};

#endif
