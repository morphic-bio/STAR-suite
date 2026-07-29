#ifndef STAR_SPATIAL_GEX_H
#define STAR_SPATIAL_GEX_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "SpatialR1Decoder.h"

namespace spatial_gex {

enum Product : std::uint8_t {
    ProductStrict = 1u << 0,
    ProductSoftExpected = 1u << 1,
    ProductHard = 1u << 2,
    ProductGatedHard = 1u << 3,
    ProductAll = ProductStrict | ProductSoftExpected | ProductHard | ProductGatedHard
};

enum Scale : std::uint8_t {
    Scale2um = 1u << 0,
    Scale8um = 1u << 1,
    Scale16um = 1u << 2,
    ScaleAll = Scale2um | Scale8um | Scale16um
};

enum class OverflowPolicy : std::uint8_t {
    Fail = 0,
    Spill = 1
};

// Production hot-path records. Keep these fixed-width and free of owning STL
// members. A source ordinal is sufficient for deterministic scheduling; read
// names are not retained.
struct ReadEvidence {
    std::uint32_t geneIndex;
    std::uint32_t rawUmi;
    std::uint32_t candidateBegin;
    std::uint32_t sourceOrdinal;
    std::uint16_t candidateCount;
    std::uint16_t flags;
};

struct CandidateEvidence {
    std::uint32_t coordinateIndex;
    std::uint32_t auditBits;
    double logSequenceLikelihood;
};

struct PolicySupport {
    std::uint32_t geneIndex;
    std::uint32_t coordinateIndex;
    std::uint32_t rawUmi;
    std::uint32_t strictCount;
    std::uint32_t hardCount;
    std::uint32_t gatedHardCount;
    double softPosteriorSupport;
    double softLogAbsent;
};

struct CorrectedSupport {
    std::uint32_t coordinateIndex;
    std::uint32_t geneIndex;
    std::uint32_t correctedUmi;
    std::uint32_t originalUmi;
    std::uint32_t correctedCount;
    std::uint32_t originalAtCorrectedCount;
    double softExpected;
};

struct FinalMolecule {
    std::uint32_t geneIndex;
    std::uint32_t coordinateIndex;
    std::uint32_t correctedUmi;
    std::uint32_t policy;
    double weight;
};

static_assert(sizeof(ReadEvidence) == 20, "spatial GEX read evidence size changed");
static_assert(sizeof(CandidateEvidence) == 16, "spatial GEX candidate evidence size changed");
static_assert(sizeof(PolicySupport) == 40, "spatial GEX policy support size changed");
static_assert(sizeof(CorrectedSupport) == 32, "spatial GEX corrected support size changed");
static_assert(sizeof(FinalMolecule) == 24, "spatial GEX final molecule size changed");

struct Capacity {
    std::uint64_t reads = 0;
    std::uint64_t candidates = 0;
    std::uint32_t threads = 1;
};

struct MemoryModel {
    std::uint64_t accumulationFixedBytes = 0;
    std::uint64_t accumulationBytes = 0;
    std::uint64_t cliqueBytes = 0;
    std::uint64_t correctionBytes = 0;
    std::uint64_t reconciliationBytes = 0;
    std::uint64_t materializationBytes = 0;
    // End-to-end Spill uses a fixed downstream workspace and coordinate
    // partitions instead of the slide-sized phase estimates above.
    std::uint64_t downstreamSpoolBytes = 0;
    // Conservative simultaneous temporary-run and uncompressed-MEX capacity.
    std::uint64_t downstreamSpoolDiskBytes = 0;
    std::uint64_t peakBytes = 0;
};

bool parseProducts(const std::string &value, std::uint8_t &mask, std::string &error);
bool parseScales(const std::string &value, std::uint8_t &mask, std::string &error);
bool parseOverflowPolicy(const std::string &value, OverflowPolicy &policy,
                         std::string &error);
bool estimateMemory(const Capacity &capacity, MemoryModel &model, std::string &error);
bool memoryFits(const MemoryModel &model, std::uint64_t availableBytes,
                double fraction, std::uint64_t &budgetBytes, std::string &error);
// Spill uses the compact accumulation runs plus a bounded downstream
// coordinate-partition workspace. This guard checks the fixed end-to-end
// resident requirement rather than the all-memory downstream phase peaks.
bool spillBudgetFits(const MemoryModel &model, std::uint64_t budgetBytes,
                     std::string &error);

// Returns the smaller of host MemAvailable and an active cgroup limit after
// subtracting current cgroup usage. Zero means that no reliable value could be
// obtained.
std::uint64_t availableMemoryBytes();

// Describes the feature decision paired with the current raw-R1 spatial
// decode. Flex cache hits and alignment fallbacks are kept distinct for
// deterministic accounting; GEX uses the final post-rescue annotation.
enum class FeatureEvidenceClass : std::uint8_t {
    Gex = 0,
    FlexH0 = 1,
    FlexH1 = 2,
    FlexHashDeny = 3,
    FlexAlignment = 4
};

struct PipelineConfig {
    std::string barcodeContractDirectory;
    std::string bc1OligosPath;
    std::string bc2OligosPath;
    std::string outputDirectory;
    std::string temporaryDirectory;
    std::string starSuiteVersion;
    std::string sourceRevision;
    std::string featureAxisPath;
    std::string featureAxisSha256;
    std::string featureCachePath;
    std::string featureCacheSha256;
    std::uint32_t featureCount = 0;
    std::uint64_t expectedReads = 0;
    std::uint64_t expectedCandidates = 0;
    std::uint32_t threads = 1;
    // Hard is the production policy. ProductAll remains available explicitly
    // for diagnostic policy comparisons.
    std::uint8_t products = ProductHard;
    std::uint8_t scales = ScaleAll;
    double memoryFraction = 0.80;
    OverflowPolicy overflowPolicy = OverflowPolicy::Fail;
    // Zero selects a budget-derived automatic limit. A nonzero value is an
    // advanced deterministic test/tuning override and is valid only with
    // OverflowPolicy::Spill.
    std::uint64_t spillHighWaterCandidatesPerThread = 0;
    // Internal downstream partitioning knobs. STAR production defaults are
    // conservative; tests may use smaller values to force boundary cases.
    std::uint32_t downstreamSpoolShards = 256;
    std::uint64_t downstreamSpoolBufferBytes = 64 * 1024;
    // Native Flex requires exactly one terminal feature decision for every
    // decoded R1. This catches early returns and thread-local state drift.
    bool requirePairedCompletion = false;
    bool flexFeatureMode = false;
    // Compatibility default target-scans barcode windows with N. False uses
    // the N-aware hash candidate path and keeps the target scan out of the run.
    bool barcodeNdpFallback = true;
};

struct PipelineSummary {
    std::uint64_t readsDecoded = 0;
    std::uint64_t readsWithCandidates = 0;
    std::uint64_t uniqueGeneReads = 0;
    std::uint64_t joinedReads = 0;
    std::uint64_t candidateRows = 0;
    std::uint64_t exactH0Reads = 0;
    std::uint64_t barcodeReadsWithN = 0;
    std::uint64_t barcodeNBases = 0;
    std::uint64_t barcodeDpRecoveredReads = 0;
    std::uint64_t barcodeDpAmbiguousReads = 0;
    std::uint64_t barcodeDpUnassignedReads = 0;
    std::uint64_t barcodeNHashRecoveredReads = 0;
    std::uint64_t barcodeNHashAmbiguousReads = 0;
    std::uint64_t barcodeNHashUnassignedReads = 0;
    std::uint64_t barcodeUnsupportedReads = 0;
    std::uint64_t umiReadsWithN = 0;
    std::uint64_t umiReadsWithInvalidBase = 0;
    std::uint64_t readCliques = 0;
    std::uint64_t strictMolecules = 0;
    double softExpectedMass = 0.0;
    std::uint64_t hardMolecules = 0;
    std::uint64_t gatedHardMolecules = 0;
    std::uint64_t spillRuns = 0;
    std::uint64_t spillBytes = 0;
    std::uint64_t spillHighWaterCandidatesPerThread = 0;
    double spillMergeSeconds = 0.0;
    std::uint64_t peakResidentReads = 0;
    std::uint64_t peakResidentCandidates = 0;
    std::uint64_t downstreamContributionRecords = 0;
    std::uint64_t downstreamContributionRuns = 0;
    std::uint64_t downstreamContributionBytes = 0;
    std::uint64_t hardAssignmentChunks = 0;
    std::uint64_t hardAssignmentOverflowCliques = 0;
    std::uint64_t hardAssignmentMaxMemberCount = 0;
    std::uint64_t downstreamLargestShardRecords = 0;
    std::uint64_t downstreamMatrixRuns = 0;
    std::uint64_t downstreamMatrixBytes = 0;
    double downstreamResolveSeconds = 0.0;
    double downstreamMaterializeSeconds = 0.0;
    std::uint64_t featureAssignedReads = 0;
    std::uint64_t featureUnassignedReads = 0;
    std::uint64_t featureAssignedReadsWithCandidates = 0;
    std::uint64_t featureAssignedReadsWithoutCandidates = 0;
    std::uint64_t featureUnassignedReadsWithCandidates = 0;
    std::uint64_t featureUnassignedReadsWithoutCandidates = 0;
    std::uint64_t flexHashH0Reads = 0;
    std::uint64_t flexHashH1Reads = 0;
    std::uint64_t flexHashDenyReads = 0;
    std::uint64_t flexAlignmentMissReads = 0;
    std::uint64_t flexAlignmentResolvedReads = 0;
    std::uint64_t flexAlignmentUnresolvedReads = 0;
};

// Run-owned, default-off spatial state. One instance is shared by all mapping
// workers; hot-path writes go only to the worker's indexed ThreadState.
class Pipeline {
  public:
    // Opaque implementation type is public only so translation-unit-local
    // numeric kernels can operate on it without exposing its fields here.
    struct Impl;
    static std::unique_ptr<Pipeline> create(const PipelineConfig &config,
                                            std::string &error);
    ~Pipeline();
    Pipeline(const Pipeline &) = delete;
    Pipeline &operator=(const Pipeline &) = delete;

    bool decode(std::uint32_t threadIndex, const char *sequence,
                std::size_t sequenceLength, const char *quality,
                std::size_t qualityLength, spatial_r1_decoder::Result &result,
                std::string &error);
    bool append(std::uint32_t threadIndex, std::uint32_t geneIndex,
                std::uint64_t sourceOrdinal,
                const spatial_r1_decoder::Result &decoded,
                std::string &error);
    void noteUniqueGene(std::uint32_t threadIndex);

    // STAR mapping-thread interface. It retains the current decoded R1 inside
    // the gated pipeline, avoiding any spatial payload or layout change in
    // ReadAlign when integrated mode is disabled.
    bool decodeCurrentThread(const char *sequence, std::size_t sequenceLength,
                             const char *quality, std::size_t qualityLength,
                             std::uint64_t sourceOrdinal,
                             std::string &error);
    bool completeCurrentThread(FeatureEvidenceClass source, bool assigned,
                               std::uint32_t geneIndex,
                               std::uint64_t sourceOrdinal,
                               std::string &error);
    bool appendCurrentThread(std::uint32_t geneIndex,
                             std::uint64_t sourceOrdinal,
                             std::string &error);

    bool finalize(const std::vector<std::string> &geneIds,
                  std::string &error);
    const PipelineSummary &summary() const;
    const MemoryModel &memoryModel() const;
    std::uint64_t memoryBudgetBytes() const;

  private:
    explicit Pipeline(std::unique_ptr<Impl> impl);
    bool currentThreadIndex(std::uint32_t &threadIndex, std::string &error);
    std::unique_ptr<Impl> impl_;
};

} // namespace spatial_gex

#endif
