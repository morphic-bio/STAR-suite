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
    std::uint64_t peakBytes = 0;
};

bool parseProducts(const std::string &value, std::uint8_t &mask, std::string &error);
bool parseScales(const std::string &value, std::uint8_t &mask, std::string &error);
bool parseOverflowPolicy(const std::string &value, OverflowPolicy &policy,
                         std::string &error);
bool estimateMemory(const Capacity &capacity, MemoryModel &model, std::string &error);
bool memoryFits(const MemoryModel &model, std::uint64_t availableBytes,
                double fraction, std::uint64_t &budgetBytes, std::string &error);
// The current compact spill boundary covers read/candidate accumulation. This
// guard prevents Spill from promising a budget that a later in-memory phase
// cannot satisfy.
bool spillBudgetFits(const MemoryModel &model, std::uint64_t budgetBytes,
                     std::string &error);

// Returns the smaller of host MemAvailable and an active cgroup limit after
// subtracting current cgroup usage. Zero means that no reliable value could be
// obtained.
std::uint64_t availableMemoryBytes();

struct PipelineConfig {
    std::string barcodeContractDirectory;
    std::string bc1OligosPath;
    std::string bc2OligosPath;
    std::string outputDirectory;
    std::string temporaryDirectory;
    std::string starSuiteVersion;
    std::string sourceRevision;
    std::uint64_t expectedReads = 0;
    std::uint64_t expectedCandidates = 0;
    std::uint32_t threads = 1;
    std::uint8_t products = ProductAll;
    std::uint8_t scales = ScaleAll;
    double memoryFraction = 0.80;
    OverflowPolicy overflowPolicy = OverflowPolicy::Fail;
    // Zero selects a budget-derived automatic limit. A nonzero value is an
    // advanced deterministic test/tuning override and is valid only with
    // OverflowPolicy::Spill.
    std::uint64_t spillHighWaterCandidatesPerThread = 0;
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
