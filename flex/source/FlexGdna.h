#ifndef H_FlexGdna
#define H_FlexGdna

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

enum FlexGdnaRegion : uint8_t {
    FlexGdnaUnknown = 0,
    FlexGdnaSpliced = 1,
    FlexGdnaUnspliced = 2,
    FlexGdnaConflicting = 3
};

static const uint32_t kFlexGdnaCountMask = 0x3FFFFFFFu;
static const uint32_t kFlexGdnaRegionShift = 30u;

inline FlexGdnaRegion flexGdnaMergeRegion(FlexGdnaRegion lhs, FlexGdnaRegion rhs)
{
    if (lhs == FlexGdnaConflicting || rhs == FlexGdnaConflicting)
        return FlexGdnaConflicting;
    if (lhs == FlexGdnaUnknown)
        return rhs;
    if (rhs == FlexGdnaUnknown)
        return lhs;
    return lhs == rhs ? lhs : FlexGdnaConflicting;
}

inline uint32_t flexGdnaValueCount(uint32_t value)
{
    return value & kFlexGdnaCountMask;
}

inline FlexGdnaRegion flexGdnaValueRegion(uint32_t value)
{
    return static_cast<FlexGdnaRegion>((value >> kFlexGdnaRegionShift) & 0x3u);
}

inline uint32_t flexGdnaPackValue(uint32_t count, FlexGdnaRegion region)
{
    if (count > kFlexGdnaCountMask)
        count = kFlexGdnaCountMask;
    return count | (static_cast<uint32_t>(region) << kFlexGdnaRegionShift);
}

inline uint32_t flexGdnaMergeValue(uint32_t lhs, uint32_t rhs)
{
    const uint64_t sum = static_cast<uint64_t>(flexGdnaValueCount(lhs))
                       + static_cast<uint64_t>(flexGdnaValueCount(rhs));
    const uint32_t count = sum > kFlexGdnaCountMask
        ? kFlexGdnaCountMask
        : static_cast<uint32_t>(sum);
    return flexGdnaPackValue(
        count,
        flexGdnaMergeRegion(flexGdnaValueRegion(lhs), flexGdnaValueRegion(rhs)));
}

inline uint32_t flexGdnaIncrementValue(uint32_t value, uint32_t count, FlexGdnaRegion region)
{
    return flexGdnaMergeValue(value, flexGdnaPackValue(count, region));
}

struct FlexGdnaGeneProbeCounts {
    uint32_t spliced = 0;
    uint32_t unspliced = 0;
};

struct FlexGdnaGeneMoleculeCounts {
    uint64_t spliced = 0;
    uint64_t unspliced = 0;
};

struct FlexGdnaEstimate {
    bool valid = false;
    std::string status;
    uint32_t controlGenes = 0;
    uint64_t totalFilteredMolecules = 0;
    uint64_t classifiedMolecules = 0;
    uint64_t unknownMolecules = 0;
    uint64_t conflictingMolecules = 0;
    uint64_t unassignedMolecules = 0;
    double estimatedGdnaPerProbe = 0.0;
    double estimatedGdnaFraction = 0.0;
    double threshold = 0.0;
    double modelConstant = 0.0;
    double modelSlope = 0.0;
    double modelCriticalPoint = 0.0;
    double modelRss = 0.0;
};

class FlexGdnaProbeMetadata {
public:
    static FlexGdnaProbeMetadata& instance();

    void reset();
    bool load(const std::string& probeCsvPath,
              const std::string& probeListPath,
              std::string* errorOut = nullptr);

    bool ready() const { return ready_; }
    const std::string& probeCsvPath() const { return probeCsvPath_; }
    uint32_t totalProbes() const { return totalProbes_; }
    uint32_t totalUnsplicedProbes() const { return totalUnsplicedProbes_; }
    uint32_t controlGeneCount() const { return controlGeneCount_; }
    const std::vector<FlexGdnaGeneProbeCounts>& geneProbeCounts() const {
        return geneProbeCounts_;
    }

    FlexGdnaRegion regionForProbeId(const std::string& probeId) const;

    static std::string discoverProbeCsv(const std::string& configuredPath,
                                        const std::string& probeListPath,
                                        const std::string& genomeDir);

private:
    FlexGdnaProbeMetadata() = default;

    bool ready_ = false;
    std::string probeCsvPath_;
    std::unordered_map<std::string, FlexGdnaRegion> probeRegionById_;
    std::vector<FlexGdnaGeneProbeCounts> geneProbeCounts_;
    uint32_t totalProbes_ = 0;
    uint32_t totalUnsplicedProbes_ = 0;
    uint32_t controlGeneCount_ = 0;
};

FlexGdnaEstimate flexGdnaEstimate(
    const FlexGdnaProbeMetadata& metadata,
    const std::vector<FlexGdnaGeneMoleculeCounts>& molecules,
    uint64_t classifiedMolecules,
    uint64_t unknownMolecules,
    uint64_t conflictingMolecules,
    uint64_t unassignedMolecules = 0);

#endif
