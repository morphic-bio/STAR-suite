#ifndef H_FlexProbeHalfIndex
#define H_FlexProbeHalfIndex

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

/**
 * Compact lookup for the experimental Flex probe seed-and-extend policy.
 *
 * The two 25-base probe halves are indexed independently. Each index contains
 * the exact half and all of its Hamming-1 neighbours. A unique half-probe seed
 * is extended against its complete 50-base probe, accepting at most ten total
 * mismatches. Probe identity is retained so halves from different probes are
 * rejected even when those probes target the same gene.
 */
class FlexProbeHalfIndex {
public:
    struct Probe {
        std::string sequence;
        uint16_t geneIdx15;
        Probe(const std::string& sequenceIn = std::string(),
              uint16_t geneIdx15In = 0)
            : sequence(sequenceIn), geneIdx15(geneIdx15In) {}
    };

    enum Status : uint8_t {
        None = 0,
        Unique = 1,
        Ambiguous = 2,
        ScoreFail = 3,
        SplitProbe = 4
    };

    struct Result {
        Status status;
        uint16_t geneIdx15;
        uint8_t hammingDistance;
        Result(Status statusIn = None, uint16_t geneIdx15In = 0,
               uint8_t hammingDistanceIn = 0xFF)
            : status(statusIn), geneIdx15(geneIdx15In),
              hammingDistance(hammingDistanceIn) {}
    };

    bool build(const std::vector<Probe>& probes, std::string* errorOut = nullptr);
    Result classify(const char* readSeq, std::size_t readLen) const;
    Result classifyCachePacked(uint64_t seqLo, uint64_t seqHi,
                               uint64_t nMask = 0) const;

    bool ready() const { return ready_; }
    std::size_t probeCount() const { return probeCount_; }
    std::size_t leftKeyCount() const { return expanded_[0].size(); }
    std::size_t rightKeyCount() const { return expanded_[1].size(); }

private:
    struct PackedProbe {
        uint64_t half[2];
        uint16_t geneIdx15;
    };

    struct Entry {
        uint64_t key;
        uint32_t probeIdxPlus1;
        Entry(uint64_t keyIn = 0, uint32_t probeIdxPlus1In = 0)
            : key(keyIn), probeIdxPlus1(probeIdxPlus1In) {}
    };

    struct HalfResult {
        Status status;
        uint32_t probeIdxPlus1;
        uint8_t hammingDistance;
        HalfResult(Status statusIn = None, uint32_t probeIdxPlus1In = 0,
                   uint8_t hammingDistanceIn = 0xFF)
            : status(statusIn), probeIdxPlus1(probeIdxPlus1In),
              hammingDistance(hammingDistanceIn) {}
    };

    static bool packHalf(const char* sequence, uint64_t& key);
    static void addExpanded(std::vector<Entry>& entries, uint64_t exactKey,
                            uint32_t probeIdxPlus1);
    void reduce(std::vector<Entry>& entries);
    static uint8_t hammingDistance(uint64_t lhs, uint64_t rhs,
                                   uint32_t nMask = 0);
    HalfResult lookup(const std::vector<Entry>& entries, uint64_t key,
                      uint32_t nMask, unsigned side) const;
    static HalfResult mergeHalf(HalfResult lhs, HalfResult rhs);
    HalfResult classifyHalf(const char* sequence, unsigned side,
                            uint64_t& key, uint32_t& nMask) const;
    HalfResult classifyPackedHalf(uint64_t key, uint32_t nMask,
                                  unsigned side) const;
    Result classifyPackedHalves(uint64_t leftKey, uint64_t rightKey,
                                uint32_t leftNMask,
                                uint32_t rightNMask) const;

    bool ready_ = false;
    std::size_t probeCount_ = 0;
    std::vector<PackedProbe> probes_;
    std::vector<Entry> exact_[2];
    std::vector<Entry> expanded_[2];
};

#endif
