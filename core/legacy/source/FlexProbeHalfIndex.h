#ifndef H_FlexProbeHalfIndex
#define H_FlexProbeHalfIndex

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

/**
 * Compact lookup for the experimental Flex residual-alignment policy.
 *
 * The two 25-base probe halves are indexed independently. Each index contains
 * the exact half and all of its Hamming-1 neighbours. Lookup combines the two
 * corresponding halves at gene level: no hit, one gene, or multiple genes.
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
        Ambiguous = 2
    };

    struct Result {
        Status status;
        uint16_t geneIdx15;
        Result(Status statusIn = None, uint16_t geneIdx15In = 0)
            : status(statusIn), geneIdx15(geneIdx15In) {}
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
    struct Entry {
        uint64_t key;
        uint16_t geneIdx15;
        Entry(uint64_t keyIn = 0, uint16_t geneIdx15In = 0)
            : key(keyIn), geneIdx15(geneIdx15In) {}
    };

    static bool packHalf(const char* sequence, uint64_t& key);
    static void addExpanded(std::vector<Entry>& entries, uint64_t exactKey,
                            uint16_t geneIdx15);
    static void reduce(std::vector<Entry>& entries);
    static Result lookup(const std::vector<Entry>& entries, uint64_t key);
    static Result merge(Result lhs, Result rhs);
    Result classifyHalf(const char* sequence, unsigned side) const;
    Result classifyPackedHalf(uint64_t key, uint32_t nMask,
                              unsigned side) const;

    bool ready_ = false;
    std::size_t probeCount_ = 0;
    std::vector<Entry> exact_[2];
    std::vector<Entry> expanded_[2];
};

#endif
