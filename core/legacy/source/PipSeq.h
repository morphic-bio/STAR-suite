#ifndef STAR_PIPSEQ_H
#define STAR_PIPSEQ_H

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

// Independent implementation from published DRAGEN PIPseq documentation.
// No vendor source or packaged references are required.
namespace pipseq {

enum class DecodeStatus { Accepted, InvalidBarcode, NoMatch, PhaseTie, Ambiguous };
// A tier with no exact entry and several entries one substitution away:
// reject the read (result does not depend on list order), or take the first
// listed entry (the published wording; depends on list order).
enum class MultiMatch { Reject, FirstListed };
struct Decoded {
    DecodeStatus status = DecodeStatus::InvalidBarcode;
    uint64_t cell = 0;
    int phase = -1;
    unsigned mismatches = 0;
    int bin = -1; // Invalid BI does not invalidate the cell assignment.
};

class Decoder {
public:
    explicit Decoder(const std::vector<std::string>& paths, MultiMatch multiMatch = MultiMatch::Reject);
    Decoded decode(const std::string& read) const;
    static std::string cellSequence(uint64_t cell);
private:
    // entry: list index chosen for this tier value, or -1. best: smallest
    // Hamming distance to any entry. ambiguous: no exact entry, more than one
    // entry at distance one, and MultiMatch::Reject.
    struct Match { int entry = -1; unsigned best = 255; bool ambiguous = false; };
    std::array<std::vector<std::string>, 4> tiers_;
    std::array<std::vector<Match>, 4> lookup_;
};

struct Group {
    uint64_t bins = 0;
    uint64_t reads = 0;
    std::unordered_set<uint64_t> fragments; // (global aligned 5' base, strand, BI)
    void add(unsigned bin, uint64_t coordinate, bool reverse);
    void merge(Group& other);
    unsigned occupiedBins() const;
};
using Groups = std::map<std::pair<uint64_t, uint32_t>, Group>; // cell, gene

struct Calibration {
    double raw = 0;
    double effective = 1;
    uint64_t groups = 0;
    std::string source;
};
double expectedMolecules(unsigned occupiedBins);
Calibration calibrate(const Groups& groups, double suppliedIPM);
uint64_t transcripts(const Group& group, double ipm);

} // namespace pipseq
#endif
