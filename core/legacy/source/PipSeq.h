#ifndef STAR_PIPSEQ_H
#define STAR_PIPSEQ_H

#include <array>
#include <cstdint>
#include <string>
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
    // As read, before correction: the four tier words at this phase (28
    // bases) and the three BI bases; "-" when the read does not reach them.
    static std::string rawCell(const std::string& read, int phase);
    static std::string rawBin(const std::string& read, int phase);
private:
    // entry: list index chosen for this tier value, or -1. best: smallest
    // Hamming distance to any entry. ambiguous: no exact entry, more than one
    // entry at distance one, and MultiMatch::Reject.
    struct Match { int entry = -1; unsigned best = 255; bool ambiguous = false; };
    std::array<std::vector<std::string>, 4> tiers_;
    std::array<std::vector<Match>, 4> lookup_;
};

// One distinct (cell, gene, BI, IMI) observation and its read count.
// key = (global aligned 5' base << 7) | (strand << 6) | BI.
struct Fragment {
    uint64_t cell, key;
    uint32_t gene, reads;
};
struct GroupSummary {
    uint64_t cell = 0;
    uint32_t gene = 0;
    uint64_t reads = 0, fragments = 0;
    unsigned bins = 0; // occupied BI bins, 0-64
};

// Flat store of fragments: 24 bytes per distinct fragment, no per-group
// containers. Each mapping thread fills its own store; the stores are merged
// after mapping, so duplicate reads seen by different threads collapse.
class Fragments {
public:
    void add(uint64_t cell, uint32_t gene, unsigned bin, uint64_t coordinate, bool reverse);
    void finalize();                // sort; collapse identical fragments, summing reads
    void absorb(Fragments& other);  // merge another store into this one and empty it
    void clear();
    size_t size() const { return records_.size(); } // distinct fragments after finalize()
    // Calls f(const GroupSummary&) for every (cell, gene), ordered by cell then
    // gene. Requires finalize() or absorb() since the last add().
    template <class F> void forEachGroup(F&& f) const {
        for (size_t i = 0; i < records_.size();) {
            GroupSummary g;
            g.cell = records_[i].cell; g.gene = records_[i].gene;
            uint64_t bins = 0;
            for (; i < records_.size() && records_[i].cell == g.cell && records_[i].gene == g.gene; ++i) {
                bins |= uint64_t(1) << (records_[i].key & 63);
                g.reads += records_[i].reads;
                ++g.fragments;
            }
            g.bins = __builtin_popcountll(bins);
            f(g);
        }
    }
private:
    void collapse();
    std::vector<Fragment> records_;
    size_t sorted_ = 0; // records_[0, sorted_) is sorted and collapsed
};

struct Calibration {
    double raw = 0;
    double effective = 1;
    uint64_t groups = 0;
    std::string source;
};
double expectedMolecules(unsigned occupiedBins);
Calibration calibrate(const Fragments& fragments, double suppliedIPM);
uint64_t transcripts(unsigned occupiedBins, uint64_t fragments, double ipm);

// cDNA read trimming as in the published pipeline description for STAR
// alignment: one leading base; the template-switch and SMART PCR adapter
// sequences at the 5' end; poly-A and poly-G at the 3' end; a minimum overlap
// of 12 bases; no indels. Sequences are numeric (A, C, G, T = 0..3); any other
// value never matches and separates the adapters of a list.
extern const char* const kitAdapters5p; // "AGAGTGAATGGG,TCAACGCAGAGT"
extern const char* const kitAdapters3p; // "polyA,polyG"
constexpr uint32_t minTrimmedLength = 20;
// Bases to remove from the 5' end. Every placement of every adapter is scored
// by matches minus mismatches; a placement needs at least minOverlap aligned
// bases and at most floor(maxMismatchFraction * aligned bases) mismatches.
// An adapter longer than minOverlap may hang off the read start; a complete
// adapter may lie anywhere in the read. The best score wins; on ties, the
// placement that removes most.
uint32_t clip5pAdapter(const char* read, uint32_t length, const char* adapters, uint32_t adaptersLength,
                       double maxMismatchFraction, uint32_t minOverlap = 12);
// Bases to remove from the 3' end: from the start of a window of `window`
// bases with at most floor(maxMismatchFraction * window) bases other than
// `base` to the read end. The window with the fewest wins, the leftmost on ties.
uint32_t clip3pHomopolymer(const char* read, uint32_t length, char base, double maxMismatchFraction, uint32_t window = 12);

} // namespace pipseq
#endif
