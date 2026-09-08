#include "FlexProbeHalfIndex.h"

#include <algorithm>
#include <cstring>

namespace {

uint8_t baseCode(char base)
{
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

} // namespace

bool FlexProbeHalfIndex::packHalf(const char* sequence, uint64_t& key)
{
    key = 0;
    for (unsigned ii = 0; ii < 25; ++ii) {
        const uint8_t code = baseCode(sequence[ii]);
        if (code > 3)
            return false;
        key = (key << 2) | code;
    }
    return true;
}

void FlexProbeHalfIndex::addExpanded(std::vector<Entry>& entries,
                                     uint64_t exactKey,
                                     uint16_t geneIdx15)
{
    entries.push_back(Entry{exactKey, geneIdx15});
    for (unsigned position = 0; position < 25; ++position) {
        const unsigned shift = 2U * (24U - position);
        const uint64_t current = (exactKey >> shift) & UINT64_C(3);
        const uint64_t cleared = exactKey & ~(UINT64_C(3) << shift);
        for (uint64_t alternative = 0; alternative < 4; ++alternative) {
            if (alternative == current)
                continue;
            entries.push_back(Entry{
                cleared | (alternative << shift), geneIdx15});
        }
    }
}

void FlexProbeHalfIndex::reduce(std::vector<Entry>& entries)
{
    std::sort(entries.begin(), entries.end(),
              [](const Entry& lhs, const Entry& rhs) {
                  if (lhs.key != rhs.key)
                      return lhs.key < rhs.key;
                  return lhs.geneIdx15 < rhs.geneIdx15;
              });
    std::size_t write = 0;
    for (std::size_t begin = 0; begin < entries.size();) {
        std::size_t end = begin + 1;
        uint16_t gene = entries[begin].geneIdx15;
        while (end < entries.size() && entries[end].key == entries[begin].key) {
            if (entries[end].geneIdx15 != gene)
                gene = 0;
            ++end;
        }
        entries[write++] = Entry{entries[begin].key, gene};
        begin = end;
    }
    entries.resize(write);
    entries.shrink_to_fit();
}

bool FlexProbeHalfIndex::build(const std::vector<Probe>& probes,
                               std::string* errorOut)
{
    ready_ = false;
    probeCount_ = 0;
    exact_[0].clear();
    exact_[1].clear();
    expanded_[0].clear();
    expanded_[1].clear();

    exact_[0].reserve(probes.size());
    exact_[1].reserve(probes.size());
    expanded_[0].reserve(probes.size() * 76U);
    expanded_[1].reserve(probes.size() * 76U);

    for (const Probe& probe : probes) {
        if (probe.sequence.size() != 50 || probe.geneIdx15 == 0) {
            if (errorOut != nullptr)
                *errorOut = "active Flex probes must have 50 A/C/G/T bases and a nonzero gene index";
            return false;
        }
        for (unsigned side = 0; side < 2; ++side) {
            uint64_t exactKey = 0;
            if (!packHalf(probe.sequence.data() + side * 25U, exactKey)) {
                if (errorOut != nullptr)
                    *errorOut = "active Flex probe contains a non-ACGT base";
                return false;
            }
            exact_[side].push_back(Entry{exactKey, probe.geneIdx15});
            addExpanded(expanded_[side], exactKey, probe.geneIdx15);
        }
        ++probeCount_;
    }
    if (probeCount_ == 0) {
        if (errorOut != nullptr)
            *errorOut = "active Flex probe set is empty";
        return false;
    }

    reduce(exact_[0]);
    reduce(exact_[1]);
    reduce(expanded_[0]);
    reduce(expanded_[1]);
    ready_ = true;
    return true;
}

FlexProbeHalfIndex::Result FlexProbeHalfIndex::lookup(
    const std::vector<Entry>& entries, uint64_t key)
{
    const auto it = std::lower_bound(
        entries.begin(), entries.end(), key,
        [](const Entry& entry, uint64_t value) { return entry.key < value; });
    if (it == entries.end() || it->key != key)
        return Result{};
    if (it->geneIdx15 == 0)
        return Result{Ambiguous, 0};
    return Result{Unique, it->geneIdx15};
}

FlexProbeHalfIndex::Result FlexProbeHalfIndex::merge(Result lhs, Result rhs)
{
    if (lhs.status == Ambiguous || rhs.status == Ambiguous)
        return Result{Ambiguous, 0};
    if (lhs.status == None)
        return rhs;
    if (rhs.status == None)
        return lhs;
    if (lhs.geneIdx15 != rhs.geneIdx15)
        return Result{Ambiguous, 0};
    return lhs;
}

FlexProbeHalfIndex::Result FlexProbeHalfIndex::classifyHalf(
    const char* sequence, unsigned side) const
{
    unsigned invalidCount = 0;
    unsigned invalidPosition = 0;
    uint64_t key = 0;
    for (unsigned position = 0; position < 25; ++position) {
        const uint8_t code = baseCode(sequence[position]);
        key <<= 2;
        if (code > 3) {
            ++invalidCount;
            invalidPosition = position;
        } else {
            key |= code;
        }
    }
    if (invalidCount == 0)
        return lookup(expanded_[side], key);
    if (invalidCount > 1)
        return Result{};

    // A single N/invalid base consumes the complete Hamming-1 allowance. Try
    // its four possible bases only against exact probe halves; do not allow an
    // additional mismatch elsewhere in the half.
    Result result;
    const unsigned shift = 2U * (24U - invalidPosition);
    const uint64_t cleared = key & ~(UINT64_C(3) << shift);
    for (uint64_t alternative = 0; alternative < 4; ++alternative) {
        result = merge(result,
                       lookup(exact_[side], cleared | (alternative << shift)));
    }
    return result;
}

FlexProbeHalfIndex::Result FlexProbeHalfIndex::classify(
    const char* readSeq, std::size_t readLen) const
{
    if (!ready_ || readSeq == nullptr || readLen < 50)
        return Result{};
    return merge(classifyHalf(readSeq, 0), classifyHalf(readSeq + 25, 1));
}

FlexProbeHalfIndex::Result FlexProbeHalfIndex::classifyPackedHalf(
    uint64_t key, uint32_t nMask, unsigned side) const
{
    if (nMask == 0)
        return lookup(expanded_[side], key);
    if ((nMask & (nMask - 1U)) != 0)
        return Result{};

    // A single N consumes the complete Hamming-1 allowance, as in the ASCII
    // path. nMask is in read-coordinate order, so bit zero is the first base
    // of this half while the packed key stores that base in bits 49..48.
    const unsigned position = static_cast<unsigned>(__builtin_ctz(nMask));
    const unsigned shift = 2U * (24U - position);
    const uint64_t cleared = key & ~(UINT64_C(3) << shift);
    Result result;
    for (uint64_t alternative = 0; alternative < 4; ++alternative) {
        result = merge(result,
                       lookup(exact_[side], cleared | (alternative << shift)));
    }
    return result;
}

FlexProbeHalfIndex::Result FlexProbeHalfIndex::classifyCachePacked(
    uint64_t seqLo, uint64_t seqHi, uint64_t nMask) const
{
    if (!ready_)
        return Result{};

    const uint64_t halfMask = (UINT64_C(1) << 50U) - 1U;
    const uint64_t leftKey = ((seqHi << 14U) | (seqLo >> 50U)) & halfMask;
    const uint64_t rightKey = seqLo & halfMask;
    const uint32_t leftNMask = static_cast<uint32_t>(nMask & ((UINT64_C(1) << 25U) - 1U));
    const uint32_t rightNMask = static_cast<uint32_t>((nMask >> 25U) & ((UINT64_C(1) << 25U) - 1U));
    return merge(classifyPackedHalf(leftKey, leftNMask, 0),
                 classifyPackedHalf(rightKey, rightNMask, 1));
}
