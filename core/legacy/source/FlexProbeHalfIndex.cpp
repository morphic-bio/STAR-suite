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

uint8_t hammingByte(uint8_t difference)
{
    uint8_t count = 0;
    for (unsigned shift = 0; shift < 8; shift += 2)
        count += ((difference >> shift) & 3U) != 0;
    return count;
}

const uint8_t* hammingLut()
{
    static uint8_t lut[256];
    static const bool initialized = []() {
        for (unsigned value = 0; value < 256; ++value)
            lut[value] = hammingByte(static_cast<uint8_t>(value));
        return true;
    }();
    (void)initialized;
    return lut;
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
                                     uint32_t probeIdxPlus1)
{
    entries.push_back(Entry{exactKey, probeIdxPlus1});
    for (unsigned position = 0; position < 25; ++position) {
        const unsigned shift = 2U * (24U - position);
        const uint64_t current = (exactKey >> shift) & UINT64_C(3);
        const uint64_t cleared = exactKey & ~(UINT64_C(3) << shift);
        for (uint64_t alternative = 0; alternative < 4; ++alternative) {
            if (alternative == current)
                continue;
            entries.push_back(Entry{
                cleared | (alternative << shift), probeIdxPlus1});
        }
    }
}

void FlexProbeHalfIndex::reduce(std::vector<Entry>& entries)
{
    std::sort(entries.begin(), entries.end(),
              [](const Entry& lhs, const Entry& rhs) {
                  if (lhs.key != rhs.key)
                      return lhs.key < rhs.key;
                  return lhs.probeIdxPlus1 < rhs.probeIdxPlus1;
              });
    std::size_t write = 0;
    for (std::size_t begin = 0; begin < entries.size();) {
        std::size_t end = begin + 1;
        const uint32_t firstProbeIdxPlus1 = entries[begin].probeIdxPlus1;
        uint32_t probeIdxPlus1 = firstProbeIdxPlus1;
        const PackedProbe& first = probes_[firstProbeIdxPlus1 - 1U];
        while (end < entries.size() && entries[end].key == entries[begin].key) {
            if (entries[end].probeIdxPlus1 != firstProbeIdxPlus1) {
                const PackedProbe& other = probes_[entries[end].probeIdxPlus1 - 1U];
                if (first.geneIdx15 != other.geneIdx15 ||
                    first.half[0] != other.half[0] ||
                    first.half[1] != other.half[1]) {
                    probeIdxPlus1 = 0;
                }
            }
            ++end;
        }
        entries[write++] = Entry{entries[begin].key, probeIdxPlus1};
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
    probes_.clear();
    exact_[0].clear();
    exact_[1].clear();
    expanded_[0].clear();
    expanded_[1].clear();

    probes_.reserve(probes.size());
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
        PackedProbe packed{{0, 0}, probe.geneIdx15};
        for (unsigned side = 0; side < 2; ++side) {
            if (!packHalf(probe.sequence.data() + side * 25U,
                          packed.half[side])) {
                if (errorOut != nullptr)
                    *errorOut = "active Flex probe contains a non-ACGT base";
                return false;
            }
        }
        probes_.push_back(packed);
        const uint32_t probeIdxPlus1 = static_cast<uint32_t>(probes_.size());
        for (unsigned side = 0; side < 2; ++side) {
            exact_[side].push_back(Entry{packed.half[side], probeIdxPlus1});
            addExpanded(expanded_[side], packed.half[side], probeIdxPlus1);
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

uint8_t FlexProbeHalfIndex::hammingDistance(uint64_t lhs, uint64_t rhs,
                                            uint32_t nMask)
{
    uint64_t difference = lhs ^ rhs;
    while (nMask != 0) {
        const unsigned position = static_cast<unsigned>(__builtin_ctz(nMask));
        difference |= UINT64_C(1) << (2U * (24U - position));
        nMask &= nMask - 1U;
    }
    const uint8_t* lut = hammingLut();
    return static_cast<uint8_t>(
        lut[(difference >> 0) & 0xFFU] +
        lut[(difference >> 8) & 0xFFU] +
        lut[(difference >> 16) & 0xFFU] +
        lut[(difference >> 24) & 0xFFU] +
        lut[(difference >> 32) & 0xFFU] +
        lut[(difference >> 40) & 0xFFU] +
        lut[(difference >> 48) & 0xFFU]);
}

FlexProbeHalfIndex::HalfResult FlexProbeHalfIndex::lookup(
    const std::vector<Entry>& entries, uint64_t key, uint32_t nMask,
    unsigned side) const
{
    const auto it = std::lower_bound(
        entries.begin(), entries.end(), key,
        [](const Entry& entry, uint64_t value) { return entry.key < value; });
    if (it == entries.end() || it->key != key)
        return HalfResult{};
    if (it->probeIdxPlus1 == 0)
        return HalfResult{Ambiguous, 0};
    const PackedProbe& probe = probes_[it->probeIdxPlus1 - 1U];
    return HalfResult{Unique, it->probeIdxPlus1,
                      hammingDistance(key, probe.half[side], nMask)};
}

FlexProbeHalfIndex::HalfResult FlexProbeHalfIndex::mergeHalf(
    HalfResult lhs, HalfResult rhs)
{
    if (lhs.status == Ambiguous || rhs.status == Ambiguous)
        return HalfResult{Ambiguous, 0};
    if (lhs.status == None)
        return rhs;
    if (rhs.status == None)
        return lhs;
    if (lhs.probeIdxPlus1 != rhs.probeIdxPlus1)
        return HalfResult{Ambiguous, 0};
    return lhs;
}

FlexProbeHalfIndex::HalfResult FlexProbeHalfIndex::classifyHalf(
    const char* sequence, unsigned side, uint64_t& key,
    uint32_t& nMask) const
{
    unsigned invalidCount = 0;
    unsigned invalidPosition = 0;
    key = 0;
    nMask = 0;
    for (unsigned position = 0; position < 25; ++position) {
        const uint8_t code = baseCode(sequence[position]);
        key <<= 2;
        if (code > 3) {
            ++invalidCount;
            invalidPosition = position;
            nMask |= UINT32_C(1) << position;
        } else {
            key |= code;
        }
    }
    if (invalidCount == 0)
        return lookup(expanded_[side], key, 0, side);
    if (invalidCount > 1)
        return HalfResult{};

    // A single N/invalid base consumes the complete Hamming-1 allowance. Try
    // its four possible bases only against exact probe halves; do not allow an
    // additional mismatch elsewhere in the half.
    HalfResult result;
    const unsigned shift = 2U * (24U - invalidPosition);
    const uint64_t cleared = key & ~(UINT64_C(3) << shift);
    for (uint64_t alternative = 0; alternative < 4; ++alternative) {
        result = mergeHalf(result,
                           lookup(exact_[side],
                                  cleared | (alternative << shift),
                                  nMask, side));
    }
    return result;
}

FlexProbeHalfIndex::Result FlexProbeHalfIndex::classifyPackedHalves(
    uint64_t leftKey, uint64_t rightKey, uint32_t leftNMask,
    uint32_t rightNMask) const
{
    const HalfResult left = classifyPackedHalf(leftKey, leftNMask, 0);
    const HalfResult right = classifyPackedHalf(rightKey, rightNMask, 1);
    if (left.status == Ambiguous || right.status == Ambiguous)
        return Result{Ambiguous, 0};
    if (left.status == None && right.status == None)
        return Result{};
    if (left.status == Unique && right.status == Unique) {
        if (left.probeIdxPlus1 != right.probeIdxPlus1)
            return Result{SplitProbe, 0};
        const PackedProbe& probe = probes_[left.probeIdxPlus1 - 1U];
        return Result{Unique, probe.geneIdx15,
                      static_cast<uint8_t>(left.hammingDistance +
                                           right.hammingDistance)};
    }

    const HalfResult& anchor = left.status == Unique ? left : right;
    const unsigned otherSide = left.status == Unique ? 1U : 0U;
    const uint64_t otherKey = otherSide == 0 ? leftKey : rightKey;
    const uint32_t otherNMask = otherSide == 0 ? leftNMask : rightNMask;
    const PackedProbe& probe = probes_[anchor.probeIdxPlus1 - 1U];
    const uint8_t otherDistance = hammingDistance(
        otherKey, probe.half[otherSide], otherNMask);
    const uint8_t totalDistance = static_cast<uint8_t>(
        anchor.hammingDistance + otherDistance);
    if (totalDistance > 10)
        return Result{ScoreFail, 0, totalDistance};
    return Result{Unique, probe.geneIdx15, totalDistance};
}

FlexProbeHalfIndex::Result FlexProbeHalfIndex::classify(
    const char* readSeq, std::size_t readLen) const
{
    if (!ready_ || readSeq == nullptr || readLen < 50)
        return Result{};
    uint64_t leftKey = 0;
    uint64_t rightKey = 0;
    uint32_t leftNMask = 0;
    uint32_t rightNMask = 0;
    const HalfResult left = classifyHalf(
        readSeq, 0, leftKey, leftNMask);
    const HalfResult right = classifyHalf(
        readSeq + 25, 1, rightKey, rightNMask);
    if (left.status == Ambiguous || right.status == Ambiguous)
        return Result{Ambiguous, 0};
    if (left.status == None && right.status == None)
        return Result{};
    if (left.status == Unique && right.status == Unique) {
        if (left.probeIdxPlus1 != right.probeIdxPlus1)
            return Result{SplitProbe, 0};
        const PackedProbe& probe = probes_[left.probeIdxPlus1 - 1U];
        return Result{Unique, probe.geneIdx15,
                      static_cast<uint8_t>(left.hammingDistance +
                                           right.hammingDistance)};
    }
    const HalfResult& anchor = left.status == Unique ? left : right;
    const unsigned otherSide = left.status == Unique ? 1U : 0U;
    const uint64_t otherKey = otherSide == 0 ? leftKey : rightKey;
    const uint32_t otherNMask = otherSide == 0 ? leftNMask : rightNMask;
    const PackedProbe& probe = probes_[anchor.probeIdxPlus1 - 1U];
    const uint8_t totalDistance = static_cast<uint8_t>(
        anchor.hammingDistance + hammingDistance(
            otherKey, probe.half[otherSide], otherNMask));
    if (totalDistance > 10)
        return Result{ScoreFail, 0, totalDistance};
    return Result{Unique, probe.geneIdx15, totalDistance};
}

FlexProbeHalfIndex::HalfResult FlexProbeHalfIndex::classifyPackedHalf(
    uint64_t key, uint32_t nMask, unsigned side) const
{
    if (nMask == 0)
        return lookup(expanded_[side], key, 0, side);
    if ((nMask & (nMask - 1U)) != 0)
        return HalfResult{};

    // A single N consumes the complete Hamming-1 allowance, as in the ASCII
    // path. nMask is in read-coordinate order, so bit zero is the first base
    // of this half while the packed key stores that base in bits 49..48.
    const unsigned position = static_cast<unsigned>(__builtin_ctz(nMask));
    const unsigned shift = 2U * (24U - position);
    const uint64_t cleared = key & ~(UINT64_C(3) << shift);
    HalfResult result;
    for (uint64_t alternative = 0; alternative < 4; ++alternative) {
        result = mergeHalf(result,
                           lookup(exact_[side],
                                  cleared | (alternative << shift),
                                  nMask, side));
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
    return classifyPackedHalves(leftKey, rightKey, leftNMask, rightNMask);
}
