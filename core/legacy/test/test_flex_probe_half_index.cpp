#include "FlexProbeHalfIndex.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

std::string mutate(std::string sequence, std::size_t position, char base)
{
    sequence[position] = base;
    return sequence;
}

void packCacheOrder(const std::string& sequence, uint64_t& lo, uint64_t& hi)
{
    lo = 0;
    hi = 0;
    for (char base : sequence) {
        const uint64_t code = base == 'A' ? 0 : base == 'C' ? 1
                              : base == 'G' ? 2 : 3;
        hi = (hi << 2) | (lo >> 62);
        lo = (lo << 2) | code;
    }
}

} // namespace

int main()
{
    const std::string probeA =
        "AACCGGTTAACCGGTTAACCGGTTA"
        "GGTTAACCGGTTAACCGGTTAACCG";
    const std::string probeB =
        "TTGGCCAATTGGCCAATTGGCCAAT"
        "CCAATTGGCCAATTGGCCAATTGGC";

    FlexProbeHalfIndex index;
    std::string error;
    std::vector<FlexProbeHalfIndex::Probe> probes = {
        {probeA, 11},
        {probeA, 11}, // repeated probes for one gene are not gene ambiguity
        {probeB, 22},
    };
    require(index.build(probes, &error), error.c_str());

    auto result = index.classify(probeA.data(), probeA.size());
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "an exact two-half match must resolve to gene A");

    const std::string oneMismatch = mutate(probeA, 3, 'T');
    result = index.classify(oneMismatch.data(), oneMismatch.size());
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "one mismatch in a half must retain its unique gene");

    uint64_t packedLo = 0;
    uint64_t packedHi = 0;
    packCacheOrder(oneMismatch, packedLo, packedHi);
    result = index.classifyCachePacked(packedLo, packedHi);
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "packed CBQ/cache-order lookup must match ASCII lookup");

    std::string oneAnchored = probeA;
    oneAnchored.replace(25, 25, "AAAAAAAAAAAAAAAAAAAAAAAAA");
    result = index.classify(oneAnchored.data(), oneAnchored.size());
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "one qualifying half is sufficient");

    std::string conflicting = probeA.substr(0, 25) + probeB.substr(25, 25);
    result = index.classify(conflicting.data(), conflicting.size());
    require(result.status == FlexProbeHalfIndex::Ambiguous,
            "qualifying halves from different genes must be ambiguous");

    std::string noAnchor(50, 'A');
    result = index.classify(noAnchor.data(), noAnchor.size());
    require(result.status == FlexProbeHalfIndex::None,
            "a read without a qualifying half must not be anchored");

    std::string oneN = probeA;
    oneN[7] = 'N';
    result = index.classify(oneN.data(), oneN.size());
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "one N otherwise exact must consume the one-mismatch allowance");

    oneN[8] = oneN[8] == 'A' ? 'C' : 'A';
    oneN.replace(25, 25, "AAAAAAAAAAAAAAAAAAAAAAAAA");
    result = index.classify(oneN.data(), oneN.size());
    require(result.status == FlexProbeHalfIndex::None,
            "an N plus another mismatch in the same half must exceed Hamming one");

    std::string packedN = probeA;
    packedN[7] = 'A';
    packCacheOrder(packedN, packedLo, packedHi);
    result = index.classifyCachePacked(packedLo, packedHi, UINT64_C(1) << 7);
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "packed N-mask positions must retain read-coordinate orientation");

    packedN = probeA;
    packedN[41] = 'A';
    packCacheOrder(packedN, packedLo, packedHi);
    result = index.classifyCachePacked(packedLo, packedHi, UINT64_C(1) << 41);
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "packed right-half N-mask positions must retain read-coordinate orientation");

    result = index.classifyCachePacked(
        packedLo, packedHi, (UINT64_C(1) << 40) | (UINT64_C(1) << 41));
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "two Ns in one packed half must not erase a unique opposite-half anchor");

    result = index.classifyCachePacked(
        packedLo, packedHi,
        (UINT64_C(1) << 6) | (UINT64_C(1) << 7) |
        (UINT64_C(1) << 40) | (UINT64_C(1) << 41));
    require(result.status == FlexProbeHalfIndex::None,
            "packed halves with more than one N each must not anchor");

    std::cout << "PASS: Flex split-half active-probe index\n";
    return 0;
}
