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

std::string mutateAway(std::string sequence, std::size_t begin,
                       std::size_t count)
{
    for (std::size_t position = begin; position < begin + count; ++position)
        sequence[position] = sequence[position] == 'A' ? 'C' : 'A';
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
    require(result.hammingDistance == 0,
            "an exact two-half match must report Hamming zero");

    const std::string oneMismatch = mutate(probeA, 3, 'T');
    result = index.classify(oneMismatch.data(), oneMismatch.size());
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "one mismatch in a half must retain its unique gene");
    require(result.hammingDistance == 1,
            "one mismatch must be reported by the seed-and-extend result");

    uint64_t packedLo = 0;
    uint64_t packedHi = 0;
    packCacheOrder(oneMismatch, packedLo, packedHi);
    result = index.classifyCachePacked(packedLo, packedHi);
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "packed CBQ/cache-order lookup must match ASCII lookup");

    const std::string tenMismatch = mutateAway(probeA, 25, 10);
    result = index.classify(tenMismatch.data(), tenMismatch.size());
    require(result.status == FlexProbeHalfIndex::Unique && result.geneIdx15 == 11,
            "an exact half seed plus ten total mismatches must pass");
    require(result.hammingDistance == 10,
            "the inclusive ten-mismatch boundary must be reported");

    const std::string elevenMismatch = mutateAway(probeA, 25, 11);
    result = index.classify(elevenMismatch.data(), elevenMismatch.size());
    require(result.status == FlexProbeHalfIndex::ScoreFail &&
                result.hammingDistance == 11,
            "an exact half seed plus eleven total mismatches must fail");

    packCacheOrder(tenMismatch, packedLo, packedHi);
    result = index.classifyCachePacked(packedLo, packedHi);
    require(result.status == FlexProbeHalfIndex::Unique &&
                result.hammingDistance == 10,
            "packed CBQ lookup must include the ten-mismatch boundary");
    packCacheOrder(elevenMismatch, packedLo, packedHi);
    result = index.classifyCachePacked(packedLo, packedHi);
    require(result.status == FlexProbeHalfIndex::ScoreFail &&
                result.hammingDistance == 11,
            "packed CBQ lookup must reject eleven total mismatches");

    std::string h1PlusNine = mutateAway(probeA, 0, 1);
    h1PlusNine = mutateAway(h1PlusNine, 25, 9);
    result = index.classify(h1PlusNine.data(), h1PlusNine.size());
    require(result.status == FlexProbeHalfIndex::Unique &&
                result.hammingDistance == 10,
            "an H1 half seed leaves nine mismatches for the other half");

    const std::string h1PlusTen = mutateAway(h1PlusNine, 34, 1);
    result = index.classify(h1PlusTen.data(), h1PlusTen.size());
    require(result.status == FlexProbeHalfIndex::ScoreFail &&
                result.hammingDistance == 11,
            "an H1 half seed plus ten other-half mismatches must fail");

    std::string conflicting = probeA.substr(0, 25) + probeB.substr(25, 25);
    result = index.classify(conflicting.data(), conflicting.size());
    require(result.status == FlexProbeHalfIndex::SplitProbe,
            "qualifying halves from different probes must be rejected");

    FlexProbeHalfIndex sameGeneSplitIndex;
    const std::string probeC =
        "CCCCCCCCCCCCCCCCCCCCCCCCC"
        "TTTTTTTTTTTTTTTTTTTTTTTTT";
    require(sameGeneSplitIndex.build({{probeA, 11}, {probeC, 11}}, &error),
            error.c_str());
    const std::string sameGeneSplit =
        probeA.substr(0, 25) + probeC.substr(25, 25);
    result = sameGeneSplitIndex.classify(
        sameGeneSplit.data(), sameGeneSplit.size());
    require(result.status == FlexProbeHalfIndex::SplitProbe,
            "two probe IDs remain split even when they target one gene");

    FlexProbeHalfIndex collisionIndex;
    const std::string sharedHalfProbe =
        probeA.substr(0, 25) + probeB.substr(25, 25);
    require(collisionIndex.build({{probeA, 11}, {sharedHalfProbe, 22}}, &error),
            error.c_str());
    result = collisionIndex.classify(
        elevenMismatch.data(), elevenMismatch.size());
    require(result.status == FlexProbeHalfIndex::Ambiguous,
            "a half key belonging to multiple probes must remain ambiguous");

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
            "two Ns in one half count toward the full-probe mismatch budget");

    result = index.classifyCachePacked(
        packedLo, packedHi,
        (UINT64_C(1) << 6) | (UINT64_C(1) << 7) |
        (UINT64_C(1) << 40) | (UINT64_C(1) << 41));
    require(result.status == FlexProbeHalfIndex::None,
            "packed halves with more than one N each must not anchor");

    std::cout << "PASS: Flex split-half active-probe index\n";
    return 0;
}
