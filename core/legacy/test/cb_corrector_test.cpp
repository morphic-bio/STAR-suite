#include "solo/CbCorrector.h"

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

uint32_t packCbq(const std::string &sequence) {
    uint32_t packed = 0;
    for (size_t i = 0; i < sequence.size(); ++i) {
        uint32_t code = 0;
        switch (sequence[i]) {
            case 'A': code = 0; break;
            case 'C': code = 1; break;
            case 'G': code = 2; break;
            case 'T': code = 3; break;
            default: return UINT32_MAX;
        }
        packed |= code << (2 * i);
    }
    return packed;
}

bool require(bool condition, const char *message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
    }
    return condition;
}

} // namespace

int main() {
    const std::vector<std::string> whitelist = {
        "AAAAAAAA", "AAAAAAAC", "CCCCCCCC"
    };
    CbCorrector legacy(whitelist, 1);
    CbCorrector cbq(whitelist, 1, true);

    bool ok = true;
    ok &= require(legacy.correct("AAAAAAAA").whitelistIdx == 1,
                  "legacy exact string lookup");
    ok &= require(cbq.correct("AAAAAAAA").whitelistIdx == 1,
                  "CBQ-order exact string lookup");
    ok &= require(cbq.correct("GCCCCCCC").whitelistIdx == 3,
                  "CBQ-order H1 string lookup");
    ok &= require(cbq.correct("AAAAAAAG").ambiguous,
                  "CBQ-order ambiguous string lookup");
    ok &= require(cbq.correct("NCCCCCCC").whitelistIdx == 3,
                  "CBQ-order N expansion");
    ok &= require(cbq.decodePackedKey(packCbq("ACGTACGT"), 8) == "ACGTACGT",
                  "CBQ-order decode");

    uint32_t index = 0;
    uint8_t distance = 255;
    ok &= require(cbq.correctPackedCbq(packCbq("AAAAAAAA"), index, distance)
                      && index == 1 && distance == 0,
                  "packed CBQ exact lookup");
    ok &= require(cbq.correctPackedCbq(packCbq("GCCCCCCC"), index, distance)
                      && index == 3 && distance == 1,
                  "packed CBQ unique H1 lookup");
    ok &= require(!cbq.correctPackedCbq(packCbq("AAAAAAAG"), index, distance),
                  "packed CBQ ambiguous lookup falls back");
    ok &= require(!cbq.correctPackedCbq(packCbq("GGGGGGGG"), index, distance),
                  "packed CBQ miss falls back");
    ok &= require(!legacy.correctPackedCbq(packCbq("AAAAAAAA"), index, distance),
                  "legacy-order table rejects native packed lookup");

    static const char bases[] = "ACGT";
    std::string sequence(8, 'A');
    for (uint32_t packed = 0; packed < (1u << 16); ++packed) {
        for (size_t pos = 0; pos < sequence.size(); ++pos) {
            sequence[pos] = bases[(packed >> (2 * pos)) & 3u];
        }
        const CbMatch legacyMatch = legacy.correct(sequence);
        const CbMatch cbqMatch = cbq.correct(sequence);
        if (legacyMatch.whitelistIdx != cbqMatch.whitelistIdx ||
            legacyMatch.hammingDist != cbqMatch.hammingDist ||
            legacyMatch.ambiguous != cbqMatch.ambiguous ||
            legacyMatch.ambiguousIdx != cbqMatch.ambiguousIdx) {
            std::cerr << "FAIL: table-order mismatch for " << sequence << '\n';
            ok = false;
            break;
        }

        uint32_t packedIndex = 0;
        uint8_t packedDistance = 255;
        const bool packedMatched = cbq.correctPackedCbq(
            packed, packedIndex, packedDistance);
        const bool uniqueStringMatch = !cbqMatch.ambiguous
            && cbqMatch.whitelistIdx != 0 && cbqMatch.hammingDist <= 1;
        if (packedMatched != uniqueStringMatch ||
            (packedMatched && (packedIndex != cbqMatch.whitelistIdx ||
                               packedDistance != cbqMatch.hammingDist))) {
            std::cerr << "FAIL: packed/string mismatch for " << sequence << '\n';
            ok = false;
            break;
        }
    }

    if (!ok) {
        return 1;
    }
    std::cout << "PASS: CbCorrector legacy and CBQ-native lookup\n";
    return 0;
}
