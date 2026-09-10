#ifndef FLEX_HASH_SCREEN_DECISION_H
#define FLEX_HASH_SCREEN_DECISION_H
#include <cstdint>
#include "FlexGdna.h"

struct FlexHashScreenDecision {
    enum Action : uint8_t {
        Disabled = 0,
        Pass = 1,
        Keep = 2,
        Deny = 3
    };

    Action action = Disabled;
    uint16_t geneIdx15 = 0;
    uint8_t cacheClass = 0;
    uint8_t negativeCode = 0;
    int8_t offset = 0;
    FlexGdnaRegion probeRegion = FlexGdnaUnknown;
    // Diagnostic provenance for the conservative exactly-one-N retry. The
    // runtime cacheClass remains H1 for compatibility; singleNCacheClass
    // identifies the underlying cache tier that supplied the unique gene.
    bool singleN = false;
    uint8_t singleNCacheClass = 0xFF;
    // Total fixed-position 50-base probe distance for an H1X2
    // seed-and-extend decision; 0xFF for other tiers or no scored probe.
    uint8_t probeHammingDistance = 0xFF;
    // Legacy diagnostic field for the former residual-alignment gate. New
    // H1X2 seed-and-extend decisions are terminal and leave this field zero.
    uint16_t residualAnchorGeneIdx15 = 0;
};

enum FlexHashScreenCacheClass : uint8_t {
    FlexHashCacheH0 = 0,
    FlexHashCacheH1 = 1,
    FlexHashCacheNegative = 2,
    FlexHashCacheH2 = 3,
    FlexHashCacheH1X2 = 4
};

// Binary negative class codes from scripts/flex_h01_pilot.py.
enum FlexHashScreenNegativeCode : uint8_t {
    FlexHashNegNone = 0,
    FlexHashNegProbeAmbig = 1,
    FlexHashNegHalfNoAnchor = 2,
    FlexHashNegHalfGeneAmbig = 3,
    FlexHashNegHalfScoreFail = 4,
    FlexHashNegHalfSplitProbe = 5
};

#endif
