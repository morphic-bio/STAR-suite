#include "FlexProbeRegion.h"
#include <cassert>
#include <cstdint>

int main()
{
    // These are STAR wire-format tests, not biological estimator fixtures.
    for (uint32_t region = 0; region != 4; ++region) {
        const uint32_t wire = 123u | (region << 30);
        assert(flexProbeValueCount(wire) == 123u);
        assert(static_cast<uint32_t>(flexProbeValueRegion(wire)) == region);
        assert(flexProbePackValue(123, static_cast<FlexProbeRegion>(region)) == wire);
        assert(flexProbeValueCount(flexProbeMergeValue(wire, 7u)) == 130u);
        assert(flexProbeValueRegion(flexProbeMergeValue(wire, 7u))
               == static_cast<FlexProbeRegion>(region));
    }
    const uint32_t spliced = flexProbePackValue(11, FlexProbeRegionSpliced);
    const uint32_t unspliced = flexProbePackValue(13, FlexProbeRegionUnspliced);
    assert(flexProbeValueCount(flexProbeMergeValue(spliced, unspliced)) == 24u);
    assert(flexProbeValueRegion(flexProbeMergeValue(spliced, unspliced))
           == FlexProbeRegionConflicting);
    assert(flexProbeValueCount(flexProbeMergeValue(UINT32_MAX, UINT32_MAX))
           == kFlexProbeCountMask);
    assert(flexProbeValueCount(flexProbePackValue(UINT32_MAX, FlexProbeRegionUnknown))
           == kFlexProbeCountMask);
}
