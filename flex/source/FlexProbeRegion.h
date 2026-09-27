#ifndef H_FlexProbeRegion
#define H_FlexProbeRegion

#include <cstdint>

// STAR's existing packed count/cache representation reserves the high two
// bits for probe region. Keep decoding old caches without treating these bits
// as read counts. No diagnostic or probe-metadata loading is implemented here.
enum FlexProbeRegion : uint8_t {
    FlexProbeRegionUnknown = 0,
    FlexProbeRegionSpliced = 1,
    FlexProbeRegionUnspliced = 2,
    FlexProbeRegionConflicting = 3
};

static const uint32_t kFlexProbeCountMask = 0x3FFFFFFFu;
static const uint32_t kFlexProbeRegionShift = 30u;

inline FlexProbeRegion flexProbeMergeRegion(FlexProbeRegion lhs, FlexProbeRegion rhs)
{
    if (lhs == FlexProbeRegionConflicting || rhs == FlexProbeRegionConflicting)
        return FlexProbeRegionConflicting;
    if (lhs == FlexProbeRegionUnknown)
        return rhs;
    if (rhs == FlexProbeRegionUnknown)
        return lhs;
    return lhs == rhs ? lhs : FlexProbeRegionConflicting;
}

inline uint32_t flexProbeValueCount(uint32_t value)
{
    return value & kFlexProbeCountMask;
}

inline FlexProbeRegion flexProbeValueRegion(uint32_t value)
{
    return static_cast<FlexProbeRegion>((value >> kFlexProbeRegionShift) & 0x3u);
}

inline uint32_t flexProbePackValue(uint32_t count, FlexProbeRegion region)
{
    if (count > kFlexProbeCountMask)
        count = kFlexProbeCountMask;
    return count | (static_cast<uint32_t>(region) << kFlexProbeRegionShift);
}

inline uint32_t flexProbeMergeValue(uint32_t lhs, uint32_t rhs)
{
    const uint64_t sum = static_cast<uint64_t>(flexProbeValueCount(lhs))
                       + static_cast<uint64_t>(flexProbeValueCount(rhs));
    const uint32_t count = sum > kFlexProbeCountMask
        ? kFlexProbeCountMask : static_cast<uint32_t>(sum);
    return flexProbePackValue(
        count, flexProbeMergeRegion(flexProbeValueRegion(lhs), flexProbeValueRegion(rhs)));
}

#endif
