#ifndef STAR_SUITE_ADAPTIVE_AMBIENT_WINDOW_H
#define STAR_SUITE_ADAPTIVE_AMBIENT_WINDOW_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

struct AdaptiveAmbientWindow {
    uint32_t start = 0;
    uint32_t end = 0;
    uint64_t umiMass = 0;
};

// Counts must be sorted in descending UMI order.  The established rank window
// is always retained; only its low-count endpoint may be extended.
inline AdaptiveAmbientWindow selectAdaptiveAmbientWindow(
    const std::vector<std::pair<uint32_t, uint32_t>>& sortedUmi,
    uint32_t requestedStart,
    uint32_t requestedBaseEnd,
    uint64_t targetUmiMass)
{
    const size_t n = sortedUmi.size();
    const size_t start = std::min<size_t>(requestedStart, n);
    size_t end = std::min<size_t>(std::max(requestedBaseEnd, requestedStart), n);
    uint64_t mass = 0;
    for (size_t rank = start; rank < end; ++rank) {
        mass += sortedUmi[rank].first;
    }
    while (end < n && mass < targetUmiMass) {
        mass += sortedUmi[end].first;
        ++end;
    }
    AdaptiveAmbientWindow result;
    result.start = static_cast<uint32_t>(start);
    result.end = static_cast<uint32_t>(end);
    result.umiMass = mass;
    return result;
}

#endif
