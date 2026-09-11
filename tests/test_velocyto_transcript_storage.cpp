// The preserved translation unit supplies the actual former merge procedure.
#include BASELINE_SOURCE
#include "PooledTranscriptMap.h"
#include <cassert>
#include <random>

using Pool = PooledTranscriptMap<trTypeStruct>;
static void compare(const std::unordered_map<uintUMI, std::vector<trTypeStruct>>& expected,
                    const Pool& actual) {
    assert(expected.size() == actual.size());
    size_t live = 0;
    for (const auto entry : actual) {
        const auto& old = expected.at(entry.first);
        assert(old.size() == entry.second.size());
        assert(actual.rejected(entry.first) == old.empty());
        for (size_t i = 0; i < old.size(); ++i) {
            assert(old[i].tr == entry.second[i].tr);
            assert(old[i].type == entry.second[i].type);
        }
        live += old.size();
    }
    assert(actual.liveSlots() == live && actual.allocatedSlots() >= live);
}
int main() {
    std::mt19937 rng(333444);
    std::unordered_map<uintUMI, std::vector<trTypeStruct>> old;
    Pool current;
    current.reserve(2000);
    for (unsigned n = 0; n < 60000; ++n) {
        uintUMI umi = n < 2000 ? n : rng() % 2000;
        if (n % 311 == 0) umi = UINT32_MAX;
        std::vector<trTypeStruct> incoming;
        for (unsigned t = 0; t < 12; ++t)
            if (rng() % 3 != 0) incoming.push_back({umi % 7 * 100 + t, static_cast<uint8_t>(rng() % 16)});
        if (n % 401 == 0) incoming.clear();
        if (n % 499 == 0 && !incoming.empty()) incoming.push_back(incoming.back());
        applyVelocytoMerge(old, umi, incoming);
        current.merge(umi, incoming);
        if (n % 997 == 0) compare(old, current);
    }
    compare(old, current);
    std::vector<Pool> moved;
    moved.push_back(std::move(current));
    for (unsigned i = 0; i < 40; ++i) moved.emplace_back();
    compare(old, moved.front());
    assert(current.size() == 0 && current.liveSlots() == 0);
    current = std::move(moved.front());
    compare(old, current);
    std::cout << "Velocyto pooled transcript parity PASS: keys=" << current.size()
              << " slots=" << current.allocatedSlots() << " live=" << current.liveSlots() << '\n';
    current.clearAndFree();
    assert(current.size() == 0 && current.allocatedSlots() == 0 && !current.rejected(0));
    bool rejected = false;
    try { current.reserve(SIZE_MAX); } catch (const std::length_error&) { rejected = true; }
    assert(rejected);
}
