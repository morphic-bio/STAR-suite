#include "BridgeReadCounts.h"
#include <cassert>
#include <map>
#include <random>
#include <vector>

static std::map<uint32_t, uint64_t> snapshot(const BridgeReadCounts& counts) {
    std::map<uint32_t, uint64_t> result;
    for (const auto& entry : counts) assert(result.emplace(entry.first, entry.second).second);
    assert(result.size() == counts.size());
    return result;
}
int main() {
    BridgeReadCounts a;
    assert(a.empty() && snapshot(a).empty());
    std::map<uint32_t, uint64_t> expected;
    std::mt19937 random(7311);
    for (unsigned i = 0; i < 200000; ++i) {
        uint32_t key = random() % 50000;
        uint64_t value = i % 2 ? 1ull : 1ull << 32;
        a[key] += value; expected[key] += value;
    }
    a[UINT32_MAX] = UINT64_MAX; expected[UINT32_MAX] = UINT64_MAX;
    a[0] += UINT64_MAX; expected[0] += UINT64_MAX;
    assert(snapshot(a) == expected);
    a.reserve(200000); assert(snapshot(a) == expected);
    bool rejected = false;
    try { a.reserve(SIZE_MAX); } catch (const std::length_error&) { rejected = true; }
    assert(rejected && snapshot(a) == expected);
    BridgeReadCounts copy(a);
    BridgeReadCounts moved(std::move(copy));
    assert(copy.empty() && snapshot(moved) == expected);
    moved.clear(); assert(moved.empty());
    moved[7] = 9; assert(moved.size() == 1 && moved[7] == 9);
    BridgeReadCounts assigned; assigned = a;
    assert(snapshot(assigned) == expected);
    for (const auto& entry : a) assigned[entry.first] += entry.second;
    for (auto& entry : expected) entry.second += entry.second;
    assert(snapshot(assigned) == expected);
    BridgeReadCounts().swap(assigned); assert(assigned.empty());
}
