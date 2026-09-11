#include "BridgeReadInfo.h"
#include <cassert>
#include <map>
#include <tuple>

int main() {
    BridgeReadInfo info, worker;
    // Same original UMI in different genes must follow each gene's decision.
    worker.reads = {{2, 100, 7, 1}, {2, 100, 8, 2}, {2, 100, 7, 3},
                    {9, 100, 7, 4}, {2, 55, UINT32_MAX, 5},
                    {2, 90, 7, 6}, {3, 100, 7, 7}};
    worker.pending = {{99, {0, 100, 7, 8}}, {98, {0, 100, 7, 9}},
                      {97, {0, 45, UINT32_MAX, 10}}};
    info.resolvedCBs = {{99, 2}, {97, 3}};
    info.append(worker);
    assert(worker.reads.empty() && worker.pending.empty());
    std::map<uint32_t, std::tuple<uint32_t, uint32_t, uint8_t>> result;
    auto sink = [&](uint32_t read, uint32_t cb, uint32_t umi, uint8_t status) {
        assert(result.emplace(read, std::make_tuple(cb, umi, status)).second);
    };
    info.prepare([](uint32_t cb) { return cb != 9; });
    std::vector<BridgeReadInfo::Decision> decisions = {
        {8, 100, UINT32_MAX}, {7, 100, 101}, {7, 90, 90}};
    info.replay(2, decisions, sink);
    decisions = {{7, 100, 100}};
    info.replay(3, decisions, sink);
    assert(result.size() == 6);
    assert(result.count(4) == 0 && result.count(9) == 0); // gated CB / unresolved CB
    for (uint32_t read : {1u, 3u, 8u})
        assert(result.at(read) == std::make_tuple(2u, 101u, uint8_t{1}));
    assert(result.at(2) == std::make_tuple(2u, UINT32_MAX, uint8_t{2}));
    assert(result.count(5) == 0); // no assigned gene, as in CountingSink
    assert(result.at(6) == std::make_tuple(2u, 90u, uint8_t{1}));
    assert(result.at(7) == std::make_tuple(3u, 100u, uint8_t{1}));
    assert(result.count(10) == 0);
}
