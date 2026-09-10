#include "OcmCountStore.h"
#include <cassert>
#include <iostream>
#include <stdexcept>

static void checkRecords(ocm::CountStore& store, uint64_t expected, uint32_t salt) {
    uint64_t i = 0;
    store.forEach([&](const ocm::CountRecord& record) {
        assert(record.gene == (i * 17 + salt) % 313);
        assert(record.cell == (i * 29 + salt) % 103);
        assert(record.count == i + 1);
        ++i;
    });
    assert(i == expected);
}

int main() {
    // Two stores compete for the same blocks. Exercise empty/partial/exact
    // block boundaries, a spill after retention, and an immediate spill.
    for (uint64_t limit : {uint64_t(0), uint64_t(ocm::CountStore::blockBytes),
                           uint64_t(8 * ocm::CountStore::blockBytes)}) {
        ocm::CountBufferBudget budget(limit);
        {
            ocm::CountStore a(budget, "/tmp"), b(budget, "/tmp"), empty(budget, "/tmp");
            const uint64_t n = 2 * ocm::CountStore::blockRecords + 19;
            for (uint64_t i = 0; i < n; ++i) {
                a.append({uint32_t(i * 17 % 313), uint32_t(i * 29 % 103), uint32_t(i + 1)});
                if (i < ocm::CountStore::blockRecords)
                    b.append({uint32_t((i * 17 + 3) % 313), uint32_t((i * 29 + 3) % 103), uint32_t(i + 1)});
                assert(budget.used <= limit);
            }
            a.finish(); b.finish(); empty.finish();
            checkRecords(a, n, 0);
            checkRecords(b, ocm::CountStore::blockRecords, 3);
            checkRecords(empty, 0, 0);
            // A callback failure must not make a later replay start mid-file.
            try {
                a.forEach([](const ocm::CountRecord&) { throw std::runtime_error("injected"); });
                assert(false);
            } catch (const std::runtime_error&) {}
            checkRecords(a, n, 0);
            if (limit == 0) assert(a.spilled() && b.spilled());
            if (limit >= 8 * ocm::CountStore::blockBytes) assert(!a.spilled() && !b.spilled());
            bool rejected = false;
            try { a.append({0, 0, 1}); } catch (const std::logic_error&) { rejected = true; }
            assert(rejected);
        }
        assert(budget.used == 0);
    }
    ocm::CountBufferBudget budget(ocm::CountStore::blockBytes);
    try {
        ocm::CountStore store(budget, "/nonexistent/ocm-test");
        for (size_t i = 0; i <= ocm::CountStore::blockRecords; ++i) store.append({0, 0, 1});
        assert(false);
    } catch (const std::runtime_error&) {}
    assert(budget.used == 0);
    std::cout << "PASS: OCM count retention, spill, replay and cleanup\n";
}
