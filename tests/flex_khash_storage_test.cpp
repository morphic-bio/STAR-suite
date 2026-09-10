#include "FlexHashCacheStorage.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstring>

static void require(bool ok, const std::string& s) { if (!ok) throw std::runtime_error(s); }
int main(int argc, char** argv) {
    require(argc == 4, "usage: test input.bin output.khash corrupt.khash");
    std::string error;
    FlexHashCacheStorage legacy;
    require(legacy.open(argv[1], &error), error);
    require(!legacy.persisted() && legacy.verify(&error), error);
    require(legacy.writeSnapshot(argv[2], &error), error);
    FlexHashCacheStorage stored;
    require(stored.open(argv[2], &error), error);
    require(stored.persisted() && stored.verify(&error), error);
    require(stored.recordCount() == legacy.recordCount(), "record counts");
    for (uint64_t i = 0; i < legacy.recordCount(); ++i) {
        const auto& r = legacy.record(i);
        require(memcmp(&r, &stored.record(i), sizeof(r)) == 0, "embedded records changed");
        for (unsigned t = 0; t < 2; ++t) {
            const auto key = FlexHashCacheStorage::cbqKey(r.key);
            const auto* a = legacy.lookup(t, key); const auto* b = stored.lookup(t, key);
            require(bool(a) == bool(b) && (!a || memcmp(a, b, sizeof(*a)) == 0), "table payload changed");
        }
        for (uint16_t sample = 0; sample < 4; ++sample) for (bool h0 : {false, true}) {
            FlexProbeRecord a {}, b {};
            require(legacy.find(r.key, sample, h0, a) == stored.find(r.key, sample, h0, b) &&
                    memcmp(&a, &b, sizeof(a)) == 0, "sample-aware lookup changed");
        }
    }
    require(!legacy.writeSnapshot(argv[2], &error), "existing output was overwritten");
    require(stored.verify(&error), error);
    std::ifstream in(argv[2], std::ios::binary);
    std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    for (unsigned offset : {8u, 16u, 24u, 40u, 72u, 80u, 112u}) {
        std::string bad = bytes; bad[offset] ^= 1;
        { std::ofstream out(argv[3], std::ios::binary); out.write(bad.data(), bad.size()); }
        FlexHashCacheStorage rejected;
        require(!rejected.open(argv[3], &error), "invalid snapshot metadata accepted");
    }
    { std::ofstream out(argv[3], std::ios::binary); out.write(bytes.data(), bytes.size()-1); }
    FlexHashCacheStorage truncated;
    require(!truncated.open(argv[3], &error), "truncated snapshot accepted");
    std::cout << "PASS: persisted/legacy payloads, sample lookup, immutable output, invalid headers and truncation\n";
}
