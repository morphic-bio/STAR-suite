#include "SampleTagCache.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

std::uint16_t encode(const std::string& sequence)
{
    if (sequence.size() != 8u) std::abort();
    std::uint16_t code = 0;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        unsigned base = 0;
        switch (sequence[i]) {
            case 'A': base = 0; break;
            case 'C': base = 1; break;
            case 'G': base = 2; break;
            case 'T': base = 3; break;
            default: std::abort();
        }
        code = static_cast<std::uint16_t>(code | (base << (2u * i)));
    }
    return code;
}

unsigned hamming(std::uint16_t lhs, std::uint16_t rhs)
{
    const std::uint16_t x = static_cast<std::uint16_t>(lhs ^ rhs);
    const std::uint16_t differing = static_cast<std::uint16_t>((x | (x >> 1u)) & 0x5555u);
    return static_cast<unsigned>(__builtin_popcount(static_cast<unsigned>(differing)));
}

void require(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

star::flex::SampleTagCache::Result oracle(
    std::uint16_t query,
    const std::vector<std::uint16_t>& codes,
    const std::vector<std::uint16_t>& owners,
    bool allowMismatch)
{
    std::set<std::uint16_t> matches;
    for (std::size_t i = 0; i < codes.size(); ++i) {
        if (codes[i] == query) matches.insert(owners[i]);
    }
    if (!matches.empty()) {
        if (matches.size() == 1u) {
            return star::flex::SampleTagCache::Result(
                *matches.begin(), star::flex::SampleTagCache::Tier::H0, false);
        }
        return star::flex::SampleTagCache::Result(
            0u, star::flex::SampleTagCache::Tier::H0, true);
    }
    if (!allowMismatch) return star::flex::SampleTagCache::Result();

    for (std::size_t i = 0; i < codes.size(); ++i) {
        if (hamming(codes[i], query) == 1u) matches.insert(owners[i]);
    }
    if (matches.empty()) return star::flex::SampleTagCache::Result();
    if (matches.size() == 1u) {
        return star::flex::SampleTagCache::Result(
            *matches.begin(), star::flex::SampleTagCache::Tier::H1, false);
    }
    return star::flex::SampleTagCache::Result(
        0u, star::flex::SampleTagCache::Tier::H1, true);
}

void compare(const star::flex::SampleTagCache::Result& actual,
             const star::flex::SampleTagCache::Result& expected,
             std::uint16_t query,
             bool allowMismatch)
{
    if (actual.owner != expected.owner || actual.tier != expected.tier ||
        actual.ambiguous != expected.ambiguous) {
        std::cerr << "FAIL: query=" << query
                  << " allowMismatch=" << allowMismatch
                  << " actual(owner,tier,ambiguous)=" << actual.owner << ','
                  << static_cast<unsigned>(actual.tier) << ',' << actual.ambiguous
                  << " expected=" << expected.owner << ','
                  << static_cast<unsigned>(expected.tier) << ',' << expected.ambiguous
                  << '\n';
        std::exit(1);
    }
}

} // namespace

int main()
{
    using star::flex::SampleTagCache;

    const std::vector<std::uint16_t> codes{
        encode("AAAAAAAA"), encode("CCAAAAAA"), encode("TTTTTTTT"),
        encode("ACGTACGT"), encode("ACGTACGT")
    };
    const std::vector<std::uint16_t> owners{1u, 2u, 3u, 4u, 4u};
    SampleTagCache cache;
    cache.build(codes, owners);

    SampleTagCache::Result result = cache.lookup(encode("AAAAAAAA"), true);
    require(result.owner == 1u && result.tier == SampleTagCache::Tier::H0 && !result.ambiguous,
            "exact H0 must return before H1");
    result = cache.lookup(encode("GAAAAAAA"), true);
    require(result.owner == 1u && result.tier == SampleTagCache::Tier::H1 && !result.ambiguous,
            "unique H1 key must resolve");
    result = cache.lookup(encode("CAAAAAAA"), true);
    require(result.owner == 0u && result.tier == SampleTagCache::Tier::H1 && result.ambiguous,
            "cross-owner H1 key must be ambiguous");
    result = cache.lookup(encode("GAAAAAAA"), false);
    require(result.owner == 0u && result.tier == SampleTagCache::Tier::Miss && !result.ambiguous,
            "H1 must not be consulted in exact-only mode");

    // Exhaustively compare every possible 8-mer with an independent slow
    // distance oracle, both with and without the H1 tier.
    for (std::uint32_t query = 0; query < 65536u; ++query) {
        const std::uint16_t packed = static_cast<std::uint16_t>(query);
        compare(cache.lookup(packed, false), oracle(packed, codes, owners, false), packed, false);
        compare(cache.lookup(packed, true), oracle(packed, codes, owners, true), packed, true);
    }

    // An exact key remains authoritative even when another H0 key generates
    // the same sequence as an H1 neighbor.
    SampleTagCache precedence;
    precedence.build(
        std::vector<std::uint16_t>{encode("AAAAAAAA"), encode("CAAAAAAA")},
        std::vector<std::uint16_t>{1u, 2u});
    result = precedence.lookup(encode("CAAAAAAA"), true);
    require(result.owner == 2u && result.tier == SampleTagCache::Tier::H0 && !result.ambiguous,
            "H0 must take precedence over a neighboring H1 generation");

    const SampleTagCache::Stats& stats = cache.stats();
    require(stats.exact == 4u, "duplicate same-owner H0 entry must be harmless");
    require(stats.mismatch1 > 0u, "test table must construct H1 entries");
    require(stats.ambiguous > 0u, "test table must construct ambiguous entries");

    std::cout << "PASS: two-tier sample-tag cache matches exhaustive oracle"
              << " exact=" << stats.exact
              << " h1=" << stats.mismatch1
              << " ambiguous=" << stats.ambiguous << '\n';
    return 0;
}
