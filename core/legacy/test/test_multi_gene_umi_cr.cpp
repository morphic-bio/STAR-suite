#include "MultiGeneUmiCr.h"

#include <cassert>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>
#include <vector>

namespace mg = multi_gene_umi_cr;

namespace {

struct Expected {
    bool accepted = false;
    std::uint32_t gene = UINT32_MAX;
    std::string reason;
    Expected(bool acceptedIn, std::uint32_t geneIn, const std::string &reasonIn)
        : accepted(acceptedIn), gene(geneIn), reason(reasonIn) {}
};

Expected legacyReference(const std::map<std::uint32_t, std::uint64_t> &corrected,
                         const std::map<std::uint32_t, std::uint64_t> &original)
{
    std::uint64_t maximum = 0;
    std::uint32_t winner = UINT32_MAX;
    for (const auto &entry : corrected) {
        if (entry.second > maximum) {
            maximum = entry.second;
            winner = entry.first;
        } else if (entry.second == maximum) {
            winner = UINT32_MAX;
        }
    }
    if (winner == UINT32_MAX) return {false, UINT32_MAX, "corrected_count_tie"};
    const std::uint64_t winnerOriginal = original.count(winner) ? original.at(winner) : 0;
    for (const auto &entry : original) {
        if (entry.second > winnerOriginal) {
            return {false, UINT32_MAX, "original_umi_dominance"};
        }
    }
    return {true, winner, "unique_corrected_count_winner"};
}

Expected bridgeReference(
    const std::vector<std::pair<std::uint32_t, std::uint64_t> > &corrected,
    const std::vector<std::pair<std::uint32_t, std::uint64_t> > &original)
{
    std::uint64_t maximum = 0;
    std::uint32_t winner = UINT32_MAX;
    for (const auto &entry : corrected) {
        if (entry.second > maximum) {
            maximum = entry.second;
            winner = entry.first;
        } else if (entry.second == maximum) {
            winner = UINT32_MAX;
        }
    }
    if (winner == UINT32_MAX) return {false, UINT32_MAX, "corrected_count_tie"};
    std::uint64_t winnerOriginal = 0;
    for (const auto &entry : original) {
        if (entry.first == winner) winnerOriginal = entry.second;
    }
    for (const auto &entry : original) {
        if (entry.second > winnerOriginal) {
            return {false, UINT32_MAX, "original_umi_dominance"};
        }
    }
    return {true, winner, "unique_corrected_count_winner"};
}

void assertSame(const mg::Result &actual, const Expected &expected)
{
    assert(actual.accepted == expected.accepted);
    assert(actual.reason == expected.reason);
    if (expected.accepted) assert(actual.gene == expected.gene);
}

void randomizedParity()
{
    std::mt19937_64 random(0x4d554c544947454eULL);
    for (unsigned iteration = 0; iteration < 20000; ++iteration) {
        const std::size_t geneCount = 1 + random() % 8;
        std::vector<std::uint32_t> genes;
        for (std::size_t index = 0; index < geneCount; ++index) {
            genes.push_back(static_cast<std::uint32_t>(100 + iteration * 11 + index));
        }
        std::shuffle(genes.begin(), genes.end(), random);

        std::map<std::uint32_t, std::uint64_t> legacyCorrected;
        std::map<std::uint32_t, std::uint64_t> legacyOriginal;
        std::vector<std::pair<std::uint32_t, std::uint64_t> > bridgeCorrected;
        std::vector<std::pair<std::uint32_t, std::uint64_t> > bridgeOriginal;
        std::vector<mg::GeneSupport> legacySupports;
        for (std::uint32_t gene : genes) {
            const std::uint64_t corrected = 1 + random() % 50;
            const std::uint64_t original = random() % (corrected + 1);
            legacyCorrected[gene] = corrected;
            if (original != 0) legacyOriginal[gene] = original;
            legacySupports.push_back({gene, corrected, original});
            bridgeCorrected.push_back({gene, corrected});
            if (original != 0) bridgeOriginal.push_back({gene, original});
        }
        std::sort(bridgeCorrected.begin(), bridgeCorrected.end());
        std::sort(bridgeOriginal.begin(), bridgeOriginal.end());

        assertSame(mg::resolve(legacySupports),
                   legacyReference(legacyCorrected, legacyOriginal));

        std::vector<mg::GeneSupport> bridgeSupports;
        std::size_t originalIndex = 0;
        for (const auto &entry : bridgeCorrected) {
            while (originalIndex < bridgeOriginal.size()
                   && bridgeOriginal[originalIndex].first < entry.first) {
                ++originalIndex;
            }
            const std::uint64_t original = originalIndex < bridgeOriginal.size()
                    && bridgeOriginal[originalIndex].first == entry.first
                ? bridgeOriginal[originalIndex].second : 0;
            bridgeSupports.push_back({entry.first, entry.second, original});
        }
        assertSame(mg::resolve(bridgeSupports),
                   bridgeReference(bridgeCorrected, bridgeOriginal));
    }
}

} // namespace

int main()
{
    mg::Result result = mg::resolve({{7, 5, 2}});
    assert(result.accepted && result.gene == 7);

    result = mg::resolve({{3, 10, 4}, {8, 9, 4}});
    assert(result.accepted && result.gene == 3);

    result = mg::resolve({{3, 10, 4}, {8, 10, 3}});
    assert(!result.accepted && result.reason == "corrected_count_tie");

    // Gene 3 wins after correction, but gene 8 has more reads whose original
    // UMI was already the corrected UMI. STAR rejects this molecule.
    result = mg::resolve({{3, 10, 2}, {8, 9, 3}});
    assert(!result.accepted && result.reason == "original_umi_dominance");

    result = mg::resolve({{3, 10, 3}, {8, 9, 3}});
    assert(result.accepted && result.gene == 3);

    bool threw = false;
    try { mg::resolve({{1, 1, 1}, {1, 2, 1}}); }
    catch (const std::invalid_argument &) { threw = true; }
    assert(threw);
    randomizedParity();
    std::cout << "MultiGeneUMI_CR helper and 20,000-case legacy/bridge parity tests passed\n";
    return 0;
}
