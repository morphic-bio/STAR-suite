#include "solo/CbBayesianResolver.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

struct OracleResult {
    bool resolved = false;
    uint32_t bestIdx = 0;
    double posterior = 0.0;
};

OracleResult legacyOracle(const std::vector<uint32_t> &candidateIdx,
                          const std::vector<uint8_t> &candidateQual,
                          const std::vector<uint32_t> &exactCounts,
                          int qsBase,
                          uint32_t qsMax,
                          double minPosterior)
{
    OracleResult out;
    double total = 0.0;
    double best = 0.0;
    for (size_t ii = 0; ii < candidateIdx.size(); ++ii) {
        const uint32_t idx0 = candidateIdx[ii] - 1u;
        if (exactCounts[idx0] == 0) {
            continue;
        }
        int quality = static_cast<int>(candidateQual[ii]) - qsBase;
        quality = quality < static_cast<int>(qsMax)
            ? quality : static_cast<int>(qsMax);
        const double weight = exactCounts[idx0]
            * std::pow(10.0, -static_cast<double>(quality) / 10.0);
        total += weight;
        if (weight > best) {
            best = weight;
            out.bestIdx = candidateIdx[ii];
        }
    }
    if (total > 0.0) {
        out.posterior = best / total;
        out.resolved = best >= minPosterior * total;
    }
    return out;
}

void require(bool condition, const char *message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

void compareWithOracle(const CbBayesianResolver &resolver,
                       const std::vector<uint32_t> &candidateIdx,
                       const std::vector<uint8_t> &candidateQual,
                       const std::vector<uint32_t> &exactCounts,
                       double minPosterior)
{
    const int qsBase = 33;
    const uint32_t qsMax = 40;
    const OracleResult expected = legacyOracle(candidateIdx, candidateQual,
                                                exactCounts, qsBase, qsMax,
                                                minPosterior);
    const BayesianResult actual = resolver.resolveStarSolo(
        candidateIdx, candidateQual, exactCounts, qsBase, qsMax, minPosterior);

    require((actual.status == BayesianResult::Resolved) == expected.resolved,
            "resolution status differs from legacy posterior");
    require(actual.bestIdx == expected.bestIdx,
            "winning whitelist index differs from legacy posterior");
    require(std::fabs(actual.posteriorBest - expected.posterior) < 1e-12,
            "posterior differs from legacy posterior");
}

} // namespace

int main()
{
    CbBayesianResolver resolver(4, nullptr);

    // Abundance prior resolves otherwise identical mismatch evidence.
    compareWithOracle(resolver, {1, 2}, {'I', 'I'}, {100, 1, 0, 0}, 0.975);

    // Lower quality at the mismatching base makes that candidate more likely.
    compareWithOracle(resolver, {1, 2}, {'+', '?'}, {10, 10, 0, 0}, 0.975);

    // Equal priors and equal mismatch qualities remain ambiguous.
    const BayesianResult tie = resolver.resolveStarSolo(
        {1, 2}, {'I', 'I'}, {1, 1, 0, 0}, 33, 40, 0.975);
    require(tie.status == BayesianResult::Ambiguous,
            "equal candidates should remain ambiguous");

    // The configured cbMinP threshold, not a separate hard-coded threshold or
    // runner-up ratio, controls acceptance.
    compareWithOracle(resolver, {1, 2}, {'I', 'I'}, {100, 1, 0, 0}, 0.995);

    // Qualities above QSmax are capped exactly as in SoloReadFeature::inputRecords.
    compareWithOracle(resolver, {1, 2}, {'I', ']'}, {7, 7, 0, 0}, 0.975);

    const BayesianResult malformed = resolver.resolveStarSolo(
        {1, 2}, {'I'}, {1, 1, 0, 0}, 33, 40, 0.975);
    require(malformed.status == BayesianResult::Unresolved,
            "malformed candidate-quality payload should not resolve");

    std::cout << "PASS: STARsolo ambiguous-CB posterior contract\n";
    return 0;
}
