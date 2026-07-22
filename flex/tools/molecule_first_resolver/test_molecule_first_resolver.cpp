#include "solo/MoleculeFirstResolver.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

molecule_first::CandidateRead read(const std::string &id, const std::string &umi,
                                   const std::vector<std::pair<std::string, double> > &scores)
{
    molecule_first::CandidateRead output;
    output.readId = id;
    output.featureId = "gene";
    output.rawUmi = umi;
    for (const auto &entry : scores) {
        molecule_first::CandidateScore score;
        score.candidate = entry.first;
        score.logLikelihood = entry.second;
        output.scores.push_back(score);
    }
    return output;
}

bool near(double left, double right, double tolerance = 1e-12)
{
    return std::fabs(left - right) <= tolerance;
}

void expectInvalid(const std::function<void()> &operation)
{
    bool threw = false;
    try {
        operation();
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    assert(threw);
}

} // namespace

int main()
{
    using molecule_first::CandidateRead;
    const molecule_first::PriorCounts priors = {{"A", 0}, {"B", 0}, {"C", 0}, {"D", 0}};
    molecule_first::Config config;

    const std::vector<CandidateRead> chain = {
        read("a", "AAAA", {{"A", 0.0}, {"B", 0.0}}),
        read("b", "AAAA", {{"B", 0.0}, {"C", 0.0}}),
        read("c", "AAAA", {{"C", 0.0}, {"D", 0.0}}),
    };
    const std::vector<molecule_first::ReadClique> chainCliques =
        molecule_first::buildReadCliques(chain, priors, config);
    assert(chainCliques.size() == 2);
    for (const auto &clique : chainCliques) {
        assert(!clique.candidates.empty());
    }

    std::vector<CandidateRead> reversed = chain;
    std::reverse(reversed.begin(), reversed.end());
    const std::vector<molecule_first::ReadClique> reversedCliques =
        molecule_first::buildReadCliques(reversed, priors, config);
    assert(chainCliques.size() == reversedCliques.size());
    for (std::size_t index = 0; index < chainCliques.size(); ++index) {
        assert(chainCliques[index].cliqueId == reversedCliques[index].cliqueId);
        assert(chainCliques[index].memberReadIds == reversedCliques[index].memberReadIds);
        assert(chainCliques[index].candidates == reversedCliques[index].candidates);
    }

    const molecule_first::PriorCounts pcrPriors = {{"A", 8}, {"B", 0}};
    const std::vector<CandidateRead> pcr = {
        read("pcr1", "AAAA", {{"A", 0.0}, {"B", 0.0}}),
        read("pcr2", "AAAA", {{"A", 0.0}, {"B", 0.0}}),
    };
    const auto pcrCliques = molecule_first::buildReadCliques(pcr, pcrPriors, config);
    assert(pcrCliques.size() == 1);
    assert(near(pcrCliques[0].posterior[0], 0.9));
    assert(near(pcrCliques[0].posterior[0] + pcrCliques[0].posterior[1], 1.0));

    const std::vector<CandidateRead> separate = {
        read("one", "AAAA", {{"A", 0.0}, {"B", 0.0}}),
        read("two", "AAAT", {{"A", 0.0}, {"C", 0.0}}),
    };
    const auto separateCliques = molecule_first::buildReadCliques(separate, priors, config);
    const auto occupancies = molecule_first::weightedOccupancies(separateCliques, "1mm_cr");
    bool sawA = false;
    for (const auto &row : occupancies) {
        if (row.candidate == "A") {
            assert(near(row.expectedCount, 0.75));
            assert(row.correctedUmi == "AAAA");
            sawA = true;
        }
        assert(row.expectedCount >= 0.0 && row.expectedCount <= 1.0);
    }
    assert(sawA);

    const auto calls = molecule_first::gatedHardCalls(pcrCliques, config);
    assert(!calls[0].assigned);
    molecule_first::Config permissive = config;
    permissive.gateMinPosterior = 0.89;
    permissive.gateMinMargin = 0.79;
    const auto accepted = molecule_first::gatedHardCalls(pcrCliques, permissive);
    assert(accepted[0].assigned && accepted[0].candidate == "A");

    const auto strict = molecule_first::policyMolecules(pcrCliques, "exact", "strict", calls);
    const auto hard = molecule_first::policyMolecules(pcrCliques, "exact", "hard", calls);
    const auto gated = molecule_first::policyMolecules(pcrCliques, "exact", "gated_hard", calls);
    assert(strict.empty());
    assert(hard.size() == 1 && hard[0].candidate == "A");
    assert(gated.empty());

    molecule_first::ReadClique nearTie;
    nearTie.cliqueId = "near_tie";
    nearTie.featureId = "gene";
    nearTie.rawUmi = "AAAA";
    nearTie.memberReadIds = {"near_tie_read"};
    nearTie.candidates = {"A", "B"};
    nearTie.posterior = {0.5, 0.5 + 5e-15};
    assert(molecule_first::topCandidate(nearTie).first == "A");

    molecule_first::ReadClique lexicalRoot;
    lexicalRoot.cliqueId = "lexical_root";
    lexicalRoot.featureId = "gene";
    lexicalRoot.rawUmi = "AAAA";
    lexicalRoot.candidates = {"A"};
    lexicalRoot.posterior = {0.5};
    molecule_first::ReadClique numericalRoot = lexicalRoot;
    numericalRoot.cliqueId = "numerical_root";
    numericalRoot.rawUmi = "TAAA";
    numericalRoot.posterior = {0.5 + 5e-15};
    const auto nearTieCorrections = molecule_first::correctedUmis(
        {lexicalRoot, numericalRoot}, "1mm_cr");
    assert(nearTieCorrections.at(std::make_tuple("gene", "A", "AAAA")) == "AAAA");
    assert(nearTieCorrections.at(std::make_tuple("gene", "A", "TAAA")) == "AAAA");

    for (const auto &clique : pcrCliques) {
        for (const std::string &candidate : clique.candidates) {
            assert(candidate == "A" || candidate == "B");
        }
    }

    const molecule_first::PriorCounts missingPrior = {{"A", 1}};
    expectInvalid([&]() {
        molecule_first::buildReadCliques(pcr, missingPrior, config);
    });
    molecule_first::Config invalidConfig = config;
    invalidConfig.temperature = 0.0;
    expectInvalid([&]() {
        molecule_first::buildReadCliques(pcr, pcrPriors, invalidConfig);
    });

    std::cout << "molecule-first native unit tests passed\n";
    return 0;
}
