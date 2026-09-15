// White-box diagnostic: exercise the production resolver without alignments,
// a reference index, or a slide-sized input. Compile with the audited source
// directory on the include path; do not also link SpatialGexPipeline.o.
#include "SpatialGexPipeline.cpp"
#include <iostream>
#include <random>

using namespace spatial_gex;
using Contribution = spatial_gex::downstream_spool::Contribution;

static void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

static bool equal(const FinalMolecule &a, const FinalMolecule &b) {
    return a.geneIndex == b.geneIndex && a.coordinateIndex == b.coordinateIndex
        && a.correctedUmi == b.correctedUmi && a.policy == b.policy
        && a.weight == b.weight;
}

static std::vector<FinalMolecule> subset(const std::vector<FinalMolecule> &all,
                                        unsigned mask) {
    std::vector<FinalMolecule> result;
    for (const auto &m : all) if (m.policy & mask) result.push_back(m);
    return result;
}

static void same(const std::vector<FinalMolecule> &a,
                 const std::vector<FinalMolecule> &b) {
    require(a.size() == b.size(), "enabled policies changed molecule count");
    for (size_t i = 0; i < a.size(); ++i)
        require(equal(a[i], b[i]), "enabled policies changed molecule identity/weight");
}

static void add(std::vector<Contribution> &rows, unsigned gene, unsigned umi,
                unsigned coordinate, double probability, unsigned count,
                unsigned flags, unsigned clique, unsigned candidate) {
    Contribution c = {};
    c.gene = gene; c.rawUmi = umi; c.coordinate = coordinate;
    c.posterior = probability; c.memberCount = count; c.flags = flags;
    c.cliqueOrdinal = clique; c.candidateOrdinal = candidate;
    rows.push_back(c);
}

static void isolation() {
    std::vector<Contribution> rows;
    std::mt19937 random(20260915);
    for (unsigned i = 0; i < 100000; ++i) {
        unsigned gene = random() % 40, coord = random() % 128;
        unsigned umi = random() % 1024, count = 1 + random() % 15;
        unsigned candidates = 1 + random() % 4;
        for (unsigned j = 0; j < candidates; ++j) {
            unsigned flags = j == 0 ? downstream_spool::ContributionHard : 0;
            if (candidates == 1) flags |= downstream_spool::ContributionStrict
                | downstream_spool::ContributionGatedHard;
            const double p = candidates == 1 ? 1.0
                : (j == 0 ? 0.97 : 0.03 / (candidates - 1));
            if (j == 0 && candidates > 1) flags |= downstream_spool::ContributionGatedHard;
            add(rows, gene, umi, coord + j * 3350, p, count, flags, i, j);
        }
    }
    const auto before = rows;
    const auto all = resolveContributionShard(rows, ProductAll);
    for (const unsigned mask : {1u, 2u, 4u, 8u, 5u}) {
        const auto begin = std::chrono::steady_clock::now();
        const auto result = resolveContributionShard(rows, mask);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - begin).count();
        same(result, subset(all, mask));
        std::cout << "policy\t" << mask << "\tmolecules\t" << result.size()
                  << "\tseconds\t" << std::setprecision(17) << seconds << '\n';
    }
    require(rows.size() == before.size(), "resolver resized source contributions");
    for (size_t i = 0; i < rows.size(); ++i) {
        const auto &a = rows[i]; const auto &b = before[i];
        require(a.gene == b.gene && a.rawUmi == b.rawUmi
                && a.coordinate == b.coordinate && a.posterior == b.posterior
                && a.memberCount == b.memberCount && a.flags == b.flags
                && a.cliqueOrdinal == b.cliqueOrdinal
                && a.candidateOrdinal == b.candidateOrdinal,
                "resolver mutated source contributions");
    }
    std::cout << "policy_isolation\tPASS\tcontributions\t" << rows.size() << '\n';
}

static bool separated(unsigned a, unsigned b) {
    return a != b && !hammingOneUmi(a, b);
}

static void softOrderProbe() {
    // Three individually valid posterior values form a comparison cycle under
    // the production 'nearly equal, then UMI' ranking rule. Each clique also
    // has its complementary candidate, so posterior sums remain one.
    const double p[] = {0.5, 0.50000000004, 0.50000000008};
    auto less = [&](unsigned a, unsigned b) {
        if (!nearlyEqual(p[a], p[b])) return p[a] > p[b];
        return a < b;
    };
    const bool cycle = less(0, 1) && less(1, 2) && less(2, 0);
    std::cout << "soft_order_comparison_cycle\t" << (cycle ? "FOUND" : "absent") << '\n';
    std::vector<unsigned> decoys;
    for (unsigned u = 3; u < (1u << 18) && decoys.size() < 64; ++u) {
        bool independent = separated(u, 0) && separated(u, 1) && separated(u, 2);
        for (auto other : decoys) independent = independent && separated(u, other);
        if (independent) decoys.push_back(u);
    }
    require(decoys.size() == 64, "could not construct disconnected UMI controls");
    for (unsigned n = 0; n <= decoys.size(); ++n) {
        std::vector<Contribution> rows;
        for (unsigned u = 0; u < 3; ++u) {
            add(rows, 0, u, 0, p[u], 1, downstream_spool::ContributionHard, u, 0);
            add(rows, 0, u, 3350, 1 - p[u], 1, 0, u, 1);
        }
        add(rows, 1, 0, 0, 0.9, 1, downstream_spool::ContributionHard, 3, 0);
        add(rows, 1, 0, 3350, 0.1, 1, 0, 3, 1);
        for (unsigned i = 0; i < n; ++i) {
            add(rows, 0, decoys[i], 0, 0.1, 1, 0, 4 + i, 0);
            add(rows, 0, decoys[i], 3350, 0.9, 1,
                downstream_spool::ContributionHard, 4 + i, 1);
        }
        const auto result = resolveContributionShard(rows, ProductSoftExpected);
        std::vector<Clique> cliques;
        std::vector<CliqueCandidate> candidates;
        for (size_t i = 0; i < rows.size(); i += 2) {
            Clique clique;
            clique.gene = rows[i].gene; clique.rawUmi = rows[i].rawUmi;
            clique.memberCount = rows[i].memberCount;
            clique.candidateBegin = candidates.size(); clique.candidateCount = 2;
            cliques.push_back(clique);
            for (size_t j = i; j < i + 2; ++j) {
                CliqueCandidate candidate;
                candidate.coordinate = rows[j].coordinate;
                candidate.posterior = rows[j].posterior;
                candidates.push_back(candidate);
            }
        }
        auto memory = resolveSoft(cliques, candidates);
        std::sort(memory.begin(), memory.end(), [](const FinalMolecule &a, const FinalMolecule &b) {
            return std::tie(a.policy, a.coordinateIndex, a.correctedUmi, a.geneIndex)
                < std::tie(b.policy, b.coordinateIndex, b.correctedUmi, b.geneIndex);
        });
        same(result, memory);
        const auto raw = buildSoftRawSupport(cliques, candidates);
        const auto &original = findSoftSupport(raw, 0, 0, 0);
        double targetMass = 0;
        unsigned root = UINT32_MAX;
        for (const auto &m : result) {
            if (m.coordinateIndex == 0 && m.geneIndex == 0 && m.correctedUmi < 3) {
                targetMass += m.weight; root = m.correctedUmi;
            }
        }
        std::cout << "soft_order_probe\tdecoys\t" << n << "\ttarget_mass\t"
                  << std::setprecision(17) << targetMass << "\troot\t" << root
                  << "\tpre_reconcile_root\t" << original.correctedUmi << '\n';
    }
}

int main() {
    isolation();
    softOrderProbe();
}
