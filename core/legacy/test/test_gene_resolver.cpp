#include "GeneResolver.h"

#include <cassert>
#include <iostream>
#include <vector>

namespace {

CandidateView candidate(bool genomic, uint16_t gene, int score)
{
    CandidateView out;
    out.isGenomic = genomic;
    out.geneIdx15 = gene;
    if (gene != 0) {
        out.zgGeneIdx15.push_back(gene);
    }
    out.asScore = score;
    out.nm = 0;
    return out;
}

} // namespace

int main()
{
    {
        std::vector<CandidateView> candidates{candidate(false, 11, 47)};
        assert(resolveGeneFromCandidates(candidates) == 11);
    }

    {
        std::vector<CandidateView> candidates{
            candidate(false, 11, 47), candidate(false, 12, 46)};
        assert(resolveGeneFromCandidates(candidates) == 11);
    }

    {
        std::vector<CandidateView> candidates{
            candidate(false, 11, 47), candidate(false, 12, 47)};
        assert(resolveGeneFromCandidates(candidates) == 0);
    }

    {
        std::vector<CandidateView> candidates{
            candidate(false, 11, 47), candidate(true, 12, 47)};
        assert(resolveGeneFromCandidates(candidates) == 0);
    }

    {
        std::vector<CandidateView> candidates{
            candidate(false, 11, 48), candidate(true, 12, 47)};
        assert(resolveGeneFromCandidates(candidates) == 11);
    }

    {
        std::vector<CandidateView> candidates{
            candidate(false, 11, 47), candidate(true, 12, 48)};
        assert(resolveGeneFromCandidates(candidates) == 12);
    }

    {
        std::vector<CandidateView> candidates{
            candidate(false, 11, 47), candidate(true, 11, 47)};
        assert(resolveGeneFromCandidates(candidates) == 11);
    }

    {
        std::vector<CandidateView> candidates{
            candidate(true, 11, 47), candidate(true, 12, 46)};
        assert(resolveGeneFromCandidates(candidates) == 11);
    }

    {
        CandidateView multigene = candidate(true, 11, 47);
        multigene.zgGeneIdx15.push_back(12);
        std::vector<CandidateView> candidates{multigene};
        assert(resolveGeneFromCandidates(candidates) == 0);
    }

    std::cout << "gene resolver tests passed\n";
    return 0;
}
