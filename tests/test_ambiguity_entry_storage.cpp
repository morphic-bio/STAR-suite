#include "SoloReadFeature.h"
#include "solo/CbBayesianResolver.h"
#include <iomanip>
#include <iostream>
#include <utility>

int main() {
    using Entry = SoloReadFeature::ExtendedAmbiguousEntry;
    std::cerr << "entry_bytes=" << sizeof(Entry) << '\n';
    Entry a;
    a.cbSeq = "ACGTACGTACGTACGN";
    a.cbQual = "I!H"; // exercises padding
    a.candidateIdx = {1, 7, 10};
    a.candidateQual = {33, 40, 55};
    a.umiCounts[10] = 7;
    for (unsigned n = 0; n < 80; ++n)
        cb_bayesian::accumulateCbQualityEvidence(a.cbSeq, a.cbQual,
            a.cbLogLikMatch, a.cbLogLikMismatch, a.cbEvidenceReads);
    Entry b = a;
    Entry c = std::move(b);
    cb_bayesian::mergeCbQualityEvidence(a.cbLogLikMatch, a.cbLogLikMismatch,
        a.cbEvidenceReads, c.cbLogLikMatch, c.cbLogLikMismatch, c.cbEvidenceReads);
    std::cout << c.cbSeq << ' ' << c.cbQual << ' ' << c.cbEvidenceReads << '\n';
    for (auto v : c.candidateIdx) std::cout << v << ' ';
    std::cout << '\n' << c.umiCounts.at(10) << '\n' << std::hexfloat;
    for (size_t p = 0; p < c.cbLogLikMatch.size(); ++p)
        std::cout << c.cbLogLikMatch[p] << ' ' << c.cbLogLikMismatch[p] << '\n';
}
