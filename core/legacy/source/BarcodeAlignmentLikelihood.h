#ifndef STAR_BARCODE_ALIGNMENT_LIKELIHOOD_H
#define STAR_BARCODE_ALIGNMENT_LIKELIHOOD_H
#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>
#include <vector>

namespace barcode_quality {
inline double errorProbability(int phred) {
    const double probability =
        std::pow(10.0, -static_cast<double>(std::max(0, phred)) / 10.0);
    return std::max(1.0e-10, std::min(0.75, probability));
}

inline unsigned char baseCode(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

inline double alignmentLogLikelihood(
    const std::string& observed, const std::string& qualities,
    const std::string& candidate, int gap_q = 30, int missing_q = 30) {
    std::vector<double> previous(candidate.size() + 1, 0.0);
    std::vector<double> current(candidate.size() + 1, 0.0);
    std::vector<unsigned char> candidate_bases(candidate.size());
    for (size_t j = 0; j < candidate.size(); ++j) {
        candidate_bases[j] = baseCode(candidate[j]);
    }
    const double deletion = std::log(errorProbability(gap_q));
    const double ambiguous = std::log(0.25);
    for (size_t j = 1; j <= candidate.size(); ++j) previous[j] = previous[j - 1] + deletion;
    for (size_t i = 1; i <= observed.size(); ++i) {
        const int q = i - 1 < qualities.size()
            ? std::max(0, static_cast<int>(static_cast<unsigned char>(qualities[i - 1])) - 33)
            : missing_q;
        // These terms depend only on this observed base and its quality, not
        // the candidate column. Keep the original operations and DP order so
        // likelihoods (including ties) remain bit-identical.
        const double error = errorProbability(q);
        const double insertion = std::log(error / 4.0);
        const double match = std::log1p(-error);
        const double mismatch = std::log(error / 3.0);
        const unsigned char obs = baseCode(observed[i - 1]);
        current[0] = previous[0] + insertion;
        for (size_t j = 1; j <= candidate.size(); ++j) {
            const unsigned char truth = candidate_bases[j - 1];
            const double pair = obs == 4 || truth == 4
                ? ambiguous : (obs == truth ? match : mismatch);
            current[j] = std::max({
                previous[j - 1] + pair,
                previous[j] + insertion,
                current[j - 1] + deletion,
            });
        }
        previous.swap(current);
    }
    return previous[candidate.size()];
}

} // namespace barcode_quality
#endif
