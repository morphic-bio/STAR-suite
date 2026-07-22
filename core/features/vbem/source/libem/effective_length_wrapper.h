#ifndef EFFECTIVE_LENGTH_WRAPPER_H
#define EFFECTIVE_LENGTH_WRAPPER_H

#include <vector>
#include <cstdint>
#include <cstddef>

namespace libem {
class Transcriptome;
}

struct DynamicGCEffectiveLengthResult {
    std::vector<double> effective_lengths;
    std::vector<double> observed_gc;
    std::vector<double> expected_gc;
    std::vector<double> gc_bias;
    size_t background_transcripts = 0;
    bool applied = false;
};

// Wrapper to avoid Transcriptome class name conflict in STAR.cpp
// Forward declaration approach: compute effective lengths from FLD PMF without requiring Transcriptome
std::vector<double> computeEffectiveLengthsFromPMFWrapper(
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths
);

std::vector<double> computeGCBiasedEffectiveLengthsWrapper(
    const libem::Transcriptome& transcriptome,
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths,
    const std::vector<double>& gc_bias
);

DynamicGCEffectiveLengthResult computeDynamicGCBiasedEffectiveLengthsWrapper(
    const libem::Transcriptome& transcriptome,
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths,
    const std::vector<double>& alpha_counts,
    const std::vector<double>& effective_lengths_in,
    const std::vector<double>& observed_gc_101
);

#endif // EFFECTIVE_LENGTH_WRAPPER_H
