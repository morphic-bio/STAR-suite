#ifndef EFFECTIVE_LENGTH_WRAPPER_H
#define EFFECTIVE_LENGTH_WRAPPER_H

#include <vector>
#include <cstdint>

// Wrapper to avoid Transcriptome class name conflict in STAR.cpp
// Forward declaration approach: compute effective lengths from FLD PMF without requiring Transcriptome
std::vector<double> computeEffectiveLengthsFromPMFWrapper(
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths
);

#endif // EFFECTIVE_LENGTH_WRAPPER_H
