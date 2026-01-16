#include "effective_length_wrapper.h"
#include "effective_length.h"

std::vector<double> computeEffectiveLengthsFromPMFWrapper(
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths)
{
    EffectiveLengthCalculator calc;
    calc.setFLDPMF(fld_pmf);
    return calc.computeEffectiveLengthsFromPMF(fld_pmf, raw_lengths);
}
