#include "effective_length_wrapper.h"
#include "effective_length.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace {

constexpr int kSalmonFragGCBins = 25;
constexpr int kBiasSpeedSamp = 5;
constexpr double kMinAlpha = 1e-8;
constexpr double kMinCDFMass = 1e-10;
constexpr double kGCPrior = 0.1;
constexpr double kMaxGCBiasRatio = 1000.0;

int salmonGCBin(int32_t gcPct) {
    if (gcPct < 0) {
        gcPct = 0;
    } else if (gcPct > 100) {
        gcPct = 100;
    }
    const double width = 100.0 / static_cast<double>(kSalmonFragGCBins);
    int bin = static_cast<int>(static_cast<double>(gcPct) / width);
    if (bin >= kSalmonFragGCBins) {
        bin = kSalmonFragGCBins - 1;
    }
    return bin;
}

std::vector<double> normalizeWithPrior(const std::vector<double>& counts) {
    std::vector<double> normalized(counts.size(), 0.0);
    double total = 0.0;
    for (double count : counts) {
        total += count + kGCPrior;
    }
    if (total <= 0.0 || !std::isfinite(total)) {
        const double uniform = counts.empty() ? 0.0 : 1.0 / static_cast<double>(counts.size());
        std::fill(normalized.begin(), normalized.end(), uniform);
        return normalized;
    }
    for (size_t i = 0; i < counts.size(); ++i) {
        normalized[i] = (counts[i] + kGCPrior) / total;
    }
    return normalized;
}

std::vector<double> computeCDF(const std::vector<double>& pmf) {
    std::vector<double> cdf(pmf.size(), 0.0);
    if (pmf.empty()) {
        return cdf;
    }
    cdf[0] = pmf[0];
    for (size_t i = 1; i < pmf.size(); ++i) {
        cdf[i] = cdf[i - 1] + pmf[i];
    }
    return cdf;
}

void computeQuantileBounds(const std::vector<double>& cdf, int32_t& fldLow, int32_t& fldHigh) {
    fldLow = 0;
    fldHigh = cdf.empty() ? 0 : static_cast<int32_t>(cdf.size() - 1);
    bool lb = false;
    for (size_t i = 0; i < cdf.size(); ++i) {
        if (!lb && cdf[i] >= 0.005) {
            fldLow = static_cast<int32_t>(i);
            lb = true;
        }
        if (cdf[i] >= 0.995) {
            fldHigh = static_cast<int32_t>(i);
            break;
        }
    }
}

} // namespace

std::vector<double> computeEffectiveLengthsFromPMFWrapper(
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths)
{
    EffectiveLengthCalculator calc;
    calc.setFLDPMF(fld_pmf);
    return calc.computeEffectiveLengthsFromPMF(fld_pmf, raw_lengths);
}

std::vector<double> computeGCBiasedEffectiveLengthsWrapper(
    const libem::Transcriptome& transcriptome,
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths,
    const std::vector<double>& gc_bias)
{
    EffectiveLengthCalculator calc;
    calc.setFLDPMF(fld_pmf);
    calc.loadGCBias(gc_bias);

    std::vector<double> raw_lengths_double;
    raw_lengths_double.reserve(raw_lengths.size());
    for (int32_t len : raw_lengths) {
        raw_lengths_double.push_back(static_cast<double>(len));
    }
    return calc.computeAllEffectiveLengths(transcriptome, raw_lengths_double);
}

DynamicGCEffectiveLengthResult computeDynamicGCBiasedEffectiveLengthsWrapper(
    const libem::Transcriptome& transcriptome,
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths,
    const std::vector<double>& alpha_counts,
    const std::vector<double>& effective_lengths_in,
    const std::vector<double>& observed_gc_101)
{
    DynamicGCEffectiveLengthResult result;
    result.effective_lengths = effective_lengths_in;
    result.observed_gc.assign(kSalmonFragGCBins, 0.0);
    result.expected_gc.assign(kSalmonFragGCBins, 0.0);
    result.gc_bias.assign(kSalmonFragGCBins, 1.0);

    if (fld_pmf.empty() ||
        raw_lengths.size() != effective_lengths_in.size() ||
        alpha_counts.size() != effective_lengths_in.size() ||
        transcriptome.size() < raw_lengths.size()) {
        return result;
    }

    std::vector<double> observedCounts(kSalmonFragGCBins, 0.0);
    for (size_t pct = 0; pct < observed_gc_101.size(); ++pct) {
        observedCounts[salmonGCBin(static_cast<int32_t>(pct))] += observed_gc_101[pct];
    }
    result.observed_gc = normalizeWithPrior(observedCounts);

    const std::vector<double> cdf = computeCDF(fld_pmf);
    if (cdf.empty()) {
        return result;
    }
    int32_t fldLow = 0;
    int32_t fldHigh = 0;
    computeQuantileBounds(cdf, fldLow, fldHigh);

    std::vector<double> expectedCounts(kSalmonFragGCBins, 0.0);
    const int32_t cdfMaxIdx = static_cast<int32_t>(cdf.size() - 1);

    for (size_t txpIdx = 0; txpIdx < raw_lengths.size(); ++txpIdx) {
        const libem::TranscriptSequence* txp =
            transcriptome.getTranscript(static_cast<uint32_t>(txpIdx));
        if (txp == nullptr) {
            continue;
        }

        const int32_t refLen = raw_lengths[txpIdx];
        const int32_t elen = static_cast<int32_t>(effective_lengths_in[txpIdx]);
        const int32_t unprocessedLen = std::max(0, refLen - elen);
        const int32_t cdfMaxArg = std::min(cdfMaxIdx, refLen);
        if (cdfMaxArg < 0) {
            continue;
        }
        const double cdfMaxVal = cdf[static_cast<size_t>(cdfMaxArg)];
        if (cdfMaxVal < kMinCDFMass ||
            alpha_counts[txpIdx] < kMinAlpha ||
            unprocessedLen <= 0 ||
            effective_lengths_in[txpIdx] <= 0.0) {
            continue;
        }

        ++result.background_transcripts;
        const double weight = alpha_counts[txpIdx] / effective_lengths_in[txpIdx];
        auto conditionalCDF = [&](int32_t x) -> double {
            if (x > cdfMaxArg) {
                return 1.0;
            }
            if (x < 0) {
                x = 0;
            }
            return cdf[static_cast<size_t>(x)] / cdfMaxVal;
        };

        const int32_t locFLDLow = (refLen < cdfMaxArg) ? 1 : fldLow;
        const int32_t locFLDHigh = (refLen < cdfMaxArg) ? cdfMaxArg : fldHigh;

        for (int32_t fragStart = 0; fragStart < refLen - 1; ++fragStart) {
            size_t sp = static_cast<size_t>((locFLDLow > 0) ? locFLDLow - 1 : 0);
            double prevFLMass = conditionalCDF(static_cast<int32_t>(sp));
            for (int32_t fl = locFLDLow; fl <= locFLDHigh; fl += kBiasSpeedSamp) {
                const int32_t fragEnd = fragStart + fl - 1;
                if (fragEnd < refLen) {
                    const double flMass = conditionalCDF(fl);
                    const double flWeight = flMass - prevFLMass;
                    prevFLMass = flMass;
                    if (flWeight > 0.0) {
                        const int32_t gcFrac = txp->gcFrac(fragStart, fragEnd);
                        expectedCounts[salmonGCBin(gcFrac)] += weight * flWeight;
                    }
                } else {
                    break;
                }
            }
        }
    }

    result.expected_gc = normalizeWithPrior(expectedCounts);
    bool haveExpected = false;
    for (double v : expectedCounts) {
        if (v > 0.0) {
            haveExpected = true;
            break;
        }
    }
    if (!haveExpected) {
        return result;
    }

    const double minRatio = 1.0 / kMaxGCBiasRatio;
    for (size_t i = 0; i < result.gc_bias.size(); ++i) {
        double ratio = 1.0;
        if (result.expected_gc[i] > 0.0) {
            ratio = result.observed_gc[i] / result.expected_gc[i];
            ratio = std::max(minRatio, std::min(kMaxGCBiasRatio, ratio));
        }
        result.gc_bias[i] = ratio;
    }

    #pragma omp parallel for schedule(dynamic, 64)
    for (ptrdiff_t idx = 0; idx < static_cast<ptrdiff_t>(raw_lengths.size()); ++idx) {
        const size_t txpIdx = static_cast<size_t>(idx);
        const libem::TranscriptSequence* txp =
            transcriptome.getTranscript(static_cast<uint32_t>(txpIdx));
        if (txp == nullptr) {
            continue;
        }

        const int32_t refLen = raw_lengths[txpIdx];
        const int32_t elen = static_cast<int32_t>(effective_lengths_in[txpIdx]);
        const int32_t unprocessedLen = std::max(0, refLen - elen);
        const int32_t cdfMaxArg = std::min(cdfMaxIdx, refLen);
        if (cdfMaxArg < 0) {
            result.effective_lengths[txpIdx] = static_cast<double>(elen);
            continue;
        }
        const double cdfMaxVal = cdf[static_cast<size_t>(cdfMaxArg)];

        double effLength = 0.0;
        bool wasProcessed = false;
        if (alpha_counts[txpIdx] >= kMinAlpha &&
            unprocessedLen > 0 &&
            cdfMaxVal > kMinCDFMass) {
            auto conditionalCDF = [&](int32_t x) -> double {
                if (x > cdfMaxArg) {
                    return 1.0;
                }
                if (x < 0) {
                    x = 0;
                }
                return cdf[static_cast<size_t>(x)] / cdfMaxVal;
            };

            const int32_t locFLDLow = (refLen < cdfMaxArg) ? 1 : fldLow;
            const int32_t locFLDHigh = (refLen < cdfMaxArg) ? cdfMaxArg : fldHigh;
            int32_t fl = locFLDLow;
            const int32_t maxLen = std::min(refLen, locFLDHigh + 1);
            bool done = (fl >= maxLen);
            size_t sp = static_cast<size_t>((fl > 0) ? fl - 1 : 0);
            double prevFLMass = conditionalCDF(static_cast<int32_t>(sp));

            while (!done) {
                if (fl >= maxLen) {
                    done = true;
                    fl = maxLen - 1;
                }
                const double flWeight = conditionalCDF(fl) - prevFLMass;
                prevFLMass = conditionalCDF(fl);

                double flMassTotal = 0.0;
                for (int32_t fragStart = 0; fragStart < refLen - fl; ++fragStart) {
                    const int32_t fragEnd = fragStart + fl - 1;
                    if (fragStart < refLen && fragEnd < refLen) {
                        const int32_t gcFrac = txp->gcFrac(fragStart, fragEnd);
                        flMassTotal += result.gc_bias[salmonGCBin(gcFrac)];
                    } else {
                        break;
                    }
                }

                effLength += flWeight * flMassTotal;
                fl += kBiasSpeedSamp;
            }
            wasProcessed = true;
        }

        if (wasProcessed) {
            const double offset = std::max(1.0, static_cast<double>(unprocessedLen));
            const double effLengthNoBias = static_cast<double>(elen);
            result.effective_lengths[txpIdx] =
                std::max(effLength, std::min(effLengthNoBias, offset));
        } else {
            result.effective_lengths[txpIdx] = static_cast<double>(elen);
        }
    }

    result.applied = true;
    return result;
}
