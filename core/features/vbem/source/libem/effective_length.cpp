#include "effective_length.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

bool EffectiveLengthCalculator::loadFLD(const std::string& fld_path) {
    std::ifstream in(fld_path);
    if (!in.is_open()) {
        std::cerr << "Error: Cannot open FLD file: " << fld_path << "\n";
        return false;
    }
    
    // Initialize FLD (max 2000)
    constexpr int MAX_FRAG_LEN = 2000;
    fld_.resize(MAX_FRAG_LEN, 0.0);
    
    std::string line;
    bool has_tab = false;
    
    // Check first line to detect format
    if (std::getline(in, line)) {
        if (line.find('\t') != std::string::npos) {
            has_tab = true;  // TSV format
        }
        in.seekg(0);  // Rewind
    }
    
    double sum = 0.0;
    
    if (has_tab) {
        // TSV format: length <tab> count <tab> probability
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            
            int len;
            unsigned long count;
            double prob;
            if (sscanf(line.c_str(), "%d %lu %lf", &len, &count, &prob) >= 2) {
                if (len >= 0 && len < MAX_FRAG_LEN) {
                    // Use probability if available, otherwise use count
                    if (prob > 0) {
                        fld_[len] = prob;
                    } else if (count > 0) {
                        fld_[len] = static_cast<double>(count);
                    }
                    sum += fld_[len];
                }
            }
        }
    } else {
        // Legacy format: one value per line
        for (int i = 0; i < MAX_FRAG_LEN; ++i) {
            if (std::getline(in, line)) {
                try {
                    fld_[i] = std::stod(line);
                    sum += fld_[i];
                } catch (...) {
                    fld_[i] = 0.0;
                }
            } else {
                break;
            }
        }
    }
    
    // Normalize if needed
    if (sum > 0 && abs(sum - 1.0) > 0.01) {
        for (size_t i = 0; i < fld_.size(); ++i) {
            fld_[i] /= sum;
        }
    }
    
    // Build CDF and compute quantile bounds
    buildCDF();
    computeQuantileBounds();
    
    return true;
}

void EffectiveLengthCalculator::setFLDPMF(const std::vector<double>& fld_pmf) {
    fld_ = fld_pmf;
    // Resize to MAX_FRAG_LEN if needed
    if (fld_.size() < static_cast<size_t>(2000 + 1)) {
        fld_.resize(2001, 0.0);
    }
    // Build CDF and compute quantile bounds
    buildCDF();
    computeQuantileBounds();
}

void EffectiveLengthCalculator::buildCDF() {
    fld_cdf_.resize(fld_.size());
    if (fld_.empty()) return;
    
    fld_cdf_[0] = fld_[0];
    for (size_t i = 1; i < fld_.size(); ++i) {
        fld_cdf_[i] = fld_cdf_[i-1] + fld_[i];
    }
    
    // Normalize CDF to [0, 1]
    double cdf_max = fld_cdf_.back();
    if (cdf_max > 0) {
        for (size_t i = 0; i < fld_cdf_.size(); ++i) {
            fld_cdf_[i] /= cdf_max;
        }
    }
}

void EffectiveLengthCalculator::computeQuantileBounds(double quantile_low, double quantile_high) {
    if (fld_cdf_.empty()) {
        fld_low_ = 0;
        fld_high_ = 1000;
        return;
    }
    
    fld_low_ = 0;
    fld_high_ = static_cast<int32_t>(fld_cdf_.size() - 1);
    
    for (size_t i = 0; i < fld_cdf_.size(); ++i) {
        if (fld_cdf_[i] >= quantile_low && fld_low_ == 0) {
            fld_low_ = static_cast<int32_t>(i);
        }
        if (fld_cdf_[i] >= quantile_high) {
            fld_high_ = static_cast<int32_t>(i);
            break;
        }
    }
}

void EffectiveLengthCalculator::loadGCBias(const std::vector<double>& bias_ratio) {
    gc_bias_ = bias_ratio;
    if (gc_bias_.size() != 101) {
        std::cerr << "Warning: GC bias vector size (" << gc_bias_.size() 
                  << ") != 101, resizing\n";
        gc_bias_.resize(101, 1.0);
    }
}

double EffectiveLengthCalculator::computeEffectiveLength(
    const libem::TranscriptSequence& txp, int32_t refLen) const 
{
    if (fld_.empty()) {
        // No FLD: fall back to raw length.
        return static_cast<double>(refLen);
    }

    std::vector<int32_t> raw_lengths{refLen};
    std::vector<double> fld_lengths = computeEffectiveLengthsFromPMF(fld_, raw_lengths);
    double fld_effective_length = fld_lengths.empty()
        ? static_cast<double>(refLen)
        : fld_lengths.front();
    return computeGCBiasedEffectiveLength(txp, refLen, fld_effective_length);
}

double EffectiveLengthCalculator::computeGCBiasedEffectiveLength(
    const libem::TranscriptSequence& txp,
    int32_t refLen,
    double fld_effective_length) const
{
    if (refLen <= 1 || fld_.empty() || fld_cdf_.empty() || gc_bias_.empty()) {
        return fld_effective_length;
    }

    constexpr int32_t biasSpeedSamp = 5; // Salmon default: --biasSpeedSamp 5
    constexpr double minCDFMass = 1e-10;

    const int32_t cdfMaxIdx = static_cast<int32_t>(fld_cdf_.size() - 1);
    const int32_t cdfMaxArg = std::min(cdfMaxIdx, refLen);
    const double cdfMaxVal = fld_cdf_[static_cast<size_t>(cdfMaxArg)];
    if (cdfMaxVal <= minCDFMass) {
        return fld_effective_length;
    }

    auto conditionalCDF = [&](int32_t x) -> double {
        if (x <= 0) {
            return fld_cdf_[0] / cdfMaxVal;
        }
        if (x > cdfMaxArg) {
            return 1.0;
        }
        return fld_cdf_[static_cast<size_t>(x)] / cdfMaxVal;
    };

    int32_t locFLDLow = (refLen < cdfMaxIdx) ? 1 : fld_low_;
    int32_t locFLDHigh = (refLen < cdfMaxIdx) ? cdfMaxArg : fld_high_;
    locFLDLow = std::max(1, locFLDLow);
    locFLDHigh = std::min(locFLDHigh, cdfMaxArg);

    int32_t fl = locFLDLow;
    const int32_t maxLen = std::min(refLen, locFLDHigh + 1);
    if (fl >= maxLen) {
        return fld_effective_length;
    }

    double effLength = 0.0;
    double prevFLMass = conditionalCDF(fl > 0 ? fl - 1 : 0);
    bool done = false;

    while (!done) {
        if (fl >= maxLen) {
            done = true;
            fl = maxLen - 1;
        }

        const double flMass = conditionalCDF(fl);
        const double flWeight = flMass - prevFLMass;
        prevFLMass = flMass;

        if (flWeight > 0.0) {
            double flMassTotal = 0.0;
            const int32_t lastStart = refLen - fl;
            for (int32_t fragStart = 0; fragStart < lastStart; ++fragStart) {
                const int32_t fragEnd = fragStart + fl - 1;
                int32_t gcFrac = txp.gcFrac(fragStart, fragEnd);
                if (gcFrac >= 0 && gcFrac < static_cast<int32_t>(gc_bias_.size())) {
                    flMassTotal += gc_bias_[static_cast<size_t>(gcFrac)];
                }
            }
            effLength += flWeight * flMassTotal;
        }

        fl += biasSpeedSamp;
    }

    // Match Salmon's bias-length barrier: if correction cannot be trusted,
    // retain the FLD-only effective length rather than over-correcting.
    const double unprocessedLen = std::max(0.0, static_cast<double>(refLen) - fld_effective_length);
    const double offset = std::max(1.0, unprocessedLen);
    return std::max(effLength, std::min(fld_effective_length, offset));
}

std::vector<double> EffectiveLengthCalculator::computeAllEffectiveLengths(
    const libem::Transcriptome& txome,
    const std::vector<double>& raw_lengths) const
{
    std::vector<int32_t> raw_lengths_int;
    raw_lengths_int.reserve(raw_lengths.size());
    for (double raw_len : raw_lengths) {
        raw_lengths_int.push_back(static_cast<int32_t>(raw_len));
    }

    std::vector<double> eff_lengths = computeEffectiveLengthsFromPMF(fld_, raw_lengths_int);
    if (eff_lengths.size() != raw_lengths.size()) {
        eff_lengths.assign(raw_lengths.begin(), raw_lengths.end());
    }

    #pragma omp parallel for schedule(dynamic, 64)
    for (ptrdiff_t idx = 0; idx < static_cast<ptrdiff_t>(raw_lengths.size()); ++idx) {
        const size_t i = static_cast<size_t>(idx);
        const libem::TranscriptSequence* txp = txome.getTranscript(static_cast<uint32_t>(i));
        if (!txp) {
            continue;
        }
        
        int32_t raw_len = static_cast<int32_t>(raw_lengths[i]);
        eff_lengths[i] = computeGCBiasedEffectiveLength(*txp, raw_len, eff_lengths[i]);
    }
    
    return eff_lengths;
}

std::vector<double> EffectiveLengthCalculator::computeEffectiveLengthsFromPMF(
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths) const
{
    std::vector<double> eff_lengths;
    eff_lengths.reserve(raw_lengths.size());

    if (fld_pmf.empty()) {
        for (int32_t raw_len : raw_lengths) {
            eff_lengths.push_back(static_cast<double>(raw_len));
        }
        return eff_lengths;
    }

    // Match Salmon's alignment-mode effective length correction:
    // build cumulative E[fragment length | fragment length <= transcript length],
    // then subtract that correction factor from each transcript's raw length.
    std::vector<double> correction_factors(fld_pmf.size(), 0.0);
    std::vector<double> cumulative_weighted_lengths(fld_pmf.size(), 0.0);
    std::vector<double> cumulative_masses(fld_pmf.size(), 0.0);

    cumulative_masses[0] = fld_pmf[0];
    for (size_t i = 1; i < fld_pmf.size(); ++i) {
        const double mass = fld_pmf[i];
        cumulative_weighted_lengths[i] =
            (mass * static_cast<double>(i)) + cumulative_weighted_lengths[i - 1];
        cumulative_masses[i] = mass + cumulative_masses[i - 1];
        if (cumulative_masses[i] > 0.0) {
            correction_factors[i] = cumulative_weighted_lengths[i] / cumulative_masses[i];
        }
    }

    for (int32_t raw_len : raw_lengths) {
        const double orig_len = static_cast<double>(raw_len);
        const size_t correction_idx =
            (static_cast<size_t>(raw_len) >= fld_pmf.size())
                ? fld_pmf.size() - 1
                : static_cast<size_t>(raw_len);

        double eff_len = orig_len - correction_factors[correction_idx];
        if (eff_len < 1.0) {
            eff_len = orig_len;
        }
        eff_lengths.push_back(eff_len);
    }
    
    return eff_lengths;
}
