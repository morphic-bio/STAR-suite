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
    if (gc_bias_.empty()) {
        // No GC bias: effective length = raw length
        return static_cast<double>(refLen);
    }
    
    if (fld_.empty()) {
        // No FLD: use uniform distribution
        return static_cast<double>(refLen);
    }
    
    double effLength = 0.0;
    
    // Iterate over all valid fragment start positions
    for (int32_t fragStart = 0; fragStart < refLen; ++fragStart) {
        double flMassTotal = 0.0;
        
        // Iterate over fragment lengths (within quantile bounds)
        for (int32_t fl = fld_low_; fl <= fld_high_ && fl < static_cast<int32_t>(fld_.size()); ++fl) {
            int32_t fragEnd = fragStart + fl - 1;
            if (fragEnd >= refLen) break;
            
            // Base fragment factor (FLD weight)
            double fragFactor = fld_[fl];
            if (fragFactor <= 0) continue;
            
            // Apply GC bias
            int32_t gcFrac = txp.gcFrac(fragStart, fragEnd);
            if (gcFrac >= 0 && gcFrac < static_cast<int32_t>(gc_bias_.size())) {
                fragFactor *= gc_bias_[gcFrac];
            }
            
            flMassTotal += fragFactor;
        }
        
        effLength += flMassTotal;
    }
    
    // Normalize by number of start positions
    if (refLen > 0) {
        effLength /= refLen;
    }
    
    return effLength;
}

std::vector<double> EffectiveLengthCalculator::computeAllEffectiveLengths(
    const libem::Transcriptome& txome,
    const std::vector<double>& raw_lengths) const
{
    std::vector<double> eff_lengths;
    eff_lengths.reserve(raw_lengths.size());
    
    for (size_t i = 0; i < raw_lengths.size(); ++i) {
        const libem::TranscriptSequence* txp = txome.getTranscript(static_cast<uint32_t>(i));
        if (!txp) {
            eff_lengths.push_back(raw_lengths[i]);  // Fall back to raw length
            continue;
        }
        
        int32_t raw_len = static_cast<int32_t>(raw_lengths[i]);
        double eff_len = computeEffectiveLength(*txp, raw_len);
        eff_lengths.push_back(eff_len);
    }
    
    return eff_lengths;
}

std::vector<double> EffectiveLengthCalculator::computeEffectiveLengthsFromPMF(
    const std::vector<double>& fld_pmf,
    const std::vector<int32_t>& raw_lengths) const
{
    std::vector<double> eff_lengths;
    eff_lengths.reserve(raw_lengths.size());
    
    for (int32_t raw_len : raw_lengths) {
        if (fld_pmf.empty()) {
            eff_lengths.push_back(static_cast<double>(raw_len));
            continue;
        }
        
        // Compute effective length using FLD PMF (Salmon formula)
        // effLen = sum_{fl} PMF[fl] * max(0, L - fl + 1)
        // where (L - fl + 1) is the number of valid start positions for fragment length fl
        double effLen = 0.0;
        int32_t refLen = raw_len;
        
        // Iterate over fragment lengths (within quantile bounds if set, otherwise all)
        int32_t fl_low = (fld_low_ > 0) ? fld_low_ : 1;  // Start at 1 (min fragment length)
        int32_t fl_high = (fld_high_ > 0 && fld_high_ < static_cast<int32_t>(fld_pmf.size())) 
                        ? fld_high_ 
                        : static_cast<int32_t>(fld_pmf.size() - 1);
        
        for (int32_t fl = fl_low; fl <= fl_high && fl < static_cast<int32_t>(fld_pmf.size()); ++fl) {
            // Number of valid start positions for this fragment length
            int32_t numStartPositions = refLen - fl + 1;
            if (numStartPositions > 0) {
                effLen += fld_pmf[fl] * static_cast<double>(numStartPositions);
            }
        }
        
        // Clamp to valid range
        if (effLen < 1.0) effLen = 1.0;
        if (effLen > raw_len) effLen = static_cast<double>(raw_len);
        
        eff_lengths.push_back(effLen);
    }
    
    return eff_lengths;
}
