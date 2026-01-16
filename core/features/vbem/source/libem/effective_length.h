#ifndef EFFECTIVE_LENGTH_H
#define EFFECTIVE_LENGTH_H

#include "alignment_model.h"
#include <vector>
#include <string>

// Effective Length Calculator
// Computes GC-corrected effective lengths matching Salmon's algorithm
class EffectiveLengthCalculator {
public:
    EffectiveLengthCalculator() 
        : fld_low_(0), fld_high_(1000) {}
    
    // Load fragment length distribution from TSV file
    // Format: length <tab> count <tab> probability
    bool loadFLD(const std::string& fld_path);
    
    // Set FLD PMF directly (for use with observed FLD)
    void setFLDPMF(const std::vector<double>& fld_pmf);
    
    // Load GC bias ratios (from computeBiasRatio)
    void loadGCBias(const std::vector<double>& bias_ratio);
    
    // Compute GC-corrected effective length for a transcript
    double computeEffectiveLength(
        const libem::TranscriptSequence& txp,
        int32_t raw_length
    ) const;
    
    // Compute for all transcripts
    std::vector<double> computeAllEffectiveLengths(
        const libem::Transcriptome& txome,
        const std::vector<double>& raw_lengths
    ) const;
    
    // Compute effective lengths using FLD PMF only (no GC bias, no TranscriptSequence required)
    // More efficient than computeEffectiveLength when GC bias is disabled
    std::vector<double> computeEffectiveLengthsFromPMF(
        const std::vector<double>& fld_pmf,
        const std::vector<int32_t>& raw_lengths
    ) const;
    
    // Set quantile bounds (default: 0.5% and 99.5%)
    void setQuantileBounds(int32_t low, int32_t high) {
        fld_low_ = low;
        fld_high_ = high;
    }
    
private:
    std::vector<double> fld_;           // Fragment length distribution (PDF)
    std::vector<double> fld_cdf_;      // CDF for quantile bounds
    std::vector<double> gc_bias_;      // Bias ratio per GC bin (101 bins)
    int32_t fld_low_;                  // Lower quantile bound
    int32_t fld_high_;                 // Upper quantile bound
    
    // Compute quantile bounds from FLD
    void computeQuantileBounds(double quantile_low = 0.005, double quantile_high = 0.995);
    
    // Build CDF from FLD
    void buildCDF();
};

#endif // EFFECTIVE_LENGTH_H
