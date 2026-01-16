#include "gc_bias.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

// Log-space addition helper
static inline double logAdd(double a, double b) {
    if (a == -std::numeric_limits<double>::infinity()) return b;
    if (b == -std::numeric_limits<double>::infinity()) return a;
    if (a > b) {
        return a + std::log(1.0 + std::exp(b - a));
    } else {
        return b + std::log(1.0 + std::exp(a - b));
    }
}

void GCFragModel::inc(int32_t gc_pct, double weight) {
    // weight is in log space
    int32_t bin = clampGC(gc_pct);
    if (counts_[bin] == -std::numeric_limits<double>::infinity()) {
        counts_[bin] = weight;
    } else {
        counts_[bin] = logAdd(counts_[bin], weight);
    }
    normalized_ = false;
}

void GCFragModel::incLinear(int32_t gc_pct, double weight) {
    // weight is in linear space
    int32_t bin = clampGC(gc_pct);
    counts_[bin] += weight;
    normalized_ = false;
}

void GCFragModel::normalize() {
    if (normalized_) return;
    
    // Sum all counts (convert from log space if needed)
    double total = 0.0;
    bool in_log_space = false;
    
    // Check if we're in log space (any negative values that are very negative)
    for (int i = 0; i < GC_BINS; ++i) {
        if (counts_[i] < -100) {
            in_log_space = true;
            break;
        }
    }
    
    if (in_log_space) {
        // Log space: use log-sum-exp
        double logTotal = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < GC_BINS; ++i) {
            logTotal = logAdd(logTotal, counts_[i]);
        }
        
        // Convert to linear space
        for (int i = 0; i < GC_BINS; ++i) {
            counts_[i] = std::exp(counts_[i] - logTotal);
        }
    } else {
        // Linear space: simple normalization
        for (int i = 0; i < GC_BINS; ++i) {
            total += counts_[i];
        }
        
        if (total > 0.0) {
            for (int i = 0; i < GC_BINS; ++i) {
                counts_[i] /= total;
            }
        }
    }
    
    normalized_ = true;
}

double GCFragModel::get(int32_t gc_pct) const {
    int32_t bin = clampGC(gc_pct);
    return counts_[bin];
}

std::vector<double> GCFragModel::computeBiasRatio(
    const GCFragModel& expected, double maxRatio) const 
{
    std::vector<double> bias(GC_BINS, 1.0);
    double minRatio = 1.0 / maxRatio;
    
    // Ensure both models are normalized
    const_cast<GCFragModel*>(this)->normalize();
    const_cast<GCFragModel*>(&expected)->normalize();
    
    for (int i = 0; i < GC_BINS; ++i) {
        double obs = counts_[i];
        double exp = expected.counts_[i];
        
        if (exp > 0.0) {
            double ratio = obs / exp;
            if (ratio > maxRatio) ratio = maxRatio;
            if (ratio < minRatio) ratio = minRatio;
            bias[i] = ratio;
        } else {
            // If expected is zero, use default ratio of 1.0
            bias[i] = 1.0;
        }
    }
    
    return bias;
}

void GCFragModel::combineCounts(const GCFragModel& other) {
    // Check if we're in log space
    bool in_log_space = false;
    for (int i = 0; i < GC_BINS; ++i) {
        if (counts_[i] < -100 || other.counts_[i] < -100) {
            in_log_space = true;
            break;
        }
    }
    
    if (in_log_space) {
        // Log space: use log-add
        for (int i = 0; i < GC_BINS; ++i) {
            counts_[i] = logAdd(counts_[i], other.counts_[i]);
        }
    } else {
        // Linear space: simple addition
        for (int i = 0; i < GC_BINS; ++i) {
            counts_[i] += other.counts_[i];
        }
    }
    
    normalized_ = false;
}

void GCFragModel::reset() {
    counts_.fill(0.0);
    normalized_ = false;
}

bool GCFragModel::writeToFile(const std::string& path) const {
    std::ofstream out(path);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open file for writing: " << path << "\n";
        return false;
    }
    
    // Ensure normalized
    const_cast<GCFragModel*>(this)->normalize();
    
    for (int i = 0; i < GC_BINS; ++i) {
        out << i << "\t" << counts_[i] << "\n";
    }
    
    return true;
}

bool GCFragModel::loadFromFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "Error: Cannot open file for reading: " << path << "\n";
        return false;
    }
    
    reset();
    
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        int gc_bin;
        double prob;
        if (sscanf(line.c_str(), "%d %lf", &gc_bin, &prob) == 2) {
            if (gc_bin >= 0 && gc_bin < GC_BINS) {
                counts_[gc_bin] = prob;
            }
        }
    }
    
    normalized_ = true;  // Assume loaded file is already normalized
    return true;
}

