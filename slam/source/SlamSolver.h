#ifndef SLAM_SOLVER_H
#define SLAM_SOLVER_H

#include <cstdint>
#include <map>

struct SlamResult {
    double ntr = 0.0;
    double sigma = 0.0;
    double log_likelihood = 0.0;
    bool converged = false;
};

// Key: (Total_Ts << 16) | Mismatch_Count
// Value: Frequency (weighted count of reads with this pattern)
using MismatchHistogramKey = uint32_t;
using MismatchHistogram = std::map<MismatchHistogramKey, double>;

inline MismatchHistogramKey slamPackMismatchKey(uint16_t nT, uint16_t tc) {
    return (static_cast<MismatchHistogramKey>(nT) << 16) |
           static_cast<MismatchHistogramKey>(tc);
}

inline uint16_t slamKeyNT(MismatchHistogramKey key) {
    return static_cast<uint16_t>(key >> 16);
}

inline uint16_t slamKeyTC(MismatchHistogramKey key) {
    return static_cast<uint16_t>(key & 0xFFFFu);
}

class SlamSolver {
public:
    SlamSolver(double error_rate = 0.001, double conversion_rate = 0.05)
        : p_error_rate_(error_rate), p_conversion_rate_(conversion_rate) {}

    SlamResult solve(const MismatchHistogram& gene_data) const;

private:
    double p_error_rate_;
    double p_conversion_rate_;

    double calc_log_likelihood(const MismatchHistogram& data, double pi) const;
    double log_binom_pmf(uint16_t n, uint16_t tc, double p) const;
};

#endif
