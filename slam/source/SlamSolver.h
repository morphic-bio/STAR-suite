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

// Key: (Total_Ts << 8) | Mismatch_Count
// Value: Frequency (weighted count of reads with this pattern)
using MismatchHistogram = std::map<uint16_t, double>;

class SlamSolver {
public:
    SlamSolver(double error_rate = 0.001, double conversion_rate = 0.05)
        : p_error_rate_(error_rate), p_conversion_rate_(conversion_rate) {}

    SlamResult solve(const MismatchHistogram& gene_data) const;

private:
    double p_error_rate_;
    double p_conversion_rate_;

    double calc_log_likelihood(const MismatchHistogram& data, double pi) const;
    double log_binom_pmf(uint16_t n, uint8_t k, double p) const;
};

#endif
