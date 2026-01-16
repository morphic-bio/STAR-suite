#ifndef FORGETTING_MASS_H
#define FORGETTING_MASS_H

#include <cmath>
#include <cstdint>
#include <mutex>
#include <vector>

// Simplified ForgettingMassCalculator ported from Salmon.
class ForgettingMassCalculator {
public:
    explicit ForgettingMassCalculator(double forgetting_factor = 0.65)
        : batch_num_(0),
          forgetting_factor_(forgetting_factor),
          log_forgetting_mass_(0.0) {}

    // Compute log forgetting mass for the next minibatch and increment timestep.
    // Matches Salmon's operator() behavior:
    //   - First increment batchNum_
    //   - If batchNum_ > 1, apply forgetting mass formula
    //   - Return the current log forgetting mass
    void getLogMassAndTimestep(double& log_forgetting_mass,
                               uint64_t& current_minibatch_timestep) {
        std::lock_guard<std::mutex> lock(mutex_);
        
        // Salmon increments batchNum_ FIRST
        ++batch_num_;
        
        // If we've already computed this timestep, return cached value
        if (batch_num_ <= log_forgetting_masses_.size()) {
            current_minibatch_timestep = batch_num_ - 1;
            log_forgetting_mass = log_forgetting_masses_[batch_num_ - 1];
            return;
        }
        
        // Compute new forgetting mass
        // Salmon: if (batchNum_ > 1) { apply formula }
        double fm = log_forgetting_mass_;  // Start with current accumulator value
        if (batch_num_ > 1) {
            // Salmon formula: forgetting_factor * log(t-1) - log(t^forgetting_factor - 1)
            double t = static_cast<double>(batch_num_);
            double denom = std::pow(t, forgetting_factor_) - 1.0;
            if (denom > 0) {
                fm += forgetting_factor_ * std::log(t - 1.0) - std::log(denom);
            }
        }
        
        log_forgetting_mass_ = fm;  // Update accumulator
        log_forgetting_masses_.push_back(fm);
        
        // Update cumulative
        if (cumulative_log_forgetting_masses_.empty()) {
            cumulative_log_forgetting_masses_.push_back(fm);
        } else {
            cumulative_log_forgetting_masses_.push_back(
                logAdd(cumulative_log_forgetting_masses_.back(), fm));
        }

        current_minibatch_timestep = batch_num_ - 1;
        log_forgetting_mass = fm;
    }

private:
    static double logAdd(double a, double b) {
        if (std::isinf(a)) return b;
        if (std::isinf(b)) return a;
        if (a > b) {
            return a + std::log(1.0 + std::exp(b - a));
        }
        return b + std::log(1.0 + std::exp(a - b));
    }

    uint64_t batch_num_;
    double forgetting_factor_;
    double log_forgetting_mass_;
    std::vector<double> log_forgetting_masses_;
    std::vector<double> cumulative_log_forgetting_masses_;
    std::mutex mutex_;
};

#endif // FORGETTING_MASS_H
