#ifndef STAR_WORKLOAD_DRAIN_ESTIMATE_H
#define STAR_WORKLOAD_DRAIN_ESTIMATE_H
#include <algorithm>
#include <cstdint>
#include <limits>

// Total and completed work are read pairs. Decoder blocks never enter this
// estimate. Exhausting an approximate count is not application completion.
class WorkloadDrainEstimate {
public:
    struct Result {
        uint64_t estimate = 0, delta = 0;
        double rate = 0, eta = std::numeric_limits<double>::infinity();
    };
    Result observe(uint64_t completed, uint64_t suppliedTotal, double seconds, bool complete) {
        Result r;
        r.delta = completed >= previous_ ? completed - previous_ : 0;
        previous_ = completed;
        const double instant = seconds > 0 ? r.delta / seconds : 0;
        rate_ = rate_ > 0 ? 0.5 * rate_ + 0.5 * instant : instant;
        estimate_ = std::max(estimate_, suppliedTotal);
        if (estimate_ && estimate_ <= completed && !complete) {
            const uint64_t bump = std::max<uint64_t>(1, completed / 10);
            const uint64_t max = std::numeric_limits<uint64_t>::max();
            estimate_ = completed > max - bump ? max : completed + bump;
        }
        r.estimate = estimate_;
        r.rate = rate_;
        if (complete) r.eta = 0;
        else if (estimate_ > completed && rate_ > 0) r.eta = (estimate_ - completed) / rate_;
        return r;
    }
    void allocationChanged() { rate_ = 0; }
private:
    uint64_t previous_ = 0, estimate_ = 0;
    double rate_ = 0;
};
#endif
