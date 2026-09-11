#ifndef STAR_BGZF_RATE_CONTROLLER_H
#define STAR_BGZF_RATE_CONTROLLER_H

#include <algorithm>
#include <cstdint>

// Inner policy only. Its objective is completed parent reads/second at a given
// parent budget, never decompressed bytes/second or a decoder completion ETA.
// Queue pressure proposes a one-permit change; measured parent throughput
// accepts it or rolls it back. The outer allocator owns the parent budget.
class BgzfRateController {
public:
    struct Observation {
        int budget = 1;
        unsigned workers = 0;
        uint64_t completedPairs = 0, inputWaitNs = 0;
        uint64_t ready = 0, capacity = 0;
        double seconds = 0.25;
        bool decodePending = false;
        bool processorReady = false;
        bool allowNewTrial = true;
    };
    struct Decision {
        int limit = 1;
        double rate = 0;
        bool changed = false;
        const char* reason = "observe";
    };

    Decision observe(const Observation& o) {
        const int ceiling = std::min(std::max(1, static_cast<int>(o.workers)),
                                     std::max(1, o.budget - 1));
        Decision d;
        if (!primed_ || o.budget != budget_ || ceiling != ceiling_) {
            const int previous = limit_;
            limit_ = primed_ ? std::min(limit_, ceiling)
                            : std::min(ceiling, std::max(1, o.budget / 4));
            primed_ = true;
            budget_ = o.budget;
            ceiling_ = ceiling;
            previousPairs_ = o.completedPairs;
            previousWait_ = o.inputWaitNs;
            rate_ = 0;
            settle_ = 2;
            pending_ = false;
            rejectedLimit_ = -1;
            pressure_ = 0;
            d.limit = limit_;
            d.changed = previous != limit_;
            d.reason = "budget_or_supply_changed";
            return d;
        }
        const uint64_t reads = o.completedPairs >= previousPairs_
            ? o.completedPairs - previousPairs_ : 0;
        const uint64_t inputWait = o.inputWaitNs >= previousWait_
            ? o.inputWaitNs - previousWait_ : 0;
        previousPairs_ = o.completedPairs;
        previousWait_ = o.inputWaitNs;
        const double instant = o.seconds > 0 ? reads / o.seconds : 0;
        rate_ = rate_ > 0 ? 0.5 * rate_ + 0.5 * instant : instant;
        d.limit = limit_;
        d.rate = rate_;
        if (settle_ > 0) {
            --settle_;
            d.reason = "settling";
            return d;
        }
        if (pending_) {
            trialReads_ += reads;
            trialSeconds_ += o.seconds;
            if (++trialWindows_ < 4) {
                d.reason = "measure_trial";
                return d;
            }
            const double trialRate = trialSeconds_ > 0 ? trialReads_ / trialSeconds_ : 0;
            pending_ = false;
            // A no-progress library transition is not a throughput observation.
            if (trialReads_ > 0 && baselineRate_ > 0 && trialRate < baselineRate_ * 0.90) {
                rejectedLimit_ = limit_;
                rejectedAtRate_ = baselineRate_;
                limit_ = previousLimit_;
                rate_ = baselineRate_;
                d.limit = limit_;
                d.changed = true;
                d.reason = "restore_parent_rate";
                settle_ = 4;
            } else {
                d.reason = "accept_parent_rate";
                settle_ = 2;
            }
            return d;
        }

        if (!o.allowNewTrial) return d;
        const double fill = o.capacity ? double(o.ready) / o.capacity : 0;
        int direction = 0;
        if (o.capacity && fill >= 0.70 && o.processorReady && limit_ > 1) {
            direction = -1;
        } else if (o.capacity && fill <= 0.25 && inputWait > 0 &&
                   o.decodePending && limit_ < ceiling_) {
            direction = 1;
        }
        if (!direction) {
            pressure_ = 0;
            return d;
        }
        if (limit_ + direction == rejectedLimit_ && rate_ >= rejectedAtRate_ * 0.75 &&
            rate_ <= rejectedAtRate_ * 1.25) {
            pressure_ = 0;
            d.reason = "retain_better_parent_rate";
            return d;
        }
        pressure_ = (pressure_ * direction > 0) ? pressure_ + direction : direction;
        if (pressure_ < 2 && pressure_ > -2) return d;
        pressure_ = 0;
        previousLimit_ = limit_;
        baselineRate_ = rate_;
        limit_ += direction;
        pending_ = true;
        trialWindows_ = 0;
        trialReads_ = 0;
        trialSeconds_ = 0;
        settle_ = 2;
        d.limit = limit_;
        d.changed = true;
        d.reason = direction > 0 ? "input_starved_trial" : "buffer_full_trial";
        return d;
    }

    bool settled() const { return primed_ && settle_ == 0 && !pending_; }
    int limit() const { return limit_; }

private:
    bool primed_ = false, pending_ = false;
    int budget_ = 0, ceiling_ = 1, limit_ = 1, previousLimit_ = 1;
    int settle_ = 0, pressure_ = 0, trialWindows_ = 0;
    int rejectedLimit_ = -1;
    uint64_t previousPairs_ = 0, previousWait_ = 0, trialReads_ = 0;
    double rate_ = 0, baselineRate_ = 0, trialSeconds_ = 0;
    double rejectedAtRate_ = 0;
};
#endif
