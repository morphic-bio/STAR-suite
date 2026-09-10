#ifndef STAR_SCRNA_BOUNDED_SAMPLE_TASKS_H
#define STAR_SCRNA_BOUNDED_SAMPLE_TASKS_H
#include "ParallelTasks.h"
#include <cstdint>
#include <limits>

namespace scrna {
struct SampleTaskStats { size_t peakWorkers = 0, permitsReturned = 0; uint64_t peakEstimatedBytes = 0; };
// Prepare only a bounded window of matrices. A sample larger than the budget
// runs alone; the estimate is an admission limit, not a process RSS guarantee.
template<class Function>
SampleTaskStats boundedSampleTasks(const std::vector<uint64_t>& bytes, size_t workers,
                                  uint64_t memoryBudget, const Function& function) {
    workers = std::max<size_t>(1, workers);
    ThreadPermitPool pool(workers);
    SampleTaskStats stats;
    {
        ThreadPermit base(pool);
        ScopedThreadPermits context(&pool);
        for (size_t begin = 0; begin < bytes.size();) {
            size_t end = begin;
            uint64_t used = 0;
            do {
                const uint64_t value = bytes[end];
                if (end > begin && (value > memoryBudget || used > memoryBudget - value)) break;
                if (value > std::numeric_limits<uint64_t>::max() - used)
                    throw std::overflow_error("Sample matrix memory estimate overflow");
                used += value; ++end;
            } while (end < bytes.size() && end - begin < workers && used < memoryBudget);
            stats.peakEstimatedBytes = std::max(stats.peakEstimatedBytes, used);
            parallelFor(end - begin, workers, [&](size_t task, size_t) { function(begin + task); });
            begin = end;
        }
    }
    stats.peakWorkers = pool.peakUsed(); stats.permitsReturned = pool.available();
    return stats;
}
}
#endif
