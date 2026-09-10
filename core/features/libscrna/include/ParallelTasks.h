#ifndef STAR_SUITE_SCRNA_PARALLEL_TASKS_H
#define STAR_SUITE_SCRNA_PARALLEL_TASKS_H

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>
#include "ThreadPermits.h"

namespace scrna {
// Bounded independent tasks. The calling thread participates, so workers is
// the total execution budget. Join every worker before propagating failures.
// Each task owns its output slot; worker IDs can select disjoint sub-budgets.
template<class Function>
void parallelFor(size_t tasks, size_t workers, const Function& function, size_t* launchedWorkers = nullptr) {
    if (!tasks) return;
    workers = std::max<size_t>(1, std::min(workers, tasks));
    ThreadPermitPool* pool = currentThreadPermitPool();
    std::atomic<size_t> next(0);
    std::atomic<bool> failed(false);
    std::exception_ptr error;
    std::mutex errorMutex;
    auto recordError = [&]() {
        std::lock_guard<std::mutex> lock(errorMutex);
        if (!error) error = std::current_exception();
        failed.store(true);
    };
    auto run = [&](size_t worker) {
        ScopedThreadPermits context(pool);
        try {
            while (!failed.load()) {
                const size_t task = next.fetch_add(1);
                if (task >= tasks) break;
                function(task, worker);
            }
        } catch (...) { recordError(); }
    };
    std::vector<std::thread> threads;
    threads.reserve(workers - 1);
    try {
        if (pool) {
            // Poll for returned permits between tasks. A running sampler can
            // expand immediately after another sample releases its workers.
            while (!failed.load()) {
                while (threads.size() + 1 < workers && next.load() < tasks && pool->tryAcquire()) {
                    const size_t worker = threads.size() + 1;
                    try {
                        threads.emplace_back([&, worker]() {
                            ThreadPermit permit(*pool, true);
                            run(worker);
                        });
                    } catch (...) { pool->release(); throw; }
                }
                const size_t task = next.fetch_add(1);
                if (task >= tasks) break;
                function(task, 0);
            }
        } else {
            for (size_t worker = 1; worker < workers; ++worker)
                threads.emplace_back(run, worker);
            run(0);
        }
    } catch (...) { recordError(); }
    for (auto& thread : threads) thread.join();
    if (launchedWorkers) *launchedWorkers = threads.size() + 1;
    if (error) std::rethrow_exception(error);
}
}
#endif
