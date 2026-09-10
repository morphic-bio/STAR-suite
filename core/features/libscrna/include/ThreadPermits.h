#ifndef STAR_SUITE_SCRNA_THREAD_PERMITS_H
#define STAR_SUITE_SCRNA_THREAD_PERMITS_H
#include <atomic>
#include <cassert>
#include <cstddef>
#include <stdexcept>

namespace scrna {
// Call-local CPU budget shared by independent sample coordinators and their
// helpers. Borrowing is nonblocking: a coordinator always does useful work
// itself, so nested callers cannot deadlock while holding their base permit.
class ThreadPermitPool {
    const size_t capacity_;
    std::atomic<size_t> free_, peak_, acquisitions_;
public:
    explicit ThreadPermitPool(size_t capacity)
        : capacity_(capacity), free_(capacity), peak_(0), acquisitions_(0) {}
    bool tryAcquire() {
        size_t available = free_.load();
        while (available) {
            if (free_.compare_exchange_weak(available, available - 1)) {
                const size_t used = capacity_ - available + 1;
                size_t peak = peak_.load();
                while (peak < used && !peak_.compare_exchange_weak(peak, used)) {}
                ++acquisitions_;
                return true;
            }
        }
        return false;
    }
    void release() {
        const size_t before = free_.fetch_add(1);
        assert(before < capacity_);
        (void)before;
    }
    size_t available() const { return free_.load(); }
    size_t peakUsed() const { return peak_.load(); }
    size_t acquisitions() const { return acquisitions_.load(); }
};
class ThreadPermit {
    ThreadPermitPool& pool_;
public:
    explicit ThreadPermit(ThreadPermitPool& pool, bool alreadyAcquired = false) : pool_(pool) {
        if (!alreadyAcquired && !pool_.tryAcquire()) throw std::logic_error("No base caller permit");
    }
    ~ThreadPermit() { pool_.release(); }
    ThreadPermit(const ThreadPermit&) = delete;
    ThreadPermit& operator=(const ThreadPermit&) = delete;
};
inline ThreadPermitPool*& currentThreadPermitPool() {
    static thread_local ThreadPermitPool* pool = nullptr;
    return pool;
}
class ScopedThreadPermits {
    ThreadPermitPool* previous_;
public:
    explicit ScopedThreadPermits(ThreadPermitPool* pool) : previous_(currentThreadPermitPool()) {
        currentThreadPermitPool() = pool;
    }
    ~ScopedThreadPermits() { currentThreadPermitPool() = previous_; }
    ScopedThreadPermits(const ScopedThreadPermits&) = delete;
    ScopedThreadPermits& operator=(const ScopedThreadPermits&) = delete;
};
}
#endif
