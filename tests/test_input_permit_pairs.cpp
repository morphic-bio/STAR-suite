#include "ThreadControl.h"
#include <atomic>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

using D = ThreadControl::PermitDomain;
using W = ThreadControl::PermitWork;

static void waitFor(const std::function<bool()>& done) {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
    while (!done()) {
        assert(std::chrono::steady_clock::now() < deadline);
        std::this_thread::yield();
    }
}

// The producer may use its parent's earlier ticket, but cannot move the
// unrelated external ticket. Repeat with a one-permit machine and a parent
// receiving only one of four permits; the remaining three belong to EXTERNAL.
static void ticketOrder(int budget, bool low) {
    ThreadControl pool;
    pool.mapPermitConfigure(true, budget, budget, true, false);
    pool.mapPermitConfigureInputPairs(true);
    int reader;
    pool.mapPermitObserveDecode(D::MAP, &reader, low ? 0 : 8, 8, 8, 1, low, 1);
    std::vector<uint64_t> external;
    for (int i = 1; i < budget; ++i) external.push_back(pool.mapPermitAcquireForDomain(D::EXTERNAL));
    const auto held = pool.mapPermitAcquireForDomain(D::MAP);
    std::atomic<int> sequence{0}, release{0};
    std::vector<std::thread> actors;
    auto queue = [&](D domain, W work, int expected) {
        actors.emplace_back([&, domain, work, expected] {
            auto token = pool.mapPermitAcquireForDomain(domain, work);
            assert(sequence.fetch_add(1) + 1 == expected);
            waitFor([&] { return release.load() >= expected; });
            pool.mapPermitReleaseForDomain(domain, token, 1, 1, 1, work);
        });
        waitFor([&] { return pool.mapPermitSnapshot().currentWaiters == actors.size(); });
    };
    queue(D::MAP, W::PROCESS, low ? 3 : 1);
    queue(D::EXTERNAL, W::PROCESS, 2);
    queue(D::MAP, W::BGZF, low ? 1 : 3);
    pool.mapPermitReleaseForDomain(D::MAP, held, 1, 1, 1);
    for (int i = 1; i <= 3; ++i) {
        waitFor([&] { return sequence.load() == i; });
        assert(pool.mapPermitSnapshot().inUsePermits == budget);
        release.store(i);
    }
    for (auto& actor : actors) actor.join();
    for (auto token : external) pool.mapPermitReleaseForDomain(D::EXTERNAL, token, 1, 1, 1);
    pool.mapPermitObserveDecode(D::MAP, &reader, 0, 0, 0, 0, 0, 0);
    assert(pool.mapPermitSnapshot().availablePermits == budget);
}

struct Pipeline {
    std::mutex mutex;
    std::condition_variable changed;
    std::deque<int> ready;
    int producers = 2, filling = 0, waiting = 0;
    std::atomic<int> next{0}, consumed{0};
};

// A buffered reader or one temporarily in I/O may have no queued acquire.
// A long mapper must still leave its parent's input slot reclaimable. Other
// parents can borrow it; unregistering the reader gives the slot to mapping.
static void protectsInputBeforeAcquire() {
    ThreadControl pool;
    pool.mapPermitConfigure(true, 4, 4, true, false);
    pool.mapPermitConfigureInputPairs(true);
    int reader;
    pool.mapPermitObserveDecode(D::MAP, &reader, 8, 8, 8, 1, 0, 1);
    std::vector<uint64_t> held;
    for (int i = 0; i < 3; ++i) held.push_back(pool.mapPermitAcquireForDomain(D::MAP));
    std::atomic<bool> mapped{false};
    std::thread consumer([&] {
        auto token = pool.mapPermitAcquireForDomain(D::MAP);
        mapped = true;
        pool.mapPermitReleaseForDomain(D::MAP, token, 1, 1, 1);
    });
    waitFor([&] { return pool.mapPermitSnapshot().mapDomain.currentWaiters == 1; });
    assert(!mapped && pool.mapPermitSnapshot().availablePermits == 1);
    auto external = pool.mapPermitAcquireForDomain(D::EXTERNAL);
    assert(!mapped && pool.mapPermitSnapshot().externalDomain.inUse == 1);
    pool.mapPermitReleaseForDomain(D::EXTERNAL, external, 1, 1, 1);
    auto input = pool.mapPermitAcquireForDomain(D::MAP, W::BGZF);
    assert(!mapped && pool.mapPermitSnapshot().mapDomain.inUse == 4);
    pool.mapPermitReleaseForDomain(D::MAP, input, 1, 1, 1, W::BGZF);
    assert(!mapped);
    pool.mapPermitObserveDecode(D::MAP, &reader, 0, 0, 0, 0, 0, 0);
    waitFor([&] { return mapped.load(); });
    consumer.join();
    for (auto token : held) pool.mapPermitReleaseForDomain(D::MAP, token, 1, 1, 1);
    assert(pool.mapPermitSnapshot().availablePermits == 4);
}

// Independent work may borrow an idle input slot, but must repay the loan
// before taking another lease when that input becomes runnable. Only that
// slot gets reclaimed; the independent job receives it again after input.
static void reclaimsBorrowedInput() {
    ThreadControl pool;
    pool.mapPermitConfigure(true, 4, 4, true, false);
    pool.mapPermitConfigureInputPairs(true);
    int reader;
    pool.mapPermitObserveDecode(D::MAP, &reader, 0, 1, 8, 2, 1, 1);
    std::vector<uint64_t> mapping;
    for (int i = 0; i < 3; ++i) mapping.push_back(pool.mapPermitAcquireForDomain(D::MAP));
    auto loan = pool.mapPermitAcquireForDomain(D::EXTERNAL);
    std::atomic<int> order{0};
    std::atomic<bool> inputDone{false};
    std::thread independent([&] {
        auto token = pool.mapPermitAcquireForDomain(D::EXTERNAL);
        assert(order.fetch_add(1) == 1);
        pool.mapPermitReleaseForDomain(D::EXTERNAL, token, 1, 1, 1);
    });
    waitFor([&] { return pool.mapPermitSnapshot().externalDomain.currentWaiters == 1; });
    std::thread input([&] {
        auto token = pool.mapPermitAcquireForDomain(D::MAP, W::BGZF);
        assert(order.fetch_add(1) == 0);
        inputDone = true;
        pool.mapPermitReleaseForDomain(D::MAP, token, 1, 1, 1, W::BGZF);
    });
    waitFor([&] { return pool.mapPermitSnapshot().mapDomain.decode.waiters == 1; });
    pool.mapPermitReleaseForDomain(D::EXTERNAL, loan, 1, 1, 1);
    input.join(); independent.join();
    assert(inputDone && order == 2);
    for (auto token : mapping) pool.mapPermitReleaseForDomain(D::MAP, token, 1, 1, 1);
    pool.mapPermitObserveDecode(D::MAP, &reader, 0, 0, 0, 0, 0, 0);
    assert(pool.mapPermitSnapshot().availablePermits == 4);
}

// Real bounded producer/consumer queues compete with independent EXTERNAL
// work. Consumers never hold a compute permit while waiting for input; input
// never holds one while waiting for buffer space. Every item must finish,
// including EOF, delayed input and a one-permit shared pool.
static void pipelineProgress(int budget, bool floors) {
    ThreadControl pool;
    pool.mapPermitConfigure(true, budget, budget, true, false);
    if (floors) pool.mapPermitConfigureDomainFloors({1, 1, 1});
    pool.mapPermitConfigureInputPairs(true);
    Pipeline pipes[2];
    std::vector<std::thread> actors;
    constexpr int reads = 200;
    for (int i = 0; i < 2; ++i) {
        const D domain = i ? D::FEATURE : D::MAP;
        auto observe = [&, i, domain] {
            auto& p = pipes[i]; // caller holds the pipeline mutex
            pool.mapPermitObserveDecode(domain, &p, p.ready.size(), p.ready.size() + p.filling,
                                        8, 2, p.waiting != 0, 1);
        };
        for (int j = 0; j < 2; ++j) actors.emplace_back([&, i, domain, observe] {
            auto& p = pipes[i];
            for (;;) {
                const int n = p.next.fetch_add(1);
                if (n >= reads) break;
                {
                    std::unique_lock<std::mutex> lock(p.mutex);
                    p.changed.wait(lock, [&] { return p.ready.size() + p.filling < 8; });
                    ++p.filling;
                    observe();
                }
                if (n % 17 == 0) std::this_thread::sleep_for(std::chrono::microseconds(200));
                auto token = pool.mapPermitAcquireForDomain(domain, W::BGZF);
                std::this_thread::sleep_for(std::chrono::microseconds(50));
                pool.mapPermitReleaseForDomain(domain, token, 1, 1, 50000, W::BGZF);
                {
                    std::lock_guard<std::mutex> lock(p.mutex);
                    --p.filling;
                    p.ready.push_back(n);
                    observe();
                }
                p.changed.notify_all();
            }
            {
                std::lock_guard<std::mutex> lock(p.mutex);
                --p.producers;
            }
            p.changed.notify_all();
        });
        for (int j = 0; j < 3; ++j) actors.emplace_back([&, i, domain, observe] {
            auto& p = pipes[i];
            for (;;) {
                {
                    std::unique_lock<std::mutex> lock(p.mutex);
                    ++p.waiting;
                    observe();
                    p.changed.wait(lock, [&] { return !p.ready.empty() || !p.producers; });
                    --p.waiting;
                    if (p.ready.empty()) { observe(); break; }
                    p.ready.pop_front();
                    observe();
                }
                p.changed.notify_all();
                auto token = pool.mapPermitAcquireForDomain(domain);
                std::this_thread::sleep_for(std::chrono::microseconds(150));
                const auto s = pool.mapPermitSnapshot();
                assert(s.inUsePermits <= budget && s.availablePermits >= 0);
                assert(s.mapDomain.inUse + s.featureDomain.inUse + s.externalDomain.inUse == s.inUsePermits);
                assert(s.mapDomain.decode.inUse <= s.mapDomain.inUse);
                assert(s.featureDomain.decode.inUse <= s.featureDomain.inUse);
                pool.mapPermitReleaseForDomain(domain, token, 1, 1, 150000);
                ++p.consumed;
            }
        });
    }
    for (int i = 0; i < 2; ++i) actors.emplace_back([&] {
        for (int n = 0; n < reads; ++n) {
            auto token = pool.mapPermitAcquireForDomain(D::EXTERNAL);
            std::this_thread::sleep_for(std::chrono::microseconds(100));
            pool.mapPermitReleaseForDomain(D::EXTERNAL, token, 1, 1, 100000);
        }
    });
    for (auto& actor : actors) actor.join();
    for (int i = 0; i < 2; ++i) {
        assert(pipes[i].consumed == reads);
        pool.mapPermitObserveDecode(i ? D::FEATURE : D::MAP, &pipes[i], 0, 0, 0, 0, 0, 0);
    }
    for (auto d : {D::FEATURE, D::EXTERNAL}) pool.mapPermitMarkDomainComplete(d);
    const auto s = pool.mapPermitSnapshot();
    assert(s.availablePermits == budget && s.inUsePermits == 0 && s.currentWaiters == 0);
    assert(s.mapDomain.completedReadPairs == reads && s.featureDomain.completedReadPairs == reads);
    assert(s.externalDomain.completedReadPairs == 2 * reads);
    for (const auto* d : {&s.mapDomain, &s.featureDomain}) {
        assert(d->decode.acquireCalls == reads && d->decode.releaseCalls == reads);
        assert(d->decode.inUse == 0 && d->decode.waiters == 0);
    }
    // The surviving parent borrows every released permit.
    std::vector<uint64_t> loans;
    for (int i = 0; i < budget; ++i) loans.push_back(pool.mapPermitAcquireForDomain(D::MAP));
    assert(pool.mapPermitSnapshot().mapDomain.inUse == budget);
    for (auto token : loans) pool.mapPermitReleaseForDomain(D::MAP, token, 1, 1, 1);
    pool.mapPermitMarkDomainComplete(D::MAP);
    assert(pool.mapPermitSnapshot().mapDomain.floor == 0);
}

int main() {
    protectsInputBeforeAcquire();
    reclaimsBorrowedInput();
    for (int budget : {1, 4}) for (bool low : {false, true}) ticketOrder(budget, low);
    for (int budget : {1, 4, 16}) pipelineProgress(budget, budget > 1);
}
