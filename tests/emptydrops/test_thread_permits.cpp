#include "ParallelTasks.h"
#include "EmptyDropsCRSampler.h"
#include <chrono>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>

static void require(bool ok, const char* text) { if (!ok) throw std::runtime_error(text); }
int main() {
    scrna::ThreadPermitPool pool(4);
    {
        scrna::ThreadPermit base(pool);
        scrna::ScopedThreadPermits context(&pool);
        std::vector<std::unique_ptr<scrna::ThreadPermit>> held;
        for (int i=0;i<3;++i) held.emplace_back(new scrna::ThreadPermit(pool));
        std::atomic<bool> started(false), released(false);
        std::thread finishingGroup([&]() {
            while (!started.load()) std::this_thread::yield();
            held.clear(); released=true;
        });
        std::vector<unsigned> visits(80,0);
        std::atomic<unsigned> active(0), peak(0);
        size_t launched=0;
        scrna::parallelFor(visits.size(), 12, [&](size_t task,size_t worker) {
            require(scrna::currentThreadPermitPool()==&pool, "helper inherits the shared budget");
            if (task==0) {
                require(worker==0, "starts on its reserved worker with no spare permits");
                started=true;
                while (!released.load()) std::this_thread::yield();
            }
            unsigned now=++active, old=peak.load();
            while (old<now && !peak.compare_exchange_weak(old,now)) {}
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            ++visits[task];--active;
        }, &launched);
        finishingGroup.join();
        require(launched==4 && peak>1 && peak<=4, "running task expands after permits return, within capacity");
        for (auto v:visits) require(v==1, "no dropped or duplicate tasks");
        require(pool.available()==3, "all helper permits returned");
        bool caught=false;
        try {
            scrna::parallelFor(40, 12, [&](size_t task,size_t) {
                if (task==3) throw std::runtime_error("deliberate worker failure");
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            });
        } catch (const std::runtime_error&) { caught=true; }
        require(caught && pool.available()==3, "failure joins helpers and returns every permit");
    }
    require(pool.available()==4 && !scrna::currentThreadPermitPool(), "coordinator releases its base permit and context");
    const vector<uint32_t> totals={20,40,70}, lengths={4,4,4};
    const vector<double> prob={-60,-53,-50,-45,-139,-126,-119,-110,-275,-254,-241,-229};
    const vector<double> ambient={.05,.1,.15,.25,.45};
    for (uint64_t seed:{1ULL,9137ULL}) {
        const auto reference=EmptyDropsCRSampler::montecarloPval(totals,lengths,prob,ambient,257,seed,1);
        bool nontrivial=false;
        for (auto n:reference) nontrivial|=n>0 && n<257;
        require(nontrivial, "fixture includes nontrivial Monte Carlo counts");
        for (uint32_t workers:{2u,7u}) {
            require(reference==EmptyDropsCRSampler::montecarloPval(totals,lengths,prob,ambient,257,seed,workers),
                "batching preserves exact per-iteration samples and integer tallies");
        }
        scrna::ThreadPermit base(pool);
        scrna::ScopedThreadPermits context(&pool);
        require(reference==EmptyDropsCRSampler::montecarloPval(totals,lengths,prob,ambient,257,seed,12),
            "permit scheduling preserves Monte Carlo results across translation units");
        require(pool.available()==3, "sampler releases helper permits");
    }
    require(pool.peakUsed()<=4 && pool.available()==4, "shared budget never exceeded or leaked");
    std::cout << "PASS: live permit borrowing, bounded concurrency, exception cleanup and exact sampler parity\n";
}
