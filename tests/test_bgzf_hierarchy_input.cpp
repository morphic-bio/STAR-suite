#include "ThreadControl.h"
#include "input/BgzfRangeReader.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <zlib.h>

using D = ThreadControl::PermitDomain;
using W = ThreadControl::PermitWork;
using namespace star::input;

static uint64_t acquire(void* p) {
    return static_cast<ThreadControl*>(p)->mapPermitAcquireForDomain(D::FEATURE, W::BGZF);
}
static void release(void* p, uint64_t wait, uint64_t units, uint64_t bytes, uint64_t ns) {
    static_cast<ThreadControl*>(p)->mapPermitReleaseForDomain(D::FEATURE, wait, units, bytes, ns, W::BGZF);
}
static void observe(void* p, const void* reader, uint64_t ready, uint64_t outstanding,
                    uint64_t capacity, unsigned workers, int waiting, int live) {
    static_cast<ThreadControl*>(p)->mapPermitObserveDecode(D::FEATURE, reader,
        ready, outstanding, capacity, workers, waiting, live);
}

int main(int argc, char** argv) {
    if (argc != 5) return 2; // R1 R2 read-limit output-prefix
    ThreadControl pool;
    pool.mapPermitConfigure(true, 1, 1, true, false);
    pool.mapPermitConfigureDomainFloors({0, 1, 0});
    pool.mapPermitStartHierarchy(std::string(argv[4]) + ".tsv", true);
    BgzfWorkPermitHooks hooks;
    hooks.context = &pool; hooks.acquire = acquire; hooks.release = release; hooks.observe = observe;
    uint64_t count = 0, limit = std::strtoull(argv[3], nullptr, 10);
    uLong crc = crc32(0, Z_NULL, 0);
    std::string error;
    {
        BgzfRangeReader readers[2];
        BgzfFastqRecord records[2];
        bool opened = true;
        for (int i = 0; i < 2; ++i) {
            if (!readers[i].open(argv[i+1], 0, UINT64_MAX, 4, true, &error, &hooks)) {
                opened = false; break;
            }
        }
        while (opened && error.empty() && (!limit || count < limit)) {
            const bool first = readers[0].next(&records[0], &error);
            if (!error.empty()) break;
            const bool second = readers[1].next(&records[1], &error);
            if (!error.empty()) break;
            if (first != second) { error = "mate count mismatch"; break; }
            if (!first) break;
            size_t n[2] = {records[0].nameLength, records[1].nameLength};
            for (int i = 0; i < 2; ++i)
                if (n[i] >= 2 && records[i].name_data()[n[i]-2] == '/') n[i] -= 2;
            if (n[0] != n[1] || std::memcmp(records[0].name_data(), records[1].name_data(), n[0])) {
                error = "mate name mismatch"; break;
            }
            const auto token = pool.mapPermitAcquireForDomain(D::FEATURE);
            for (int i = 0; i < 2; ++i) {
                crc = crc32(crc, reinterpret_cast<const Bytef*>(records[i].sequence_data()), records[i].sequenceLength);
                crc = crc32(crc, reinterpret_cast<const Bytef*>(records[i].quality_data()), records[i].qualityLength);
            }
            ++count;
            pool.mapPermitReleaseForDomain(D::FEATURE, token, 1, 0, 1);
        }
    } // Early limit and malformed-input paths must join all pending inflaters.
    pool.mapPermitMarkDomainComplete(D::FEATURE);
    pool.mapPermitStopHierarchy();
    const auto snap = pool.mapPermitSnapshot();
    const auto& decode = snap.featureDomain.decode;
    const bool balanced = snap.inUsePermits == 0 && snap.currentWaiters == 0 &&
        snap.availablePermits == 1 && decode.inUse == 0 && decode.waiters == 0 &&
        decode.acquireCalls == decode.releaseCalls && decode.workers == 0;
    assert(snap.featureDomain.completedReadPairs == count);
    std::printf("pairs=%llu crc=%08lx balanced=%d\n", static_cast<unsigned long long>(count), crc, balanced);
    if (!error.empty()) std::fprintf(stderr, "%s\n", error.c_str());
    return !balanced || !error.empty();
}
