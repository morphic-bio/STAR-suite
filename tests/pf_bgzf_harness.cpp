#include "pf_api.h"
#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <stdexcept>
#include <vector>
#include <zlib.h>

// A single bounded pool used by BOTH PF consumers and native inflate workers.
struct Permits {
    std::mutex mutex;
    std::condition_variable cv;
    int available = 1, active = 0, peak = 0;
    unsigned long long acquired = 0, released = 0;
};
static uint64_t acquire(void *v) {
    auto& p = *static_cast<Permits*>(v);
    std::unique_lock<std::mutex> lock(p.mutex);
    p.cv.wait(lock, [&] { return p.available > 0; });
    --p.available; ++p.active; ++p.acquired;
    if (p.active > p.peak) p.peak = p.active;
    return 0;
}
static uint64_t throw_acquire(void *) { throw std::runtime_error("injected inflate failure"); }
static void release(void *v, uint64_t, uint64_t, uint64_t, uint64_t) {
    auto& p = *static_cast<Permits*>(v);
    std::lock_guard<std::mutex> lock(p.mutex);
    ++p.available; --p.active; ++p.released;
    p.cv.notify_one();
}
int main(int argc, char **argv) {
    // decode THREADS PERMITS LIMIT R1 R2 [R3]
    // assign MODE THREADS PERMITS WHITELIST FEATURES OUT R1 R2 [R1 R2 ...]
    if (argc < 7) return 2;
    const bool decode = !strcmp(argv[1], "decode");
    Permits permits;
    permits.available = atoi(argv[decode ? 3 : 4]);
#ifndef PF_BGZF_BASELINE
    if (decode) {
        int streams = argc - 5;
        if (streams < 2 || streams > 3) return 2;
        char error[1024];
        pf_bgzf_permits hooks = {&permits, getenv("PF_TEST_THROW_INFLATE") ? throw_acquire : acquire, release};
        auto *input = pf_bgzf_open(const_cast<const char **>(argv + 5), streams,
            atoi(argv[2]), 1, permits.available ? &hooks : nullptr, error, sizeof(error));
        if (!input) { fprintf(stderr, "%s\n", error); return 1; }
        unsigned long count = 0, limit = strtoul(argv[4], nullptr, 10);
        uLong digest = crc32(0, Z_NULL, 0);
        int rc = 0;
        while (!limit || count < limit) {
            pf_bgzf_record records[3];
            rc = pf_bgzf_next(input, records, error, sizeof(error));
            if (rc <= 0) break;
            for (int i = 0; i < streams; ++i) {
                auto& r = records[i];
                for (auto span : {std::pair<const char*, size_t>{r.name, r.name_length},
                                  {r.sequence, r.sequence_length}, {r.quality, r.quality_length}}) {
                    digest = crc32(digest, (const Bytef*)span.first, span.second);
                    digest = crc32(digest, (const Bytef*)"\n", 1);
                }
            }
            ++count;
        }
        pf_bgzf_close(input);
        if (rc < 0) { fprintf(stderr, "%s\n", error); return 1; }
        printf("records=%lu crc32=%08lx permits_peak=%d balanced=%d\n", count, digest,
               permits.peak, permits.acquired == permits.released);
        return permits.active != 0;
    }
#endif
    if (decode || strcmp(argv[1], "assign") || argc < 10 || (argc-8)%2) return 2;
    auto *cfg = pf_config_create();
    pf_config_set_search_threads(cfg, 1);
    pf_config_set_consumer_threads(cfg, getenv("PF_TEST_CONSUMERS") ? atoi(getenv("PF_TEST_CONSUMERS")) : 2);
    pf_config_set_read_buffer_lines(cfg, getenv("PF_TEST_BUFFER_LINES") ? atoi(getenv("PF_TEST_BUFFER_LINES")) : 12); // sustained backpressure, < 64 record permit batch
    pf_config_set_skip_emptydrops(cfg, 1);
    pf_config_set_skip_qc_outputs(cfg, 1);
    pf_config_set_min_counts(cfg, 0);
    pf_config_set_feature_offset(cfg, 0);
    pf_config_set_limit_search(cfg, -1);
    if (permits.available) pf_config_set_permit_hooks(cfg, acquire, release, &permits);
#ifndef PF_BGZF_BASELINE
    pf_bgzf_mode mode = !strcmp(argv[2], "off") ? PF_BGZF_OFF :
                        !strcmp(argv[2], "range") ? PF_BGZF_RANGE : PF_BGZF_AUTO;
    if (pf_config_set_bgzf_input(cfg, mode, atoi(argv[3]), 1)) return 2;
#endif
    auto *ctx = pf_init(cfg);
    pf_config_destroy(cfg);
    if (!ctx) return 1;
    pf_error err = pf_load_whitelist(ctx, argv[5]);
    if (!err) err = pf_load_feature_ref(ctx, argv[6]);
    std::vector<const char*> r1, r2;
    for (int i = 8; i < argc; i += 2) { r1.push_back(argv[i]); r2.push_back(argv[i+1]); }
    pf_stats stats = {};
    if (!err) err = pf_process_fastqs(ctx, r1.data(), r2.data(), r1.size(), argv[7], "sample", &stats);
    if (err) fprintf(stderr, "PF error %d: %s\n", err, pf_get_error(ctx));
    pf_destroy(ctx);
    printf("status=%d permits_peak=%d balanced=%d\n", err, permits.peak,
           permits.acquired == permits.released);
    return err || permits.active != 0;
}
