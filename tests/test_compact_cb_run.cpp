#include "CbBucketStore.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unistd.h>

using star::solo::PackedCbRecord;
using star::solo::EncodedCbRun;
using star::solo::CbBucketStore;

static void same(std::vector<PackedCbRecord> a, std::vector<PackedCbRecord> b) {
    const auto less = [](const PackedCbRecord& x, const PackedCbRecord& y) {
        return x.key != y.key ? x.key < y.key : x.value < y.value;
    };
    std::sort(a.begin(), a.end(), less);
    std::sort(b.begin(), b.end(), less);
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i)
        assert(a[i].key == b[i].key && a[i].value == b[i].value);
}

static void fields() {
    const uint32_t counts[] = {0, 1, 62, 63, 64, 1000000, (1u << 30) - 1};
    std::vector<PackedCbRecord> expected;
    for (uint32_t cb : {737280u, 741375u})
        for (uint32_t umi : {0u, 0xFFFFFFu})
            for (uint16_t gene : {uint16_t(0), uint16_t(32767)})
                for (uint8_t tag : {uint8_t(0), uint8_t(31)})
                    for (uint8_t flags = 0; flags < 4; ++flags)
                        for (uint32_t count : counts)
                            expected.push_back(PackedCbRecord::make(cb, umi, gene, tag, count, flags));
    EncodedCbRun run;
    run.compact = true; run.cbBase = 737280;
    run.bytes.resize(expected.size() * 8);
    size_t escaped = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        run.set_record(i, expected[i]);
        escaped += expected[i].count30() >= 63;
        const auto decoded = run.record_at(i);
        assert(decoded.key == expected[i].key && decoded.value == expected[i].value);
        assert(run.sort_key(i) == expected[i].group_sort_key());
    }
    assert(run.overflow.size() == escaped && run.record_count() == expected.size());
    const auto legacy = run.legacy_bytes();
    assert(legacy.size() == expected.size() * 12);
    for (size_t i = 0; i < expected.size(); ++i) {
        PackedCbRecord decoded;
        assert(PackedCbRecord::decode(legacy.data() + 12 * i, &decoded));
        assert(decoded.key == expected[i].key && decoded.value == expected[i].value);
    }
    const size_t missing = run.overflow.back().record;
    run.overflow.pop_back();
    bool threw = false;
    try { run.record_at(missing); } catch (const std::runtime_error&) { threw = true; }
    assert(threw);
    threw = false;
    try { run.set_record(0, PackedCbRecord::make(737279, 0, 0, 0, 1, 0)); }
    catch (const std::out_of_range&) { threw = true; }
    assert(threw);
}

static void store(CbBucketStore::Mode mode) {
    char pattern[] = "/tmp/star_compact_bucket_XXXXXX";
    const char* directory = mkdtemp(pattern);
    assert(directory);
    CbBucketStore::Config config;
    config.mode = mode;
    config.whitelistSize = 8195; // deliberately not divisible by bucket count
    config.bucketCount = 4;
    config.scratchDirectory = directory;
    config.memoryBudgetBytes = 10000;
    config.mergeWorkerCount = mode == CbBucketStore::Mode::Spill ? 0 : 2;
    config.mergeFanIn = 2;
    {
        CbBucketStore buckets(config);
        std::vector<std::vector<PackedCbRecord>> expected(config.bucketCount);
        std::string error;
        for (uint32_t batch = 0; batch < 50; ++batch) {
            for (uint32_t b = 0; b < config.bucketCount; ++b) {
                const uint32_t begin = (uint64_t(b) * config.whitelistSize + config.bucketCount - 1) / config.bucketCount;
                const uint32_t end = (uint64_t(b + 1) * config.whitelistSize + config.bucketCount - 1) / config.bucketCount;
                assert(buckets.bucket_for_cb(begin) == b && buckets.bucket_for_cb(end - 1) == b);
                std::vector<PackedCbRecord> input;
                for (uint32_t i = 0; i < 32; ++i) {
                    const uint32_t count = i % 4 == 0 ? (1u << 30) - 1
                        : i % 4 == 1 ? 63 : i % 4 == 2 ? 62 : 1;
                    input.push_back(PackedCbRecord::make(i % 2 ? begin : end - 1,
                        (batch * 7 + i) % 19, i % 7, i % 5, count, (batch + i) % 4));
                }
                expected[b].insert(expected[b].end(), input.begin(), input.end());
                assert(buckets.append_segment(batch % 4, b, std::move(input), &error));
            }
        }
        assert(buckets.finalize(&error));
        assert(buckets.using_spill() == (mode != CbBucketStore::Mode::Ram));
        assert(buckets.transitioned_to_spill() == (mode == CbBucketStore::Mode::Auto));
        assert(buckets.payload_bytes() == 50 * config.bucketCount * 32 * 12);
        for (uint32_t b = 0; b < config.bucketCount; ++b) {
            // Reusable legacy interfaces must reconstruct exact wide values.
            std::vector<PackedCbRecord> decoded;
            assert(buckets.load_bucket(b, &decoded, &error));
            same(decoded, expected[b]);
            std::vector<std::vector<PackedCbRecord>> oldRuns;
            assert(buckets.load_sorted_segments(b, &oldRuns, &error));
            decoded.clear();
            for (const auto& run : oldRuns) decoded.insert(decoded.end(), run.begin(), run.end());
            same(decoded, expected[b]);
            std::vector<EncodedCbRun> runs;
            assert(buckets.consume_compact_segments(b, &runs, &error));
            decoded.clear(); size_t escaped = 0;
            for (const auto& run : runs) {
                assert(run.compact == (mode == CbBucketStore::Mode::Ram));
                escaped += run.overflow.size();
                for (size_t i = 0; i < run.record_count(); ++i) {
                    const auto record = run.record_at(i);
                    assert(run.sort_key(i) == record.group_sort_key());
                    if (i) assert(run.sort_key(i - 1) <= run.sort_key(i));
                    decoded.push_back(record);
                }
            }
            same(decoded, expected[b]);
            if (mode == CbBucketStore::Mode::Ram) {
                assert(escaped == 50 * 16);
                assert(!buckets.consume_compact_segments(b, &runs, &error));
            } else {
                assert(buckets.consume_compact_segments(b, &runs, &error));
            }
        }
        if (mode == CbBucketStore::Mode::Ram) assert(buckets.async_merge_count() > 0);
    }
    // The store destructor owns removal of its spill files.
    assert(rmdir(directory) == 0);
}

static void wideBucket() {
    CbBucketStore::Config config;
    config.whitelistSize = 1u << 20;
    config.bucketCount = 1;
    CbBucketStore buckets(config);
    std::string error;
    const auto record = PackedCbRecord::make((1u << 20) - 1, 0xFFFFFF, 32767, 31,
                                            (1u << 30) - 1, 3);
    assert(buckets.append_segment(0, 0, {record}, &error));
    assert(buckets.finalize(&error));
    std::vector<EncodedCbRun> runs;
    assert(buckets.consume_compact_segments(0, &runs, &error));
    assert(runs.size() == 1 && !runs[0].compact && runs[0].bytes.size() == 12);
    same({runs[0].record_at(0)}, {record});
}

int main() {
    fields();
    store(CbBucketStore::Mode::Ram);
    store(CbBucketStore::Mode::Spill);
    store(CbBucketStore::Mode::Auto);
    wideBucket();
    std::cout << "PASS: compact fields, wide-count escapes, merged overflow positions, bucket boundaries, reusable legacy readers, RAM/spill transition and wide-bucket fallback\n";
}
