#include "pf_bgzf_input.h"
#include "input/BgzfRangeReader.h"

#include <array>
#include <cstdio>
#include <cstring>
#include <exception>
#include <sys/stat.h>

using namespace star::input;

struct pf_bgzf_input {
    int streams = 0;
    std::array<BgzfRangeReader, 3> readers;
    std::array<BgzfFastqRecord, 3> records;
    std::array<BgzfBatchLease, 3> leases;
};

static void report(char *error, size_t capacity, const std::string& message) {
    if (error && capacity) std::snprintf(error, capacity, "%s", message.c_str());
}

extern "C" int pf_bgzf_detect(const char *path, char *error, size_t capacity) {
    try {
        report(error, capacity, "");
        if (!path) { report(error, capacity, "null FASTQ path"); return -1; }
        struct stat st;
        if (stat(path, &st) != 0) {
            report(error, capacity, std::string("cannot stat FASTQ: ") + path);
            return -1;
        }
        if (!S_ISREG(st.st_mode)) return 0;
        BgzfDetection detection;
        std::string message;
        if (!detect_bgzf(path, &detection, &message)) {
            report(error, capacity, message); return -1;
        }
        return detection.isBgzf ? 1 : 0;
    } catch (const std::exception& e) {
        report(error, capacity, e.what()); return -1;
    } catch (...) {
        report(error, capacity, "unknown BGZF detection exception"); return -1;
    }
}

extern "C" pf_bgzf_input *pf_bgzf_open(const char *const *paths, int streams,
        unsigned threads, int crc, const pf_bgzf_permits *permits,
        char *error, size_t capacity) {
    try {
        report(error, capacity, "");
        if (!paths || streams < 2 || streams > 3) {
            report(error, capacity, "BGZF PF input requires two or three streams");
            return nullptr;
        }
        std::unique_ptr<pf_bgzf_input> input(new pf_bgzf_input);
        input->streams = streams;
        BgzfWorkPermitHooks hooks;
        if (permits) {
            hooks.context = permits->context;
            hooks.acquire = permits->acquire;
            hooks.release = permits->release;
        }
        for (int i = 0; i < streams; ++i) {
            std::string message;
            if (!paths[i] || !input->readers[i].open(paths[i], 0, UINT64_MAX,
                    threads / streams + (unsigned(i) < threads % streams), crc != 0,
                    &message, &hooks, true, BgzfNameMode::Token)) {
                report(error, capacity, message.empty() ? "null BGZF path" : message);
                return nullptr;
            }
        }
        return input.release();
    } catch (const std::exception& e) {
        report(error, capacity, e.what()); return nullptr;
    } catch (...) {
        report(error, capacity, "unknown BGZF open exception"); return nullptr;
    }
}

static size_t normalized_length(const BgzfFastqRecord& r) {
    size_t n = r.nameLength;
    if (n >= 2 && r.name_data()[n - 2] == '/' &&
        r.name_data()[n - 1] >= '1' && r.name_data()[n - 1] <= '3') n -= 2;
    return n;
}

extern "C" int pf_bgzf_next(pf_bgzf_input *input, pf_bgzf_record *out,
                             char *error, size_t capacity) {
    try {
        report(error, capacity, "");
        if (!input || !out) { report(error, capacity, "null BGZF record output"); return -1; }
        int present = 0;
        for (int i = 0; i < input->streams; ++i) {
            input->leases[i].clear();
            std::string message;
            if (input->readers[i].next(&input->records[i], &message, &input->leases[i])) ++present;
            else if (!message.empty()) { report(error, capacity, message); return -1; }
        }
        if (!present) return 0;
        if (present != input->streams) {
            report(error, capacity, "BGZF mate record counts differ"); return -1;
        }
        const auto& first = input->records[0];
        const size_t n = normalized_length(first);
        for (int i = 0; i < input->streams; ++i) {
            const auto& r = input->records[i];
            if (r.ordinal != first.ordinal || normalized_length(r) != n ||
                std::memcmp(first.name_data(), r.name_data(), n) != 0) {
                report(error, capacity, "BGZF mate names differ at record " + std::to_string(first.ordinal));
                return -1;
            }
            out[i] = {r.name_data(), r.sequence_data(), r.quality_data(),
                      r.nameLength, r.sequenceLength, r.qualityLength, r.ordinal};
        }
        return 1;
    } catch (const std::exception& e) {
        report(error, capacity, e.what()); return -1;
    } catch (...) {
        report(error, capacity, "unknown BGZF read exception"); return -1;
    }
}

extern "C" void pf_bgzf_close(pf_bgzf_input *input) { delete input; }

#include <algorithm>
#include <atomic>
#include <chrono>
#include <deque>
#include <stdexcept>

extern "C" int pf_bgzf_process_batches(const char *const *paths, unsigned lanes, int streams,
    unsigned workers, unsigned inflater_threads, int crc, uint64_t max_reads,
    const pf_bgzf_permits *permits, pf_bgzf_batch_consumer consume, void *context,
    char *error, size_t error_size) {
    if (!paths || !lanes || !workers || streams < 2 || streams > 3 || !consume) {
        report(error, error_size, "invalid BGZF batch dispatch arguments"); return 0;
    }
    constexpr size_t batch_size = 512;
    struct Batch {
        std::vector<BgzfFastqRecord> owned;
        std::array<BgzfBatchLease, 3> leases;
        std::vector<pf_bgzf_record> views;
        size_t count = 0;
        explicit Batch(int streams) : owned(batch_size * streams), views(batch_size * streams) {}
    };
    std::mutex mutex;
    std::condition_variable ready, space;
    std::deque<std::unique_ptr<Batch>> free, queue;
    std::vector<std::thread> producers, consumers;
    std::atomic<bool> failed(false);
    unsigned remaining = lanes;
    size_t peak = 0;
    uint64_t records = 0, decode_ns = 0;
    std::string message;
    auto fail = [&](const std::string& msg) {
        std::lock_guard<std::mutex> lock(mutex);
        if (!failed.exchange(true)) message = msg;
        ready.notify_all(); space.notify_all();
    };
    try {
        for (size_t i = 0; i < workers * 2 + lanes; ++i) free.emplace_back(new Batch(streams));
        for (unsigned worker = 0; worker < workers; ++worker) consumers.emplace_back([&, worker] {
            try {
                while (true) {
                    std::unique_ptr<Batch> batch;
                    {
                        std::unique_lock<std::mutex> lock(mutex);
                        ready.wait(lock, [&] { return failed || !queue.empty() || !remaining; });
                        if (failed || queue.empty()) break;
                        batch = std::move(queue.front()); queue.pop_front();
                    }
                    if (!consume(context, worker, batch->views.data(), batch->count, streams))
                        throw std::runtime_error("PF direct batch consumer failed");
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        free.push_back(std::move(batch));
                    }
                    space.notify_one();
                }
            } catch (const std::exception& e) { fail(e.what()); }
              catch (...) { fail("unknown PF batch consumer exception"); }
        });
        for (unsigned lane = 0; lane < lanes; ++lane) producers.emplace_back([&, lane] {
            try {
                char local_error[1024] = "";
                std::unique_ptr<pf_bgzf_input, decltype(&pf_bgzf_close)> input(
                    pf_bgzf_open(paths + lane * streams, streams,
                        inflater_threads / lanes + (lane < inflater_threads % lanes),
                        crc, permits, local_error, sizeof(local_error)), pf_bgzf_close);
                if (!input) throw std::runtime_error(local_error);
                uint64_t count = 0, elapsed = 0;
                bool end = false;
                while (!failed && !end && (!max_reads || count < max_reads)) {
                    std::unique_ptr<Batch> batch;
                    {
                        std::unique_lock<std::mutex> lock(mutex);
                        space.wait(lock, [&] { return failed || !free.empty(); });
                        if (failed) break;
                        batch = std::move(free.front()); free.pop_front();
                    }
                    for (auto& lease : batch->leases) lease.clear();
                    batch->count = 0;
                    auto start = std::chrono::steady_clock::now();
                    while (!failed && batch->count < batch_size && (!max_reads || count < max_reads)) {
                        const size_t offset = batch->count * streams;
                        int present = 0;
                        for (int i = 0; i < streams; ++i) {
                            std::string err;
                            if (input->readers[i].next(&batch->owned[offset+i], &err, &batch->leases[i])) ++present;
                            else if (!err.empty()) throw std::runtime_error(err);
                        }
                        if (!present) { end = true; break; }
                        if (present != streams) throw std::runtime_error("BGZF mate record counts differ");
                        const auto& first = batch->owned[offset];
                        const size_t name_length = normalized_length(first);
                        for (int i = 0; i < streams; ++i) {
                            const auto& rec = batch->owned[offset+i];
                            if (rec.ordinal != first.ordinal || normalized_length(rec) != name_length ||
                                std::memcmp(first.name_data(), rec.name_data(), name_length))
                                throw std::runtime_error("BGZF mate names differ at record " + std::to_string(first.ordinal));
                            batch->views[offset+i] = {rec.name_data(), rec.sequence_data(), rec.quality_data(),
                                rec.nameLength, rec.sequenceLength, rec.qualityLength, rec.ordinal};
                        }
                        ++batch->count; ++count;
                    }
                    elapsed += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - start).count();
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        if (batch->count) { queue.push_back(std::move(batch)); peak = std::max(peak, queue.size()); }
                        else free.push_back(std::move(batch));
                    }
                    ready.notify_one();
                }
                std::lock_guard<std::mutex> lock(mutex);
                records += count; decode_ns += elapsed;
            } catch (const std::exception& e) { fail(e.what()); }
              catch (...) { fail("unknown BGZF batch producer exception"); }
            {
                std::lock_guard<std::mutex> lock(mutex);
                --remaining;
            }
            ready.notify_all();
        });
    } catch (const std::exception& e) { fail(e.what()); }
      catch (...) { fail("BGZF dispatcher startup failed"); }
    for (auto& thread : producers) thread.join();
    for (auto& thread : consumers) thread.join();
    report(error, error_size, message);
    std::fprintf(stderr, "[pf-bgzf] handoff=direct records=%llu batch_records=%zu queue_peak_batches=%zu batch_capacity=%u decode_pair_seconds_sum=%.6f\n",
        (unsigned long long)records, batch_size, peak, workers*2+lanes, double(decode_ns)/1e9);
    return failed ? 0 : 1;
}
