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
