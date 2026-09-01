#include "input/BgzfIndex.h"
#include "input/BgzfRangeReader.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

struct Options {
    std::string mode;
    std::string input;
    uint32_t threads = 1;
    uint32_t workers = 1;
    bool crcCheck = true;
};

Options parse_args(int argc, char* argv[]) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--mode" && index + 1 < argc) {
            options.mode = argv[++index];
        } else if (arg == "--input" && index + 1 < argc) {
            options.input = argv[++index];
        } else if (arg == "--threads" && index + 1 < argc) {
            options.threads = static_cast<uint32_t>(std::stoul(argv[++index]));
        } else if (arg == "--workers" && index + 1 < argc) {
            options.workers = static_cast<uint32_t>(std::stoul(argv[++index]));
        } else if (arg == "--crc-check" && index + 1 < argc) {
            const std::string value = argv[++index];
            if (value != "0" && value != "1") {
                throw std::runtime_error("--crc-check must be 0 or 1");
            }
            options.crcCheck = value == "1";
        } else if (arg == "--workers" || arg == "--crc-check") {
            if (index + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value");
            }
        } else {
            throw std::runtime_error("unknown or incomplete option: " + arg);
        }
    }
    if (options.mode.empty() || options.input.empty()) {
        throw std::runtime_error("--mode and --input are required");
    }
    if (options.threads == 0 || options.workers == 0) {
        throw std::runtime_error("--threads and --workers must be positive");
    }
    return options;
}

uint64_t fnv1a_update(uint64_t value, const std::string& field) {
    const uint64_t prime = 1099511628211ULL;
    for (size_t index = 0; index < field.size(); ++index) {
        value ^= static_cast<unsigned char>(field[index]);
        value *= prime;
    }
    return value;
}

uint64_t record_checksum(const star::input::BgzfFastqRecord& record) {
    uint64_t value = 14695981039346656037ULL;
    value = fnv1a_update(value, record.name);
    value ^= 0;
    value *= 1099511628211ULL;
    value = fnv1a_update(value, record.sequence);
    value ^= 0;
    value *= 1099511628211ULL;
    return fnv1a_update(value, record.quality);
}

uint64_t partition_start(uint64_t total, uint32_t part, uint32_t parts) {
    const uint64_t quotient = total / parts;
    const uint64_t remainder = total % parts;
    return quotient * part + std::min<uint64_t>(part, remainder);
}

bool scan_records(const Options& options,
                  uint64_t* record_count,
                  uint64_t* checksum,
                  std::string* error) {
    star::input::BgzfIndex index;
    if (!index.open(options.input, options.workers, error)) {
        return false;
    }
    *record_count = index.record_count();
    *checksum = 0;
    if (index.record_count() == 0) {
        return true;
    }
    const uint32_t workers = static_cast<uint32_t>(
        std::min<uint64_t>(options.workers, index.record_count()));
    std::vector<uint64_t> counts(workers, 0);
    std::vector<uint64_t> checksums(workers, 0);
    std::vector<std::thread> threads;
    std::atomic<bool> failed(false);
    std::mutex error_mutex;
    std::string worker_error;
    threads.reserve(workers);
    for (uint32_t worker = 0; worker < workers; ++worker) {
        threads.emplace_back([&, worker]() {
            const uint64_t begin = partition_start(index.record_count(), worker, workers);
            const uint64_t end = partition_start(index.record_count(), worker + 1, workers);
            star::input::BgzfRangeReader reader;
            std::string local_error;
            if (!reader.open(&index, begin, end - begin, options.crcCheck, &local_error)) {
                failed.store(true);
            } else {
                star::input::BgzfFastqRecord record;
                while (!failed.load() && reader.next(&record, &local_error)) {
                    checksums[worker] += record_checksum(record);
                    ++counts[worker];
                }
                if (local_error.empty() && counts[worker] != end - begin) {
                    local_error = "BGZF range reader returned an incomplete record range";
                }
                if (!local_error.empty()) {
                    failed.store(true);
                }
            }
            if (!local_error.empty()) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (worker_error.empty()) {
                    worker_error = local_error;
                }
            }
        });
    }
    for (size_t worker = 0; worker < threads.size(); ++worker) {
        threads[worker].join();
    }
    if (failed.load()) {
        *error = worker_error.empty() ? "BGZF record worker failed" : worker_error;
        return false;
    }
    uint64_t observed_count = 0;
    for (uint32_t worker = 0; worker < workers; ++worker) {
        observed_count += counts[worker];
        *checksum += checksums[worker];
    }
    if (observed_count != index.record_count()) {
        *error = "BGZF record workers did not cover the complete input";
        return false;
    }
    return true;
}

void print_index(const star::input::BgzfIndex& index) {
    std::cout << "{\"blocks\":[";
    const std::vector<star::input::BgzfBlock>& blocks = index.blocks();
    for (size_t block_index = 0; block_index < blocks.size(); ++block_index) {
        if (block_index != 0) {
            std::cout << ',';
        }
        const star::input::BgzfBlock& block = blocks[block_index];
        std::cout << "{\"compressed_offset\":" << block.compressedOffset
                  << ",\"compressed_size\":" << block.compressedSize
                  << ",\"first_record_offset\":";
        if (block.hasFirstRecordOffset) {
            std::cout << block.firstRecordOffset;
        } else {
            std::cout << "null";
        }
        std::cout << ",\"first_record_ordinal\":" << block.firstRecordOrdinal
                  << ",\"isize\":" << block.isize
                  << ",\"records_starting\":" << block.recordsStarting << '}';
    }
    std::cout << "],\"cache_status\":\""
              << star::input::bgzf_cache_status_name(index.cache_status())
              << "\",\"record_count\":" << index.record_count() << "}\n";
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        const Options options = parse_args(argc, argv);
        std::string error;
        if (options.mode == "detect") {
            star::input::BgzfDetection detection;
            if (!star::input::BgzfIndex::detect(options.input, &detection, &error)) {
                std::cerr << "ERROR: " << error << '\n';
                return 1;
            }
            std::cout << "{\"bgzf\":" << (detection.isBgzf ? "true" : "false")
                      << ",\"has_eof\":" << (detection.hasEofMarker ? "true" : "false")
                      << "}\n";
            return 0;
        }
        if (options.mode == "index") {
            star::input::BgzfIndex index;
            if (!index.open(options.input, options.threads, &error)) {
                std::cerr << "ERROR: " << error << '\n';
                return 1;
            }
            print_index(index);
            return 0;
        }
        if (options.mode == "records") {
            uint64_t record_count = 0;
            uint64_t checksum = 0;
            if (!scan_records(options, &record_count, &checksum, &error)) {
                std::cerr << "ERROR: " << error << '\n';
                return 1;
            }
            std::cout << "{\"checksum64\":\"" << std::hex << std::setw(16)
                      << std::setfill('0') << checksum << "\",\"record_count\":"
                      << std::dec << record_count << "}\n";
            return 0;
        }
        throw std::runtime_error("unsupported --mode: " + options.mode);
    } catch (const std::exception& exc) {
        std::cerr << "ERROR: " << exc.what() << '\n';
        return 2;
    }
}
