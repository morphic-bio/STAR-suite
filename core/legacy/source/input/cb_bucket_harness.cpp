#include "CbBucketStore.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <vector>

namespace {

struct Options {
    std::string mode;
    std::string fixture;
    std::string scratch;
    std::uint32_t producers = 1;
    std::uint32_t records = 8192;
    std::uint32_t bucketCount = 256;
    std::uint32_t whitelistSize = 16;
};

Options parse_args(int argc, char **argv)
{
    Options result;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--mode" && index + 1 < argc)
            result.mode = argv[++index];
        else if (arg == "--fixture" && index + 1 < argc)
            result.fixture = argv[++index];
        else if (arg == "--scratch" && index + 1 < argc)
            result.scratch = argv[++index];
        else if (arg == "--producers" && index + 1 < argc)
            result.producers = std::stoul(argv[++index]);
        else if (arg == "--records" && index + 1 < argc)
            result.records = std::stoul(argv[++index]);
        else if (arg == "--bucket-count" && index + 1 < argc)
            result.bucketCount = std::stoul(argv[++index]);
        else if (arg == "--whitelist-size" && index + 1 < argc)
            result.whitelistSize = std::stoul(argv[++index]);
        else
            throw std::runtime_error("unknown or incomplete option: " + arg);
    }
    if (result.mode.empty())
        throw std::runtime_error("--mode is required");
    return result;
}

std::uint32_t pack_umi(const std::string &sequence)
{
    std::uint32_t result = 0;
    for (char base : sequence) {
        const std::string alphabet = "ACGT";
        const std::size_t value = alphabet.find(base);
        if (value == std::string::npos)
            throw std::runtime_error("invalid UMI base in fixture");
        result = (result << 2) | static_cast<std::uint32_t>(value);
    }
    return result;
}

std::vector<star::solo::PackedCbRecord> load_fixture(const std::string &path)
{
    std::ifstream input(path.c_str());
    if (!input)
        throw std::runtime_error("cannot open fixture " + path);
    std::vector<star::solo::PackedCbRecord> result;
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        std::istringstream row(line);
        std::uint64_t readId = 0;
        std::uint32_t cb = 0, gene = 0, tag = 0, region = 0;
        std::string umi;
        if (!(row >> readId >> cb >> umi >> gene >> tag >> region))
            throw std::runtime_error("malformed fixture row");
        (void)readId;
        result.push_back(star::solo::PackedCbRecord::make(
            cb, pack_umi(umi), static_cast<std::uint16_t>(gene),
            static_cast<std::uint8_t>(tag), 1,
            static_cast<std::uint8_t>(region)));
    }
    return result;
}

void mkdir_if_needed(const std::string &path)
{
    if (!path.empty() && ::mkdir(path.c_str(), 0700) != 0 && errno != EEXIST)
        throw std::runtime_error("cannot create scratch directory " + path);
}

void test_roundtrip()
{
    const std::uint32_t maxima[][6] = {
        {0, 0, 0, 0, 0, 0},
        {(1u << 20) - 1, (1u << 24) - 1, (1u << 15) - 1,
         (1u << 5) - 1, (1u << 30) - 1, 3},
        {737279, 0xABCDEF, 18128, 16, 1234567, 2},
    };
    for (const auto &values : maxima) {
        const star::solo::PackedCbRecord expected = star::solo::PackedCbRecord::make(
            values[0], values[1], values[2], values[3], values[4], values[5]);
        std::uint8_t bytes[star::solo::PackedCbRecord::kSerializedBytes];
        expected.encode(bytes);
        star::solo::PackedCbRecord observed;
        if (!star::solo::PackedCbRecord::decode(bytes, &observed)
            || observed.key != expected.key || observed.value != expected.value
            || observed.cb_index() != values[0] || observed.umi24() != values[1]
            || observed.gene15() != values[2] || observed.tag5() != values[3]
            || observed.count30() != values[4] || observed.flags2() != values[5])
            throw std::runtime_error("packed CB record round-trip failed");
    }
}

void print_partition(const Options &options)
{
    star::solo::CbBucketStore::Config config;
    config.bucketCount = options.bucketCount;
    config.whitelistSize = options.whitelistSize;
    star::solo::CbBucketStore store(config);
    const std::vector<star::solo::PackedCbRecord> records = load_fixture(options.fixture);
    std::vector<std::uint64_t> counts(options.bucketCount, 0);
    for (const auto &record : records)
        ++counts[store.bucket_for_cb(record.cb_index())];
    std::cout << "{\"buckets\":[";
    bool first = true;
    for (std::uint32_t bucket = 0; bucket < options.bucketCount; ++bucket) {
        if (counts[bucket] == 0)
            continue;
        if (!first)
            std::cout << ',';
        first = false;
        std::cout << "{\"bucket\":" << bucket << ",\"records\":"
                  << counts[bucket] << '}';
    }
    std::cout << "]}\n";
}

void test_claims(const Options &options)
{
    star::solo::CbBucketStore::Config config;
    config.bucketCount = options.bucketCount;
    config.whitelistSize = 737280;
    config.mergeWorkerCount = 2;
    config.mergeFanIn = 2;
    star::solo::CbBucketStore store(config);
    std::vector<std::thread> threads;
    std::atomic<bool> failed(false);
    std::string firstError;
    std::mutex errorMutex;
    for (std::uint32_t worker = 0; worker < options.producers; ++worker) {
        threads.emplace_back([&, worker]() {
            std::vector<std::vector<star::solo::PackedCbRecord> > segments(options.bucketCount);
            for (std::uint32_t ordinal = worker; ordinal < options.records;
                 ordinal += options.producers) {
                const std::uint32_t cb = (ordinal * 7919u) % config.whitelistSize;
                const auto record = star::solo::PackedCbRecord::make(
                    cb, ordinal, ordinal % 30000u, (ordinal % 31u) + 1u, 1, ordinal % 4u);
                segments[store.bucket_for_cb(cb)].push_back(record);
            }
            for (std::uint32_t bucket = 0; bucket < options.bucketCount; ++bucket) {
                if (segments[bucket].empty())
                    continue;
                std::string error;
                if (!store.append_segment(worker, bucket, std::move(segments[bucket]), &error)) {
                    failed.store(true);
                    std::lock_guard<std::mutex> lock(errorMutex);
                    if (firstError.empty()) firstError = error;
                    return;
                }
            }
        });
    }
    for (std::thread &thread : threads)
        thread.join();
    std::string error;
    if (failed.load() || !store.finalize(&error))
        throw std::runtime_error(firstError.empty() ? error : firstError);
    std::vector<std::uint8_t> seen(options.records, 0);
    std::uint64_t observed = 0;
    for (std::uint32_t bucket = 0; bucket < options.bucketCount; ++bucket) {
        std::vector<star::solo::PackedCbRecord> records;
        if (!store.load_bucket(bucket, &records, &error))
            throw std::runtime_error(error);
        for (const auto &record : records) {
            if (record.umi24() >= seen.size() || seen[record.umi24()]++)
                throw std::runtime_error("atomic segment claim duplicated or corrupted a record");
            ++observed;
        }
    }
    if (observed != options.records
        || std::find(seen.begin(), seen.end(), 0) != seen.end())
        throw std::runtime_error("atomic segment claim lost records");
    std::cout << "{\"producers\":" << options.producers
              << ",\"records\":" << observed << "}\n";
}

void test_store(const Options &options)
{
    mkdir_if_needed(options.scratch);
    const std::vector<star::solo::PackedCbRecord> input = load_fixture(options.fixture);
    star::solo::CbBucketStore::Config ramConfig;
    ramConfig.bucketCount = options.bucketCount;
    ramConfig.whitelistSize = options.whitelistSize;
    ramConfig.mode = star::solo::CbBucketStore::Mode::Ram;
    star::solo::CbBucketStore::Config spillConfig = ramConfig;
    spillConfig.mode = star::solo::CbBucketStore::Mode::Spill;
    spillConfig.scratchDirectory = options.scratch;
    spillConfig.filePrefix = "roundtrip";
    star::solo::CbBucketStore::Config asyncConfig = ramConfig;
    asyncConfig.mergeWorkerCount = 2;
    asyncConfig.mergeFanIn = 2;
    star::solo::CbBucketStore ram(ramConfig);
    star::solo::CbBucketStore spill(spillConfig);
    star::solo::CbBucketStore asyncRam(asyncConfig);
    std::vector<std::vector<star::solo::PackedCbRecord> > pending(options.bucketCount);
    std::string error;
    std::uint32_t worker = 0;
    for (const auto &record : input) {
        const std::uint32_t bucket = ram.bucket_for_cb(record.cb_index());
        pending[bucket].push_back(record);
        if (pending[bucket].size() == 7) {
            std::vector<star::solo::PackedCbRecord> copy = pending[bucket];
            std::vector<star::solo::PackedCbRecord> asyncCopy = pending[bucket];
            if (!ram.append_segment(worker, bucket, std::move(copy), &error)
                || !asyncRam.append_segment(worker, bucket,
                                            std::move(asyncCopy), &error)
                || !spill.append_segment(worker, bucket, std::move(pending[bucket]), &error))
                throw std::runtime_error(error);
            pending[bucket].clear();
            ++worker;
        }
    }
    for (std::uint32_t bucket = 0; bucket < options.bucketCount; ++bucket) {
        if (pending[bucket].empty())
            continue;
        std::vector<star::solo::PackedCbRecord> copy = pending[bucket];
        std::vector<star::solo::PackedCbRecord> asyncCopy = pending[bucket];
        if (!ram.append_segment(worker, bucket, std::move(copy), &error)
            || !asyncRam.append_segment(worker, bucket,
                                        std::move(asyncCopy), &error)
            || !spill.append_segment(worker, bucket, std::move(pending[bucket]), &error))
            throw std::runtime_error(error);
        ++worker;
    }
    if (!ram.finalize(&error) || !spill.finalize(&error)
        || !asyncRam.finalize(&error))
        throw std::runtime_error(error);
    for (std::uint32_t bucket = 0; bucket < options.bucketCount; ++bucket) {
        std::vector<std::uint8_t> ramBytes, spillBytes;
        if (!ram.load_bucket_bytes(bucket, &ramBytes, &error)
            || !spill.load_bucket_bytes(bucket, &spillBytes, &error))
            throw std::runtime_error(error);
        if (ramBytes != spillBytes)
            throw std::runtime_error("RAM/spill bucket byte mismatch");

        std::vector<std::vector<star::solo::PackedCbRecord> > ramSegments;
        std::vector<std::vector<star::solo::PackedCbRecord> > spillSegments;
        if (!ram.load_sorted_segments(bucket, &ramSegments, &error)
            || !spill.load_sorted_segments(bucket, &spillSegments, &error))
            throw std::runtime_error(error);
        if (ramSegments.size() != spillSegments.size())
            throw std::runtime_error("RAM/spill segment count mismatch");
        for (std::size_t segment = 0; segment < ramSegments.size(); ++segment) {
            const auto &ramRun = ramSegments[segment];
            const auto &spillRun = spillSegments[segment];
            if (ramRun.size() != spillRun.size())
                throw std::runtime_error("RAM/spill segment length mismatch");
            if (!std::is_sorted(
                    ramRun.begin(), ramRun.end(),
                    [](const star::solo::PackedCbRecord &left,
                       const star::solo::PackedCbRecord &right) {
                        return left.group_sort_key() < right.group_sort_key();
                    }))
                throw std::runtime_error("producer-local CB bucket run is not sorted");
            for (std::size_t index = 0; index < ramRun.size(); ++index) {
                if (ramRun[index].key != spillRun[index].key
                    || ramRun[index].value != spillRun[index].value)
                    throw std::runtime_error("RAM/spill decoded segment mismatch");
            }
        }

        std::vector<star::solo::PackedCbRecord> originalRecords;
        std::vector<star::solo::PackedCbRecord> asyncRecords;
        if (!ram.load_bucket(bucket, &originalRecords, &error)
            || !asyncRam.load_bucket(bucket, &asyncRecords, &error))
            throw std::runtime_error(error);
        const auto order = [](const star::solo::PackedCbRecord &left,
                              const star::solo::PackedCbRecord &right) {
            if (left.group_sort_key() != right.group_sort_key())
                return left.group_sort_key() < right.group_sort_key();
            return left.value < right.value;
        };
        std::sort(originalRecords.begin(), originalRecords.end(), order);
        std::sort(asyncRecords.begin(), asyncRecords.end(), order);
        if (originalRecords.size() != asyncRecords.size())
            throw std::runtime_error("asynchronous CB bucket merge lost records");
        for (std::size_t index = 0; index < originalRecords.size(); ++index) {
            if (originalRecords[index].key != asyncRecords[index].key
                || originalRecords[index].value != asyncRecords[index].value)
                throw std::runtime_error("asynchronous CB bucket merge changed records");
        }
    }
    if (asyncRam.async_merge_count() == 0
        || asyncRam.async_merged_records() == 0)
        throw std::runtime_error("asynchronous CB bucket merge was not exercised");
}

} // namespace

int main(int argc, char **argv)
{
    try {
        const Options options = parse_args(argc, argv);
        if (options.mode == "roundtrip")
            test_roundtrip();
        else if (options.mode == "partition")
            print_partition(options);
        else if (options.mode == "claims")
            test_claims(options);
        else if (options.mode == "store")
            test_store(options);
        else
            throw std::runtime_error("unsupported --mode: " + options.mode);
        return 0;
    } catch (const std::exception &exception) {
        std::cerr << "ERROR: " << exception.what() << '\n';
        return 1;
    }
}
