#include "input/BgzfBlockReader.h"
#include "input/BgzfRangeReader.h"
#include "input/BgzfStarAdapter.h"

#include <cstdlib>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace {

struct Options {
    std::string mode;
    std::string input;
    std::string input2;  // second mate for --mode pair
    uint32_t threads = 1;
    uint32_t workers = 1;
    bool crcCheck = true;
    bool validateReadNames = true;
};

Options parse_args(int argc, char* argv[]) {
    Options options;
    for (int index = 1; index < argc; ++index) {
        const std::string arg = argv[index];
        if (arg == "--mode" && index + 1 < argc) {
            options.mode = argv[++index];
        } else if (arg == "--input" && index + 1 < argc) {
            options.input = argv[++index];
        } else if (arg == "--input2" && index + 1 < argc) {
            options.input2 = argv[++index];
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
        } else if (arg == "--validate-read-names" && index + 1 < argc) {
            const std::string value = argv[++index];
            if (value != "0" && value != "1") {
                throw std::runtime_error("--validate-read-names must be 0 or 1");
            }
            options.validateReadNames = value == "1";
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

uint64_t fnv1a_update(uint64_t value, const char* field, size_t field_size) {
    const uint64_t prime = 1099511628211ULL;
    for (size_t index = 0; index < field_size; ++index) {
        value ^= static_cast<unsigned char>(field[index]);
        value *= prime;
    }
    return value;
}

uint64_t record_checksum(const star::input::BgzfFastqRecord& record) {
    uint64_t value = 14695981039346656037ULL;
    value = fnv1a_update(value, record.name_data(), record.nameLength);
    value ^= 0;
    value *= 1099511628211ULL;
    value = fnv1a_update(value, record.sequence_data(), record.sequenceLength);
    value ^= 0;
    value *= 1099511628211ULL;
    return fnv1a_update(value, record.quality_data(), record.qualityLength);
}

bool scan_blocks(const std::string& path,
                 std::vector<star::input::BgzfBlock>* blocks,
                 std::string* error) {
    struct stat info;
    if (::stat(path.c_str(), &info) != 0 || info.st_size < 0) {
        *error = "could not stat BGZF input";
        return false;
    }
    const uint64_t file_size = static_cast<uint64_t>(info.st_size);
    const int input_fd = ::open(path.c_str(), O_RDONLY);
    if (input_fd < 0) {
        *error = "could not open BGZF input";
        return false;
    }
    uint64_t offset = 0;
    while (offset < file_size) {
        star::input::BgzfBlock block;
        bool is_eof = false;
        if (!star::input::read_bgzf_block_header_fd(
                input_fd, offset, file_size, &block, &is_eof, error)) {
            ::close(input_fd);
            return false;
        }
        if (is_eof) {
            offset = file_size;
            break;
        }
        blocks->push_back(block);
        offset += block.compressedSize;
    }
    ::close(input_fd);
    return true;
}

void print_blocks(const std::vector<star::input::BgzfBlock>& blocks) {
    std::cout << "{\"blocks\":[";
    for (size_t index = 0; index < blocks.size(); ++index) {
        if (index != 0) {
            std::cout << ',';
        }
        std::cout << "{\"compressed_offset\":" << blocks[index].compressedOffset
                  << ",\"compressed_size\":" << blocks[index].compressedSize
                  << ",\"isize\":" << blocks[index].isize << '}';
    }
    std::cout << "]}\n";
}

bool scan_records(const Options& options,
                  uint64_t* record_count,
                  uint64_t* checksum,
                  std::string* error) {
    star::input::BgzfRangeReader reader;
    if (!reader.open(options.input, 0, std::numeric_limits<uint64_t>::max(),
                     options.workers, options.crcCheck, error)) {
        return false;
    }
    *record_count = 0;
    *checksum = 0;
    star::input::BgzfFastqRecord record;
    while (reader.next(&record, error)) {
        *checksum += record_checksum(record);
        ++*record_count;
    }
    return error->empty();
}

// Pair two mate files through the STAR adapter exactly as a fused Flex lane
// does, so ordinal, record-count, and read-name validation can be tested
// without a genome.
bool scan_pairs(const Options& options,
                uint64_t* record_count,
                uint64_t* mate1_name_bytes,
                uint64_t* mate0_name_view_records,
                uint64_t* sequence_view_records,
                uint64_t* quality_view_records,
                std::string* error) {
    star::input::BgzfStarAdapter adapter;
    star::input::BgzfStarAdapterOptions adapter_options;
    adapter_options.mate0_reader_threads = options.workers;
    adapter_options.mate1_reader_threads = options.workers;
    adapter_options.crc_check = options.crcCheck;
    adapter_options.validate_read_names = options.validateReadNames;
    if (!adapter.open(options.input, options.input2, adapter_options, error)) {
        return false;
    }
    *record_count = 0;
    *mate1_name_bytes = 0;
    *mate0_name_view_records = 0;
    *sequence_view_records = 0;
    *quality_view_records = 0;
    const size_t batch_capacity = 17;
    std::vector<star::input::BgzfStarRecord> records(batch_capacity);
    star::input::BgzfStarBatchLease lease;
    while (true) {
        size_t returned = 0;
        const star::input::InputStatus status = adapter.next_records(
            records.data(), records.size(), &returned, error, &lease);
        if (status == star::input::InputStatus::Record) {
            for (size_t index = 0; index < returned; ++index) {
                const star::input::BgzfStarRecord& record = records[index];
                ++*record_count;
                *mate1_name_bytes += record.mates[1].nameLength;
                *mate0_name_view_records += record.mates[0].name_is_view() ? 1 : 0;
                *sequence_view_records += record.mates[0].sequence_is_view() ? 1 : 0;
                *sequence_view_records += record.mates[1].sequence_is_view() ? 1 : 0;
                *quality_view_records += record.mates[0].quality_is_view() ? 1 : 0;
                *quality_view_records += record.mates[1].quality_is_view() ? 1 : 0;
                // Dereference every leased field after next_records() has
                // released the adapter lock. This is the lifetime required by
                // fused Flex consumers.
                (void)record.mates[0].name_data()[0];
                (void)record.mates[0].sequence_data()[0];
                (void)record.mates[1].sequence_data()[0];
                (void)record.mates[0].quality_data()[0];
                (void)record.mates[1].quality_data()[0];
            }
            continue;
        }
        return status == star::input::InputStatus::End;
    }
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        const Options options = parse_args(argc, argv);
        std::string error;
        if (options.mode == "detect") {
            star::input::BgzfDetection detection;
            if (!star::input::detect_bgzf(options.input, &detection, &error)) {
                std::cerr << "ERROR: " << error << '\n';
                return 1;
            }
            std::cout << "{\"bgzf\":" << (detection.isBgzf ? "true" : "false")
                      << ",\"has_eof\":" << (detection.hasEofMarker ? "true" : "false")
                      << "}\n";
            return 0;
        }
        if (options.mode == "blocks") {
            std::vector<star::input::BgzfBlock> blocks;
            if (!scan_blocks(options.input, &blocks, &error)) {
                std::cerr << "ERROR: " << error << '\n';
                return 1;
            }
            print_blocks(blocks);
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
        if (options.mode == "pair") {
            if (options.input2.empty()) {
                throw std::runtime_error("--mode pair requires --input2");
            }
            uint64_t record_count = 0;
            uint64_t mate1_name_bytes = 0;
            uint64_t mate0_name_view_records = 0;
            uint64_t sequence_view_records = 0;
            uint64_t quality_view_records = 0;
            if (!scan_pairs(options, &record_count, &mate1_name_bytes,
                            &mate0_name_view_records, &sequence_view_records,
                            &quality_view_records, &error)) {
                std::cerr << "ERROR: " << error << '\n';
                return 1;
            }
            std::cout << "{\"mate1_name_bytes\":" << mate1_name_bytes
                      << ",\"mate0_name_view_records\":"
                      << mate0_name_view_records
                      << ",\"quality_view_records\":"
                      << quality_view_records
                      << ",\"sequence_view_records\":"
                      << sequence_view_records
                      << ",\"fastq_record_bytes\":"
                      << sizeof(star::input::BgzfFastqRecord)
                      << ",\"record_count\":" << record_count << "}\n";
            return 0;
        }
        throw std::runtime_error("unsupported --mode: " + options.mode);
    } catch (const std::exception& exc) {
        std::cerr << "ERROR: " << exc.what() << '\n';
        return 2;
    }
}
