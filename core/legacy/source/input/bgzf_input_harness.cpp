#include "input/BgzfIndex.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

struct Options {
    std::string mode;
    std::string input;
    uint32_t threads = 1;
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
        } else if (arg == "--workers" || arg == "--crc-check") {
            if (index + 1 >= argc) {
                throw std::runtime_error(arg + " requires a value");
            }
            ++index;
        } else {
            throw std::runtime_error("unknown or incomplete option: " + arg);
        }
    }
    if (options.mode.empty() || options.input.empty()) {
        throw std::runtime_error("--mode and --input are required");
    }
    return options;
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
            std::cerr << "ERROR: records mode requires Phase 2 BgzfRangeReader\n";
            return 1;
        }
        throw std::runtime_error("unsupported --mode: " + options.mode);
    } catch (const std::exception& exc) {
        std::cerr << "ERROR: " << exc.what() << '\n';
        return 2;
    }
}
