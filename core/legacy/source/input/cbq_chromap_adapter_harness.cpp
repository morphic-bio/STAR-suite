#include "input/CbqChromapAdapter.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

struct Args {
    std::string read_pair_cbq;
    std::string barcode_cbq;
    std::string output_dir;
    std::string sample_name = "sample";
    bool require_name_match = true;
    bool allow_missing_qualities = false;
};

void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " --readPairCbq reads.cbq --barcodeCbq barcodes.cbq --outputDir OUT"
        << " [--sampleName NAME] [--allowNameMismatch] [--allowMissingQualities]\n";
}

bool parse_args(int argc, char** argv, Args* args) {
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto need_value = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                return nullptr;
            }
            return argv[++i];
        };

        if (key == "--readPairCbq") {
            const char* value = need_value("--readPairCbq");
            if (!value) return false;
            args->read_pair_cbq = value;
        } else if (key == "--barcodeCbq") {
            const char* value = need_value("--barcodeCbq");
            if (!value) return false;
            args->barcode_cbq = value;
        } else if (key == "--outputDir") {
            const char* value = need_value("--outputDir");
            if (!value) return false;
            args->output_dir = value;
        } else if (key == "--sampleName") {
            const char* value = need_value("--sampleName");
            if (!value) return false;
            args->sample_name = value;
        } else if (key == "--allowNameMismatch") {
            args->require_name_match = false;
        } else if (key == "--allowMissingQualities") {
            args->allow_missing_qualities = true;
        } else if (key == "--help" || key == "-h") {
            return false;
        } else {
            std::cerr << "Unknown argument: " << key << "\n";
            return false;
        }
    }

    if (args->read_pair_cbq.empty() || args->barcode_cbq.empty() || args->output_dir.empty()) {
        std::cerr << "--readPairCbq, --barcodeCbq, and --outputDir are required\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, &args)) {
        usage(argv[0]);
        return 2;
    }

    star::input::CbqChromapAdapterOptions options;
    options.read_pair_cbq = args.read_pair_cbq;
    options.barcode_cbq = args.barcode_cbq;
    options.output_dir = args.output_dir;
    options.sample_name = args.sample_name;
    options.require_name_match = args.require_name_match;
    options.allow_missing_qualities = args.allow_missing_qualities;

    star::input::CbqChromapFastqPaths paths;
    std::string error;
    if (!star::input::materialize_cbq_chromap_fastqs(options, &paths, &error)) {
        std::cerr << "CBQ Chromap adapter failed";
        if (!error.empty()) {
            std::cerr << ": " << error;
        }
        std::cerr << "\n";
        return 1;
    }

    std::cout << "read1_fastq\t" << paths.read1_fastq << "\n";
    std::cout << "read2_fastq\t" << paths.read2_fastq << "\n";
    std::cout << "barcode_fastq\t" << paths.barcode_fastq << "\n";
    return 0;
}
