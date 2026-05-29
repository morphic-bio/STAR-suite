#include "input/FastxInputModule.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <openssl/sha.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using star::input::FastxInputModule;
using star::input::InputRecord;
using star::input::InputSourcePlan;
using star::input::InputStatus;

namespace {

bool starts_with_option(const std::string& value) {
    return value.size() >= 2 && value[0] == '-' && value[1] == '-';
}

void split_string(const std::string& value, char delimiter, std::vector<std::string>& out) {
    out.clear();
    std::string item;
    std::istringstream in(value);
    while (std::getline(in, item, delimiter)) {
        if (!item.empty()) {
            out.push_back(item);
        }
    }
}

bool ends_with_case_insensitive(const std::string& value, const std::string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    const size_t offset = value.size() - suffix.size();
    for (size_t ii = 0; ii < suffix.size(); ++ii) {
        const unsigned char c1 = static_cast<unsigned char>(value[offset + ii]);
        const unsigned char c2 = static_cast<unsigned char>(suffix[ii]);
        if (std::tolower(c1) != std::tolower(c2)) {
            return false;
        }
    }
    return true;
}

std::string join_command(const std::vector<std::string>& tokens) {
    std::string command;
    for (const std::string& token : tokens) {
        if (!command.empty()) {
            command += " ";
        }
        command += token;
    }
    return command;
}

std::string sha256_hex(const std::string& value) {
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char*>(value.data()), value.size(), digest);

    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (unsigned char byte : digest) {
        out << std::setw(2) << static_cast<unsigned int>(byte);
    }
    return out.str();
}

std::vector<std::string> split_tab_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    std::istringstream in(line);
    while (std::getline(in, field, '\t')) {
        fields.push_back(field);
    }
    return fields;
}

std::string read_group_id(const std::string& rg_line) {
    std::string line = rg_line;
    if (line.substr(0, 3) != "ID:") {
        line = "ID:" + line;
    }
    const size_t tab = line.find('\t');
    return line.substr(3, tab == std::string::npos ? std::string::npos : tab - 3);
}

struct HarnessOptions {
    std::vector<std::string> read_files_in;
    std::vector<std::string> read_files_command;
    std::string read_files_manifest = "-";
    std::string read_files_prefix;
    std::vector<char> read_name_separator_chars{' '};
    bool legacy_zcat = false;
    bool dump_fastq = false;
};

void usage(std::ostream& out) {
    out << "Usage: fastx_input_harness [options]\n"
        << "  --readFilesIn R1[,R1_L2] [R2[,R2_L2] ...]\n"
        << "  --readFilesManifest manifest.tsv\n"
        << "  --readFilesCommand COMMAND [ARGS...]\n"
        << "  --readFilesPrefix PREFIX\n"
        << "  --readNameSeparator space|none|CHAR[,CHAR...]\n"
        << "  --readFilesLegacyZcat Yes|No\n"
        << "  --dump-fastq\n";
}

std::vector<char> parse_separators(const std::string& value) {
    std::vector<char> separators;
    std::vector<std::string> tokens;
    split_string(value, ',', tokens);
    for (const std::string& token : tokens) {
        if (token == "space") {
            separators.push_back(' ');
        } else if (token == "none") {
            continue;
        } else if (token.size() == 1) {
            separators.push_back(token[0]);
        } else {
            throw std::runtime_error("unsupported --readNameSeparator value: " + token);
        }
    }
    return separators;
}

HarnessOptions parse_args(int argc, char* argv[]) {
    HarnessOptions opts;
    for (int ii = 1; ii < argc; ++ii) {
        const std::string arg = argv[ii];
        if (arg == "--help" || arg == "-h") {
            usage(std::cout);
            std::exit(0);
        } else if (arg == "--readFilesIn") {
            while (ii + 1 < argc && !starts_with_option(argv[ii + 1])) {
                opts.read_files_in.push_back(argv[++ii]);
            }
        } else if (arg == "--readFilesCommand") {
            while (ii + 1 < argc && !starts_with_option(argv[ii + 1])) {
                opts.read_files_command.push_back(argv[++ii]);
            }
        } else if (arg == "--readFilesManifest") {
            if (++ii >= argc) {
                throw std::runtime_error("--readFilesManifest requires a value");
            }
            opts.read_files_manifest = argv[ii];
        } else if (arg == "--readFilesPrefix") {
            if (++ii >= argc) {
                throw std::runtime_error("--readFilesPrefix requires a value");
            }
            opts.read_files_prefix = argv[ii] == std::string("-") ? "" : argv[ii];
        } else if (arg == "--readNameSeparator") {
            if (++ii >= argc) {
                throw std::runtime_error("--readNameSeparator requires a value");
            }
            opts.read_name_separator_chars = parse_separators(argv[ii]);
        } else if (arg == "--readFilesLegacyZcat") {
            if (++ii >= argc) {
                throw std::runtime_error("--readFilesLegacyZcat requires a value");
            }
            std::string value = argv[ii];
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (value == "yes") {
                opts.legacy_zcat = true;
            } else if (value == "no" || value == "-") {
                opts.legacy_zcat = false;
            } else {
                throw std::runtime_error("--readFilesLegacyZcat must be Yes or No");
            }
        } else if (arg == "--dump-fastq") {
            opts.dump_fastq = true;
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    return opts;
}

InputSourcePlan build_plan(const HarnessOptions& opts) {
    std::vector<std::vector<std::string>> read_files_names;
    std::vector<std::string> read_groups;

    if (opts.read_files_manifest != "-") {
        std::ifstream manifest(opts.read_files_manifest.c_str());
        if (!manifest.good()) {
            throw std::runtime_error("could not open manifest: " + opts.read_files_manifest);
        }

        std::vector<std::string> mate1;
        std::vector<std::string> mate2;
        std::string line;
        while (std::getline(manifest, line)) {
            if (line.find_first_not_of(" \t") == std::string::npos) {
                continue;
            }
            const std::vector<std::string> fields = split_tab_line(line);
            if (fields.size() < 3) {
                throw std::runtime_error("manifest rows require Read1<TAB>Read2<TAB>ReadGroup");
            }
            mate1.push_back(opts.read_files_prefix + fields[0]);
            if (fields[1] != "-") {
                mate2.push_back(opts.read_files_prefix + fields[1]);
            }
            read_groups.push_back(read_group_id(fields[2]));
        }

        if (mate1.empty()) {
            throw std::runtime_error("manifest did not contain any reads");
        }
        read_files_names.push_back(mate1);
        if (!mate2.empty()) {
            if (mate2.size() != mate1.size()) {
                throw std::runtime_error("manifest mixes single-end and paired-end rows");
            }
            read_files_names.push_back(mate2);
        }
    } else {
        if (opts.read_files_in.empty()) {
            throw std::runtime_error("--readFilesIn or --readFilesManifest is required");
        }

        read_files_names.resize(opts.read_files_in.size());
        for (size_t imate = 0; imate < opts.read_files_in.size(); ++imate) {
            split_string(opts.read_files_in[imate], ',', read_files_names[imate]);
            if (read_files_names[imate].empty()) {
                throw std::runtime_error("empty mate file list in --readFilesIn");
            }
            if (imate > 0 && read_files_names[imate].size() != read_files_names[imate - 1].size()) {
                throw std::runtime_error("all mates must have the same number of lane files");
            }
            for (std::string& path : read_files_names[imate]) {
                path = opts.read_files_prefix + path;
            }
        }
    }

    bool all_gz = true;
    for (const auto& mate_files : read_files_names) {
        for (const std::string& path : mate_files) {
            all_gz = all_gz && ends_with_case_insensitive(path, ".gz");
        }
    }

    const bool has_command = !opts.read_files_command.empty();
    bool uses_internal_gzip = !has_command && all_gz && !opts.legacy_zcat;
    std::string command_string;
    if (has_command) {
        command_string = join_command(opts.read_files_command);
    } else if (uses_internal_gzip) {
        command_string = "INTERNAL_GZIP";
    } else if (all_gz && opts.legacy_zcat) {
        command_string = "zcat";
    } else if (!read_files_names.empty() && read_files_names.front().size() > 1) {
        command_string = "cat";
    }

    InputSourcePlan plan = star::input::make_fastx_input_source_plan(
        read_files_names,
        read_groups,
        command_string,
        opts.read_files_prefix,
        uses_internal_gzip);
    plan.read_name_separator_chars = opts.read_name_separator_chars;
    return plan;
}

void emit_tsv(const InputRecord& record, std::ostream& out) {
    for (uint32_t imate = 0; imate < record.mate_count; ++imate) {
        const auto& mate = record.mates[imate];
        out << record.lane_index << '\t'
            << record.read_ordinal << '\t'
            << record.read_filter << '\t'
            << record.read_name << '\t'
            << (imate + 1) << '\t'
            << mate.original_length << '\t'
            << sha256_hex(mate.sequence) << '\t'
            << sha256_hex(mate.quality) << '\n';
    }
}

void emit_fastq(const InputRecord& record, std::ostream& out) {
    for (uint32_t imate = 0; imate < record.mate_count; ++imate) {
        const auto& mate = record.mates[imate];
        out << '@' << record.read_name
            << " lane:" << record.lane_index
            << " ordinal:" << record.read_ordinal
            << " mate:" << (imate + 1)
            << " filter:" << record.read_filter << '\n'
            << mate.sequence << "\n+\n"
            << mate.quality << '\n';
    }
}

int run(const HarnessOptions& opts) {
    InputSourcePlan plan = build_plan(opts);
    FastxInputModule module;
    std::string error;
    if (!module.configure(plan, &error)) {
        std::cerr << "configure failed: " << error << "\n";
        return 2;
    }
    if (!module.open(&error)) {
        std::cerr << "open failed: " << error << "\n";
        return 2;
    }

    InputRecord record;
    while (true) {
        error.clear();
        const InputStatus status = module.next_record(&record, &error);
        if (status == InputStatus::End) {
            break;
        }
        if (status == InputStatus::Error) {
            std::cerr << "read failed: " << error << "\n";
            return 3;
        }
        if (opts.dump_fastq) {
            emit_fastq(record, std::cout);
        } else {
            emit_tsv(record, std::cout);
        }
    }
    module.close();
    return 0;
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        return run(parse_args(argc, argv));
    } catch (const std::exception& ex) {
        std::cerr << "ERROR: " << ex.what() << "\n\n";
        usage(std::cerr);
        return 1;
    }
}
