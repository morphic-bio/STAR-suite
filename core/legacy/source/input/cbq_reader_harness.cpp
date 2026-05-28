#include "input/CbqInputModule.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <openssl/sha.h>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

using star::input::CbqInputModule;
using star::input::InputRecord;
using star::input::InputSourcePlan;
using star::input::InputStatus;
using star::input::SourceFormat;

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

struct HarnessOptions {
    std::vector<std::string> read_files_in;
    std::string read_files_manifest = "-";
    std::string read_files_prefix;
    std::vector<char> read_name_separator_chars{' '};
    uint32_t mate_count = 2;
    bool dump_fastq = false;
};

void usage(std::ostream& out) {
    out << "Usage: cbq_reader_harness [options]\n"
        << "  --readFilesIn sample.cbq[,lane2.cbq]\n"
        << "  --readFilesManifest manifest.tsv  # CBQ<TAB>-<TAB>ReadGroup\n"
        << "  --mateCount 1|2\n"
        << "  --readFilesPrefix PREFIX\n"
        << "  --readNameSeparator space|none|CHAR[,CHAR...]\n"
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
            if (++ii >= argc || starts_with_option(argv[ii])) {
                throw std::runtime_error("--readFilesIn requires one comma-separated CBQ lane list");
            }
            opts.read_files_in.push_back(argv[ii]);
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
        } else if (arg == "--mateCount") {
            if (++ii >= argc) {
                throw std::runtime_error("--mateCount requires a value");
            }
            opts.mate_count = static_cast<uint32_t>(std::stoul(argv[ii]));
            if (opts.mate_count != 1 && opts.mate_count != 2) {
                throw std::runtime_error("--mateCount must be 1 or 2");
            }
        } else if (arg == "--readNameSeparator") {
            if (++ii >= argc) {
                throw std::runtime_error("--readNameSeparator requires a value");
            }
            opts.read_name_separator_chars = parse_separators(argv[ii]);
        } else if (arg == "--dump-fastq") {
            opts.dump_fastq = true;
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    return opts;
}

InputSourcePlan build_plan_from_manifest(const HarnessOptions& opts) {
    std::ifstream manifest(opts.read_files_manifest.c_str());
    if (!manifest.good()) {
        throw std::runtime_error("could not open manifest: " + opts.read_files_manifest);
    }

    std::vector<std::string> cbq_files;
    std::vector<std::string> read_groups;
    std::string line;
    while (std::getline(manifest, line)) {
        if (line.find_first_not_of(" \t") == std::string::npos) {
            continue;
        }
        const std::vector<std::string> fields = split_tab_line(line);
        if (fields.size() < 3) {
            throw std::runtime_error("CBQ manifest rows require CBQ<TAB>-<TAB>ReadGroup");
        }
        if (fields[1] != "-") {
            throw std::runtime_error("CBQ manifest second column must be '-' for records packed in one CBQ");
        }
        cbq_files.push_back(opts.read_files_prefix + fields[0]);
        read_groups.push_back(read_group_id(fields[2]));
    }
    if (cbq_files.empty()) {
        throw std::runtime_error("manifest did not contain any CBQ files");
    }

    InputSourcePlan plan;
    plan.format = SourceFormat::Binseq;
    plan.module_name = "Cbq";
    plan.mate_files.push_back(cbq_files);
    plan.read_groups = read_groups;
    plan.read_files_n = static_cast<uint32_t>(cbq_files.size());
    plan.mate_count = opts.mate_count;
    plan.read_name_separator_chars = opts.read_name_separator_chars;
    return plan;
}

InputSourcePlan build_plan(const HarnessOptions& opts) {
    if (opts.read_files_manifest != "-") {
        return build_plan_from_manifest(opts);
    }
    if (opts.read_files_in.size() != 1) {
        throw std::runtime_error("--readFilesIn requires one comma-separated CBQ lane list");
    }

    std::vector<std::string> cbq_files;
    split_string(opts.read_files_in.front(), ',', cbq_files);
    if (cbq_files.empty()) {
        throw std::runtime_error("--readFilesIn did not contain any CBQ files");
    }
    for (std::string& path : cbq_files) {
        path = opts.read_files_prefix + path;
    }

    InputSourcePlan plan;
    plan.format = SourceFormat::Binseq;
    plan.module_name = "Cbq";
    plan.mate_files.push_back(cbq_files);
    plan.read_files_n = static_cast<uint32_t>(cbq_files.size());
    plan.mate_count = opts.mate_count;
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
    CbqInputModule module;
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
