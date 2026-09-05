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
using star::input::CbqByteSpan;
using star::input::CbqReadBatchView;
using star::input::CbqReadView;
using star::input::CbqSegmentView;
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

std::string sha256_hex(CbqByteSpan value) {
    unsigned char digest[SHA256_DIGEST_LENGTH];
    const unsigned char* data = reinterpret_cast<const unsigned char*>("");
    if (value.data != nullptr && value.size != 0) {
        data = reinterpret_cast<const unsigned char*>(value.data);
    }
    SHA256(data, value.size, digest);

    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (unsigned char byte : digest) {
        out << std::setw(2) << static_cast<unsigned int>(byte);
    }
    return out.str();
}

std::string sha256_hex_string(const std::string& value) {
    CbqByteSpan span;
    span.data = value.data();
    span.size = value.size();
    return sha256_hex(span);
}

struct HarnessOptions {
    std::vector<std::string> read_files_in;
    std::string read_files_manifest = "-";
    std::string read_files_prefix;
    std::vector<char> read_name_separator_chars{' '};
    uint32_t mate_count = 2;
    bool dump_fastq = false;
    bool verify_packed_windows = false;
};

void usage(std::ostream& out) {
    out << "Usage: cbq_reader_harness [options]\n"
        << "  --readFilesIn sample.cbq[,lane2.cbq]\n"
        << "  --readFilesManifest manifest.tsv  # CBQ<TAB>-<TAB>ReadGroup\n"
        << "  --mateCount 1|2\n"
        << "  --readFilesPrefix PREFIX\n"
        << "  --readNameSeparator space|none|CHAR[,CHAR...]\n"
        << "  --verify-packed-windows\n"
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
        } else if (arg == "--verify-packed-windows") {
            opts.verify_packed_windows = true;
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    return opts;
}

bool verify_packed_windows(const CbqSegmentView& segment, std::string* error) {
    std::string sequence;
    materialize_cbq_segment_sequence(segment, &sequence);
    const size_t widths[] = {1, 8, 16, 32, 50, 64};
    for (size_t width : widths) {
        if (width > sequence.size()) {
            continue;
        }
        const size_t offsets[] = {0, (sequence.size() - width) / 2,
                                  sequence.size() - width};
        for (size_t offset : offsets) {
            if (width <= 32) {
                uint64_t packed = 0;
                uint32_t n_mask = 0;
                if (!star::input::cbq_pack_segment_window_lsb(
                        segment, offset, width, &packed, &n_mask)) {
                    *error = "LSB packed-window extraction failed";
                    return false;
                }
                uint64_t expected = 0;
                uint32_t expected_n_mask = 0;
                for (size_t ii = 0; ii < width; ++ii) {
                    const unsigned char code =
                        star::input::cbq_segment_base_star_number(segment, offset + ii);
                    if (code == 4) {
                        expected_n_mask |= static_cast<uint32_t>(1U << ii);
                    } else {
                        expected |= static_cast<uint64_t>(code) << (2U * ii);
                    }
                }
                if (packed != expected || n_mask != expected_n_mask) {
                    *error = "LSB packed-window result differs from decoded sequence";
                    return false;
                }
            }

            uint64_t lo = 0;
            uint64_t hi = 0;
            uint64_t n_mask = 0;
            const bool encoded = star::input::cbq_pack_segment_window_lsb_pair(
                segment, offset, width, &lo, &hi, &n_mask);
            uint64_t expected_lo = 0;
            uint64_t expected_hi = 0;
            uint64_t expected_n_mask = 0;
            for (size_t ii = 0; ii < width; ++ii) {
                const unsigned char code =
                    star::input::cbq_segment_base_star_number(segment, offset + ii);
                if (code == 4) {
                    expected_n_mask |= UINT64_C(1) << ii;
                } else if (ii < 32U) {
                    expected_lo |= static_cast<uint64_t>(code) << (2U * ii);
                } else {
                    expected_hi |= static_cast<uint64_t>(code) << (2U * (ii - 32U));
                }
            }
            if (!encoded || lo != expected_lo || hi != expected_hi ||
                n_mask != expected_n_mask) {
                *error = "paired LSB packed-window result differs from decoded sequence";
                return false;
            }
        }
    }
    return true;
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

void write_span(std::ostream& out, CbqByteSpan span) {
    if (span.data != nullptr && span.size != 0) {
        out.write(span.data, static_cast<std::streamsize>(span.size));
    }
}

void emit_tsv(const CbqReadView& record, std::ostream& out) {
    for (uint32_t isegment = 0; isegment < record.segment_count; ++isegment) {
        const CbqSegmentView& segment = record.segments[isegment];
        std::string sequence;
        materialize_cbq_segment_sequence(segment, &sequence);
        out << record.lane_index << '\t'
            << record.read_ordinal << '\t'
            << record.read_filter << '\t';
        write_span(out, record.read_name);
        out << '\t'
            << (segment.source_index + 1) << '\t'
            << cbq_segment_sequence_length(segment) << '\t'
            << sha256_hex_string(sequence) << '\t'
            << sha256_hex(segment.quality) << '\n';
    }
}

void emit_fastq(const CbqReadView& record, std::ostream& out) {
    for (uint32_t isegment = 0; isegment < record.segment_count; ++isegment) {
        const CbqSegmentView& segment = record.segments[isegment];
        std::string sequence;
        materialize_cbq_segment_sequence(segment, &sequence);
        out << '@';
        write_span(out, record.read_name);
        out
            << " lane:" << record.lane_index
            << " ordinal:" << record.read_ordinal
            << " mate:" << (segment.source_index + 1)
            << " filter:" << record.read_filter << '\n';
        out.write(sequence.data(), static_cast<std::streamsize>(sequence.size()));
        out << "\n+\n";
        write_span(out, segment.quality);
        out << '\n';
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

    CbqReadBatchView batch;
    while (true) {
        error.clear();
        const InputStatus status = module.next_batch(&batch, &error);
        if (status == InputStatus::End) {
            break;
        }
        if (status == InputStatus::Error) {
            std::cerr << "read failed: " << error << "\n";
            return 3;
        }
        for (uint32_t irecord = 0; irecord < batch.record_count; ++irecord) {
            const CbqReadView& record = batch.records[irecord];
            if (opts.verify_packed_windows) {
                for (uint32_t isegment = 0; isegment < record.segment_count; ++isegment) {
                    if (!verify_packed_windows(record.segments[isegment], &error)) {
                        std::cerr << "packed-window verification failed: " << error << "\n";
                        return 4;
                    }
                }
            }
            if (opts.dump_fastq) {
                emit_fastq(record, std::cout);
            } else {
                emit_tsv(record, std::cout);
            }
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
