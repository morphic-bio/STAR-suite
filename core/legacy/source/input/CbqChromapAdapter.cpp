#include "input/CbqChromapAdapter.h"

#include "input/CbqInputModule.h"

#include <cerrno>
#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

namespace star {
namespace input {
namespace {

bool set_error(std::string* error, const std::string& message) {
    if (error != nullptr) {
        *error = message;
    }
    return false;
}

bool ensure_directory(const std::string& path, std::string* error) {
    struct stat st;
    if (stat(path.c_str(), &st) == 0) {
        if (S_ISDIR(st.st_mode)) {
            return true;
        }
        return set_error(error, path + " exists but is not a directory");
    }
    if (mkdir(path.c_str(), 0755) == 0) {
        return true;
    }
    std::ostringstream msg;
    msg << "failed to create directory " << path << ": " << std::strerror(errno);
    return set_error(error, msg.str());
}

std::string sanitize_sample_name(const std::string& sample_name) {
    if (sample_name.empty()) {
        return "sample";
    }
    std::string out;
    out.reserve(sample_name.size());
    for (char value : sample_name) {
        const unsigned char ch = static_cast<unsigned char>(value);
        if (std::isalnum(ch) || value == '.' || value == '_' || value == '-') {
            out.push_back(value);
        } else {
            out.push_back('_');
        }
    }
    return out.empty() ? "sample" : out;
}

std::string join_path(const std::string& dir, const std::string& name) {
    if (dir.empty() || dir[dir.size() - 1] == '/') {
        return dir + name;
    }
    return dir + "/" + name;
}

InputSourcePlan make_plan(const std::string& cbq_path, uint32_t mate_count) {
    std::vector<std::vector<std::string>> read_files(1);
    read_files[0].push_back(cbq_path);
    InputSourcePlan plan = make_cbq_input_source_plan(read_files, std::vector<std::string>(), mate_count);
    plan.read_name_separator_chars.clear();
    plan.read_name_separator_chars.push_back('/');
    plan.read_name_separator_chars.push_back(' ');
    return plan;
}

bool open_module(CbqInputModule* module,
                 const std::string& cbq_path,
                 uint32_t mate_count,
                 std::string* error) {
    InputSourcePlan plan = make_plan(cbq_path, mate_count);
    return module->configure(plan, error) && module->open(error);
}

std::string display_read_name(const InputRecord& record) {
    if (!record.read_name.empty()) {
        return record.read_name;
    }
    std::ostringstream name;
    name << "read_" << record.read_ordinal;
    return name.str();
}

bool write_fastq_record(std::ofstream* out,
                        const std::string& read_name,
                        const InputMateRecord& mate,
                        bool allow_missing_qualities,
                        std::string* error) {
    if (mate.sequence.empty()) {
        return set_error(error, "CBQ Chromap adapter received an empty sequence");
    }
    std::string quality = mate.quality;
    if (!mate.has_quality || quality.empty()) {
        if (!allow_missing_qualities) {
            return set_error(error, "CBQ Chromap adapter requires quality scores for FASTQ output");
        }
        quality.assign(mate.sequence.size(), 'I');
    }
    if (quality.size() != mate.sequence.size()) {
        std::ostringstream msg;
        msg << "CBQ Chromap adapter sequence/quality length mismatch for " << read_name
            << ": sequence=" << mate.sequence.size() << " quality=" << quality.size();
        return set_error(error, msg.str());
    }

    (*out) << '@' << read_name << '\n'
           << mate.sequence << "\n+\n"
           << quality << '\n';
    if (!out->good()) {
        return set_error(error, "failed while writing Chromap FASTQ record");
    }
    return true;
}

} // namespace

bool materialize_cbq_chromap_fastqs(const CbqChromapAdapterOptions& options,
                                    CbqChromapFastqPaths* paths,
                                    std::string* error) {
    if (options.read_pair_cbq.empty()) {
        return set_error(error, "read_pair_cbq is required");
    }
    if (options.barcode_cbq.empty()) {
        return set_error(error, "barcode_cbq is required");
    }
    if (options.output_dir.empty()) {
        return set_error(error, "output_dir is required");
    }
    if (paths == nullptr) {
        return set_error(error, "paths output pointer is required");
    }

    if (!ensure_directory(options.output_dir, error)) {
        return false;
    }

    const std::string sample = sanitize_sample_name(options.sample_name);
    paths->read1_fastq = join_path(options.output_dir, sample + ".chromap.R1.fastq");
    paths->read2_fastq = join_path(options.output_dir, sample + ".chromap.R2.fastq");
    paths->barcode_fastq = join_path(options.output_dir, sample + ".chromap.barcode.fastq");

    CbqInputModule reads;
    CbqInputModule barcodes;
    if (!open_module(&reads, options.read_pair_cbq, 2, error)) {
        return false;
    }
    if (!open_module(&barcodes, options.barcode_cbq, 1, error)) {
        reads.close();
        return false;
    }

    std::ofstream read1_out(paths->read1_fastq.c_str());
    std::ofstream read2_out(paths->read2_fastq.c_str());
    std::ofstream barcode_out(paths->barcode_fastq.c_str());
    if (!read1_out.good() || !read2_out.good() || !barcode_out.good()) {
        reads.close();
        barcodes.close();
        return set_error(error, "failed to open one or more Chromap FASTQ outputs");
    }

    uint64_t n_records = 0;
    for (;;) {
        InputRecord read_record;
        InputRecord barcode_record;
        const InputStatus read_status = reads.next_record(&read_record, error);
        if (read_status == InputStatus::Error) {
            reads.close();
            barcodes.close();
            return false;
        }
        const InputStatus barcode_status = barcodes.next_record(&barcode_record, error);
        if (barcode_status == InputStatus::Error) {
            reads.close();
            barcodes.close();
            return false;
        }
        if (read_status == InputStatus::End || barcode_status == InputStatus::End) {
            reads.close();
            barcodes.close();
            if (read_status == barcode_status) {
                return true;
            }
            return set_error(error, "paired-read CBQ and barcode CBQ have different record counts");
        }

        if (read_record.mates.size() < 2 || barcode_record.mates.empty()) {
            reads.close();
            barcodes.close();
            return set_error(error, "CBQ Chromap adapter expects paired reads plus one barcode mate");
        }
        if (options.require_name_match && read_record.read_name != barcode_record.read_name) {
            std::ostringstream msg;
            msg << "CBQ Chromap adapter read/barcode name mismatch at record "
                << (n_records + 1) << ": reads=" << read_record.read_name
                << " barcode=" << barcode_record.read_name;
            reads.close();
            barcodes.close();
            return set_error(error, msg.str());
        }

        const std::string name = display_read_name(read_record);
        if (!write_fastq_record(&read1_out, name, read_record.mates[0],
                                options.allow_missing_qualities, error) ||
            !write_fastq_record(&read2_out, name, read_record.mates[1],
                                options.allow_missing_qualities, error) ||
            !write_fastq_record(&barcode_out, name, barcode_record.mates[0],
                                options.allow_missing_qualities, error)) {
            reads.close();
            barcodes.close();
            return false;
        }
        ++n_records;
    }
}

} // namespace input
} // namespace star
