#include "input/BinseqInputProbeModule.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <sys/stat.h>
#include <sys/wait.h>
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

std::string shell_quote(const std::string& value) {
    std::string quoted = "'";
    for (char c : value) {
        if (c == '\'') {
            quoted += "'\\''";
        } else {
            quoted += c;
        }
    }
    quoted += "'";
    return quoted;
}

bool path_exists(const std::string& path) {
    struct stat st;
    return ::stat(path.c_str(), &st) == 0;
}

bool ensure_dir(const std::string& path, std::string* error) {
    if (path.empty()) {
        return set_error(error, "BINSEQ probe work directory is empty");
    }

    struct stat st;
    if (::stat(path.c_str(), &st) == 0) {
        if (S_ISDIR(st.st_mode)) {
            return true;
        }
        return set_error(error, "BINSEQ probe work path exists but is not a directory: " + path);
    }

    if (::mkdir(path.c_str(), 0775) != 0) {
        return set_error(error, "could not create BINSEQ probe work directory " + path + ": " +
                                std::strerror(errno));
    }
    return true;
}

std::string trim_trailing_slashes(std::string path) {
    while (path.size() > 1 && path.back() == '/') {
        path.pop_back();
    }
    return path;
}

std::string lane_prefix(const std::string& work_dir, uint32_t lane_index) {
    std::ostringstream out;
    out << trim_trailing_slashes(work_dir) << "/lane" << lane_index;
    return out.str();
}

std::string first_existing(const std::vector<std::string>& candidates) {
    for (const std::string& candidate : candidates) {
        if (path_exists(candidate)) {
            return candidate;
        }
    }
    return "";
}

std::string decode_output_for_mate(const std::string& prefix, uint32_t mate_index) {
    const uint32_t mate_number = mate_index + 1;
    std::ostringstream r;
    r << "R" << mate_number;
    std::ostringstream n;
    n << mate_number;

    return first_existing({
        prefix + "_" + r.str() + ".fastq",
        prefix + "." + r.str() + ".fastq",
        prefix + "_" + n.str() + ".fastq",
        prefix + "." + n.str() + ".fastq",
        prefix + "_" + r.str() + ".fq",
        prefix + "." + r.str() + ".fq",
        prefix + "_" + n.str() + ".fq",
        prefix + "." + n.str() + ".fq"
    });
}

std::string describe_system_status(int status) {
    std::ostringstream out;
    if (status == -1) {
        out << "could not launch command: " << std::strerror(errno);
    } else if (WIFEXITED(status)) {
        out << "exit " << WEXITSTATUS(status);
    } else if (WIFSIGNALED(status)) {
        out << "signal " << WTERMSIG(status);
    } else {
        out << "status " << status;
    }
    return out.str();
}

bool run_decode_command(const std::string& command,
                        const std::string& stdout_path,
                        const std::string& stderr_path,
                        std::string* error) {
    const std::string full_command = command + " > " + shell_quote(stdout_path) +
                                     " 2> " + shell_quote(stderr_path);
    const int status = std::system(full_command.c_str());
    if (status != 0) {
        std::ostringstream msg;
        msg << "bqtools decode failed (" << describe_system_status(status) << "): "
            << command << "\n"
            << "stdout: " << stdout_path << "\n"
            << "stderr: " << stderr_path;
        return set_error(error, msg.str());
    }
    return true;
}

} // namespace

BinseqInputProbeModule::BinseqInputProbeModule() {
    plan_.format = SourceFormat::Binseq;
    plan_.module_name = "BinseqProbe";
    configured_ = false;
}

BinseqInputProbeModule::~BinseqInputProbeModule() = default;

const char* BinseqInputProbeModule::name() const {
    return "BinseqInputProbeModule";
}

bool BinseqInputProbeModule::configure(const InputSourcePlan& plan, std::string* error) {
    if (plan.format != SourceFormat::Binseq) {
        return set_error(error, "BinseqInputProbeModule received a non-BINSEQ input source plan");
    }
    if (plan.mate_files.size() != 1) {
        return set_error(error, "BINSEQ probe currently expects one CBQ path per lane, not split external mate files");
    }
    if (plan.mate_files.front().empty()) {
        return set_error(error, "BINSEQ probe requires at least one CBQ input path");
    }
    if (plan.mate_count != 1 && plan.mate_count != 2) {
        return set_error(error, "BINSEQ probe currently supports mate_count 1 or 2");
    }
    if (plan.read_files_n != 0 && plan.read_files_n != plan.mate_files.front().size()) {
        std::ostringstream msg;
        msg << "BINSEQ probe read_files_n=" << plan.read_files_n
            << " does not match CBQ lane count=" << plan.mate_files.front().size();
        return set_error(error, msg.str());
    }
    if (!plan.read_groups.empty() && plan.read_groups.size() != plan.mate_files.front().size()) {
        std::ostringstream msg;
        msg << "BINSEQ probe has " << plan.mate_files.front().size()
            << " lanes but " << plan.read_groups.size() << " read groups";
        return set_error(error, msg.str());
    }
    if (plan.command_string.empty()) {
        return set_error(error, "BINSEQ probe requires plan.command_string to name the bqtools executable");
    }
    if (plan.prefix.empty()) {
        return set_error(error, "BINSEQ probe requires plan.prefix to name a temporary work directory");
    }

    plan_ = plan;
    plan_.format = SourceFormat::Binseq;
    plan_.module_name = "BinseqProbe";
    plan_.read_files_n = static_cast<uint32_t>(plan.mate_files.front().size());
    configured_ = true;
    return true;
}

const InputSourcePlan& BinseqInputProbeModule::plan() const {
    return plan_;
}

bool BinseqInputProbeModule::supports_record_iteration() const {
    return true;
}

bool BinseqInputProbeModule::open(std::string* error) {
    if (!configured_) {
        return set_error(error, "BinseqInputProbeModule must be configured before open()");
    }

    close();
    if (!ensure_dir(plan_.prefix, error)) {
        return false;
    }

    std::vector<std::vector<std::string>> decoded_mate_files(plan_.mate_count);
    for (uint32_t lane = 0; lane < plan_.read_files_n; ++lane) {
        if (!decode_lane(lane, decoded_mate_files, error)) {
            return false;
        }
    }

    InputSourcePlan fastx_plan = make_fastx_input_source_plan(
        decoded_mate_files,
        plan_.read_groups,
        "",
        "",
        false);
    fastx_plan.read_name_separator_chars = plan_.read_name_separator_chars;
    fastx_plan.out_qs_conversion_add = plan_.out_qs_conversion_add;

    fastx_.reset(new FastxInputModule());
    if (!fastx_->configure(fastx_plan, error)) {
        return false;
    }
    if (!fastx_->open(error)) {
        return false;
    }
    return true;
}

InputStatus BinseqInputProbeModule::next_record(InputRecord* record, std::string* error) {
    if (!fastx_) {
        if (error != nullptr) {
            *error = "BINSEQ probe module is not open";
        }
        return InputStatus::Error;
    }
    return fastx_->next_record(record, error);
}

void BinseqInputProbeModule::close() {
    if (fastx_) {
        fastx_->close();
        fastx_.reset();
    }
}

bool BinseqInputProbeModule::decode_lane(uint32_t lane_index,
                                         std::vector<std::vector<std::string>>& decoded_mate_files,
                                         std::string* error) {
    const std::string prefix = lane_prefix(plan_.prefix, lane_index);
    const std::string input_path = plan_.mate_files.front().at(lane_index);
    const std::string stdout_path = prefix + ".decode.stdout";
    const std::string stderr_path = prefix + ".decode.stderr";

    for (uint32_t imate = 0; imate < plan_.mate_count; ++imate) {
        const std::string existing = decode_output_for_mate(prefix, imate);
        if (!existing.empty()) {
            std::remove(existing.c_str());
        }
    }

    std::string command = shell_quote(plan_.command_string) + " decode " + shell_quote(input_path) + " -f q";
    if (plan_.mate_count == 1) {
        command += " -o " + shell_quote(prefix + "_R1.fastq");
    } else {
        command += " --prefix " + shell_quote(prefix);
    }

    if (!run_decode_command(command, stdout_path, stderr_path, error)) {
        return false;
    }

    for (uint32_t imate = 0; imate < plan_.mate_count; ++imate) {
        const std::string decoded_path = decode_output_for_mate(prefix, imate);
        if (decoded_path.empty()) {
            std::ostringstream msg;
            msg << "bqtools decode succeeded but expected mate " << (imate + 1)
                << " FASTQ output was not found for prefix " << prefix;
            return set_error(error, msg.str());
        }
        decoded_mate_files[imate].push_back(decoded_path);
    }
    return true;
}

} // namespace input
} // namespace star
