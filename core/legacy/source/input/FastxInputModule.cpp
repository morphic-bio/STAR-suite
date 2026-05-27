#include "input/FastxInputModule.h"

#include "IncludeDefine.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <sys/wait.h>
#include <zlib.h>

namespace star {
namespace input {
namespace {

bool set_error(std::string* error, const std::string& message) {
    if (error != nullptr) {
        *error = message;
    }
    return false;
}

std::string strip_line_ending(std::string line) {
    if (!line.empty() && line.back() == '\n') {
        line.pop_back();
    }
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    return line;
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

std::string trim_command(const std::string& command) {
    const size_t first = command.find_first_not_of(" \t");
    if (first == std::string::npos) {
        return "";
    }
    const size_t last = command.find_last_not_of(" \t");
    return command.substr(first, last - first + 1);
}

std::string first_field(const std::string& value) {
    std::istringstream in(value);
    std::string field;
    in >> field;
    return field;
}

bool starts_with_header(const std::string& line) {
    return !line.empty() && (line[0] == '@' || line[0] == '>');
}

class LineReader {
public:
    virtual ~LineReader() {}
    virtual bool open(const std::string& source, std::string* error) = 0;
    virtual bool getline(std::string& line, std::string* error) = 0;
    virtual void close() = 0;
};

class PlainLineReader : public LineReader {
public:
    bool open(const std::string& source, std::string* error) override {
        stream_.open(source.c_str());
        if (!stream_.good()) {
            return set_error(error, "could not open FASTX file: " + source);
        }
        return true;
    }

    bool getline(std::string& line, std::string* error) override {
        if (std::getline(stream_, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            return true;
        }
        if (stream_.bad()) {
            return set_error(error, "error reading FASTX file stream");
        }
        return false;
    }

    void close() override {
        if (stream_.is_open()) {
            stream_.close();
        }
    }

private:
    std::ifstream stream_;
};

class GzipLineReader : public LineReader {
public:
    GzipLineReader() : file_(nullptr) {}

    ~GzipLineReader() override {
        close();
    }

    bool open(const std::string& source, std::string* error) override {
        file_ = gzopen(source.c_str(), "rb");
        if (file_ == nullptr) {
            return set_error(error, "could not open gzip FASTX file: " + source);
        }
        source_ = source;
        return true;
    }

    bool getline(std::string& line, std::string* error) override {
        line.clear();
        int c = 0;
        while ((c = gzgetc(file_)) != -1) {
            if (c == '\n') {
                return true;
            }
            if (c != '\r') {
                line.push_back(static_cast<char>(c));
            }
        }

        if (!line.empty()) {
            return true;
        }

        if (!gzeof(file_)) {
            int zerr = Z_OK;
            const char* zmsg = gzerror(file_, &zerr);
            if (zerr != Z_OK && zerr != Z_STREAM_END) {
                return set_error(error, "error reading gzip FASTX file " + source_ + ": " +
                                      (zmsg != nullptr ? std::string(zmsg) : std::string("unknown zlib error")));
            }
        }
        return false;
    }

    void close() override {
        if (file_ != nullptr) {
            gzclose(file_);
            file_ = nullptr;
        }
    }

private:
    gzFile file_;
    std::string source_;
};

class CommandLineReader : public LineReader {
public:
    CommandLineReader() : pipe_(nullptr), buffer_(nullptr), buffer_size_(0) {}

    ~CommandLineReader() override {
        close();
        std::free(buffer_);
    }

    bool open(const std::string& source, std::string* error) override {
        command_ = source;
        pipe_ = popen(command_.c_str(), "r");
        if (pipe_ == nullptr) {
            return set_error(error, "could not execute FASTX read command: " + command_ +
                                  ": " + std::strerror(errno));
        }
        return true;
    }

    bool getline(std::string& line, std::string* error) override {
        errno = 0;
        const ssize_t nread = ::getline(&buffer_, &buffer_size_, pipe_);
        if (nread >= 0) {
            line.assign(buffer_, static_cast<size_t>(nread));
            line = strip_line_ending(line);
            return true;
        }
        if (errno != 0) {
            return set_error(error, "error reading FASTX command output: " + command_ +
                                  ": " + std::strerror(errno));
        }
        const int status = pclose(pipe_);
        pipe_ = nullptr;
        if (status != 0) {
            std::ostringstream msg;
            msg << "FASTX command failed: " << command_;
            if (WIFEXITED(status)) {
                msg << " (exit " << WEXITSTATUS(status) << ")";
            } else if (WIFSIGNALED(status)) {
                msg << " (signal " << WTERMSIG(status) << ")";
            }
            return set_error(error, msg.str());
        }
        return false;
    }

    void close() override {
        if (pipe_ != nullptr) {
            pclose(pipe_);
            pipe_ = nullptr;
        }
    }

private:
    FILE* pipe_;
    char* buffer_;
    size_t buffer_size_;
    std::string command_;
};

struct ParsedMate {
    std::string read_name;
    std::string read_name_extra;
    char read_filter = 'N';
    InputMateRecord mate;
};

std::string canonical_read_name(const std::string& raw_name, const std::vector<char>& separators) {
    std::string name = raw_name;
    for (char separator : separators) {
        const size_t pos = name.find(separator);
        if (pos != std::string::npos) {
            name.resize(pos);
        }
    }
    return name;
}

char parse_read_filter(char record_type, const std::string& read_name_extra) {
    if (record_type != '@') {
        return 'N';
    }

    const std::string field2 = first_field(read_name_extra);
    if (field2.size() > 3 && field2[1] == ':' && field2[2] == 'Y' && field2[3] == ':') {
        return 'Y';
    }
    return 'N';
}

void apply_quality_conversion(std::string& quality, int add) {
    if (add == 0) {
        return;
    }
    for (char& qchar : quality) {
        int qs = static_cast<int>(qchar) + add;
        if (qs < 33) {
            qs = 33;
        } else if (qs > 126) {
            qs = 126;
        }
        qchar = static_cast<char>(qs);
    }
}

InputStatus parse_one_mate(
    LineReader& reader,
    std::string& pending_header,
    const std::vector<char>& read_name_separator_chars,
    int out_qs_conversion_add,
    ParsedMate& parsed,
    std::string* error) {

    std::string header;
    if (!pending_header.empty()) {
        header = pending_header;
        pending_header.clear();
    } else if (!reader.getline(header, error)) {
        return (error != nullptr && !error->empty()) ? InputStatus::Error : InputStatus::End;
    }

    if (!starts_with_header(header)) {
        return set_error(error, "FASTX record header does not start with @ or >: " + header)
            ? InputStatus::Record : InputStatus::Error;
    }

    const char record_type = header[0];
    std::string payload = header.substr(1);
    std::istringstream payload_stream(payload);
    std::string raw_name;
    payload_stream >> raw_name;
    payload_stream >> std::ws;
    std::string extra;
    std::getline(payload_stream, extra);

    std::string sequence;
    if (!reader.getline(sequence, error)) {
        return set_error(error, "truncated FASTX record after header: " + header)
            ? InputStatus::Record : InputStatus::Error;
    }

    std::string quality;
    bool has_quality = false;
    if (record_type == '@') {
        std::string plus_line;
        if (!reader.getline(plus_line, error) || plus_line.empty() || plus_line[0] != '+') {
            return set_error(error, "FASTQ record is missing + line for read: " + raw_name)
                ? InputStatus::Record : InputStatus::Error;
        }
        if (!reader.getline(quality, error)) {
            return set_error(error, "FASTQ record is missing quality line for read: " + raw_name)
                ? InputStatus::Record : InputStatus::Error;
        }
        if (quality.size() != sequence.size()) {
            std::ostringstream msg;
            msg << "FASTQ sequence/quality length mismatch for read " << raw_name
                << ": sequence=" << sequence.size() << " quality=" << quality.size();
            return set_error(error, msg.str()) ? InputStatus::Record : InputStatus::Error;
        }
        has_quality = true;
        apply_quality_conversion(quality, out_qs_conversion_add);
    } else {
        std::string line;
        while (reader.getline(line, error)) {
            if (starts_with_header(line)) {
                pending_header = line;
                break;
            }
            sequence += line;
        }
        if (error != nullptr && !error->empty()) {
            return InputStatus::Error;
        }
        quality.assign(sequence.size(), 'A');
    }

    parsed.read_name = canonical_read_name(raw_name, read_name_separator_chars);
    parsed.read_name_extra = extra;
    parsed.read_filter = parse_read_filter(record_type, extra);
    parsed.mate.sequence = sequence;
    parsed.mate.quality = quality;
    parsed.mate.original_length = static_cast<uint32_t>(sequence.size());
    parsed.mate.has_quality = has_quality;
    return InputStatus::Record;
}

} // namespace

struct FastxInputModule::Impl {
    std::vector<std::unique_ptr<LineReader>> readers;
    std::vector<std::string> pending_headers;
    uint32_t lane_index = 0;
    uint64_t read_ordinal = 0;
    bool opened = false;
};

InputSourcePlan make_fastx_input_source_plan(
    const std::vector<std::vector<std::string>>& read_files_names,
    const std::vector<std::string>& read_groups,
    const std::string& command_string,
    const std::string& prefix,
    bool uses_internal_gzip) {

    InputSourcePlan plan;
    plan.format = SourceFormat::Fastx;
    plan.module_name = "Fastx";
    plan.mate_files = read_files_names;
    plan.read_groups = read_groups;
    plan.read_files_n = read_files_names.empty() ? 0 : static_cast<uint32_t>(read_files_names.front().size());
    plan.mate_count = static_cast<uint32_t>(read_files_names.size());
    plan.uses_helper_stream = !command_string.empty();
    plan.uses_internal_gzip = uses_internal_gzip;
    plan.command_string = command_string;
    plan.prefix = prefix;
    return plan;
}

FastxInputModule::FastxInputModule() {
    plan_.format = SourceFormat::Fastx;
    plan_.module_name = "Fastx";
    configured_ = false;
    impl_.reset(new Impl());
}

FastxInputModule::~FastxInputModule() = default;

const char* FastxInputModule::name() const {
    return "FastxInputModule";
}

bool FastxInputModule::configure(const InputSourcePlan& plan, std::string* error) {
    if (plan.format != SourceFormat::Fastx) {
        return set_error(error, "FastxInputModule received a non-Fastx input source plan");
    }

    if (plan.mate_files.empty()) {
        return set_error(error, "Fastx input requires at least one read mate");
    }

    if (plan.mate_files.size() > MAX_N_MATES) {
        std::ostringstream msg;
        msg << "Fastx input has " << plan.mate_files.size()
            << " read mates, which exceeds MAX_N_MATES=" << MAX_N_MATES;
        return set_error(error, msg.str());
    }

    const size_t read_files_n = plan.mate_files.front().size();
    if (read_files_n == 0) {
        return set_error(error, "Fastx input requires at least one file per mate");
    }

    for (size_t imate = 0; imate < plan.mate_files.size(); ++imate) {
        if (plan.mate_files[imate].size() != read_files_n) {
            std::ostringstream msg;
            msg << "Fastx mate " << (imate + 1) << " has "
                << plan.mate_files[imate].size() << " files, expected "
                << read_files_n;
            return set_error(error, msg.str());
        }

        for (size_t ifile = 0; ifile < plan.mate_files[imate].size(); ++ifile) {
            if (plan.mate_files[imate][ifile].empty()) {
                std::ostringstream msg;
                msg << "Fastx mate " << (imate + 1)
                    << " has an empty file path at lane " << (ifile + 1);
                return set_error(error, msg.str());
            }
        }
    }

    if (plan.read_files_n != 0 && plan.read_files_n != read_files_n) {
        std::ostringstream msg;
        msg << "Fastx input source plan read_files_n=" << plan.read_files_n
            << " does not match mate file count=" << read_files_n;
        return set_error(error, msg.str());
    }

    if (!plan.read_groups.empty() && plan.read_groups.size() != read_files_n) {
        std::ostringstream msg;
        msg << "Fastx input has " << read_files_n
            << " files per mate but " << plan.read_groups.size()
            << " read groups";
        return set_error(error, msg.str());
    }

    plan_ = plan;
    plan_.module_name = "Fastx";
    plan_.read_files_n = static_cast<uint32_t>(read_files_n);
    plan_.mate_count = static_cast<uint32_t>(plan.mate_files.size());
    if (plan_.read_name_separator_chars.empty()) {
        plan_.read_name_separator_chars.push_back(' ');
    }
    configured_ = true;
    return true;
}

const InputSourcePlan& FastxInputModule::plan() const {
    return plan_;
}

bool FastxInputModule::supports_record_iteration() const {
    return true;
}

bool FastxInputModule::open(std::string* error) {
    if (!configured_) {
        return set_error(error, "FastxInputModule must be configured before open()");
    }
    close();
    impl_->lane_index = 0;
    impl_->read_ordinal = 0;
    impl_->opened = false;

    if (plan_.read_files_n == 0) {
        return set_error(error, "Fastx input source plan has no lanes");
    }

    return open_lane(impl_->lane_index, error);
}

InputStatus FastxInputModule::next_record(InputRecord* record, std::string* error) {
    if (record == nullptr) {
        if (error != nullptr) {
            *error = "next_record() requires a non-null InputRecord";
        }
        return InputStatus::Error;
    }
    if (!impl_->opened) {
        if (error != nullptr) {
            *error = "FastxInputModule is not open";
        }
        return InputStatus::Error;
    }

    while (impl_->lane_index < plan_.read_files_n) {
        std::vector<ParsedMate> parsed(plan_.mate_count);
        std::vector<InputStatus> statuses(plan_.mate_count, InputStatus::End);

        bool saw_record = false;
        bool saw_end = false;
        for (uint32_t imate = 0; imate < plan_.mate_count; ++imate) {
            statuses[imate] = parse_one_mate(
                *impl_->readers[imate],
                impl_->pending_headers[imate],
                plan_.read_name_separator_chars,
                plan_.out_qs_conversion_add,
                parsed[imate],
                error);

            if (statuses[imate] == InputStatus::Error) {
                return InputStatus::Error;
            }
            saw_record = saw_record || statuses[imate] == InputStatus::Record;
            saw_end = saw_end || statuses[imate] == InputStatus::End;
        }

        if (!saw_record) {
            close_lane();
            ++impl_->lane_index;
            if (impl_->lane_index >= plan_.read_files_n) {
                impl_->opened = false;
                return InputStatus::End;
            }
            if (!open_lane(impl_->lane_index, error)) {
                return InputStatus::Error;
            }
            continue;
        }

        if (saw_end) {
            if (error != nullptr) {
                *error = "FASTX mate streams ended at different records";
            }
            return InputStatus::Error;
        }

        record->read_name = parsed[0].read_name;
        record->read_name_extra = parsed[0].read_name_extra;
        record->lane_index = impl_->lane_index;
        record->read_ordinal = ++impl_->read_ordinal;
        record->read_filter = parsed[0].read_filter;
        record->mate_count = plan_.mate_count;
        record->mates.clear();
        record->mates.reserve(plan_.mate_count);
        for (const ParsedMate& mate : parsed) {
            record->mates.push_back(mate.mate);
        }
        return InputStatus::Record;
    }

    impl_->opened = false;
    return InputStatus::End;
}

void FastxInputModule::close() {
    close_lane();
    if (impl_) {
        impl_->lane_index = 0;
        impl_->opened = false;
    }
}

bool FastxInputModule::open_lane(uint32_t lane_index, std::string* error) {
    close_lane();
    impl_->readers.resize(plan_.mate_count);
    impl_->pending_headers.assign(plan_.mate_count, std::string());

    for (uint32_t imate = 0; imate < plan_.mate_count; ++imate) {
        const std::string& path = plan_.mate_files[imate][lane_index];
        if (plan_.uses_internal_gzip) {
            impl_->readers[imate].reset(new GzipLineReader());
            if (!impl_->readers[imate]->open(path, error)) {
                close_lane();
                return false;
            }
        } else if (!plan_.command_string.empty()) {
            const std::string command = trim_command(plan_.command_string);
            if (command.empty()) {
                return set_error(error, "FASTX command string is empty");
            }
            impl_->readers[imate].reset(new CommandLineReader());
            if (!impl_->readers[imate]->open(command + " " + shell_quote(path), error)) {
                close_lane();
                return false;
            }
        } else {
            impl_->readers[imate].reset(new PlainLineReader());
            if (!impl_->readers[imate]->open(path, error)) {
                close_lane();
                return false;
            }
        }
    }

    impl_->opened = true;
    return true;
}

void FastxInputModule::close_lane() {
    if (!impl_) {
        return;
    }
    for (auto& reader : impl_->readers) {
        if (reader) {
            reader->close();
        }
    }
    impl_->readers.clear();
    impl_->pending_headers.clear();
}

} // namespace input
} // namespace star
