#include "input/BgzfRangeReader.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <sstream>
#include <unistd.h>

namespace star {
namespace input {
namespace {

bool set_error(std::string* error, const std::string& message) {
    if (error != nullptr) {
        *error = message;
    }
    return false;
}

void strip_carriage_return(std::string* line) {
    if (!line->empty() && line->back() == '\r') {
        line->pop_back();
    }
}

} // namespace

BgzfRangeReader::BgzfRangeReader() = default;

BgzfRangeReader::~BgzfRangeReader() {
    close_input();
}

void BgzfRangeReader::close_input() {
    if (inputFd_ >= 0) {
        ::close(inputFd_);
        inputFd_ = -1;
    }
}

bool BgzfRangeReader::open(const BgzfIndex* index,
                           uint64_t first_record,
                           uint64_t record_count,
                           bool check_crc,
                           std::string* error) {
    close_input();
    index_ = nullptr;
    buffer_.clear();
    cursor_ = 0;
    recordsRead_ = 0;
    firstRecord_ = first_record;
    currentOrdinal_ = first_record;
    if (error != nullptr) {
        error->clear();
    }
    if (index == nullptr) {
        return set_error(error, "BGZF range reader requires a non-null index");
    }
    if (first_record > index->record_count() ||
        record_count > index->record_count() - first_record) {
        std::ostringstream message;
        message << "BGZF record range [" << first_record << ',';
        if (record_count > std::numeric_limits<uint64_t>::max() - first_record) {
            message << "overflow";
        } else {
            message << (first_record + record_count);
        }
        message << ") exceeds record count " << index->record_count();
        return set_error(error, message.str());
    }
    endRecord_ = first_record + record_count;
    index_ = index;
    checkCrc_ = check_crc;
    if (record_count == 0) {
        return true;
    }

    uint32_t decompressed_offset = 0;
    if (!index_->locate_record(first_record, &nextBlockIndex_,
                               &decompressed_offset, error)) {
        index_ = nullptr;
        return false;
    }
    currentBlockIndex_ = nextBlockIndex_;
    currentOrdinal_ = index_->blocks()[nextBlockIndex_].firstRecordOrdinal;
    inputFd_ = ::open(index_->path().c_str(), O_RDONLY);
    if (inputFd_ < 0) {
        const std::string message = "could not open BGZF input " + index_->path() +
                                    ": " + std::strerror(errno);
        index_ = nullptr;
        return set_error(error, message);
    }
    if (!append_block(error)) {
        close_input();
        index_ = nullptr;
        return false;
    }
    if (decompressed_offset > buffer_.size()) {
        close_input();
        index_ = nullptr;
        return set_error(error, "BGZF index decompressed record offset is outside its block");
    }
    cursor_ = decompressed_offset;

    BgzfFastqRecord skipped;
    while (currentOrdinal_ < firstRecord_) {
        if (!parse_record(&skipped, error)) {
            close_input();
            index_ = nullptr;
            return false;
        }
        ++currentOrdinal_;
    }
    return true;
}

bool BgzfRangeReader::append_block(std::string* error) {
    if (index_ == nullptr || inputFd_ < 0) {
        return set_error(error, "BGZF range reader is not open");
    }
    if (nextBlockIndex_ >= index_->blocks().size()) {
        return set_error(error, "unexpected end of BGZF blocks while completing a FASTQ record");
    }
    if (cursor_ != 0) {
        buffer_.erase(buffer_.begin(), buffer_.begin() + static_cast<std::ptrdiff_t>(cursor_));
        cursor_ = 0;
    }
    std::vector<unsigned char> inflated;
    currentBlockIndex_ = nextBlockIndex_;
    if (!inflate_bgzf_block_fd(inputFd_, index_->blocks()[nextBlockIndex_],
                               checkCrc_, &inflated, error)) {
        return false;
    }
    ++nextBlockIndex_;
    buffer_.insert(buffer_.end(), inflated.begin(), inflated.end());
    return true;
}

bool BgzfRangeReader::read_line(std::string* line, std::string* error) {
    if (line == nullptr) {
        return set_error(error, "BGZF FASTQ line destination is null");
    }
    while (true) {
        const std::vector<unsigned char>::iterator begin =
            buffer_.begin() + static_cast<std::ptrdiff_t>(cursor_);
        const std::vector<unsigned char>::iterator newline =
            std::find(begin, buffer_.end(), static_cast<unsigned char>('\n'));
        if (newline != buffer_.end()) {
            line->assign(reinterpret_cast<const char*>(&*begin),
                         static_cast<size_t>(newline - begin));
            cursor_ = static_cast<size_t>(newline - buffer_.begin()) + 1;
            strip_carriage_return(line);
            return true;
        }
        if (nextBlockIndex_ >= index_->blocks().size()) {
            if (cursor_ == buffer_.size()) {
                return set_error(error, "unexpected end of BGZF FASTQ record");
            }
            line->assign(reinterpret_cast<const char*>(buffer_.data() + cursor_),
                         buffer_.size() - cursor_);
            cursor_ = buffer_.size();
            strip_carriage_return(line);
            return true;
        }
        if (!append_block(error)) {
            return false;
        }
    }
}

bool BgzfRangeReader::parse_record(BgzfFastqRecord* record, std::string* error) {
    if (record == nullptr) {
        return set_error(error, "BGZF FASTQ record destination is null");
    }
    std::string header;
    if (!read_line(&header, error) ||
        !read_line(&record->sequence, error) ||
        !read_line(&record->plus, error) ||
        !read_line(&record->quality, error)) {
        return false;
    }
    if (header.empty() || header[0] != '@') {
        std::ostringstream message;
        message << "BGZF FASTQ record " << currentOrdinal_
                << " does not start with @ (block offset "
                << index_->blocks()[currentBlockIndex_].compressedOffset << ')';
        return set_error(error, message.str());
    }
    if (record->plus.empty() || record->plus[0] != '+') {
        std::ostringstream message;
        message << "BGZF FASTQ record " << currentOrdinal_
                << " has an invalid plus line (block offset "
                << index_->blocks()[currentBlockIndex_].compressedOffset << ')';
        return set_error(error, message.str());
    }
    if (record->sequence.size() != record->quality.size()) {
        std::ostringstream message;
        message << "BGZF FASTQ record " << currentOrdinal_
                << " has sequence length " << record->sequence.size()
                << " but quality length " << record->quality.size()
                << " (block offset "
                << index_->blocks()[currentBlockIndex_].compressedOffset << ')';
        return set_error(error, message.str());
    }
    record->name.assign(header.data() + 1, header.size() - 1);
    record->ordinal = currentOrdinal_;
    return true;
}

bool BgzfRangeReader::next(BgzfFastqRecord* record, std::string* error) {
    if (error != nullptr) {
        error->clear();
    }
    if (index_ == nullptr) {
        return set_error(error, "BGZF range reader is not open");
    }
    if (currentOrdinal_ >= endRecord_) {
        return false;
    }
    if (!parse_record(record, error)) {
        return false;
    }
    ++currentOrdinal_;
    ++recordsRead_;
    return true;
}

uint64_t BgzfRangeReader::first_record() const {
    return firstRecord_;
}

uint64_t BgzfRangeReader::end_record() const {
    return endRecord_;
}

uint64_t BgzfRangeReader::records_read() const {
    return recordsRead_;
}

} // namespace input
} // namespace star
