#include "input/BgzfRangeReader.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

namespace star {
namespace input {
namespace {

const size_t kMaxBlocksPerWork = 64;

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
    {
        std::lock_guard<std::mutex> lock(mutex_);
        stopping_ = true;
    }
    readyCv_.notify_all();
    spaceCv_.notify_all();
    for (size_t index = 0; index < workers_.size(); ++index) {
        if (workers_[index].joinable()) {
            workers_[index].join();
        }
    }
    workers_.clear();
    if (inputFd_ >= 0) {
        ::close(inputFd_);
        inputFd_ = -1;
    }
}

void BgzfRangeReader::fail_locked(const std::string& message) {
    failed_ = true;
    if (workerError_.empty()) {
        workerError_ = message;
    }
    readyCv_.notify_all();
    spaceCv_.notify_all();
}

bool BgzfRangeReader::open(const std::string& path,
                           uint64_t range_start,
                           uint64_t range_end,
                           uint32_t worker_threads,
                           bool check_crc,
                           std::string* error) {
    close_input();
    if (error != nullptr) {
        error->clear();
    }
    struct stat info;
    if (::stat(path.c_str(), &info) != 0) {
        return set_error(error, "could not stat BGZF input " + path + ": " +
                                std::strerror(errno));
    }
    if (info.st_size < 0) {
        return set_error(error, "BGZF input has a negative file size: " + path);
    }
    const uint64_t physical_end = static_cast<uint64_t>(info.st_size);
    if (range_end == std::numeric_limits<uint64_t>::max()) {
        range_end = physical_end;
    }
    if (range_start > range_end || range_end > physical_end) {
        std::ostringstream message;
        message << "BGZF compressed range [" << range_start << ',' << range_end
                << ") is outside file size " << physical_end;
        return set_error(error, message.str());
    }

    inputFd_ = ::open(path.c_str(), O_RDONLY);
    if (inputFd_ < 0) {
        return set_error(error, "could not open BGZF input " + path + ": " +
                                std::strerror(errno));
    }
    path_ = path;
    checkCrc_ = check_crc;
    rangeStart_ = range_start;
    rangeEnd_ = range_end;
    claimedOffset_ = range_start;
    nextClaimSequence_ = 0;
    nextConsumeSequence_ = 0;
    currentBlockOffset_ = range_start;
    recordsRead_ = 0;
    workerCount_ = worker_threads;
    const uint64_t range_bytes = range_end - range_start;
    const uint64_t planning_threads = std::max<uint32_t>(worker_threads, 1);
    const uint64_t planned = range_bytes / (planning_threads * 128U);
    targetCompressedBytes_ = std::max<uint64_t>(64U * 1024U,
        std::min<uint64_t>(1024U * 1024U, planned));
    maxOutstandingWork_ = std::max<size_t>(4, static_cast<size_t>(worker_threads) * 2);
    buffer_.clear();
    cursor_ = 0;
    completed_.clear();
    claimsFinished_ = range_start == range_end;
    failed_ = false;
    stopping_ = false;
    workerError_.clear();

    workers_.reserve(worker_threads);
    try {
        for (uint32_t worker = 0; worker < worker_threads; ++worker) {
            workers_.emplace_back(&BgzfRangeReader::worker_loop, this);
        }
    } catch (const std::exception& exc) {
        close_input();
        return set_error(error, std::string("could not start BGZF inflate worker: ") +
                                exc.what());
    }
    return true;
}

void BgzfRangeReader::worker_loop() {
    while (true) {
        std::vector<BgzfBlock> blocks;
        uint64_t sequence = 0;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            spaceCv_.wait(lock, [&]() {
                return stopping_ || failed_ || claimsFinished_ ||
                       nextClaimSequence_ - nextConsumeSequence_ < maxOutstandingWork_;
            });
            if (stopping_ || failed_ || claimsFinished_) {
                return;
            }
            std::string header_error;
            bool at_end = false;
            if (!claim_work(&blocks, &sequence, &at_end, &header_error)) {
                fail_locked(header_error);
                return;
            }
            if (at_end) {
                readyCv_.notify_all();
                spaceCv_.notify_all();
                return;
            }
        }

        InflatedBlock result;
        result.compressedOffset = blocks.front().compressedOffset;
        size_t output_bytes = 0;
        for (size_t index = 0; index < blocks.size(); ++index) {
            output_bytes += blocks[index].isize;
        }
        result.bytes.reserve(output_bytes);
        for (size_t index = 0; index < blocks.size(); ++index) {
            std::vector<unsigned char> inflated;
            std::string inflate_error;
            if (!inflate_bgzf_block_fd(inputFd_, blocks[index], checkCrc_, &inflated,
                                       &inflate_error)) {
                std::lock_guard<std::mutex> lock(mutex_);
                fail_locked(inflate_error);
                return;
            }
            result.bytes.insert(result.bytes.end(), inflated.begin(), inflated.end());
        }
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (stopping_ || failed_) {
                return;
            }
            completed_.emplace(sequence, std::move(result));
        }
        readyCv_.notify_all();
    }
}

bool BgzfRangeReader::claim_work(std::vector<BgzfBlock>* blocks,
                                 uint64_t* sequence,
                                 bool* at_end,
                                 std::string* error) {
    blocks->clear();
    *at_end = false;
    uint64_t compressed_bytes = 0;
    while (claimedOffset_ < rangeEnd_ &&
           blocks->size() < kMaxBlocksPerWork &&
           (blocks->empty() || compressed_bytes < targetCompressedBytes_)) {
        BgzfBlock block;
        bool is_eof_marker = false;
        if (!read_bgzf_block_header_fd(inputFd_, claimedOffset_, rangeEnd_, &block,
                                       &is_eof_marker, error)) {
            return false;
        }
        if (is_eof_marker) {
            claimedOffset_ = rangeEnd_;
            claimsFinished_ = true;
            break;
        }
        claimedOffset_ += block.compressedSize;
        compressed_bytes += block.compressedSize;
        blocks->push_back(block);
    }
    if (claimedOffset_ == rangeEnd_) {
        claimsFinished_ = true;
    }
    if (blocks->empty()) {
        *at_end = true;
        return true;
    }
    *sequence = nextClaimSequence_++;
    return true;
}

bool BgzfRangeReader::claim_and_inflate_sync(InflatedBlock* result,
                                              bool* at_end,
                                              std::string* error) {
    std::vector<BgzfBlock> blocks;
    uint64_t sequence = 0;
    if (!claim_work(&blocks, &sequence, at_end, error)) {
        return false;
    }
    if (*at_end) {
        return true;
    }
    result->compressedOffset = blocks.front().compressedOffset;
    result->bytes.clear();
    size_t output_bytes = 0;
    for (size_t index = 0; index < blocks.size(); ++index) {
        output_bytes += blocks[index].isize;
    }
    result->bytes.reserve(output_bytes);
    for (size_t index = 0; index < blocks.size(); ++index) {
        std::vector<unsigned char> inflated;
        if (!inflate_bgzf_block_fd(inputFd_, blocks[index], checkCrc_, &inflated, error)) {
            return false;
        }
        result->bytes.insert(result->bytes.end(), inflated.begin(), inflated.end());
    }
    return true;
}

bool BgzfRangeReader::append_next_block(std::string* error) {
    InflatedBlock block;
    if (workerCount_ == 0) {
        bool at_end = false;
        if (!claim_and_inflate_sync(&block, &at_end, error)) {
            return false;
        }
        if (at_end) {
            return false;
        }
    } else {
        std::unique_lock<std::mutex> lock(mutex_);
        readyCv_.wait(lock, [&]() {
            return stopping_ || failed_ ||
                   completed_.find(nextConsumeSequence_) != completed_.end() ||
                   (claimsFinished_ && nextConsumeSequence_ == nextClaimSequence_);
        });
        if (failed_) {
            return set_error(error, workerError_.empty()
                ? "BGZF inflate worker failed" : workerError_);
        }
        const std::map<uint64_t, InflatedBlock>::iterator found =
            completed_.find(nextConsumeSequence_);
        if (found == completed_.end()) {
            return false;
        }
        block = std::move(found->second);
        completed_.erase(found);
        ++nextConsumeSequence_;
        lock.unlock();
        spaceCv_.notify_all();
    }

    if (cursor_ != 0) {
        buffer_.erase(buffer_.begin(),
                      buffer_.begin() + static_cast<std::ptrdiff_t>(cursor_));
        cursor_ = 0;
    }
    currentBlockOffset_ = block.compressedOffset;
    buffer_.insert(buffer_.end(), block.bytes.begin(), block.bytes.end());
    return true;
}

bool BgzfRangeReader::read_line(std::string* line,
                                bool allow_clean_end,
                                std::string* error) {
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

        std::string append_error;
        if (append_next_block(&append_error)) {
            continue;
        }
        if (!append_error.empty()) {
            return set_error(error, append_error);
        }
        if (cursor_ != buffer_.size()) {
            line->assign(reinterpret_cast<const char*>(buffer_.data() + cursor_),
                         buffer_.size() - cursor_);
            cursor_ = buffer_.size();
            strip_carriage_return(line);
            return true;
        }
        if (allow_clean_end) {
            return false;
        }
        std::ostringstream message;
        message << "unexpected end of BGZF FASTQ record after block offset "
                << currentBlockOffset_;
        return set_error(error, message.str());
    }
}

bool BgzfRangeReader::parse_record(BgzfFastqRecord* record, std::string* error) {
    if (record == nullptr) {
        return set_error(error, "BGZF FASTQ record destination is null");
    }
    std::string header;
    if (!read_line(&header, true, error)) {
        return false;
    }
    if (!read_line(&record->sequence, false, error) ||
        !read_line(&record->plus, false, error) ||
        !read_line(&record->quality, false, error)) {
        return false;
    }
    if (header.empty() || header[0] != '@') {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_
                << " does not start with @ (block offset " << currentBlockOffset_ << ')';
        return set_error(error, message.str());
    }
    if (record->plus.empty() || record->plus[0] != '+') {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_
                << " has an invalid plus line (block offset " << currentBlockOffset_ << ')';
        return set_error(error, message.str());
    }
    if (record->sequence.size() != record->quality.size()) {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_ << " has sequence length "
                << record->sequence.size() << " but quality length "
                << record->quality.size() << " (block offset "
                << currentBlockOffset_ << ')';
        return set_error(error, message.str());
    }
    record->name.assign(header.data() + 1, header.size() - 1);
    record->ordinal = recordsRead_;
    return true;
}

bool BgzfRangeReader::next(BgzfFastqRecord* record, std::string* error) {
    if (error != nullptr) {
        error->clear();
    }
    if (inputFd_ < 0) {
        return set_error(error, "BGZF range reader is not open");
    }
    if (!parse_record(record, error)) {
        return false;
    }
    ++recordsRead_;
    return true;
}

uint64_t BgzfRangeReader::records_read() const {
    return recordsRead_;
}

uint64_t BgzfRangeReader::range_start() const {
    return rangeStart_;
}

uint64_t BgzfRangeReader::range_end() const {
    return rangeEnd_;
}

} // namespace input
} // namespace star
