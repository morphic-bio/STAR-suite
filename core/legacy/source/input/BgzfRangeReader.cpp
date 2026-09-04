#include "input/BgzfRangeReader.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
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
                           std::string* error,
                           const BgzfWorkPermitHooks* permit_hooks,
                           bool store_quality) {
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
    if (permit_hooks != nullptr &&
        ((permit_hooks->acquire == nullptr) !=
         (permit_hooks->release == nullptr))) {
        return set_error(error,
                         "BGZF inflate permit hooks require both acquire and release");
    }

    inputFd_ = ::open(path.c_str(), O_RDONLY);
    if (inputFd_ < 0) {
        return set_error(error, "could not open BGZF input " + path + ": " +
                                std::strerror(errno));
    }
    path_ = path;
    checkCrc_ = check_crc;
    storeQuality_ = store_quality;
    rangeStart_ = range_start;
    rangeEnd_ = range_end;
    claimedOffset_ = range_start;
    nextClaimSequence_ = 0;
    nextConsumeSequence_ = 0;
    currentBlockOffset_ = range_start;
    recordsRead_ = 0;
    workerCount_ = worker_threads;
    permitHooks_ = permit_hooks == nullptr
        ? BgzfWorkPermitHooks() : *permit_hooks;
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
    BgzfInflater inflater;
    CompressedWork work;
    while (true) {
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
            if (!claim_work(&work, &sequence, &at_end, &header_error)) {
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
        std::string inflate_error;
        if (!inflate_work_permitted(&inflater, work, &result, &inflate_error)) {
            std::lock_guard<std::mutex> lock(mutex_);
            fail_locked(inflate_error);
            return;
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

bool BgzfRangeReader::claim_work(CompressedWork* work,
                                 uint64_t* sequence,
                                 bool* at_end,
                                 std::string* error) {
    work->compressedOffset = claimedOffset_;
    work->blocks.clear();
    work->bytes.clear();
    *at_end = false;
    bool reached_end = false;
    if (!read_bgzf_work_fd(inputFd_, claimedOffset_, rangeEnd_,
                           targetCompressedBytes_, kMaxBlocksPerWork,
                           &work->blocks, &work->bytes, &reached_end, error)) {
        return false;
    }
    if (reached_end) {
        claimedOffset_ = rangeEnd_;
        claimsFinished_ = true;
    } else {
        claimedOffset_ += work->bytes.size();
    }
    if (work->blocks.empty()) {
        *at_end = true;
        return true;
    }
    *sequence = nextClaimSequence_++;
    return true;
}

bool BgzfRangeReader::inflate_work(BgzfInflater* inflater,
                                   const CompressedWork& work,
                                   InflatedBlock* result,
                                   std::string* error) {
    if (inflater == nullptr || result == nullptr || work.blocks.empty()) {
        return set_error(error, "invalid BGZF compressed work item");
    }
    result->compressedOffset = work.blocks.front().compressedOffset;
    size_t output_bytes = 0;
    for (size_t index = 0; index < work.blocks.size(); ++index) {
        if (work.blocks[index].isize >
            std::numeric_limits<size_t>::max() - output_bytes) {
            return set_error(error, "BGZF work output exceeds platform memory limits");
        }
        output_bytes += work.blocks[index].isize;
    }
    result->bytes.resize(output_bytes);
    size_t output_offset = 0;
    for (size_t index = 0; index < work.blocks.size(); ++index) {
        const BgzfBlock& block = work.blocks[index];
        if (block.compressedOffset < work.compressedOffset) {
            return set_error(error, "BGZF member precedes its claimed work range");
        }
        const uint64_t relative64 = block.compressedOffset - work.compressedOffset;
        if (relative64 > work.bytes.size() ||
            block.compressedSize > work.bytes.size() - static_cast<size_t>(relative64)) {
            return set_error(error, "BGZF member exceeds its claimed work range");
        }
        unsigned char* destination = block.isize == 0
            ? nullptr : result->bytes.data() + output_offset;
        if (!inflater->inflate_block(
                work.bytes.data() + static_cast<size_t>(relative64),
                block.compressedSize, block.compressedOffset, checkCrc_,
                destination, block.isize, error)) {
            return false;
        }
        output_offset += block.isize;
    }
    return true;
}

bool BgzfRangeReader::inflate_work_permitted(BgzfInflater* inflater,
                                             const CompressedWork& work,
                                             InflatedBlock* result,
                                             std::string* error) {
    if (!permitHooks_.enabled()) {
        return inflate_work(inflater, work, result, error);
    }

    const uint64_t wait_ns = permitHooks_.acquire(permitHooks_.context);
    const std::chrono::steady_clock::time_point work_start =
        std::chrono::steady_clock::now();
    const bool ok = inflate_work(inflater, work, result, error);
    const uint64_t work_ns = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - work_start).count());
    permitHooks_.release(permitHooks_.context,
                         wait_ns,
                         static_cast<uint64_t>(work.blocks.size()),
                         static_cast<uint64_t>(work.bytes.size()),
                         work_ns);
    return ok;
}

bool BgzfRangeReader::claim_and_inflate_sync(InflatedBlock* result,
                                              bool* at_end,
                                              std::string* error) {
    CompressedWork work;
    uint64_t sequence = 0;
    if (!claim_work(&work, &sequence, at_end, error)) {
        return false;
    }
    if (*at_end) {
        return true;
    }
    return inflate_work_permitted(&syncInflater_, work, result, error);
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

    if (cursor_ == buffer_.size()) {
        buffer_ = std::move(block.bytes);
        cursor_ = 0;
        currentBlockOffset_ = block.compressedOffset;
        return true;
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

bool BgzfRangeReader::read_line_view(const unsigned char** line,
                                     size_t* line_size,
                                     bool allow_clean_end,
                                     std::string* error) {
    if (line == nullptr || line_size == nullptr) {
        return set_error(error, "BGZF FASTQ line view destination is null");
    }
    *line = nullptr;
    *line_size = 0;
    while (true) {
        const size_t available = buffer_.size() - cursor_;
        const unsigned char *begin = available == 0
            ? nullptr : buffer_.data() + cursor_;
        const void *found = available == 0
            ? nullptr : std::memchr(begin, '\n', available);
        if (found != nullptr) {
            const unsigned char *newline =
                static_cast<const unsigned char *>(found);
            size_t size = static_cast<size_t>(newline - begin);
            cursor_ += size + 1;
            if (size != 0 && begin[size - 1] == '\r') {
                --size;
            }
            *line = begin;
            *line_size = size;
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
            const unsigned char *tail = buffer_.data() + cursor_;
            size_t size = buffer_.size() - cursor_;
            cursor_ = buffer_.size();
            if (size != 0 && tail[size - 1] == '\r') {
                --size;
            }
            *line = tail;
            *line_size = size;
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
    const unsigned char *line = nullptr;
    size_t line_size = 0;
    if (!read_line_view(&line, &line_size, true, error)) {
        return false;
    }
    if (line_size == 0 || line[0] != '@') {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_
                << " does not start with @ (block offset " << currentBlockOffset_ << ')';
        return set_error(error, message.str());
    }
    if (line_size - 1 > kBgzfFastqNameCapacity) {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_ << " has read name length "
                << (line_size - 1) << " exceeding fixed capacity "
                << kBgzfFastqNameCapacity;
        return set_error(error, message.str());
    }
    record->nameLength = static_cast<uint16_t>(line_size - 1);
    if (record->nameLength != 0) {
        std::memcpy(record->name, line + 1, record->nameLength);
    }

    if (!read_line_view(&line, &line_size, false, error)) {
        return false;
    }
    if (line_size > kBgzfFastqSequenceCapacity) {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_ << " has sequence length "
                << line_size << " exceeding fixed capacity "
                << kBgzfFastqSequenceCapacity;
        return set_error(error, message.str());
    }
    record->sequenceLength = static_cast<uint16_t>(line_size);
    if (line_size != 0) {
        std::memcpy(record->sequence, line, line_size);
    }

    if (!read_line_view(&line, &line_size, false, error)) {
        return false;
    }
    if (line_size == 0 || line[0] != '+') {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_
                << " has an invalid plus line (block offset " << currentBlockOffset_ << ')';
        return set_error(error, message.str());
    }

    if (!read_line_view(&line, &line_size, false, error)) {
        return false;
    }
    if (line_size > kBgzfFastqSequenceCapacity) {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_ << " has quality length "
                << line_size << " exceeding fixed capacity "
                << kBgzfFastqSequenceCapacity;
        return set_error(error, message.str());
    }
    record->qualityLength = static_cast<uint16_t>(line_size);
    if (storeQuality_ && line_size != 0) {
        std::memcpy(record->quality, line, line_size);
    }
    if (record->sequenceLength != record->qualityLength) {
        std::ostringstream message;
        message << "BGZF FASTQ record " << recordsRead_ << " has sequence length "
                << record->sequenceLength << " but quality length "
                << record->qualityLength << " (block offset "
                << currentBlockOffset_ << ')';
        return set_error(error, message.str());
    }
    record->ordinal = recordsRead_;
    return true;
}

bool BgzfRangeReader::next(BgzfFastqRecord* record, std::string* error) {
    if (error != nullptr && !error->empty()) {
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
