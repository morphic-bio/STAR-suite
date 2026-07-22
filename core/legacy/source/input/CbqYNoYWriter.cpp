#include "input/CbqYNoYWriter.h"

#include "GlobalVariables.h"
#include "Parameters.h"

#include <sstream>

namespace {

void set_failure(bool* failed, std::string* failure, const std::string& message) {
    if (failed != nullptr) {
        *failed = true;
    }
    if (failure != nullptr) {
        *failure = message;
    }
}

} // namespace

namespace star {
namespace input {

CbqYNoYOrderedWriter::CbqYNoYOrderedWriter()
    : next_chunk_id_(0),
      failed_(false),
      opened_(false) {}

bool CbqYNoYOrderedWriter::open(const std::string& y_path,
                                const std::string& noy_path,
                                uint32_t mate_count,
                                int compression_level,
                                uint64_t block_size,
                                std::string* error) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (opened_) {
        if (error != nullptr) {
            *error = "CBQ Y/noY writer is already open";
        }
        return false;
    }
    next_chunk_id_ = 0;
    failed_ = false;
    failure_.clear();

    if (!y_writer_.open(y_path, mate_count, compression_level, block_size, error)) {
        return false;
    }
    if (!noy_writer_.open(noy_path, mate_count, compression_level, block_size, error)) {
        return false;
    }
    opened_ = true;
    return true;
}

bool CbqYNoYOrderedWriter::submit_chunk(uint64_t chunk_id,
                                        const std::vector<CbqYNoYRecord>& records,
                                        std::string* error) {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [&]() { return failed_ || chunk_id == next_chunk_id_; });
    if (failed_) {
        if (error != nullptr) {
            *error = failure_;
        }
        return false;
    }
    if (!opened_) {
        set_failure(&failed_, &failure_, "CBQ Y/noY writer is not open");
        cv_.notify_all();
        if (error != nullptr) {
            *error = failure_;
        }
        return false;
    }
    for (size_t iread = 0; iread < records.size(); ++iread) {
        CbqWriter& writer = records[iread].is_y ? y_writer_ : noy_writer_;
        std::string local_error;
        if (!writer.add_record(records[iread].segments, &local_error)) {
            set_failure(&failed_, &failure_, local_error);
            cv_.notify_all();
            if (error != nullptr) {
                *error = failure_;
            }
            return false;
        }
    }

    ++next_chunk_id_;
    cv_.notify_all();
    return true;
}

bool CbqYNoYOrderedWriter::finish(std::string* error) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (failed_) {
        if (error != nullptr) {
            *error = failure_;
        }
        return false;
    }
    if (!opened_) {
        return true;
    }
    if (!y_writer_.finish(error)) {
        return false;
    }
    if (!noy_writer_.finish(error)) {
        return false;
    }
    opened_ = false;
    return true;
}

bool openGlobalCbqYNoYWriter(const Parameters& P, std::string* error) {
    if (g_cbqYNoYWriter != nullptr) {
        delete g_cbqYNoYWriter;
        g_cbqYNoYWriter = nullptr;
    }
    g_cbqYNoYWriter = new CbqYNoYOrderedWriter();
    if (!g_cbqYNoYWriter->open(P.outYCbqFile,
                               P.outNoYCbqFile,
                               P.yFastqEmitReadCount(),
                               P.emitYNoYCbqCompressionLevel,
                               P.emitYNoYCbqBlockSize,
                               error)) {
        delete g_cbqYNoYWriter;
        g_cbqYNoYWriter = nullptr;
        return false;
    }
    return true;
}

bool finishGlobalCbqYNoYWriter(std::string* error) {
    if (g_cbqYNoYWriter == nullptr) {
        return true;
    }
    const bool ok = g_cbqYNoYWriter->finish(error);
    delete g_cbqYNoYWriter;
    g_cbqYNoYWriter = nullptr;
    return ok;
}

} // namespace input
} // namespace star
