#include "input/BgzfStarAdapter.h"

#include <limits>
#include <sstream>

namespace star {
namespace input {
namespace {

bool set_error(std::string* error, const std::string& message) {
    if (error != nullptr) {
        *error = message;
    }
    return false;
}

void split_read_name(const std::string& raw,
                     std::string* read_name,
                     std::string* read_name_extra) {
    const size_t separator = raw.find_first_of(" \t");
    if (separator == std::string::npos) {
        *read_name = raw;
        read_name_extra->clear();
        return;
    }
    read_name->assign(raw.data(), separator);
    size_t extra_begin = separator;
    while (extra_begin < raw.size() &&
           (raw[extra_begin] == ' ' || raw[extra_begin] == '\t')) {
        ++extra_begin;
    }
    read_name_extra->assign(raw.data() + extra_begin, raw.size() - extra_begin);
}

} // namespace

bool BgzfStarAdapter::open(const std::string& mate0_path,
                           const std::string& mate1_path,
                           const BgzfStarAdapterOptions& options,
                           std::string* error) {
    std::lock_guard<std::mutex> lock(nextMutex_);
    open_ = false;
    ended_ = false;
    recordsRead_ = 0;
    if (error != nullptr) {
        error->clear();
    }
    const uint64_t physical_end = std::numeric_limits<uint64_t>::max();
    if (!readers_[0].open(mate0_path, 0, physical_end,
                          options.mate0_reader_threads, options.crc_check, error)) {
        return false;
    }
    if (!readers_[1].open(mate1_path, 0, physical_end,
                          options.mate1_reader_threads, options.crc_check, error)) {
        return false;
    }
    options_ = options;
    open_ = true;
    return true;
}

InputStatus BgzfStarAdapter::next_record(BgzfStarRecord* record,
                                         std::string* error) {
    std::lock_guard<std::mutex> lock(nextMutex_);
    if (error != nullptr) {
        error->clear();
    }
    return next_record_locked(record, error);
}

InputStatus BgzfStarAdapter::next_records(BgzfStarRecord* records,
                                          size_t capacity,
                                          size_t* records_returned,
                                          std::string* error) {
    std::lock_guard<std::mutex> lock(nextMutex_);
    if (error != nullptr) {
        error->clear();
    }
    if (records_returned == nullptr) {
        set_error(error, "BGZF STAR adapter batch count destination is null");
        return InputStatus::Error;
    }
    *records_returned = 0;
    if (records == nullptr || capacity == 0) {
        set_error(error, "BGZF STAR adapter batch destination is empty");
        return InputStatus::Error;
    }
    while (*records_returned < capacity) {
        const InputStatus status =
            next_record_locked(records + *records_returned, error);
        if (status == InputStatus::Record) {
            ++*records_returned;
            continue;
        }
        if (status == InputStatus::End && *records_returned != 0) {
            return InputStatus::Record;
        }
        if (status == InputStatus::Error) {
            *records_returned = 0;
        }
        return status;
    }
    return InputStatus::Record;
}

InputStatus BgzfStarAdapter::next_record_locked(BgzfStarRecord* record,
                                                std::string* error) {
    if (!open_) {
        set_error(error, "BGZF STAR adapter is not open");
        return InputStatus::Error;
    }
    if (ended_) {
        return InputStatus::End;
    }
    if (record == nullptr) {
        set_error(error, "BGZF STAR adapter record destination is null");
        return InputStatus::Error;
    }

    std::string mate_error[2];
    const bool have_mate0 = readers_[0].next(&record->mates[0], &mate_error[0]);
    if (!have_mate0 && !mate_error[0].empty()) {
        set_error(error, "BGZF mate 0: " + mate_error[0]);
        return InputStatus::Error;
    }
    const bool have_mate1 = readers_[1].next(&record->mates[1], &mate_error[1]);
    if (!have_mate1 && !mate_error[1].empty()) {
        set_error(error, "BGZF mate 1: " + mate_error[1]);
        return InputStatus::Error;
    }
    if (have_mate0 != have_mate1) {
        std::ostringstream message;
        message << "BGZF mate record-count mismatch at record " << recordsRead_
                << ": mate 0 " << (have_mate0 ? "has a record" : "ended")
                << " while mate 1 " << (have_mate1 ? "has a record" : "ended");
        set_error(error, message.str());
        return InputStatus::Error;
    }
    if (!have_mate0) {
        ended_ = true;
        return InputStatus::End;
    }
    if (record->mates[0].ordinal != recordsRead_ ||
        record->mates[1].ordinal != recordsRead_) {
        std::ostringstream message;
        message << "BGZF mate ordinal mismatch at expected record " << recordsRead_
                << ": mate ordinals are " << record->mates[0].ordinal
                << " and " << record->mates[1].ordinal;
        set_error(error, message.str());
        return InputStatus::Error;
    }

    split_read_name(record->mates[0].name,
                    &record->read_name, &record->read_name_extra);
    record->lane_index = options_.lane_index;
    record->read_ordinal = recordsRead_;
    record->read_filter = 'Y';
    ++recordsRead_;
    return InputStatus::Record;
}

uint64_t BgzfStarAdapter::records_read() const {
    std::lock_guard<std::mutex> lock(nextMutex_);
    return recordsRead_;
}

} // namespace input
} // namespace star
