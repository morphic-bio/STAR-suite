#include "input/BgzfStarAdapter.h"

#include <cstring>
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

// The range reader has already retained only the whitespace-delimited token.
// Strip a legacy mate suffix when diagnostic pairing validation is enabled.
size_t read_name_stem_length(const BgzfFastqRecord& raw) {
    size_t end = raw.nameLength;
    const char* name = raw.name_data();
    if (end >= 2 && name[end - 2] == '/' &&
        (name[end - 1] == '1' || name[end - 1] == '2' ||
         name[end - 1] == '3')) {
        end -= 2;
    }
    return end;
}

std::string printable_read_name(const BgzfFastqRecord& record) {
    return std::string(record.name_data(), record.nameLength);
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
                          options.mate0_reader_threads, options.crc_check, error,
                          &options.inflate_permit_hooks,
                          options.store_mate0_quality,
                          options.parse_illumina_filter
                              ? BgzfNameMode::TokenAndIlluminaFilter
                              : BgzfNameMode::Token)) {
        return false;
    }
    if (!readers_[1].open(mate1_path, 0, physical_end,
                          options.mate1_reader_threads, options.crc_check, error,
                          &options.inflate_permit_hooks,
                          options.store_mate1_quality,
                          options.validate_read_names
                              ? BgzfNameMode::Token : BgzfNameMode::Skip)) {
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
    return next_record_locked(record, error, nullptr);
}

InputStatus BgzfStarAdapter::next_records(BgzfStarRecord* records,
                                          size_t capacity,
                                          size_t* records_returned,
                                          std::string* error,
                                          BgzfStarBatchLease* lease) {
    std::lock_guard<std::mutex> lock(nextMutex_);
    if (error != nullptr) {
        error->clear();
    }
    if (records_returned == nullptr) {
        set_error(error, "BGZF STAR adapter batch count destination is null");
        return InputStatus::Error;
    }
    *records_returned = 0;
    if (lease != nullptr) {
        lease->clear();
    }
    if (records == nullptr || capacity == 0) {
        set_error(error, "BGZF STAR adapter batch destination is empty");
        return InputStatus::Error;
    }
    while (*records_returned < capacity) {
        const InputStatus status =
            next_record_locked(records + *records_returned, error, lease);
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
                                                std::string* error,
                                                BgzfStarBatchLease* lease) {
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

    BgzfBatchLease* mate0_lease = lease == nullptr ? nullptr : &lease->mates[0];
    BgzfBatchLease* mate1_lease = lease == nullptr ? nullptr : &lease->mates[1];
    const bool have_mate0 = readers_[0].next(
        &record->mates[0], &readerError_, mate0_lease);
    if (!have_mate0 && !readerError_.empty()) {
        set_error(error, "BGZF mate 0: " + readerError_);
        return InputStatus::Error;
    }
    const bool have_mate1 = readers_[1].next(
        &record->mates[1], &readerError_, mate1_lease);
    if (!have_mate1 && !readerError_.empty()) {
        set_error(error, "BGZF mate 1: " + readerError_);
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
    if (options_.validate_read_names) {
        const size_t stem0 = read_name_stem_length(record->mates[0]);
        const size_t stem1 = read_name_stem_length(record->mates[1]);
        if (stem0 != stem1 ||
            std::memcmp(record->mates[0].name_data(),
                        record->mates[1].name_data(), stem0) != 0) {
            std::ostringstream message;
            message << "BGZF mate read-name mismatch at record " << recordsRead_
                    << ": mate 0 '" << printable_read_name(record->mates[0])
                    << "' does not pair with mate 1 '"
                    << printable_read_name(record->mates[1])
                    << "'";
            set_error(error, message.str());
            return InputStatus::Error;
        }
    }

    record->read_name_length = record->mates[0].nameLength;
    record->lane_index = options_.lane_index;
    record->read_ordinal = recordsRead_;
    record->read_filter = options_.parse_illumina_filter
        ? record->mates[0].readFilter : 'Y';
    ++recordsRead_;
    return InputStatus::Record;
}

uint64_t BgzfStarAdapter::records_read() const {
    std::lock_guard<std::mutex> lock(nextMutex_);
    return recordsRead_;
}

} // namespace input
} // namespace star
