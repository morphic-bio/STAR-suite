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

void split_read_name(const BgzfFastqRecord& raw,
                     uint16_t* read_name_length,
                     uint16_t* read_name_extra_offset,
                     uint16_t* read_name_extra_length) {
    size_t separator = 0;
    while (separator < raw.nameLength && raw.name[separator] != ' ' &&
           raw.name[separator] != '\t') {
        ++separator;
    }
    *read_name_length = static_cast<uint16_t>(separator);
    size_t extra_begin = separator;
    while (extra_begin < raw.nameLength &&
           (raw.name[extra_begin] == ' ' || raw.name[extra_begin] == '\t')) {
        ++extra_begin;
    }
    *read_name_extra_offset = static_cast<uint16_t>(extra_begin);
    *read_name_extra_length = static_cast<uint16_t>(raw.nameLength - extra_begin);
}

// The read-name stem is the name up to the first whitespace, minus a legacy
// "/1" or "/2" mate suffix. Illumina writes the same stem into both mates.
struct NameStem {
    const char* data;
    size_t size;
};

NameStem read_name_stem(const BgzfFastqRecord& raw) {
    size_t end = 0;
    while (end < raw.nameLength && raw.name[end] != ' ' && raw.name[end] != '\t') {
        ++end;
    }
    if (end >= 2 && raw.name[end - 2] == '/' &&
        (raw.name[end - 1] == '1' || raw.name[end - 1] == '2' ||
         raw.name[end - 1] == '3')) {
        end -= 2;
    }
    return NameStem{raw.name, end};
}

std::string printable_read_name(const BgzfFastqRecord& record) {
    return std::string(record.name, record.nameLength);
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
                          &options.inflate_permit_hooks)) {
        return false;
    }
    if (!readers_[1].open(mate1_path, 0, physical_end,
                          options.mate1_reader_threads, options.crc_check, error,
                          &options.inflate_permit_hooks)) {
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

    const bool have_mate0 = readers_[0].next(&record->mates[0], &readerError_);
    if (!have_mate0 && !readerError_.empty()) {
        set_error(error, "BGZF mate 0: " + readerError_);
        return InputStatus::Error;
    }
    const bool have_mate1 = readers_[1].next(&record->mates[1], &readerError_);
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
    // Pairing is by ordinal; the read names are the only evidence inside the
    // files that record i of R1 belongs with record i of R2. Without this
    // check a mate file that was truncated, filtered, or re-sorted on its
    // own would pair silently, and the end-of-stream count check catches it
    // only if the counts happen to differ.
    if (options_.validate_read_names) {
        const NameStem stem0 = read_name_stem(record->mates[0]);
        const NameStem stem1 = read_name_stem(record->mates[1]);
        if (stem0.size != stem1.size ||
            std::memcmp(stem0.data, stem1.data, stem0.size) != 0) {
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

    split_read_name(record->mates[0], &record->read_name_length,
                    &record->read_name_extra_offset,
                    &record->read_name_extra_length);
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
