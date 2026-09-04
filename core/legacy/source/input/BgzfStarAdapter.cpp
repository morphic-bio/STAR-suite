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

// The read-name stem is the name up to the first whitespace, minus a legacy
// "/1" or "/2" mate suffix. Illumina writes the same stem into both mates.
struct NameParts {
    const char* data;
    uint16_t stemLength;
    uint16_t readNameLength;
    uint16_t extraOffset;
    uint16_t extraLength;
};

size_t read_name_token_length(const BgzfFastqRecord& raw) {
    size_t end = raw.nameLength;
    const void* space = std::memchr(raw.name, ' ', end);
    if (space != nullptr) {
        end = static_cast<const char*>(space) - raw.name;
    }
    const void* tab = std::memchr(raw.name, '\t', end);
    if (tab != nullptr) {
        end = static_cast<const char*>(tab) - raw.name;
    }
    return end;
}

size_t read_name_stem_length(const BgzfFastqRecord& raw,
                             size_t token_length) {
    if (token_length >= 2 && raw.name[token_length - 2] == '/' &&
        (raw.name[token_length - 1] == '1' ||
         raw.name[token_length - 1] == '2' ||
         raw.name[token_length - 1] == '3')) {
        return token_length - 2;
    }
    return token_length;
}

NameParts read_name_parts(const BgzfFastqRecord& raw) {
    const size_t read_name_length = read_name_token_length(raw);
    const size_t stem_length = read_name_stem_length(raw, read_name_length);
    size_t extra_begin = read_name_length;
    while (extra_begin < raw.nameLength &&
           (raw.name[extra_begin] == ' ' || raw.name[extra_begin] == '\t')) {
        ++extra_begin;
    }
    return NameParts{raw.name,
                     static_cast<uint16_t>(stem_length),
                     static_cast<uint16_t>(read_name_length),
                     static_cast<uint16_t>(extra_begin),
                     static_cast<uint16_t>(raw.nameLength - extra_begin)};
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
                          &options.inflate_permit_hooks,
                          options.store_mate0_quality)) {
        return false;
    }
    if (!readers_[1].open(mate1_path, 0, physical_end,
                          options.mate1_reader_threads, options.crc_check, error,
                          &options.inflate_permit_hooks,
                          options.store_mate1_quality)) {
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
    const NameParts name0 = read_name_parts(record->mates[0]);
    if (options_.validate_read_names) {
        const size_t name1_token_length =
            read_name_token_length(record->mates[1]);
        const size_t name1_stem_length =
            read_name_stem_length(record->mates[1], name1_token_length);
        if (name0.stemLength != name1_stem_length ||
            std::memcmp(name0.data, record->mates[1].name,
                        name0.stemLength) != 0) {
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

    record->read_name_length = name0.readNameLength;
    record->read_name_extra_offset = name0.extraOffset;
    record->read_name_extra_length = name0.extraLength;
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
