#include "input/BgzfStarAdapter.h"

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

bool BgzfStarAdapter::open(const BgzfIndex* mate0_index,
                           const BgzfIndex* mate1_index,
                           const BgzfStarAdapterOptions& options,
                           std::string* error) {
    open_ = false;
    recordsRead_ = 0;
    if (error != nullptr) {
        error->clear();
    }
    if (mate0_index == nullptr || mate1_index == nullptr) {
        return set_error(error, "BGZF STAR adapter requires two mate indexes");
    }
    if (mate0_index->record_count() != mate1_index->record_count()) {
        std::ostringstream message;
        message << "BGZF mate record-count mismatch: mate 0 has "
                << mate0_index->record_count() << " records but mate 1 has "
                << mate1_index->record_count();
        return set_error(error, message.str());
    }
    if (options.first_record > mate0_index->record_count() ||
        options.record_count > mate0_index->record_count() - options.first_record) {
        return set_error(error, "BGZF STAR adapter range exceeds the paired record count");
    }
    if (!readers_[0].open(mate0_index, options.first_record, options.record_count,
                          options.crc_check, error)) {
        return false;
    }
    if (!readers_[1].open(mate1_index, options.first_record, options.record_count,
                          options.crc_check, error)) {
        return false;
    }
    options_ = options;
    open_ = true;
    return true;
}

InputStatus BgzfStarAdapter::next_record(BgzfStarRecord* record,
                                         std::string* error) {
    if (error != nullptr) {
        error->clear();
    }
    if (!open_) {
        set_error(error, "BGZF STAR adapter is not open");
        return InputStatus::Error;
    }
    if (recordsRead_ >= options_.record_count) {
        return InputStatus::End;
    }
    if (record == nullptr) {
        set_error(error, "BGZF STAR adapter record destination is null");
        return InputStatus::Error;
    }
    std::string mate_error;
    if (!readers_[0].next(&record->mates[0], &mate_error)) {
        set_error(error, mate_error.empty()
            ? "BGZF mate 0 ended before its owned record range"
            : "BGZF mate 0: " + mate_error);
        return InputStatus::Error;
    }
    if (!readers_[1].next(&record->mates[1], &mate_error)) {
        set_error(error, mate_error.empty()
            ? "BGZF mate 1 ended before its owned record range"
            : "BGZF mate 1: " + mate_error);
        return InputStatus::Error;
    }
    const uint64_t expected = options_.first_record + recordsRead_;
    if (record->mates[0].ordinal != expected ||
        record->mates[1].ordinal != expected) {
        std::ostringstream message;
        message << "BGZF mate ordinal mismatch at expected record " << expected
                << ": mate ordinals are " << record->mates[0].ordinal
                << " and " << record->mates[1].ordinal;
        set_error(error, message.str());
        return InputStatus::Error;
    }
    split_read_name(record->mates[0].name,
                    &record->read_name, &record->read_name_extra);
    record->lane_index = options_.lane_index;
    record->read_ordinal = expected;
    record->read_filter = 'Y';
    ++recordsRead_;
    return InputStatus::Record;
}

uint64_t BgzfStarAdapter::records_read() const {
    return recordsRead_;
}

uint64_t BgzfStarAdapter::record_count() const {
    return options_.record_count;
}

} // namespace input
} // namespace star
