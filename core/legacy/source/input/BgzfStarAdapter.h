#ifndef CODE_input_BgzfStarAdapter
#define CODE_input_BgzfStarAdapter

#include "input/BgzfRangeReader.h"
#include "input/InputContract.h"

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>

namespace star {
namespace input {

struct BgzfStarAdapterOptions {
    uint32_t lane_index = 0;
    uint32_t mate0_reader_threads = 0;
    uint32_t mate1_reader_threads = 0;
    bool crc_check = true;
    bool store_mate0_quality = true;
    bool store_mate1_quality = true;
    // Require the read-name stem of record i in both mates to agree.
    bool validate_read_names = true;
    BgzfWorkPermitHooks inflate_permit_hooks;
};

struct BgzfStarRecord {
    BgzfFastqRecord mates[2];
    uint16_t read_name_length = 0;
    uint32_t lane_index = 0;
    uint64_t read_ordinal = 0;
    char read_filter = 'Y';
};

struct BgzfStarBatchLease {
    BgzfBatchLease mates[2];

    void clear() {
        mates[0].clear();
        mates[1].clear();
    }
};

// Pairs two ordered BGZF streams by source ordinal. Member inflation is
// parallel inside each stream; ordered parse/pair work is claimed in batches
// so multiple fused Flex consumers can safely share one lane adapter without
// a mutex handoff for every record.
class BgzfStarAdapter {
public:
    bool open(const std::string& mate0_path,
              const std::string& mate1_path,
              const BgzfStarAdapterOptions& options,
              std::string* error);

    InputStatus next_record(BgzfStarRecord* record, std::string* error);
    InputStatus next_records(BgzfStarRecord* records,
                             size_t capacity,
                             size_t* records_returned,
                             std::string* error,
                             BgzfStarBatchLease* lease = nullptr);

    uint64_t records_read() const;

private:
    InputStatus next_record_locked(BgzfStarRecord* record, std::string* error,
                                   BgzfStarBatchLease* lease);

    BgzfRangeReader readers_[2];
    BgzfStarAdapterOptions options_;
    mutable std::mutex nextMutex_;
    std::string readerError_;
    uint64_t recordsRead_ = 0;
    bool open_ = false;
    bool ended_ = false;
};

} // namespace input
} // namespace star

#endif
