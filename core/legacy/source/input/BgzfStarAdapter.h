#ifndef CODE_input_BgzfStarAdapter
#define CODE_input_BgzfStarAdapter

#include "input/BgzfRangeReader.h"
#include "input/InputContract.h"

#include <cstdint>
#include <string>

namespace star {
namespace input {

struct BgzfStarAdapterOptions {
    uint32_t lane_index = 0;
    uint64_t first_record = 0;
    uint64_t record_count = 0;
    bool crc_check = true;
};

struct BgzfStarRecord {
    BgzfFastqRecord mates[2];
    std::string read_name;
    std::string read_name_extra;
    uint32_t lane_index = 0;
    uint64_t read_ordinal = 0;
    char read_filter = 'Y';
};

// Pairs two independently indexed BGZF mates by source record ordinal and
// exposes one record-at-a-time producer semantics like the other STAR input
// adapters. The indexes remain owned by the caller and must outlive the adapter.
class BgzfStarAdapter {
public:
    bool open(const BgzfIndex* mate0_index,
              const BgzfIndex* mate1_index,
              const BgzfStarAdapterOptions& options,
              std::string* error);

    InputStatus next_record(BgzfStarRecord* record, std::string* error);

    uint64_t records_read() const;
    uint64_t record_count() const;

private:
    BgzfRangeReader readers_[2];
    BgzfStarAdapterOptions options_;
    uint64_t recordsRead_ = 0;
    bool open_ = false;
};

} // namespace input
} // namespace star

#endif
