#ifndef CODE_input_BgzfStarAdapter
#define CODE_input_BgzfStarAdapter

#include "input/BgzfRangeReader.h"
#include "input/InputContract.h"

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
};

struct BgzfStarRecord {
    BgzfFastqRecord mates[2];
    std::string read_name;
    std::string read_name_extra;
    uint32_t lane_index = 0;
    uint64_t read_ordinal = 0;
    char read_filter = 'Y';
};

// Pairs two ordered BGZF streams by source ordinal. Member inflation is
// parallel inside each stream; the short parse/pair section is serialized so
// multiple fused Flex consumers can safely share one lane adapter.
class BgzfStarAdapter {
public:
    bool open(const std::string& mate0_path,
              const std::string& mate1_path,
              const BgzfStarAdapterOptions& options,
              std::string* error);

    InputStatus next_record(BgzfStarRecord* record, std::string* error);

    uint64_t records_read() const;

private:
    BgzfRangeReader readers_[2];
    BgzfStarAdapterOptions options_;
    mutable std::mutex nextMutex_;
    uint64_t recordsRead_ = 0;
    bool open_ = false;
    bool ended_ = false;
};

} // namespace input
} // namespace star

#endif
