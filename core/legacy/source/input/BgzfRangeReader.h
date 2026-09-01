#ifndef CODE_input_BgzfRangeReader
#define CODE_input_BgzfRangeReader

#include "input/BgzfIndex.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace star {
namespace input {

struct BgzfFastqRecord {
    std::string name;
    std::string sequence;
    std::string plus;
    std::string quality;
    uint64_t ordinal = 0;
};

class BgzfRangeReader {
public:
    BgzfRangeReader();
    ~BgzfRangeReader();

    BgzfRangeReader(const BgzfRangeReader&) = delete;
    BgzfRangeReader& operator=(const BgzfRangeReader&) = delete;

    bool open(const BgzfIndex* index,
              uint64_t first_record,
              uint64_t record_count,
              bool check_crc,
              std::string* error);

    // Returns false with an empty error at the end of the owned record range.
    bool next(BgzfFastqRecord* record, std::string* error);

    uint64_t first_record() const;
    uint64_t end_record() const;
    uint64_t records_read() const;

private:
    bool append_block(std::string* error);
    bool read_line(std::string* line, std::string* error);
    bool parse_record(BgzfFastqRecord* record, std::string* error);
    void close_input();

    const BgzfIndex* index_ = nullptr;
    int inputFd_ = -1;
    bool checkCrc_ = true;
    uint64_t firstRecord_ = 0;
    uint64_t endRecord_ = 0;
    uint64_t currentOrdinal_ = 0;
    uint64_t recordsRead_ = 0;
    size_t nextBlockIndex_ = 0;
    size_t currentBlockIndex_ = 0;
    std::vector<unsigned char> buffer_;
    size_t cursor_ = 0;
};

} // namespace input
} // namespace star

#endif
