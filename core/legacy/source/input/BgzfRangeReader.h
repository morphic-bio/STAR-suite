#ifndef CODE_input_BgzfRangeReader
#define CODE_input_BgzfRangeReader

#include "input/BgzfBlockReader.h"

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>
#include <thread>
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

// Inflate workers claim one member at a time by reading BC/BSIZE at the
// current compressed offset. Results are assembled in member order before
// FASTQ parsing, so no block or record index is needed.
class BgzfRangeReader {
public:
    BgzfRangeReader();
    ~BgzfRangeReader();

    BgzfRangeReader(const BgzfRangeReader&) = delete;
    BgzfRangeReader& operator=(const BgzfRangeReader&) = delete;

    // range_end == UINT64_MAX selects physical EOF. range_start must be a BGZF
    // member boundary; STAR currently opens complete files as [0, EOF).
    bool open(const std::string& path,
              uint64_t range_start,
              uint64_t range_end,
              uint32_t worker_threads,
              bool check_crc,
              std::string* error);

    // Returns false with an empty error at clean stream end.
    bool next(BgzfFastqRecord* record, std::string* error);

    uint64_t records_read() const;
    uint64_t range_start() const;
    uint64_t range_end() const;

private:
    struct InflatedBlock {
        uint64_t compressedOffset = 0;
        std::vector<unsigned char> bytes;
    };

    void worker_loop();
    void close_input();
    void fail_locked(const std::string& message);
    bool claim_work(std::vector<BgzfBlock>* blocks,
                    uint64_t* sequence,
                    bool* at_end,
                    std::string* error);
    bool claim_and_inflate_sync(InflatedBlock* result, bool* at_end,
                                std::string* error);
    bool append_next_block(std::string* error);
    bool read_line(std::string* line, bool allow_clean_end, std::string* error);
    bool parse_record(BgzfFastqRecord* record, std::string* error);

    std::string path_;
    int inputFd_ = -1;
    bool checkCrc_ = true;
    uint64_t rangeStart_ = 0;
    uint64_t rangeEnd_ = 0;
    uint64_t claimedOffset_ = 0;
    uint64_t nextClaimSequence_ = 0;
    uint64_t nextConsumeSequence_ = 0;
    uint64_t currentBlockOffset_ = 0;
    uint64_t recordsRead_ = 0;
    uint64_t targetCompressedBytes_ = 64 * 1024;
    size_t maxOutstandingWork_ = 4;
    uint32_t workerCount_ = 0;

    std::vector<unsigned char> buffer_;
    size_t cursor_ = 0;
    std::vector<std::thread> workers_;
    std::map<uint64_t, InflatedBlock> completed_;
    std::mutex mutex_;
    std::condition_variable readyCv_;
    std::condition_variable spaceCv_;
    bool claimsFinished_ = false;
    bool failed_ = false;
    bool stopping_ = false;
    std::string workerError_;
};

} // namespace input
} // namespace star

#endif
