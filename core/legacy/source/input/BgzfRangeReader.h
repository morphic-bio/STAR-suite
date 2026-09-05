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

static constexpr size_t kBgzfFastqNameCapacity = 512;
static constexpr size_t kBgzfFastqSequenceCapacity = 650;

enum class BgzfNameMode {
    Full,
    Token,
    Skip
};

struct BgzfFastqRecord {
    char name[kBgzfFastqNameCapacity];
    char sequence[kBgzfFastqSequenceCapacity];
    char quality[kBgzfFastqSequenceCapacity];
    uint16_t nameLength = 0;
    uint16_t sequenceLength = 0;
    uint16_t qualityLength = 0;
    uint64_t ordinal = 0;
};

// Optional scheduler hooks around one bounded inflate work item. The BGZF
// reader deliberately knows nothing about STAR's thread controller; callers
// can supply any acquire/release pair with the same telemetry contract.
struct BgzfWorkPermitHooks {
    using Acquire = uint64_t (*)(void* context);
    using Release = void (*)(void* context,
                             uint64_t wait_ns,
                             uint64_t work_units,
                             uint64_t work_bytes,
                             uint64_t work_ns);

    void* context = nullptr;
    Acquire acquire = nullptr;
    Release release = nullptr;

    bool enabled() const {
        return acquire != nullptr && release != nullptr;
    }
};

// Inflate workers claim bounded contiguous work by reading every member's
// BC/BSIZE at the current compressed frontier. Results are assembled in claim
// order before FASTQ parsing, so no block or record index is needed.
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
              std::string* error,
              const BgzfWorkPermitHooks* permit_hooks = nullptr,
              bool store_quality = true,
              BgzfNameMode name_mode = BgzfNameMode::Full);

    // Returns false with an empty error at clean stream end.
    bool next(BgzfFastqRecord* record, std::string* error);

    uint64_t records_read() const;
    uint64_t range_start() const;
    uint64_t range_end() const;

private:
    struct CompressedWork {
        uint64_t compressedOffset = 0;
        std::vector<BgzfBlock> blocks;
        std::vector<unsigned char> bytes;
    };

    struct InflatedBlock {
        uint64_t compressedOffset = 0;
        std::vector<unsigned char> bytes;
    };

    void worker_loop();
    void close_input();
    void fail_locked(const std::string& message);
    bool claim_work(CompressedWork* work,
                    uint64_t* sequence,
                    bool* at_end,
                    std::string* error);
    bool inflate_work(BgzfInflater* inflater,
                      const CompressedWork& work,
                      InflatedBlock* result,
                      std::string* error);
    bool inflate_work_permitted(BgzfInflater* inflater,
                                const CompressedWork& work,
                                InflatedBlock* result,
                                std::string* error);
    bool claim_and_inflate_sync(InflatedBlock* result, bool* at_end,
                                std::string* error);
    bool append_next_block(std::string* error);
    bool read_line_view(const unsigned char** line,
                        size_t* line_size,
                        bool allow_clean_end,
                        std::string* error);
    bool read_name_token(BgzfFastqRecord* record,
                         bool allow_clean_end,
                         std::string* error);
    bool parse_record(BgzfFastqRecord* record, std::string* error);

    std::string path_;
    int inputFd_ = -1;
    bool checkCrc_ = true;
    bool storeQuality_ = true;
    BgzfNameMode nameMode_ = BgzfNameMode::Full;
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
    BgzfWorkPermitHooks permitHooks_;
    BgzfInflater syncInflater_;

    std::vector<unsigned char> buffer_;
    // One extra byte permits CRLF at the logical FASTQ line capacity.
    unsigned char lineScratch_[kBgzfFastqSequenceCapacity + 1];
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
