#ifndef CODE_input_BgzfRangeReader
#define CODE_input_BgzfRangeReader

#include "input/BgzfBlockReader.h"

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <memory>
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
    TokenAndIlluminaFilter,
    Skip
};

// Keeps decompressed work buffers alive while a caller consumes one batch of
// zero-copy record views. Clear or destroy the lease only after every record
// returned with it has been consumed.
class BgzfBatchLease {
public:
    void clear();

private:
    friend class BgzfRangeReader;
    void retain(const std::shared_ptr<std::vector<unsigned char>>& buffer);

    std::vector<std::shared_ptr<std::vector<unsigned char>>> buffers_;
};

struct BgzfFastqOwnedFields {
    char name[kBgzfFastqNameCapacity];
    char sequence[kBgzfFastqSequenceCapacity];
    char quality[kBgzfFastqSequenceCapacity];
};

struct BgzfFastqRecord {
    uint16_t nameLength = 0;
    uint16_t sequenceLength = 0;
    uint16_t qualityLength = 0;
    uint64_t ordinal = 0;
    char readFilter = 'N';
    const char* nameView = nullptr;
    const char* sequenceView = nullptr;
    const char* qualityView = nullptr;
    std::unique_ptr<BgzfFastqOwnedFields> ownedFields;

    const char* name_data() const {
        return nameView != nullptr ? nameView
                                   : ownedFields == nullptr ? nullptr
                                                            : ownedFields->name;
    }
    const char* sequence_data() const {
        return sequenceView != nullptr ? sequenceView
                                       : ownedFields == nullptr ? nullptr
                                                                : ownedFields->sequence;
    }
    const char* quality_data() const {
        return qualityView != nullptr ? qualityView
                                      : ownedFields == nullptr ? nullptr
                                                               : ownedFields->quality;
    }
    bool name_is_view() const {
        return nameView != nullptr;
    }
    bool sequence_is_view() const {
        return sequenceView != nullptr;
    }
    bool quality_is_view() const {
        return qualityView != nullptr;
    }

    char* name_storage();
    char* sequence_storage();
    char* quality_storage();

private:
    void ensure_owned_fields();
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
    bool next(BgzfFastqRecord* record, std::string* error,
              BgzfBatchLease* lease = nullptr);

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
        std::shared_ptr<std::vector<unsigned char>> bytes;
    };

    struct CompletedSlot {
        uint64_t sequence = 0;
        InflatedBlock block;
        bool ready = false;
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
                        std::string* error,
                        bool* line_is_buffer_view = nullptr);
    bool read_name_token(BgzfFastqRecord* record,
                         bool allow_clean_end,
                         std::string* error,
                         BgzfBatchLease* lease);
    bool parse_record(BgzfFastqRecord* record, std::string* error,
                      BgzfBatchLease* lease);

    std::string path_;
    int inputFd_ = -1;
    bool checkCrc_ = true;
    bool storeQuality_ = true;
    BgzfNameMode nameMode_ = BgzfNameMode::Full;
    uint64_t rangeStart_ = 0;
    uint64_t rangeEnd_ = 0;
    uint64_t claimedOffset_ = 0;
    uint64_t nextClaimSequence_ = 0;
    uint64_t claimedWorkCount_ = 0;
    uint64_t nextConsumeSequence_ = 0;
    uint64_t currentBlockOffset_ = 0;
    uint64_t recordsRead_ = 0;
    uint64_t targetCompressedBytes_ = 64 * 1024;
    size_t maxOutstandingWork_ = 4;
    size_t outstandingWork_ = 0;
    uint32_t workerCount_ = 0;
    BgzfWorkPermitHooks permitHooks_;
    BgzfInflater syncInflater_;

    std::shared_ptr<std::vector<unsigned char>> buffer_;
    // One extra byte permits CRLF at the logical FASTQ line capacity.
    unsigned char lineScratch_[kBgzfFastqSequenceCapacity + 1];
    size_t cursor_ = 0;
    std::vector<std::thread> workers_;
    std::vector<CompletedSlot> completed_;
    std::mutex claimMutex_;
    std::mutex mutex_;
    std::condition_variable readyCv_;
    std::condition_variable spaceCv_;
    bool claimsFinished_ = false;
    bool claimExhausted_ = false;
    bool failed_ = false;
    bool stopping_ = false;
    std::string workerError_;
};

} // namespace input
} // namespace star

#endif
