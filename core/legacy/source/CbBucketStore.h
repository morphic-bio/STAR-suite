#ifndef H_CbBucketStore
#define H_CbBucketStore

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace star {
namespace solo {

// Defined on-disk record layout (12 bytes, little-endian):
//
//   bytes  0..7  key   [CB20][UMI24][GENE15][TAG5], most to least significant
//   bytes 8..11  value [FLAGS2][COUNT30]
//
// FLAGS2 currently carries the Flex probe-region status. The count saturates at
// 2^30-1, matching the existing fused Flex hash value. Reserved flag values and
// future spill schemas must use a new schema number rather than reinterpret v1.
struct PackedCbRecord {
    std::uint64_t key = 0;
    std::uint32_t value = 0;

    static const std::size_t kSerializedBytes = 12;

    static PackedCbRecord make(std::uint32_t cbIndex, std::uint32_t umi24,
                               std::uint16_t gene15, std::uint8_t tag5,
                               std::uint32_t count30, std::uint8_t flags2);
    // Defined here rather than in the .cpp so they inline. The bucket collapse
    // sorts 1.6 billion records by a key built from four of these, so a
    // comparison sort calls them tens of billions of times; out of line, in
    // another translation unit and with no LTO, each shift-and-mask became a
    // real function call and the four together accounted for ~15% of the Flex
    // tail profile.
    std::uint32_t cb_index() const { return (key >> 44) & 0xFFFFFu; }
    std::uint32_t umi24() const { return (key >> 20) & 0xFFFFFFu; }
    std::uint16_t gene15() const { return static_cast<std::uint16_t>((key >> 5) & 0x7FFFu); }
    std::uint8_t tag5() const { return static_cast<std::uint8_t>(key & 0x1Fu); }
    // Reorders the packed fields from [CB][UMI][GENE][TAG] to the collapse
    // order [CB][TAG][GENE][UMI]. This is a permutation of every key bit, so
    // equality of sort keys is equality of packed keys.
    std::uint64_t group_sort_key() const {
        return (static_cast<std::uint64_t>(cb_index()) << 44)
             | (static_cast<std::uint64_t>(tag5()) << 39)
             | (static_cast<std::uint64_t>(gene15()) << 24)
             | umi24();
    }
    std::uint32_t count30() const;
    std::uint8_t flags2() const;

    void encode(std::uint8_t output[kSerializedBytes]) const;
    static bool decode(const std::uint8_t input[kSerializedBytes],
                       PackedCbRecord *record);
};

class CbBucketStore {
  public:
    enum class Mode { Ram, Spill, Auto };

    struct Config {
        Mode mode = Mode::Ram;
        std::uint32_t bucketCount = 256;
        std::uint32_t whitelistSize = 0;
        std::uint64_t memoryBudgetBytes = 0;
        std::string scratchDirectory;
        std::string filePrefix = "cb_bucket";
        // RAM/auto modes may consolidate sorted producer runs in the
        // background. Zero leaves every run for the final k-way merge.
        std::uint32_t mergeWorkerCount = 0;
        std::uint32_t mergeFanIn = 64;
    };

    explicit CbBucketStore(const Config &config);
    ~CbBucketStore();
    CbBucketStore(const CbBucketStore &) = delete;
    CbBucketStore &operator=(const CbBucketStore &) = delete;

    std::uint32_t bucket_for_cb(std::uint32_t cbIndex) const;
    std::uint32_t bucket_count() const { return config_.bucketCount; }

    // The caller seals and gives up ownership of a homogeneous bucket segment.
    // Offset/sequence ownership is claimed atomically; serialization and pwrite
    // happen outside the claim lock.
    bool append_segment(std::uint32_t workerIndex, std::uint32_t bucketIndex,
                        std::vector<PackedCbRecord> records,
                        std::string *error);
    bool finalize(std::string *error);

    bool load_bucket(std::uint32_t bucketIndex,
                     std::vector<PackedCbRecord> *records,
                     std::string *error) const;
    // Returns the producer-local runs separately. append_segment orders each
    // run before publishing it, allowing consumers to k-way merge instead of
    // sorting a whole bucket again. Runs are returned in publication order.
    bool load_sorted_segments(
        std::uint32_t bucketIndex,
        std::vector<std::vector<PackedCbRecord> > *segments,
        std::string *error) const;
    bool load_bucket_bytes(std::uint32_t bucketIndex,
                           std::vector<std::uint8_t> *bytes,
                           std::string *error) const;

    void reset_bucket_claims();
    bool claim_bucket(std::uint32_t *bucketIndex);

    bool using_spill() const;
    bool transitioned_to_spill() const { return transitioned_.load(); }
    std::uint64_t payload_bytes() const { return payloadBytes_.load(); }
    std::uint64_t async_merge_count() const { return asyncMergeCount_.load(); }
    std::uint64_t async_merged_records() const { return asyncMergedRecords_.load(); }

  private:
    struct RamSegment {
        std::uint64_t sequence = 0;
        std::uint32_t worker = 0;
        std::uint32_t level = 0;
        std::vector<std::uint8_t> bytes;
    };
    struct RamBucket {
        mutable std::mutex mutex;
        std::vector<RamSegment> segments;
    };
    struct SpillSegment {
        SpillSegment(std::uint64_t offsetIn, std::uint64_t bytesIn)
            : offset(offsetIn), bytes(bytesIn) {}
        std::uint64_t offset;
        std::uint64_t bytes;
    };
    struct RamMergeTask {
        std::uint32_t bucket = 0;
        std::uint32_t level = 0;
        std::vector<RamSegment> runs;
    };
    enum class Backend { Ram, Transitioning, Spill };

    bool validate_config(std::string *error) const;
    bool ensure_spill_files(std::string *error);
    bool append_ram(std::uint32_t workerIndex, std::uint32_t bucketIndex,
                    std::vector<std::uint8_t> bytes, std::string *error);
    bool append_spill(std::uint32_t bucketIndex,
                      const std::vector<std::uint8_t> &bytes,
                      std::string *error);
    bool transition_to_spill(std::string *error);
    bool finalize_spill(std::string *error);
    bool schedule_ram_merge_locked(std::uint32_t bucketIndex);
    void ram_merge_worker();
    RamSegment merge_ram_runs(RamMergeTask *task) const;
    bool wait_for_ram_merges(std::string *error);
    void stop_ram_merge_workers();
    std::string spill_path(std::uint32_t bucketIndex) const;

    Config config_;
    std::vector<RamBucket> ramBuckets_;
    std::unique_ptr<std::atomic<std::uint64_t>[]> ramSequences_;
    std::unique_ptr<std::atomic<std::uint64_t>[]> spillOffsets_;
    std::unique_ptr<std::atomic<std::uint64_t>[]> spillRecordCounts_;
    std::unique_ptr<std::mutex[]> spillClaimMutexes_;
    std::unique_ptr<std::uint64_t[]> spillChecksums_;
    std::vector<std::vector<SpillSegment> > spillSegments_;
    std::vector<int> spillFds_;
    mutable std::mutex spillInitMutex_;
    std::atomic<bool> spillFilesReady_{false};

    mutable std::mutex stateMutex_;
    mutable std::condition_variable stateCv_;
    Backend backend_ = Backend::Ram;
    std::uint32_t activeAppends_ = 0;
    bool finalizing_ = false;
    bool finalized_ = false;
    bool failed_ = false;
    std::string failureMessage_;
    std::atomic<bool> transitioned_{false};
    std::atomic<std::uint64_t> payloadBytes_{0};
    std::atomic<std::uint32_t> nextBucketClaim_{0};

    std::vector<std::uint8_t> ramMergeInFlight_;
    std::mutex ramMergeMutex_;
    std::condition_variable ramMergeCv_;
    std::deque<RamMergeTask> ramMergeQueue_;
    std::vector<std::thread> ramMergeWorkers_;
    std::uint64_t ramMergeOutstanding_ = 0;
    bool ramMergeStopping_ = false;
    bool ramMergeFailed_ = false;
    std::string ramMergeFailureMessage_;
    std::atomic<std::uint64_t> asyncMergeCount_{0};
    std::atomic<std::uint64_t> asyncMergedRecords_{0};
};

} // namespace solo
} // namespace star

#endif
