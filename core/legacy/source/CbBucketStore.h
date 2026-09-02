#ifndef H_CbBucketStore
#define H_CbBucketStore

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
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
    std::uint32_t cb_index() const;
    std::uint32_t umi24() const;
    std::uint16_t gene15() const;
    std::uint8_t tag5() const;
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
    bool load_bucket_bytes(std::uint32_t bucketIndex,
                           std::vector<std::uint8_t> *bytes,
                           std::string *error) const;

    void reset_bucket_claims();
    bool claim_bucket(std::uint32_t *bucketIndex);

    bool using_spill() const;
    bool transitioned_to_spill() const { return transitioned_.load(); }
    std::uint64_t payload_bytes() const { return payloadBytes_.load(); }

  private:
    struct RamSegment {
        std::uint64_t sequence = 0;
        std::uint32_t worker = 0;
        std::vector<std::uint8_t> bytes;
    };
    struct RamBucket {
        mutable std::mutex mutex;
        std::vector<RamSegment> segments;
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
    std::string spill_path(std::uint32_t bucketIndex) const;

    Config config_;
    std::vector<RamBucket> ramBuckets_;
    std::unique_ptr<std::atomic<std::uint64_t>[]> ramSequences_;
    std::unique_ptr<std::atomic<std::uint64_t>[]> spillOffsets_;
    std::unique_ptr<std::atomic<std::uint64_t>[]> spillRecordCounts_;
    std::vector<int> spillFds_;
    mutable std::mutex spillInitMutex_;
    bool spillFilesReady_ = false;

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
};

} // namespace solo
} // namespace star

#endif
