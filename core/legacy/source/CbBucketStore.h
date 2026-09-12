#ifndef H_CbBucketStore
#define H_CbBucketStore

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <memory>
#include <algorithm>
#include <stdexcept>
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

    // The enclosing run has already been checked for 12-byte alignment.
    // memcpy permits unaligned source addresses on every target architecture.
    static PackedCbRecord from_encoded(const std::uint8_t *input) {
        PackedCbRecord record;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        std::memcpy(&record.key, input, sizeof(record.key));
        std::memcpy(&record.value, input + 8, sizeof(record.value));
#else
        for (unsigned byte = 0; byte < 8; ++byte)
            record.key |= static_cast<std::uint64_t>(input[byte]) << (byte * 8);
        for (unsigned byte = 0; byte < 4; ++byte)
            record.value |= static_cast<std::uint32_t>(input[8 + byte]) << (byte * 8);
#endif
        return record;
    }
    void encode(std::uint8_t output[kSerializedBytes]) const;
    static bool decode(const std::uint8_t input[kSerializedBytes],
                       PackedCbRecord *record);
};

// RAM-only run encoding. Compact words are
// [local CB12][TAG5][GENE15][UMI24][REGION2][COUNT6]. Count 63 escapes to
// a sorted, contiguous overflow vector; counts are never narrowed or capped.
// Spill files and the legacy byte API keep their original 12-byte schema.
struct EncodedCbRun {
    struct OverflowCount { std::uint64_t record; std::uint32_t count; };
    std::vector<std::uint8_t> bytes;
    std::vector<OverflowCount> overflow;
    std::uint32_t cbBase = 0;
    bool compact = false;

    std::size_t record_bytes() const { return compact ? 8 : PackedCbRecord::kSerializedBytes; }
    std::size_t record_count() const { return bytes.size() / record_bytes(); }
    static std::uint64_t word_at(const std::uint8_t* source) {
        std::uint64_t word = 0;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        std::memcpy(&word, source, sizeof(word));
#else
        for (unsigned i = 0; i < 8; ++i) word |= std::uint64_t(source[i]) << (8 * i);
#endif
        return word;
    }
    static void put_word(std::uint8_t* target, std::uint64_t word) {
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        std::memcpy(target, &word, sizeof(word));
#else
        for (unsigned i = 0; i < 8; ++i) target[i] = word >> (8 * i);
#endif
    }
    std::uint32_t overflow_count(std::size_t index) const {
        const auto found = std::lower_bound(overflow.begin(), overflow.end(), index,
            [](const OverflowCount& item, std::size_t at) { return item.record < at; });
        if (found == overflow.end() || found->record != index)
            throw std::runtime_error("Missing compact CB record count overflow");
        return found->count;
    }
    std::uint64_t sort_key(std::size_t index) const {
        if (compact) return (std::uint64_t(cbBase) << 44) + (word_at(bytes.data() + index * 8) >> 8);
        return PackedCbRecord::from_encoded(bytes.data() + index * 12).group_sort_key();
    }
    PackedCbRecord record_at(std::size_t index) const {
        if (!compact) return PackedCbRecord::from_encoded(bytes.data() + index * 12);
        const std::uint64_t word = word_at(bytes.data() + index * 8);
        const std::uint32_t count = (word & 63) == 63 ? overflow_count(index) : word & 63;
        PackedCbRecord record;
        record.key = (std::uint64_t(cbBase + (word >> 52)) << 44)
                   | (((word >> 8) & 0xFFFFFFu) << 20)
                   | (((word >> 32) & 0x7FFFu) << 5) | ((word >> 47) & 31);
        record.value = count | (std::uint32_t((word >> 6) & 3) << 30);
        return record;
    }
    // Called in increasing index order after bytes has been sized once.
    void set_record(std::size_t index, const PackedCbRecord& record) {
        if (!compact) { record.encode(bytes.data() + index * 12); return; }
        const std::uint32_t cb = record.cb_index();
        if (cb < cbBase || cb - cbBase >= 4096)
            throw std::out_of_range("Compact CB record exceeds its bucket-relative range");
        const std::uint32_t count = record.value & 0x3FFFFFFFu;
        const std::uint64_t word = (std::uint64_t(cb - cbBase) << 52)
            | (std::uint64_t(record.tag5()) << 47) | (std::uint64_t(record.gene15()) << 32)
            | (std::uint64_t(record.umi24()) << 8) | (std::uint64_t(record.value >> 30) << 6)
            | std::min<std::uint32_t>(count, 63);
        put_word(bytes.data() + index * 8, word);
        if (count >= 63) overflow.push_back({index, count});
    }
    std::vector<std::uint8_t> legacy_bytes() const {
        if (!compact) return bytes;
        std::vector<std::uint8_t> result(record_count() * PackedCbRecord::kSerializedBytes);
        for (std::size_t i = 0; i < record_count(); ++i)
            record_at(i).encode(result.data() + i * PackedCbRecord::kSerializedBytes);
        return result;
    }
    void release() {
        std::vector<std::uint8_t>().swap(bytes);
        std::vector<OverflowCount>().swap(overflow);
    }
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
    // Single-consumer tail API: transfer finalized RAM runs out of the store
    // and release each encoded run after decoding. Spill reads stay reusable.
    bool consume_sorted_segments(
        std::uint32_t bucketIndex,
        std::vector<std::vector<PackedCbRecord> > *segments,
        std::string *error);
    // Transfer the existing packed RAM storage without a decoded copy. Each
    // returned run contains a whole number of sorted 12-byte records. Spill
    // files are validated by the same checksummed reader and stay reusable.
    bool consume_encoded_segments(
        std::uint32_t bucketIndex,
        std::vector<std::vector<std::uint8_t>> *segments,
        std::string *error);
    // Production consumer: preserve compact RAM words and count overflows.
    bool consume_compact_segments(std::uint32_t bucketIndex,
                                 std::vector<EncodedCbRun>* segments,
                                 std::string* error);
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
        EncodedCbRun records;
    };
    struct RamBucket {
        mutable std::mutex mutex;
        std::vector<RamSegment> segments;
        bool consumed = false;
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
                    EncodedCbRun records, std::string *error);
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
