#include "CbBucketStore.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace star {
namespace solo {
namespace {

const std::uint8_t kSpillMagic[8] = {'S', 'C', 'B', 'B', 'K', 'T', '1', 0};
const std::uint32_t kSpillSchema = 1;
const std::size_t kSpillHeaderBytes = 56;
const std::uint64_t kIncomplete = std::numeric_limits<std::uint64_t>::max();
const std::uint64_t kFnvOffset = 1469598103934665603ULL;
const std::uint64_t kFnvPrime = 1099511628211ULL;

void put_u32(std::uint8_t *output, std::uint32_t value)
{
    for (unsigned shift = 0; shift < 32; shift += 8)
        *output++ = static_cast<std::uint8_t>(value >> shift);
}

void put_u64(std::uint8_t *output, std::uint64_t value)
{
    for (unsigned shift = 0; shift < 64; shift += 8)
        *output++ = static_cast<std::uint8_t>(value >> shift);
}

std::uint32_t get_u32(const std::uint8_t *input)
{
    std::uint32_t value = 0;
    for (unsigned byte = 0; byte < 4; ++byte)
        value |= static_cast<std::uint32_t>(input[byte]) << (byte * 8);
    return value;
}

std::uint64_t get_u64(const std::uint8_t *input)
{
    std::uint64_t value = 0;
    for (unsigned byte = 0; byte < 8; ++byte)
        value |= static_cast<std::uint64_t>(input[byte]) << (byte * 8);
    return value;
}

std::uint64_t encoded_group_sort_key(const std::uint8_t *input)
{
    const std::uint64_t key = get_u64(input);
    return (key & 0xFFFFF00000000000ULL)
         | ((key & 0x1FULL) << 39)
         | (((key >> 5) & 0x7FFFULL) << 24)
         | ((key >> 20) & 0xFFFFFFULL);
}

std::uint64_t fnv_update(std::uint64_t hash, const std::uint8_t *data,
                         std::size_t bytes)
{
    for (std::size_t index = 0; index < bytes; ++index) {
        hash ^= data[index];
        hash *= kFnvPrime;
    }
    return hash;
}

bool set_error(std::string *error, const std::string &message)
{
    if (error != nullptr)
        *error = message;
    return false;
}

bool pwrite_all(int fd, const std::uint8_t *data, std::size_t bytes,
                std::uint64_t offset, std::string *error)
{
    std::size_t written = 0;
    while (written < bytes) {
        const ssize_t result = ::pwrite(fd, data + written, bytes - written,
                                        static_cast<off_t>(offset + written));
        if (result < 0) {
            if (errno == EINTR)
                continue;
            return set_error(error, std::string("CB bucket pwrite failed: ")
                                      + std::strerror(errno));
        }
        if (result == 0)
            return set_error(error, "CB bucket pwrite made no progress");
        written += static_cast<std::size_t>(result);
    }
    return true;
}

bool pread_all(int fd, std::uint8_t *data, std::size_t bytes,
               std::uint64_t offset, std::string *error)
{
    std::size_t read = 0;
    while (read < bytes) {
        const ssize_t result = ::pread(fd, data + read, bytes - read,
                                       static_cast<off_t>(offset + read));
        if (result < 0) {
            if (errno == EINTR)
                continue;
            return set_error(error, std::string("CB bucket pread failed: ")
                                      + std::strerror(errno));
        }
        if (result == 0)
            return set_error(error, "truncated CB bucket spill file");
        read += static_cast<std::size_t>(result);
    }
    return true;
}

std::vector<std::uint8_t> make_header(std::uint32_t bucketIndex,
                                      std::uint32_t bucketCount,
                                      std::uint32_t whitelistSize,
                                      std::uint64_t recordCount,
                                      std::uint64_t payloadBytes,
                                      std::uint64_t checksum)
{
    std::vector<std::uint8_t> header(kSpillHeaderBytes, 0);
    std::copy(kSpillMagic, kSpillMagic + sizeof(kSpillMagic), header.begin());
    put_u32(&header[8], kSpillSchema);
    put_u32(&header[12], static_cast<std::uint32_t>(kSpillHeaderBytes));
    put_u32(&header[16], static_cast<std::uint32_t>(PackedCbRecord::kSerializedBytes));
    put_u32(&header[20], bucketIndex);
    put_u32(&header[24], bucketCount);
    put_u32(&header[28], whitelistSize);
    put_u64(&header[32], recordCount);
    put_u64(&header[40], payloadBytes);
    put_u64(&header[48], checksum);
    return header;
}

} // namespace

PackedCbRecord PackedCbRecord::make(std::uint32_t cbIndex,
                                    std::uint32_t umiValue,
                                    std::uint16_t geneValue,
                                    std::uint8_t tagValue,
                                    std::uint32_t countValue,
                                    std::uint8_t flagsValue)
{
    if (cbIndex >= (1u << 20) || umiValue >= (1u << 24)
        || geneValue >= (1u << 15) || tagValue >= (1u << 5)
        || countValue >= (1u << 30) || flagsValue >= (1u << 2)) {
        throw std::out_of_range("CB bucket record field exceeds its packed width");
    }
    PackedCbRecord result;
    result.key = (static_cast<std::uint64_t>(cbIndex) << 44)
               | (static_cast<std::uint64_t>(umiValue) << 20)
               | (static_cast<std::uint64_t>(geneValue) << 5)
               | tagValue;
    result.value = countValue | (static_cast<std::uint32_t>(flagsValue) << 30);
    return result;
}

// cb_index, umi24, gene15 and tag5 are defined inline in the header: they sit
// inside the bucket sort's comparator and must not be out-of-line calls.
std::uint32_t PackedCbRecord::count30() const { return value & 0x3FFFFFFFu; }
std::uint8_t PackedCbRecord::flags2() const { return (value >> 30) & 0x3u; }

void PackedCbRecord::encode(std::uint8_t output[kSerializedBytes]) const
{
    put_u64(output, key);
    put_u32(output + 8, value);
}

bool PackedCbRecord::decode(const std::uint8_t input[kSerializedBytes],
                            PackedCbRecord *record)
{
    if (record == nullptr)
        return false;
    record->key = get_u64(input);
    record->value = get_u32(input + 8);
    return true;
}

CbBucketStore::CbBucketStore(const Config &config)
    : config_(config), ramBuckets_(config.bucketCount),
      ramSequences_(new std::atomic<std::uint64_t>[config.bucketCount]),
      spillOffsets_(new std::atomic<std::uint64_t>[config.bucketCount]),
      spillRecordCounts_(new std::atomic<std::uint64_t>[config.bucketCount]),
      spillClaimMutexes_(new std::mutex[config.bucketCount]),
      spillChecksums_(new std::uint64_t[config.bucketCount]),
      spillSegments_(config.bucketCount),
      spillFds_(config.bucketCount, -1),
      backend_(config.mode == Mode::Spill ? Backend::Spill : Backend::Ram),
      ramMergeInFlight_(config.bucketCount, 0)
{
    std::string error;
    if (!validate_config(&error))
        throw std::invalid_argument(error);
    for (std::uint32_t bucket = 0; bucket < config_.bucketCount; ++bucket) {
        ramSequences_[bucket].store(0);
        spillOffsets_[bucket].store(0);
        spillRecordCounts_[bucket].store(0);
        spillChecksums_[bucket] = kFnvOffset;
    }
    if (config_.mode == Mode::Spill && !ensure_spill_files(&error))
        throw std::runtime_error(error);
    if (config_.mode != Mode::Spill) {
        try {
            for (std::uint32_t worker = 0; worker < config_.mergeWorkerCount; ++worker)
                ramMergeWorkers_.emplace_back(&CbBucketStore::ram_merge_worker, this);
        } catch (...) {
            stop_ram_merge_workers();
            throw;
        }
    }
}

CbBucketStore::~CbBucketStore()
{
    stop_ram_merge_workers();
    for (std::uint32_t bucket = 0; bucket < spillFds_.size(); ++bucket) {
        if (spillFds_[bucket] >= 0) {
            ::close(spillFds_[bucket]);
            ::unlink(spill_path(bucket).c_str());
        }
    }
}

bool CbBucketStore::validate_config(std::string *error) const
{
    if (config_.bucketCount == 0
        || (config_.bucketCount & (config_.bucketCount - 1)) != 0)
        return set_error(error, "CB bucket count must be a power of two");
    if (config_.whitelistSize == 0 || config_.whitelistSize > (1u << 20))
        return set_error(error, "CB bucket whitelist size must be in [1, 2^20]");
    if ((config_.mode == Mode::Spill || config_.mode == Mode::Auto)
        && config_.scratchDirectory.empty())
        return set_error(error, "CB bucket spill/auto mode requires a scratch directory");
    if (config_.filePrefix.empty())
        return set_error(error, "CB bucket spill file prefix is empty");
    if (config_.mergeWorkerCount > 0 && config_.mergeFanIn < 2)
        return set_error(error, "CB bucket asynchronous merge fan-in must be at least two");
    return true;
}

std::uint32_t CbBucketStore::bucket_for_cb(std::uint32_t cbIndex) const
{
    if (cbIndex >= config_.whitelistSize)
        throw std::out_of_range("CB index is outside the bucket whitelist domain");
    const std::uint64_t scaled = static_cast<std::uint64_t>(cbIndex)
                               * config_.bucketCount;
    const std::uint32_t bucket = static_cast<std::uint32_t>(scaled / config_.whitelistSize);
    return std::min(config_.bucketCount - 1, bucket);
}

std::string CbBucketStore::spill_path(std::uint32_t bucketIndex) const
{
    std::ostringstream path;
    path << config_.scratchDirectory << '/' << config_.filePrefix << '.'
         << std::setw(6) << std::setfill('0') << bucketIndex << ".cbb";
    return path.str();
}

bool CbBucketStore::ensure_spill_files(std::string *error)
{
    if (spillFilesReady_.load(std::memory_order_acquire))
        return true;
    std::lock_guard<std::mutex> lock(spillInitMutex_);
    if (spillFilesReady_.load(std::memory_order_relaxed))
        return true;
    if (::mkdir(config_.scratchDirectory.c_str(), 0700) != 0 && errno != EEXIST)
        return set_error(error, std::string("cannot create CB bucket spill directory: ")
                                  + std::strerror(errno));
    for (std::uint32_t bucket = 0; bucket < config_.bucketCount; ++bucket) {
        const std::string path = spill_path(bucket);
        const int fd = ::open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0600);
        if (fd < 0)
            return set_error(error, "cannot create CB bucket spill file " + path
                                      + ": " + std::strerror(errno));
        spillFds_[bucket] = fd;
        const std::vector<std::uint8_t> header = make_header(
            bucket, config_.bucketCount, config_.whitelistSize,
            kIncomplete, kIncomplete, 0);
        if (!pwrite_all(fd, header.data(), header.size(), 0, error))
            return false;
    }
    spillFilesReady_.store(true, std::memory_order_release);
    return true;
}

bool CbBucketStore::append_ram(std::uint32_t workerIndex,
                               std::uint32_t bucketIndex,
                               std::vector<std::uint8_t> bytes,
                               std::string *)
{
    RamSegment segment;
    segment.sequence = ramSequences_[bucketIndex].fetch_add(1);
    segment.worker = workerIndex;
    segment.bytes = std::move(bytes);
    std::lock_guard<std::mutex> lock(ramBuckets_[bucketIndex].mutex);
    ramBuckets_[bucketIndex].segments.push_back(std::move(segment));
    schedule_ram_merge_locked(bucketIndex);
    return true;
}

bool CbBucketStore::schedule_ram_merge_locked(std::uint32_t bucketIndex)
{
    if (ramMergeWorkers_.empty() || ramMergeInFlight_[bucketIndex])
        return false;

    RamBucket &bucket = ramBuckets_[bucketIndex];
    std::uint32_t selectedLevel = 0;
    bool found = false;
    for (std::uint32_t level = 0; level < 64 && !found; ++level) {
        std::uint32_t count = 0;
        for (const RamSegment &segment : bucket.segments) {
            if (segment.level == level && ++count == config_.mergeFanIn) {
                selectedLevel = level;
                found = true;
                break;
            }
        }
    }
    if (!found)
        return false;

    RamMergeTask task;
    task.bucket = bucketIndex;
    task.level = selectedLevel;
    task.runs.reserve(config_.mergeFanIn);
    std::vector<RamSegment> retained;
    retained.reserve(bucket.segments.size() - config_.mergeFanIn + 1);
    for (RamSegment &segment : bucket.segments) {
        if (segment.level == selectedLevel
            && task.runs.size() < config_.mergeFanIn) {
            task.runs.push_back(std::move(segment));
        } else {
            retained.push_back(std::move(segment));
        }
    }
    bucket.segments.swap(retained);
    ramMergeInFlight_[bucketIndex] = 1;
    {
        std::lock_guard<std::mutex> mergeLock(ramMergeMutex_);
        if (ramMergeStopping_ || ramMergeFailed_) {
            for (RamSegment &segment : task.runs)
                bucket.segments.push_back(std::move(segment));
            ramMergeInFlight_[bucketIndex] = 0;
            return false;
        }
        ramMergeQueue_.push_back(std::move(task));
        ++ramMergeOutstanding_;
    }
    ramMergeCv_.notify_one();
    return true;
}

CbBucketStore::RamSegment CbBucketStore::merge_ram_runs(
    RamMergeTask *task) const
{
    if (task == nullptr || task->runs.empty())
        throw std::runtime_error("CB bucket asynchronous merge received no runs");
    RamSegment output;
    output.sequence = task->runs.front().sequence;
    output.worker = task->runs.front().worker;
    output.level = task->level + 1;
    std::size_t totalBytes = 0;
    for (const RamSegment &run : task->runs) {
        if (run.bytes.size() % PackedCbRecord::kSerializedBytes != 0
            || run.bytes.size() > std::numeric_limits<std::size_t>::max() - totalBytes)
            throw std::runtime_error("invalid CB bucket run in asynchronous merge");
        totalBytes += run.bytes.size();
        output.sequence = std::min(output.sequence, run.sequence);
    }
    output.bytes.resize(totalBytes);

    const std::size_t runCount = task->runs.size();
    const std::size_t sentinel = runCount;
    std::vector<std::size_t> next(runCount, 0);
    std::vector<std::uint64_t> currentKey(runCount, 0);
    std::size_t leafCount = 1;
    while (leafCount < runCount)
        leafCount <<= 1;
    std::vector<std::size_t> tournament(leafCount * 2, sentinel);
    for (std::size_t run = 0; run < runCount; ++run) {
        if (!task->runs[run].bytes.empty()) {
            currentKey[run] = encoded_group_sort_key(task->runs[run].bytes.data());
            tournament[leafCount + run] = run;
        }
    }
    const auto winner = [&currentKey, sentinel](std::size_t left,
                                                std::size_t right) {
        if (left == sentinel) return right;
        if (right == sentinel) return left;
        if (currentKey[left] != currentKey[right])
            return currentKey[left] < currentKey[right] ? left : right;
        return std::min(left, right);
    };
    for (std::size_t node = leafCount; node-- > 1;)
        tournament[node] = winner(tournament[node * 2],
                                  tournament[node * 2 + 1]);

    std::size_t outputOffset = 0;
    while (tournament[1] != sentinel) {
        const std::size_t run = tournament[1];
        const std::size_t inputOffset =
            next[run] * PackedCbRecord::kSerializedBytes;
        std::memcpy(output.bytes.data() + outputOffset,
                    task->runs[run].bytes.data() + inputOffset,
                    PackedCbRecord::kSerializedBytes);
        outputOffset += PackedCbRecord::kSerializedBytes;
        ++next[run];
        const std::size_t leaf = leafCount + run;
        if (next[run] * PackedCbRecord::kSerializedBytes
            == task->runs[run].bytes.size()) {
            tournament[leaf] = sentinel;
        } else {
            currentKey[run] = encoded_group_sort_key(
                task->runs[run].bytes.data()
                + next[run] * PackedCbRecord::kSerializedBytes);
        }
        for (std::size_t node = leaf / 2; node > 0; node /= 2)
            tournament[node] = winner(tournament[node * 2],
                                      tournament[node * 2 + 1]);
    }
    return output;
}

void CbBucketStore::ram_merge_worker()
{
    for (;;) {
        RamMergeTask task;
        {
            std::unique_lock<std::mutex> lock(ramMergeMutex_);
            ramMergeCv_.wait(lock, [&]() {
                return ramMergeStopping_ || !ramMergeQueue_.empty();
            });
            if (ramMergeQueue_.empty()) {
                if (ramMergeStopping_)
                    return;
                continue;
            }
            task = std::move(ramMergeQueue_.front());
            ramMergeQueue_.pop_front();
        }

        bool ok = true;
        std::string failure;
        RamSegment merged;
        try {
            merged = merge_ram_runs(&task);
        } catch (const std::exception &error) {
            ok = false;
            failure = error.what();
        } catch (...) {
            ok = false;
            failure = "unknown CB bucket asynchronous merge failure";
        }

        {
            std::lock_guard<std::mutex> bucketLock(
                ramBuckets_[task.bucket].mutex);
            if (ok) {
                const std::uint64_t records =
                    merged.bytes.size() / PackedCbRecord::kSerializedBytes;
                ramBuckets_[task.bucket].segments.push_back(std::move(merged));
                asyncMergeCount_.fetch_add(1, std::memory_order_relaxed);
                asyncMergedRecords_.fetch_add(records, std::memory_order_relaxed);
            } else {
                for (RamSegment &run : task.runs)
                    ramBuckets_[task.bucket].segments.push_back(std::move(run));
            }
            ramMergeInFlight_[task.bucket] = 0;
            if (ok)
                schedule_ram_merge_locked(task.bucket);
        }
        {
            std::lock_guard<std::mutex> lock(ramMergeMutex_);
            if (!ok && !ramMergeFailed_) {
                ramMergeFailed_ = true;
                ramMergeFailureMessage_ = failure;
            }
            --ramMergeOutstanding_;
        }
        ramMergeCv_.notify_all();
    }
}

bool CbBucketStore::wait_for_ram_merges(std::string *error)
{
    if (ramMergeWorkers_.empty())
        return true;
    std::unique_lock<std::mutex> lock(ramMergeMutex_);
    ramMergeCv_.wait(lock, [&]() { return ramMergeOutstanding_ == 0; });
    if (ramMergeFailed_)
        return set_error(error, ramMergeFailureMessage_);
    return true;
}

void CbBucketStore::stop_ram_merge_workers()
{
    if (ramMergeWorkers_.empty())
        return;
    {
        std::lock_guard<std::mutex> lock(ramMergeMutex_);
        ramMergeStopping_ = true;
    }
    ramMergeCv_.notify_all();
    for (std::thread &worker : ramMergeWorkers_) {
        if (worker.joinable())
            worker.join();
    }
    ramMergeWorkers_.clear();
}

bool CbBucketStore::append_spill(std::uint32_t bucketIndex,
                                 const std::vector<std::uint8_t> &bytes,
                                 std::string *error)
{
    if (!ensure_spill_files(error))
        return false;
    std::uint64_t offset = 0;
    {
        // Claim the physical byte range and advance the checksum in that same
        // order. The pwrite remains outside the lock, so independent segments
        // and buckets can be written concurrently.
        std::lock_guard<std::mutex> lock(spillClaimMutexes_[bucketIndex]);
        offset = spillOffsets_[bucketIndex].load(std::memory_order_relaxed);
        spillOffsets_[bucketIndex].store(offset + bytes.size(),
                                         std::memory_order_relaxed);
        spillRecordCounts_[bucketIndex].fetch_add(
            bytes.size() / PackedCbRecord::kSerializedBytes,
            std::memory_order_relaxed);
        spillChecksums_[bucketIndex] = fnv_update(
            spillChecksums_[bucketIndex], bytes.data(), bytes.size());
        spillSegments_[bucketIndex].push_back(
            SpillSegment{offset, static_cast<std::uint64_t>(bytes.size())});
    }
    return pwrite_all(spillFds_[bucketIndex], bytes.data(), bytes.size(),
                      kSpillHeaderBytes + offset, error);
}

bool CbBucketStore::append_segment(std::uint32_t workerIndex,
                                   std::uint32_t bucketIndex,
                                   std::vector<PackedCbRecord> records,
                                   std::string *error)
{
    if (bucketIndex >= config_.bucketCount)
        return set_error(error, "CB bucket segment has an invalid bucket index");
    if (records.empty())
        return true;
    for (const PackedCbRecord &record : records) {
        if (record.cb_index() >= config_.whitelistSize
            || bucket_for_cb(record.cb_index()) != bucketIndex)
            return set_error(error, "CB bucket segment contains a record for another bucket");
    }
    // The producer owns this sealed vector, so ordering it here requires no
    // synchronization. Different workers sort their runs independently while
    // other workers continue producing; the collapse can later merge these
    // runs instead of re-sorting the complete bucket.
    std::sort(records.begin(), records.end(),
              [](const PackedCbRecord &left, const PackedCbRecord &right) {
                  return left.group_sort_key() < right.group_sort_key();
              });
    if (records.size() > std::numeric_limits<std::size_t>::max()
                           / PackedCbRecord::kSerializedBytes)
        return set_error(error, "CB bucket segment size overflow");
    std::vector<std::uint8_t> bytes(records.size() * PackedCbRecord::kSerializedBytes);
    for (std::size_t index = 0; index < records.size(); ++index)
        records[index].encode(bytes.data() + index * PackedCbRecord::kSerializedBytes);

    bool becomeTransitionOwner = false;
    Backend appendBackend = Backend::Ram;
    {
        std::unique_lock<std::mutex> lock(stateMutex_);
        stateCv_.wait(lock, [&]() { return backend_ != Backend::Transitioning; });
        if (failed_)
            return set_error(error, failureMessage_);
        if (finalizing_ || finalized_)
            return set_error(error, "cannot append to a finalized CB bucket store");
        const std::uint64_t bytesAfterAppend = payloadBytes_.fetch_add(bytes.size())
                                                   + bytes.size();
        if (config_.mode == Mode::Auto && backend_ == Backend::Ram
            && config_.memoryBudgetBytes > 0
            && bytesAfterAppend > config_.memoryBudgetBytes) {
            backend_ = Backend::Transitioning;
            stateCv_.wait(lock, [&]() { return activeAppends_ == 0; });
            becomeTransitionOwner = true;
        } else {
            appendBackend = backend_;
            ++activeAppends_;
        }
    }

    if (becomeTransitionOwner) {
        std::string transitionError;
        const bool ok = transition_to_spill(&transitionError);
        {
            std::lock_guard<std::mutex> lock(stateMutex_);
            if (ok) {
                backend_ = Backend::Spill;
                transitioned_.store(true);
                ++activeAppends_;
            } else {
                failed_ = true;
                failureMessage_ = transitionError;
            }
        }
        stateCv_.notify_all();
        if (!ok)
            return set_error(error, transitionError);
        appendBackend = Backend::Spill;
    }

    const bool result = appendBackend == Backend::Spill
                      ? append_spill(bucketIndex, bytes, error)
                      : append_ram(workerIndex, bucketIndex, std::move(bytes), error);
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        --activeAppends_;
        if (!result && !failed_) {
            failed_ = true;
            failureMessage_ = error == nullptr || error->empty()
                            ? "CB bucket append failed" : *error;
        }
    }
    stateCv_.notify_all();
    return result;
}

bool CbBucketStore::transition_to_spill(std::string *error)
{
    // No producer append remains when auto-transition elects its owner. Drain
    // any consolidation tasks before moving the RAM runs to their spill files.
    if (!wait_for_ram_merges(error))
        return false;
    if (!ensure_spill_files(error))
        return false;
    for (std::uint32_t bucket = 0; bucket < config_.bucketCount; ++bucket) {
        std::vector<RamSegment> segments;
        {
            std::lock_guard<std::mutex> lock(ramBuckets_[bucket].mutex);
            segments.swap(ramBuckets_[bucket].segments);
        }
        std::sort(segments.begin(), segments.end(),
                  [](const RamSegment &left, const RamSegment &right) {
                      return left.sequence < right.sequence;
                  });
        for (const RamSegment &segment : segments) {
            if (!append_spill(bucket, segment.bytes, error))
                return false;
        }
    }
    return true;
}

bool CbBucketStore::finalize_spill(std::string *error)
{
    if (!ensure_spill_files(error))
        return false;
    for (std::uint32_t bucket = 0; bucket < config_.bucketCount; ++bucket) {
        const std::uint64_t bytes = spillOffsets_[bucket].load();
        const std::vector<std::uint8_t> header = make_header(
            bucket, config_.bucketCount, config_.whitelistSize,
            spillRecordCounts_[bucket].load(), bytes,
            spillChecksums_[bucket]);
        if (!pwrite_all(spillFds_[bucket], header.data(), header.size(), 0, error))
            return false;
    }
    return true;
}

bool CbBucketStore::finalize(std::string *error)
{
    Backend backend;
    {
        std::unique_lock<std::mutex> lock(stateMutex_);
        stateCv_.wait(lock, [&]() {
            return backend_ != Backend::Transitioning && activeAppends_ == 0
                   && !finalizing_;
        });
        if (failed_)
            return set_error(error, failureMessage_);
        if (finalized_)
            return true;
        finalizing_ = true;
        backend = backend_;
    }
    bool result = true;
    if (backend == Backend::Ram)
        result = wait_for_ram_merges(error);
    if (result && backend == Backend::Spill)
        result = finalize_spill(error);
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        finalizing_ = false;
        if (result) {
            finalized_ = true;
        } else {
            failed_ = true;
            failureMessage_ = error == nullptr || error->empty()
                            ? "CB bucket finalization failed" : *error;
        }
    }
    stateCv_.notify_all();
    return result;
}

bool CbBucketStore::load_bucket_bytes(std::uint32_t bucketIndex,
                                      std::vector<std::uint8_t> *bytes,
                                      std::string *error) const
{
    if (bytes == nullptr || bucketIndex >= config_.bucketCount)
        return set_error(error, "invalid CB bucket load request");
    Backend backend;
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        if (!finalized_)
            return set_error(error, "CB bucket store must be finalized before loading");
        backend = backend_;
    }
    bytes->clear();
    if (backend == Backend::Ram) {
        std::vector<RamSegment> segments;
        {
            std::lock_guard<std::mutex> lock(ramBuckets_[bucketIndex].mutex);
            segments = ramBuckets_[bucketIndex].segments;
        }
        std::sort(segments.begin(), segments.end(),
                  [](const RamSegment &left, const RamSegment &right) {
                      return left.sequence < right.sequence;
                  });
        std::size_t total = 0;
        for (const RamSegment &segment : segments)
            total += segment.bytes.size();
        bytes->reserve(total);
        for (const RamSegment &segment : segments)
            bytes->insert(bytes->end(), segment.bytes.begin(), segment.bytes.end());
        return true;
    }

    const int fd = spillFds_[bucketIndex];
    if (fd < 0)
        return set_error(error, "CB bucket spill file is not open");
    std::vector<std::uint8_t> header(kSpillHeaderBytes);
    if (!pread_all(fd, header.data(), header.size(), 0, error))
        return false;
    if (!std::equal(kSpillMagic, kSpillMagic + sizeof(kSpillMagic), header.begin())
        || get_u32(&header[8]) != kSpillSchema
        || get_u32(&header[12]) != kSpillHeaderBytes
        || get_u32(&header[16]) != PackedCbRecord::kSerializedBytes
        || get_u32(&header[20]) != bucketIndex
        || get_u32(&header[24]) != config_.bucketCount
        || get_u32(&header[28]) != config_.whitelistSize)
        return set_error(error, "CB bucket spill header does not match the active store");
    const std::uint64_t records = get_u64(&header[32]);
    const std::uint64_t payload = get_u64(&header[40]);
    const std::uint64_t expectedChecksum = get_u64(&header[48]);
    if (records == kIncomplete || payload == kIncomplete
        || records > std::numeric_limits<std::size_t>::max()
                         / PackedCbRecord::kSerializedBytes
        || payload != records * PackedCbRecord::kSerializedBytes)
        return set_error(error, "CB bucket spill file is incomplete or has invalid counts");
    struct stat statBuffer;
    if (::fstat(fd, &statBuffer) != 0 || statBuffer.st_size < 0
        || static_cast<std::uint64_t>(statBuffer.st_size) != kSpillHeaderBytes + payload)
        return set_error(error, "CB bucket spill file size does not match its header");
    bytes->resize(static_cast<std::size_t>(payload));
    if (payload > 0 && !pread_all(fd, bytes->data(), bytes->size(),
                                  kSpillHeaderBytes, error))
        return false;
    if (fnv_update(kFnvOffset, bytes->data(), bytes->size()) != expectedChecksum)
        return set_error(error, "CB bucket spill checksum mismatch");
    return true;
}

bool CbBucketStore::load_bucket(std::uint32_t bucketIndex,
                                std::vector<PackedCbRecord> *records,
                                std::string *error) const
{
    if (records == nullptr)
        return set_error(error, "null CB bucket record output");
    std::vector<std::uint8_t> bytes;
    if (!load_bucket_bytes(bucketIndex, &bytes, error))
        return false;
    if (bytes.size() % PackedCbRecord::kSerializedBytes != 0)
        return set_error(error, "CB bucket payload is not record-aligned");
    records->resize(bytes.size() / PackedCbRecord::kSerializedBytes);
    for (std::size_t index = 0; index < records->size(); ++index) {
        if (!PackedCbRecord::decode(
                bytes.data() + index * PackedCbRecord::kSerializedBytes,
                &(*records)[index]))
            return set_error(error, "cannot decode CB bucket record");
    }
    return true;
}

bool CbBucketStore::load_sorted_segments(
    std::uint32_t bucketIndex,
    std::vector<std::vector<PackedCbRecord> > *segments,
    std::string *error) const
{
    if (segments == nullptr || bucketIndex >= config_.bucketCount)
        return set_error(error, "invalid CB bucket segment load request");
    Backend backend;
    {
        std::lock_guard<std::mutex> lock(stateMutex_);
        if (!finalized_)
            return set_error(error, "CB bucket store must be finalized before loading");
        backend = backend_;
    }
    segments->clear();

    const auto decodeRun = [error](const std::uint8_t *bytes,
                                   std::size_t byteCount,
                                   std::vector<PackedCbRecord> *run) {
        if (byteCount % PackedCbRecord::kSerializedBytes != 0)
            return set_error(error, "CB bucket segment is not record-aligned");
        run->resize(byteCount / PackedCbRecord::kSerializedBytes);
        for (std::size_t index = 0; index < run->size(); ++index) {
            if (!PackedCbRecord::decode(
                    bytes + index * PackedCbRecord::kSerializedBytes,
                    &(*run)[index]))
                return set_error(error, "cannot decode CB bucket segment record");
        }
        return true;
    };

    if (backend == Backend::Ram) {
        std::vector<RamSegment> encoded;
        {
            std::lock_guard<std::mutex> lock(ramBuckets_[bucketIndex].mutex);
            encoded = ramBuckets_[bucketIndex].segments;
        }
        std::sort(encoded.begin(), encoded.end(),
                  [](const RamSegment &left, const RamSegment &right) {
                      return left.sequence < right.sequence;
                  });
        segments->resize(encoded.size());
        for (std::size_t index = 0; index < encoded.size(); ++index) {
            if (!decodeRun(encoded[index].bytes.data(), encoded[index].bytes.size(),
                           &(*segments)[index]))
                return false;
        }
        return true;
    }

    // The spill payload remains one checksummed file per bucket. Segment
    // boundaries are lightweight in-memory metadata, so load/validation keeps
    // the existing schema while still recovering the independently sorted
    // runs for the k-way merge.
    std::vector<std::uint8_t> payload;
    if (!load_bucket_bytes(bucketIndex, &payload, error))
        return false;
    std::vector<SpillSegment> encoded;
    {
        std::lock_guard<std::mutex> lock(spillClaimMutexes_[bucketIndex]);
        encoded = spillSegments_[bucketIndex];
    }
    std::sort(encoded.begin(), encoded.end(),
              [](const SpillSegment &left, const SpillSegment &right) {
                  return left.offset < right.offset;
              });
    std::uint64_t expectedOffset = 0;
    segments->resize(encoded.size());
    for (std::size_t index = 0; index < encoded.size(); ++index) {
        const SpillSegment &segment = encoded[index];
        if (segment.offset != expectedOffset
            || expectedOffset > payload.size()
            || segment.bytes > payload.size() - expectedOffset)
            return set_error(error, "CB bucket spill segment metadata is inconsistent");
        if (!decodeRun(payload.data() + segment.offset,
                       static_cast<std::size_t>(segment.bytes),
                       &(*segments)[index]))
            return false;
        expectedOffset += segment.bytes;
    }
    if (expectedOffset != payload.size())
        return set_error(error, "CB bucket spill segment metadata is incomplete");
    return true;
}

void CbBucketStore::reset_bucket_claims()
{
    nextBucketClaim_.store(0);
}

bool CbBucketStore::claim_bucket(std::uint32_t *bucketIndex)
{
    if (bucketIndex == nullptr)
        return false;
    const std::uint32_t claimed = nextBucketClaim_.fetch_add(1);
    if (claimed >= config_.bucketCount)
        return false;
    *bucketIndex = claimed;
    return true;
}

bool CbBucketStore::using_spill() const
{
    std::lock_guard<std::mutex> lock(stateMutex_);
    return backend_ == Backend::Spill;
}

} // namespace solo
} // namespace star
