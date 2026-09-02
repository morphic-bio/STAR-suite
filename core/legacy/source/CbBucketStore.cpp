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

std::uint32_t PackedCbRecord::cb_index() const { return (key >> 44) & 0xFFFFFu; }
std::uint32_t PackedCbRecord::umi24() const { return (key >> 20) & 0xFFFFFFu; }
std::uint16_t PackedCbRecord::gene15() const { return (key >> 5) & 0x7FFFu; }
std::uint8_t PackedCbRecord::tag5() const { return key & 0x1Fu; }
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
      spillFds_(config.bucketCount, -1),
      backend_(config.mode == Mode::Spill ? Backend::Spill : Backend::Ram)
{
    std::string error;
    if (!validate_config(&error))
        throw std::invalid_argument(error);
    for (std::uint32_t bucket = 0; bucket < config_.bucketCount; ++bucket) {
        ramSequences_[bucket].store(0);
        spillOffsets_[bucket].store(0);
        spillRecordCounts_[bucket].store(0);
    }
    if (config_.mode == Mode::Spill && !ensure_spill_files(&error))
        throw std::runtime_error(error);
}

CbBucketStore::~CbBucketStore()
{
    for (int fd : spillFds_) {
        if (fd >= 0)
            ::close(fd);
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
    std::lock_guard<std::mutex> lock(spillInitMutex_);
    if (spillFilesReady_)
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
    spillFilesReady_ = true;
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
    return true;
}

bool CbBucketStore::append_spill(std::uint32_t bucketIndex,
                                 const std::vector<std::uint8_t> &bytes,
                                 std::string *error)
{
    if (!ensure_spill_files(error))
        return false;
    const std::uint64_t offset = spillOffsets_[bucketIndex].fetch_add(bytes.size());
    spillRecordCounts_[bucketIndex].fetch_add(
        bytes.size() / PackedCbRecord::kSerializedBytes);
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
    std::vector<std::uint8_t> buffer(1u << 20);
    for (std::uint32_t bucket = 0; bucket < config_.bucketCount; ++bucket) {
        const std::uint64_t bytes = spillOffsets_[bucket].load();
        std::uint64_t checksum = kFnvOffset;
        std::uint64_t offset = 0;
        while (offset < bytes) {
            const std::size_t chunk = static_cast<std::size_t>(
                std::min<std::uint64_t>(buffer.size(), bytes - offset));
            if (!pread_all(spillFds_[bucket], buffer.data(), chunk,
                           kSpillHeaderBytes + offset, error))
                return false;
            checksum = fnv_update(checksum, buffer.data(), chunk);
            offset += chunk;
        }
        const std::vector<std::uint8_t> header = make_header(
            bucket, config_.bucketCount, config_.whitelistSize,
            spillRecordCounts_[bucket].load(), bytes, checksum);
        if (!pwrite_all(spillFds_[bucket], header.data(), header.size(), 0, error))
            return false;
        if (::fsync(spillFds_[bucket]) != 0)
            return set_error(error, std::string("cannot finalize CB bucket spill file: ")
                                      + std::strerror(errno));
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
    const bool result = backend == Backend::Spill ? finalize_spill(error) : true;
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
