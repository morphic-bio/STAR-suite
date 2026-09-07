#include "FlexDecisionSidecar.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace flex_decision_sidecar {
namespace {

const unsigned char kMagic[8] = {'F', 'L', 'X', 'D', 'E', 'C', '1', 0};
const std::uint32_t kEndianMarker = 0x01020304u;

void put16(unsigned char *p, std::uint16_t value)
{
    p[0] = static_cast<unsigned char>(value);
    p[1] = static_cast<unsigned char>(value >> 8);
}

void put32(unsigned char *p, std::uint32_t value)
{
    for (int i = 0; i < 4; ++i) p[i] = static_cast<unsigned char>(value >> (8 * i));
}

void put64(unsigned char *p, std::uint64_t value)
{
    for (int i = 0; i < 8; ++i) p[i] = static_cast<unsigned char>(value >> (8 * i));
}

std::uint16_t get16(const unsigned char *p)
{
    return static_cast<std::uint16_t>(p[0])
        | static_cast<std::uint16_t>(p[1]) << 8;
}

std::uint32_t get32(const unsigned char *p)
{
    std::uint32_t value = 0;
    for (int i = 3; i >= 0; --i) value = (value << 8) | p[i];
    return value;
}

std::uint64_t get64(const unsigned char *p)
{
    std::uint64_t value = 0;
    for (int i = 7; i >= 0; --i) value = (value << 8) | p[i];
    return value;
}

void putText(unsigned char *header, std::size_t offset, std::size_t size,
             const std::string &value)
{
    const std::size_t copied = std::min(size - 1, value.size());
    std::memcpy(header + offset, value.data(), copied);
}

std::string getText(const unsigned char *header, std::size_t offset, std::size_t size)
{
    std::size_t length = 0;
    while (length < size && header[offset + length] != 0) ++length;
    return std::string(reinterpret_cast<const char *>(header + offset), length);
}

bool writeAt(int fd, const unsigned char *data, std::size_t size, off_t offset,
             std::string &error)
{
    std::size_t done = 0;
    while (done < size) {
        const ssize_t result = ::pwrite(fd, data + done, size - done,
                                        offset + static_cast<off_t>(done));
        if (result < 0 && errno == EINTR) continue;
        if (result <= 0) {
            error = std::string("pwrite failed: ") + std::strerror(errno);
            return false;
        }
        done += static_cast<std::size_t>(result);
    }
    return true;
}

bool readAtAllowMissing(int fd, unsigned char *data, std::size_t size,
                        off_t offset, std::string &error)
{
    std::memset(data, 0, size);
    std::size_t done = 0;
    while (done < size) {
        const ssize_t result = ::pread(fd, data + done, size - done,
                                       offset + static_cast<off_t>(done));
        if (result < 0 && errno == EINTR) continue;
        if (result < 0) {
            error = std::string("pread failed: ") + std::strerror(errno);
            return false;
        }
        if (result == 0) return true;
        done += static_cast<std::size_t>(result);
    }
    return true;
}

bool readAt(int fd, unsigned char *data, std::size_t size, off_t offset,
            std::string &error)
{
    std::size_t done = 0;
    while (done < size) {
        const ssize_t result = ::pread(fd, data + done, size - done,
                                       offset + static_cast<off_t>(done));
        if (result < 0 && errno == EINTR) continue;
        if (result <= 0) {
            error = result == 0 ? "unexpected end of file"
                                : std::string("pread failed: ") + std::strerror(errno);
            return false;
        }
        done += static_cast<std::size_t>(result);
    }
    return true;
}

bool decodeHeader(const unsigned char *bytes, std::size_t size, Header &header,
                  std::string &error)
{
    if (size < kHeaderBytes || std::memcmp(bytes, kMagic, sizeof(kMagic)) != 0) {
        error = "invalid Flex decision sidecar magic or short header";
        return false;
    }
    header.schemaVersion = get16(bytes + 8);
    header.headerBytes = get16(bytes + 10);
    if (get32(bytes + 12) != kEndianMarker) {
        error = "invalid Flex decision sidecar byte-order marker";
        return false;
    }
    header.recordBytes = get16(bytes + 16);
    header.complete = get16(bytes + 18) == 1;
    header.totalReads = get64(bytes + 24);
    header.recordsWritten = get64(bytes + 32);
    header.starSuiteVersion = getText(bytes, 64, 64);
    header.sourceRevision = getText(bytes, 128, 64);
    header.cachePath = getText(bytes, 192, 256);
    if (header.schemaVersion != kSchemaVersion || header.headerBytes != kHeaderBytes
        || header.recordBytes != kRecordBytes) {
        error = "unsupported Flex decision sidecar schema/header/record size";
        return false;
    }
    return true;
}

} // namespace

std::uint64_t normalizedReadNameHash(const char *qname, std::size_t length)
{
    if (qname == nullptr) return 0;
    std::size_t begin = length > 0 && qname[0] == '@' ? 1 : 0;
    std::size_t end = begin;
    while (end < length && qname[end] != '\0' && qname[end] != ' '
           && qname[end] != '\t' && qname[end] != '\r' && qname[end] != '\n') {
        ++end;
    }
    if (end >= begin + 2 && qname[end - 2] == '/'
        && (qname[end - 1] == '1' || qname[end - 1] == '2')) {
        end -= 2;
    }
    if (end == begin) return 0;
    // Standard FNV-1a-64 offset basis and prime. Keeping the published
    // constants makes name hashes interoperable with independent audit tools.
    std::uint64_t hash = 14695981039346656037ULL;
    for (std::size_t i = begin; i < end; ++i) {
        hash ^= static_cast<unsigned char>(qname[i]);
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::vector<unsigned char> encodeRecord(const Record &record)
{
    std::vector<unsigned char> bytes(kRecordBytes, 0);
    put64(bytes.data(), record.qnameHash);
    put64(bytes.data() + 8, record.laneOrdinal);
    put32(bytes.data() + 16, record.laneIndex);
    put32(bytes.data() + 20, record.statusFlags);
    put32(bytes.data() + 24, record.reserved32);
    put16(bytes.data() + 28, record.geneIdx15);
    bytes[30] = record.cacheAction;
    bytes[31] = record.cacheClass;
    bytes[32] = record.matchedCacheClass;
    bytes[33] = record.negativeCode;
    bytes[34] = record.sampleToken;
    bytes[35] = static_cast<unsigned char>(record.hashOffset);
    bytes[36] = record.probeRegion;
    bytes[37] = record.finalReason;
    return bytes;
}

bool decodeRecord(const unsigned char *bytes, std::size_t size,
                  Record &record, std::string &error)
{
    if (size != kRecordBytes) {
        error = "Flex decision sidecar record has wrong size";
        return false;
    }
    record.qnameHash = get64(bytes);
    record.laneOrdinal = get64(bytes + 8);
    record.laneIndex = get32(bytes + 16);
    record.statusFlags = get32(bytes + 20);
    record.reserved32 = get32(bytes + 24);
    record.geneIdx15 = get16(bytes + 28);
    record.cacheAction = bytes[30];
    record.cacheClass = bytes[31];
    record.matchedCacheClass = bytes[32];
    record.negativeCode = bytes[33];
    record.sampleToken = bytes[34];
    record.hashOffset = static_cast<std::int8_t>(bytes[35]);
    record.probeRegion = bytes[36];
    record.finalReason = bytes[37];
    for (std::size_t i = 38; i < size; ++i) {
        if (bytes[i] != 0) {
            error = "Flex decision sidecar record has nonzero reserved bytes";
            return false;
        }
    }
    return true;
}

bool validateRecord(const Record &record, std::string &error)
{
    if (!(record.statusFlags & kRecordPresent)) {
        error = "Flex decision sidecar record is missing";
        return false;
    }
    if (record.reserved32 != 0) {
        error = "Flex decision sidecar record has nonzero reserved field";
        return false;
    }
    if (((record.statusFlags & kNameHashPresent) != 0) != (record.qnameHash != 0)) {
        error = "Flex decision sidecar name-hash flag is inconsistent";
        return false;
    }
    if ((record.statusFlags & kSampleMatched)
        && !(record.statusFlags & kSampleChecked)) {
        error = "Flex decision sidecar sample match was not checked";
        return false;
    }
    if ((record.statusFlags & kSampleRejected)
        && !(record.statusFlags & kSampleChecked)) {
        error = "Flex decision sidecar sample reject was not checked";
        return false;
    }
    if ((record.statusFlags & kAlignmentResolved)
        && !(record.statusFlags & kAlignmentRan)) {
        error = "Flex decision sidecar alignment resolution lacks alignment run";
        return false;
    }
    if ((record.statusFlags & kAlignmentRejected)
        && !(record.statusFlags & kAlignmentRan)) {
        error = "Flex decision sidecar alignment rejection lacks alignment run";
        return false;
    }
    return true;
}

Writer::Writer()
    : fd_(-1), recordsWritten_(0), maximumOrdinalExclusive_(0), failed_(false)
{}

Writer::~Writer()
{
    if (fd_ >= 0) ::close(fd_);
}

bool Writer::fail(const std::string &message, std::string &error)
{
    {
        std::lock_guard<std::mutex> lock(errorMutex_);
        if (!failed_.load(std::memory_order_relaxed)) writeError_ = message;
    }
    failed_.store(true, std::memory_order_release);
    error = message;
    return false;
}

bool Writer::writeHeader(bool complete, std::uint64_t totalReads,
                         std::uint64_t recordsWritten, std::string &error)
{
    unsigned char bytes[kHeaderBytes] = {};
    std::memcpy(bytes, kMagic, sizeof(kMagic));
    put16(bytes + 8, kSchemaVersion);
    put16(bytes + 10, kHeaderBytes);
    put32(bytes + 12, kEndianMarker);
    put16(bytes + 16, kRecordBytes);
    put16(bytes + 18, complete ? 1 : 0);
    put64(bytes + 24, totalReads);
    put64(bytes + 32, recordsWritten);
    putText(bytes, 64, 64, config_.starSuiteVersion);
    putText(bytes, 128, 64, config_.sourceRevision);
    putText(bytes, 192, 256, config_.cachePath);
    return writeAt(fd_, bytes, sizeof(bytes), 0, error);
}

bool Writer::open(const WriterConfig &config, std::string &error)
{
    if (fd_ >= 0) return fail("Flex decision sidecar writer is already open", error);
    if (config.path.empty() || config.path == "-") {
        return fail("Flex decision sidecar path is empty or disabled", error);
    }
    config_ = config;
    temporaryPath_ = config.path + ".tmp";
    fd_ = ::open(temporaryPath_.c_str(), O_RDWR | O_CREAT | O_EXCL, 0644);
    if (fd_ < 0) {
        return fail("cannot create " + temporaryPath_ + ": " + std::strerror(errno), error);
    }
    if (!writeHeader(false, 0, 0, error)) return fail(error, error);
    return true;
}

bool Writer::readRecord(std::uint64_t ordinal, Record &record, std::string &error)
{
    unsigned char bytes[kRecordBytes];
    const off_t offset = static_cast<off_t>(kHeaderBytes)
        + static_cast<off_t>(ordinal) * kRecordBytes;
    if (!readAtAllowMissing(fd_, bytes, sizeof(bytes), offset, error)) return false;
    if (get32(bytes + 20) == 0) {
        record = Record();
        return true;
    }
    return decodeRecord(bytes, sizeof(bytes), record, error);
}

bool Writer::writeRecord(std::uint64_t ordinal, const Record &record,
                         bool wasPresent, std::string &error)
{
    const std::vector<unsigned char> bytes = encodeRecord(record);
    const off_t offset = static_cast<off_t>(kHeaderBytes)
        + static_cast<off_t>(ordinal) * kRecordBytes;
    if (!writeAt(fd_, bytes.data(), bytes.size(), offset, error)) return false;
    if (!wasPresent) recordsWritten_.fetch_add(1, std::memory_order_relaxed);
    const std::uint64_t ordinalExclusive = ordinal + 1;
    std::uint64_t current = maximumOrdinalExclusive_.load(std::memory_order_relaxed);
    while (current < ordinalExclusive
           && !maximumOrdinalExclusive_.compare_exchange_weak(
               current, ordinalExclusive, std::memory_order_relaxed)) {}
    return true;
}

bool Writer::recordTriage(std::uint64_t ordinal, std::uint32_t laneIndex,
                          std::uint64_t laneOrdinal, const char *qname,
                          std::size_t qnameLength,
                          const FlexHashScreenDecision &decision,
                          bool sampleChecked, bool sampleMatched,
                          std::uint8_t sampleToken, bool alignmentHandoff,
                          bool noAlignDropped, std::string &error)
{
    if (failed_.load(std::memory_order_acquire)) {
        std::lock_guard<std::mutex> lock(errorMutex_);
        error = writeError_;
        return false;
    }
    std::lock_guard<std::mutex> lock(locks_[ordinal % 64]);
    Record record;
    if (!readRecord(ordinal, record, error)) return fail(error, error);
    const bool wasPresent = (record.statusFlags & kRecordPresent) != 0;
    record.qnameHash = normalizedReadNameHash(qname, qnameLength);
    record.laneIndex = laneIndex;
    record.laneOrdinal = laneOrdinal;
    record.statusFlags = kRecordPresent;
    if (record.qnameHash != 0) record.statusFlags |= kNameHashPresent;
    if (decision.singleN) record.statusFlags |= kSingleNAttempted;
    if (decision.singleN && decision.action == FlexHashScreenDecision::Keep)
        record.statusFlags |= kSingleNResolved;
    if (sampleChecked) {
        record.statusFlags |= kSampleChecked;
        record.statusFlags |= sampleMatched ? kSampleMatched : kSampleRejected;
    }
    if (alignmentHandoff) record.statusFlags |= kAlignmentHandoff;
    if (noAlignDropped) record.statusFlags |= kNoAlignDropped;
    const bool sampleRejected = sampleChecked && !sampleMatched;
    if (!sampleRejected && (decision.action == FlexHashScreenDecision::Keep
                            || decision.action == FlexHashScreenDecision::Deny)) {
        record.statusFlags |= kCacheTerminal;
    }
    record.geneIdx15 = sampleRejected ? 0 : decision.geneIdx15;
    record.cacheAction = decision.action;
    const bool cacheRecord = !sampleRejected
        && (decision.action == FlexHashScreenDecision::Keep
            || decision.action == FlexHashScreenDecision::Deny);
    record.cacheClass = cacheRecord ? decision.cacheClass : 0xFF;
    record.matchedCacheClass = sampleRejected ? 0xFF : decision.singleN
        ? decision.singleNCacheClass : record.cacheClass;
    record.negativeCode = sampleRejected ? 0 : decision.negativeCode;
    record.sampleToken = sampleToken;
    record.hashOffset = decision.offset;
    record.probeRegion = static_cast<std::uint8_t>(decision.probeRegion);
    if (sampleChecked && !sampleMatched) record.finalReason = kReasonSampleTagReject;
    else if (noAlignDropped) record.finalReason = kReasonCacheMissNoAlign;
    else if (decision.action == FlexHashScreenDecision::Keep) record.finalReason = kReasonCacheKeep;
    else if (decision.action == FlexHashScreenDecision::Deny) record.finalReason = kReasonCacheDeny;
    else record.finalReason = kReasonNone;
    if (!validateRecord(record, error)) return fail(error, error);
    if (!writeRecord(ordinal, record, wasPresent, error)) return fail(error, error);
    return true;
}

bool Writer::recordAlignment(std::uint64_t ordinal, bool resolved,
                             bool genomic, std::uint16_t geneIdx15,
                             FinalReason reason, std::string &error)
{
    if (failed_.load(std::memory_order_acquire)) {
        std::lock_guard<std::mutex> lock(errorMutex_);
        error = writeError_;
        return false;
    }
    std::lock_guard<std::mutex> lock(locks_[ordinal % 64]);
    Record record;
    if (!readRecord(ordinal, record, error)) return fail(error, error);
    const bool wasPresent = (record.statusFlags & kRecordPresent) != 0;
    if (!wasPresent) return fail("alignment update has no triage record", error);
    record.statusFlags |= kAlignmentRan;
    record.statusFlags &= ~(kAlignmentResolved | kAlignmentRejected
                            | kAlignmentProbe | kAlignmentGenomic);
    if (resolved) {
        record.statusFlags |= kAlignmentResolved;
        record.statusFlags |= genomic ? kAlignmentGenomic : kAlignmentProbe;
        record.geneIdx15 = geneIdx15;
    } else {
        record.statusFlags |= kAlignmentRejected;
        record.geneIdx15 = 0;
    }
    record.finalReason = reason;
    if (!validateRecord(record, error)) return fail(error, error);
    if (!writeRecord(ordinal, record, wasPresent, error)) return fail(error, error);
    return true;
}

bool Writer::finalize(std::uint64_t totalReads, std::string &error)
{
    if (fd_ < 0) return fail("Flex decision sidecar writer is not open", error);
    if (failed_.load(std::memory_order_acquire)) {
        std::lock_guard<std::mutex> lock(errorMutex_);
        error = writeError_;
        return false;
    }
    const std::uint64_t written = recordsWritten_.load(std::memory_order_relaxed);
    const std::uint64_t extent = maximumOrdinalExclusive_.load(
        std::memory_order_relaxed);
    if (written != totalReads || extent != totalReads) {
        std::ostringstream message;
        message << "Flex decision sidecar record count " << written
                << " and ordinal extent " << extent
                << " differ from input count " << totalReads;
        return fail(message.str(), error);
    }
    if (!writeHeader(true, totalReads, written, error)) return fail(error, error);
    if (::fsync(fd_) != 0) return fail("fsync failed: " + std::string(std::strerror(errno)), error);
    if (::close(fd_) != 0) {
        fd_ = -1;
        return fail("close failed: " + std::string(std::strerror(errno)), error);
    }
    fd_ = -1;
    if (::rename(temporaryPath_.c_str(), config_.path.c_str()) != 0) {
        return fail("cannot commit " + config_.path + ": " + std::strerror(errno), error);
    }
    return true;
}

Reader::Reader() : fd_(-1) {}
Reader::~Reader() { if (fd_ >= 0) ::close(fd_); }

bool Reader::open(const std::string &path, std::string &error)
{
    if (fd_ >= 0) ::close(fd_);
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        error = "cannot open " + path + ": " + std::strerror(errno);
        return false;
    }
    unsigned char bytes[kHeaderBytes];
    if (!readAt(fd_, bytes, sizeof(bytes), 0, error)
        || !decodeHeader(bytes, sizeof(bytes), header_, error)) return false;
    if (!header_.complete) {
        error = "Flex decision sidecar is incomplete";
        return false;
    }
    struct stat info;
    if (::fstat(fd_, &info) != 0) {
        error = "fstat failed: " + std::string(std::strerror(errno));
        return false;
    }
    const std::uint64_t expected = kHeaderBytes + header_.totalReads * kRecordBytes;
    if (static_cast<std::uint64_t>(info.st_size) != expected) {
        error = "Flex decision sidecar file size does not match header";
        return false;
    }
    return true;
}

bool Reader::read(std::uint64_t ordinal, Record &record, std::string &error) const
{
    if (fd_ < 0 || ordinal >= header_.totalReads) {
        error = "Flex decision sidecar ordinal is out of range";
        return false;
    }
    unsigned char bytes[kRecordBytes];
    const off_t offset = static_cast<off_t>(kHeaderBytes)
        + static_cast<off_t>(ordinal) * kRecordBytes;
    return readAt(fd_, bytes, sizeof(bytes), offset, error)
        && decodeRecord(bytes, sizeof(bytes), record, error);
}

bool Reader::validateAll(std::string &error) const
{
    if (header_.recordsWritten != header_.totalReads) {
        error = "Flex decision sidecar header reports missing records";
        return false;
    }
    for (std::uint64_t ordinal = 0; ordinal < header_.totalReads; ++ordinal) {
        Record record;
        if (!read(ordinal, record, error) || !validateRecord(record, error)) return false;
    }
    return true;
}

const char *cacheActionName(std::uint8_t action)
{
    switch (action) {
        case FlexHashScreenDecision::Disabled: return "DISABLED";
        case FlexHashScreenDecision::Pass: return "MISS";
        case FlexHashScreenDecision::Keep: return "KEEP";
        case FlexHashScreenDecision::Deny: return "DENY";
        default: return "UNKNOWN";
    }
}

const char *cacheClassName(std::uint8_t cacheClass)
{
    switch (cacheClass) {
        case FlexHashCacheH0: return "H0";
        case FlexHashCacheH1: return "H1";
        case FlexHashCacheNegative: return "NEGATIVE";
        case FlexHashCacheH2: return "H2";
        case FlexHashCacheH1X2: return "H1X2";
        case 0xFE: return "MIXED";
        case 0xFF: return ".";
        default: return "UNKNOWN";
    }
}

const char *finalReasonName(std::uint8_t reason)
{
    switch (reason) {
        case kReasonNone: return ".";
        case kReasonCacheKeep: return "CACHE_KEEP";
        case kReasonCacheDeny: return "CACHE_DENY";
        case kReasonSampleTagReject: return "SAMPLE_TAG_REJECT";
        case kReasonCacheMissNoAlign: return "CACHE_MISS_NO_ALIGN";
        case kReasonAlignmentNoCandidates: return "ALIGN_NO_CANDIDATES";
        case kReasonAlignmentConflict: return "ALIGN_CONFLICT";
        case kReasonAlignmentProbe: return "ALIGN_PROBE";
        case kReasonAlignmentGenomic: return "ALIGN_GENOMIC";
        default: return "UNKNOWN";
    }
}

} // namespace flex_decision_sidecar
