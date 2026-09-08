#ifndef STAR_FLEX_DECISION_SIDECAR_H
#define STAR_FLEX_DECISION_SIDECAR_H

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

#include "FlexHashScreen.h"

namespace flex_decision_sidecar {

static const std::uint16_t kSchemaVersion = 1;
static const std::uint16_t kHeaderBytes = 512;
static const std::uint16_t kRecordBytes = 48;
static const std::uint64_t kMissingLaneOrdinal = UINT64_MAX;
static const std::uint32_t kMissingLane = UINT32_MAX;

enum StatusFlag : std::uint32_t {
    kRecordPresent       = 1u << 0,
    kNameHashPresent     = 1u << 1,
    kSingleNAttempted    = 1u << 2,
    kSingleNResolved     = 1u << 3,
    kSampleChecked       = 1u << 4,
    kSampleMatched       = 1u << 5,
    kSampleRejected      = 1u << 6,
    kAlignmentHandoff    = 1u << 7,
    kAlignmentRan        = 1u << 8,
    kAlignmentResolved   = 1u << 9,
    kAlignmentRejected   = 1u << 10,
    kAlignmentProbe      = 1u << 11,
    kAlignmentGenomic    = 1u << 12,
    kNoAlignDropped      = 1u << 13,
    kCacheTerminal       = 1u << 14
};

enum FinalReason : std::uint8_t {
    kReasonNone = 0,
    kReasonCacheKeep = 1,
    kReasonCacheDeny = 2,
    kReasonSampleTagReject = 3,
    kReasonCacheMissNoAlign = 4,
    kReasonAlignmentNoCandidates = 5,
    kReasonAlignmentConflict = 6,
    kReasonAlignmentProbe = 7,
    kReasonAlignmentGenomic = 8
};

struct Record {
    std::uint64_t qnameHash = 0;
    std::uint64_t laneOrdinal = kMissingLaneOrdinal;
    std::uint32_t laneIndex = kMissingLane;
    std::uint32_t statusFlags = 0;
    std::uint32_t reserved32 = 0;
    std::uint16_t geneIdx15 = 0;
    std::uint8_t cacheAction = FlexHashScreenDecision::Disabled;
    std::uint8_t cacheClass = 0xFF;
    std::uint8_t matchedCacheClass = 0xFF;
    std::uint8_t negativeCode = 0;
    std::uint8_t sampleToken = 0xFF;
    std::int8_t hashOffset = 0;
    std::uint8_t probeRegion = 0;
    std::uint8_t finalReason = kReasonNone;
};

struct Header {
    std::uint16_t schemaVersion = 0;
    std::uint16_t headerBytes = 0;
    std::uint16_t recordBytes = 0;
    bool complete = false;
    std::uint64_t totalReads = 0;
    std::uint64_t recordsWritten = 0;
    std::string starSuiteVersion;
    std::string sourceRevision;
    std::string cachePath;
};

struct WriterConfig {
    std::string path;
    std::string starSuiteVersion;
    std::string sourceRevision;
    std::string cachePath;
};

class Writer {
  public:
    Writer();
    ~Writer();
    Writer(const Writer &) = delete;
    Writer &operator=(const Writer &) = delete;

    bool open(const WriterConfig &config, std::string &error);
    bool recordTriage(std::uint64_t ordinal, std::uint32_t laneIndex,
                      std::uint64_t laneOrdinal, const char *qname,
                      std::size_t qnameLength,
                      const FlexHashScreenDecision &decision,
                      bool sampleChecked, bool sampleMatched,
                      std::uint8_t sampleToken, bool alignmentHandoff,
                      bool noAlignDropped, std::string &error);
    bool recordAlignment(std::uint64_t ordinal, bool resolved,
                         bool genomic, std::uint16_t geneIdx15,
                         FinalReason reason, std::string &error);
    bool finalize(std::uint64_t totalReads, std::string &error);
    bool isOpen() const { return fd_ >= 0; }

  private:
    bool fail(const std::string &message, std::string &error);
    bool readRecord(std::uint64_t ordinal, Record &record, std::string &error);
    bool writeRecord(std::uint64_t ordinal, const Record &record,
                     bool wasPresent, std::string &error);
    bool writeHeader(bool complete, std::uint64_t totalReads,
                     std::uint64_t recordsWritten, std::string &error);

    WriterConfig config_;
    int fd_;
    std::string temporaryPath_;
    std::atomic<std::uint64_t> recordsWritten_;
    std::atomic<std::uint64_t> maximumOrdinalExclusive_;
    std::atomic<bool> failed_;
    std::mutex locks_[64];
    std::mutex errorMutex_;
    std::string writeError_;
};

class Reader {
  public:
    Reader();
    ~Reader();
    Reader(const Reader &) = delete;
    Reader &operator=(const Reader &) = delete;

    bool open(const std::string &path, std::string &error);
    bool read(std::uint64_t ordinal, Record &record, std::string &error) const;
    bool validateAll(std::string &error) const;
    const Header &header() const { return header_; }

  private:
    int fd_;
    Header header_;
};

std::uint64_t normalizedReadNameHash(const char *qname, std::size_t length);
std::vector<unsigned char> encodeRecord(const Record &record);
bool decodeRecord(const unsigned char *bytes, std::size_t size,
                  Record &record, std::string &error);
bool validateRecord(const Record &record, std::string &error);
const char *cacheActionName(std::uint8_t action);
const char *cacheClassName(std::uint8_t cacheClass);
const char *finalReasonName(std::uint8_t reason);

} // namespace flex_decision_sidecar

#endif
