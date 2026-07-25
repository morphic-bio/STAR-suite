#ifndef STAR_SPATIAL_FEATURE_SIDECAR_H
#define STAR_SPATIAL_FEATURE_SIDECAR_H

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace spatial_feature_sidecar {

static const std::uint32_t kMissingGene = UINT32_MAX;
static const std::uint16_t kRecordPresent = 1u << 0;
static const std::uint16_t kMapped = 1u << 1;
static const std::uint16_t kUnmappedOrFiltered = 1u << 2;
static const std::uint16_t kUniqueGene = 1u << 3;
static const std::uint16_t kNoGene = 1u << 4;
static const std::uint16_t kMultiGeneRejected = 1u << 5;
static const std::uint16_t kSameGeneGenomicMultimapper = 1u << 6;
static const std::uint16_t kCrExonicRescue = 1u << 7;
static const std::uint16_t kCrIntronicFallback = 1u << 8;
static const std::uint16_t kOverlapShift = 9;
static const std::uint16_t kOverlapMask = 7u << kOverlapShift;
static const std::uint16_t kAlignmentEvidenceRejected = 1u << 12;
static const std::uint16_t kBestScoreNaDecoy = 1u << 13;
static const std::uint16_t kConflictingBestGenes = 1u << 14;
static const std::uint16_t kMultiGeneBestAlignment = 1u << 15;

static const std::uint16_t kSchemaVersion = 1;
static const std::uint16_t kHeaderBytes = 1024;
static const std::uint16_t kRecordBytes = 8;
static const std::uint32_t kNameDigestBlockReads = 16384;

struct Record {
    std::uint32_t geneIndex = kMissingGene;
    std::uint16_t statusFlags = 0;
    std::uint16_t reserved = 0;
};

struct Header {
    std::uint16_t schemaVersion = 0;
    std::uint16_t headerBytes = 0;
    std::uint16_t recordBytes = 0;
    bool complete = false;
    std::uint64_t totalReads = 0;
    std::uint64_t recordsWritten = 0;
    std::uint32_t featureCount = 0;
    std::int32_t strand = -1;
    bool crMultimapRescue = false;
    bool crIntronicFallback = false;
    std::string inputManifestSha256;
    std::string geneDictionarySha256;
    std::string starSuiteVersion;
    std::string sourceRevision;
    std::string featureType;
    std::string policy;
};

struct WriterConfig {
    std::string prefix;
    std::string starSuiteVersion;
    std::string sourceRevision;
    std::string featureType = "GeneFull";
    std::string policy;
    std::string inputManifest;
    std::int32_t strand = -1;
    bool crMultimapRescue = false;
    bool crIntronicFallback = false;
};

class Writer {
  public:
    Writer();
    ~Writer();
    Writer(const Writer &) = delete;
    Writer &operator=(const Writer &) = delete;

    bool open(const WriterConfig &config,
              const std::vector<std::string> &geneIds,
              const std::vector<std::string> &geneNames,
              std::string &error);
    bool write(std::uint64_t ordinal, std::uint32_t laneIndex,
               const std::string &normalizedReadName, const Record &record,
               std::string &error);
    bool finalize(std::uint64_t totalReads, std::string &error);
    bool isOpen() const { return binaryFd_ >= 0; }

  private:
    bool fail(const std::string &message, std::string &error);
    bool writeHeader(bool complete, std::uint64_t totalReads,
                     std::uint64_t recordsWritten, std::string &error);

    WriterConfig config_;
    Header header_;
    int binaryFd_;
    int namesFd_;
    std::string binaryTempPath_;
    std::string namesTempPath_;
    std::string featuresPath_;
    std::string digestPath_;
    std::string summaryPath_;
    std::atomic<std::uint64_t> recordsWritten_;
    std::atomic<std::uint64_t> firstOrdinal_;
    std::atomic<std::uint64_t> lastOrdinal_;
    std::atomic<bool> writeFailed_;
    std::string writeError_;
    std::map<std::uint16_t, std::uint64_t> flagCounts_;
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

std::string normalizeReadName(const std::string &value);
std::string sha256Hex(const std::string &value);
std::string sha256File(const std::string &path, std::string &error);
bool validateRecord(const Record &record, std::uint32_t featureCount,
                    std::string &error);
std::vector<unsigned char> encodeRecord(const Record &record);
bool decodeRecord(const unsigned char *bytes, std::size_t size,
                  Record &record, std::string &error);

} // namespace spatial_feature_sidecar

#endif
