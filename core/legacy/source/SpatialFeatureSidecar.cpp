#include "SpatialFeatureSidecar.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <openssl/sha.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace spatial_feature_sidecar {
namespace {

const unsigned char kMagic[8] = {'S', 'F', 'S', 'I', 'D', 'E', '1', 0};
const std::uint32_t kEndianMarker = 0x01020304u;
const std::uint32_t kNamePresent = 0x4e414d45u;
std::mutex gWriteErrorMutex;

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

std::string hexDigest(const unsigned char digest[SHA256_DIGEST_LENGTH])
{
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (std::size_t i = 0; i < SHA256_DIGEST_LENGTH; ++i) {
        output << std::setw(2) << static_cast<unsigned int>(digest[i]);
    }
    return output.str();
}

bool hexToBytes(const std::string &text, unsigned char *out, std::size_t size)
{
    if (text.size() != size * 2) return false;
    for (std::size_t i = 0; i < size; ++i) {
        unsigned int value = 0;
        std::istringstream input(text.substr(i * 2, 2));
        input >> std::hex >> value;
        if (!input || !input.eof()) return false;
        out[i] = static_cast<unsigned char>(value);
    }
    return true;
}

std::string bytesToHex(const unsigned char *bytes, std::size_t size)
{
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (std::size_t i = 0; i < size; ++i) {
        output << std::setw(2) << static_cast<unsigned int>(bytes[i]);
    }
    return output.str();
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

std::uint64_t fnv1a64(const std::string &value)
{
    std::uint64_t hash = 1469598103934665603ULL;
    for (unsigned char c : value) {
        hash ^= c;
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::string jsonEscape(const std::string &value)
{
    std::ostringstream output;
    for (unsigned char c : value) {
        switch (c) {
            case '"': output << "\\\""; break;
            case '\\': output << "\\\\"; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default:
                if (c < 0x20) output << '?';
                else output << static_cast<char>(c);
        }
    }
    return output.str();
}

bool writeAtomicText(const std::string &path, const std::string &text, std::string &error)
{
    const std::string temporary = path + ".tmp";
    const int fd = ::open(temporary.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0644);
    if (fd < 0) {
        error = "cannot create " + temporary + ": " + std::strerror(errno);
        return false;
    }
    bool ok = writeAt(fd, reinterpret_cast<const unsigned char *>(text.data()), text.size(), 0, error);
    if (ok && ::fsync(fd) != 0) {
        error = "fsync failed for " + temporary + ": " + std::strerror(errno);
        ok = false;
    }
    if (::close(fd) != 0 && ok) {
        error = "close failed for " + temporary + ": " + std::strerror(errno);
        ok = false;
    }
    if (ok && ::rename(temporary.c_str(), path.c_str()) != 0) {
        error = "cannot commit " + path + ": " + std::strerror(errno);
        ok = false;
    }
    if (!ok) ::unlink(temporary.c_str());
    return ok;
}

bool decodeHeader(const unsigned char *bytes, std::size_t size, Header &header,
                  std::string &error)
{
    if (size < kHeaderBytes || std::memcmp(bytes, kMagic, sizeof(kMagic)) != 0) {
        error = "invalid spatial feature sidecar magic or short header";
        return false;
    }
    header.schemaVersion = get16(bytes + 8);
    header.headerBytes = get16(bytes + 10);
    if (get32(bytes + 12) != kEndianMarker) {
        error = "invalid spatial feature sidecar byte-order marker";
        return false;
    }
    header.recordBytes = get16(bytes + 16);
    header.complete = get16(bytes + 18) == 1;
    header.totalReads = get64(bytes + 24);
    header.recordsWritten = get64(bytes + 32);
    header.featureCount = get32(bytes + 40);
    header.strand = static_cast<std::int32_t>(get32(bytes + 44));
    header.crMultimapRescue = bytes[48] != 0;
    header.crIntronicFallback = bytes[49] != 0;
    header.inputManifestSha256 = bytesToHex(bytes + 64, SHA256_DIGEST_LENGTH);
    header.geneDictionarySha256 = bytesToHex(bytes + 96, SHA256_DIGEST_LENGTH);
    header.starSuiteVersion = getText(bytes, 128, 64);
    header.sourceRevision = getText(bytes, 192, 384);
    header.featureType = getText(bytes, 576, 64);
    header.policy = getText(bytes, 640, 256);
    if (header.schemaVersion != kSchemaVersion || header.headerBytes != kHeaderBytes
        || header.recordBytes != kRecordBytes) {
        error = "unsupported spatial feature sidecar schema/header/record size";
        return false;
    }
    return true;
}

} // namespace

std::string normalizeReadName(const std::string &input)
{
    std::string value = input;
    const std::size_t blank = value.find_first_of(" \t\r\n");
    if (blank != std::string::npos) value.resize(blank);
    if (!value.empty() && value[0] == '@') value.erase(0, 1);
    if (value.size() > 2
        && (value.compare(value.size() - 2, 2, "/1") == 0
            || value.compare(value.size() - 2, 2, "/2") == 0)) {
        value.resize(value.size() - 2);
    }
    return value;
}

std::string sha256Hex(const std::string &value)
{
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char *>(value.data()), value.size(), digest);
    return hexDigest(digest);
}

std::string sha256File(const std::string &path, std::string &error)
{
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input) {
        error = "cannot open for SHA-256: " + path;
        return std::string();
    }
    SHA256_CTX context;
    SHA256_Init(&context);
    char buffer[1 << 20];
    while (input) {
        input.read(buffer, sizeof(buffer));
        const std::streamsize count = input.gcount();
        if (count > 0) SHA256_Update(&context, buffer, static_cast<std::size_t>(count));
    }
    if (!input.eof()) {
        error = "failed while hashing: " + path;
        return std::string();
    }
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256_Final(digest, &context);
    return hexDigest(digest);
}

std::vector<unsigned char> encodeRecord(const Record &record)
{
    std::vector<unsigned char> bytes(kRecordBytes, 0);
    put32(bytes.data(), record.geneIndex);
    put16(bytes.data() + 4, record.statusFlags);
    put16(bytes.data() + 6, record.reserved);
    return bytes;
}

bool decodeRecord(const unsigned char *bytes, std::size_t size, Record &record,
                  std::string &error)
{
    if (size != kRecordBytes) {
        error = "spatial feature record has wrong size";
        return false;
    }
    record.geneIndex = get32(bytes);
    record.statusFlags = get16(bytes + 4);
    record.reserved = get16(bytes + 6);
    return true;
}

bool validateRecord(const Record &record, std::uint32_t featureCount,
                    std::string &error)
{
    if (record.reserved != 0 || !(record.statusFlags & kRecordPresent)) {
        error = "record is missing or has nonzero reserved bits";
        return false;
    }
    const bool mapped = record.statusFlags & kMapped;
    const bool unmapped = record.statusFlags & kUnmappedOrFiltered;
    if (mapped == unmapped) {
        error = "record must be exactly one of mapped or unmapped/filtered";
        return false;
    }
    unsigned outcomes = 0;
    outcomes += (record.statusFlags & kUniqueGene) != 0;
    outcomes += (record.statusFlags & kNoGene) != 0;
    outcomes += (record.statusFlags & kMultiGeneRejected) != 0;
    if (outcomes != 1) {
        error = "record must have exactly one GeneFull outcome";
        return false;
    }
    if (record.statusFlags & kUniqueGene) {
        if (record.geneIndex >= featureCount) {
            error = "unique-gene record references an out-of-range feature";
            return false;
        }
    } else if (record.geneIndex != kMissingGene) {
        error = "non-unique record must carry UINT32_MAX gene index";
        return false;
    }
    if ((record.statusFlags & kCrExonicRescue)
        && (record.statusFlags & kCrIntronicFallback)) {
        error = "record cannot be both exonic and intronic-fallback rescued";
        return false;
    }
    return true;
}

Writer::Writer()
    : binaryFd_(-1), namesFd_(-1), recordsWritten_(0),
      firstOrdinal_(std::numeric_limits<std::uint64_t>::max()), lastOrdinal_(0),
      writeFailed_(false)
{}

Writer::~Writer()
{
    if (binaryFd_ >= 0) ::close(binaryFd_);
    if (namesFd_ >= 0) ::close(namesFd_);
}

bool Writer::fail(const std::string &message, std::string &error)
{
    writeFailed_.store(true);
    {
        std::lock_guard<std::mutex> lock(gWriteErrorMutex);
        if (writeError_.empty()) writeError_ = message;
    }
    error = message;
    return false;
}

bool Writer::open(const WriterConfig &config,
                  const std::vector<std::string> &geneIds,
                  const std::vector<std::string> &geneNames,
                  std::string &error)
{
    if (isOpen() || config.prefix.empty() || config.prefix == "-") {
        error = "invalid or already-open spatial feature sidecar prefix";
        return false;
    }
    if (geneIds.empty() || geneIds.size() != geneNames.size()
        || geneIds.size() > UINT32_MAX) {
        error = "invalid spatial feature gene dictionary";
        return false;
    }
    config_ = config;
    header_.schemaVersion = kSchemaVersion;
    header_.headerBytes = kHeaderBytes;
    header_.recordBytes = kRecordBytes;
    header_.featureCount = static_cast<std::uint32_t>(geneIds.size());
    header_.strand = config.strand;
    header_.crMultimapRescue = config.crMultimapRescue;
    header_.crIntronicFallback = config.crIntronicFallback;
    header_.starSuiteVersion = config.starSuiteVersion;
    header_.sourceRevision = config.sourceRevision;
    header_.featureType = config.featureType;
    header_.policy = config.policy;
    header_.inputManifestSha256 = sha256Hex(config.inputManifest);

    binaryTempPath_ = config.prefix + ".bin.tmp";
    namesTempPath_ = config.prefix + ".read_names.tmp";
    featuresPath_ = config.prefix + ".features.tsv";
    digestPath_ = config.prefix + ".read_name_digests.tsv";
    summaryPath_ = config.prefix + ".summary.json";
    for (const std::string &path : {config.prefix + ".bin", binaryTempPath_, namesTempPath_,
                                    featuresPath_, digestPath_, summaryPath_}) {
        if (::access(path.c_str(), F_OK) == 0) {
            error = "refusing to overwrite spatial feature sidecar path: " + path;
            return false;
        }
    }

    std::ostringstream dictionary;
    dictionary << "gene_index\tgene_id\tgene_name\n";
    for (std::size_t i = 0; i < geneIds.size(); ++i) {
        if (geneIds[i].find_first_of("\t\r\n") != std::string::npos
            || geneNames[i].find_first_of("\t\r\n") != std::string::npos) {
            error = "gene dictionary contains a tab or newline";
            return false;
        }
        dictionary << i << '\t' << geneIds[i] << '\t' << geneNames[i] << '\n';
    }
    const std::string dictionaryText = dictionary.str();
    header_.geneDictionarySha256 = sha256Hex(dictionaryText);
    if (!writeAtomicText(featuresPath_, dictionaryText, error)) return false;

    binaryFd_ = ::open(binaryTempPath_.c_str(), O_RDWR | O_CREAT | O_EXCL, 0644);
    if (binaryFd_ < 0) {
        error = "cannot create " + binaryTempPath_ + ": " + std::strerror(errno);
        return false;
    }
    namesFd_ = ::open(namesTempPath_.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
    if (namesFd_ < 0) {
        error = "cannot create " + namesTempPath_ + ": " + std::strerror(errno);
        return false;
    }
    return writeHeader(false, 0, 0, error);
}

bool Writer::writeHeader(bool complete, std::uint64_t totalReads,
                         std::uint64_t recordsWritten, std::string &error)
{
    std::vector<unsigned char> bytes(kHeaderBytes, 0);
    std::memcpy(bytes.data(), kMagic, sizeof(kMagic));
    put16(bytes.data() + 8, kSchemaVersion);
    put16(bytes.data() + 10, kHeaderBytes);
    put32(bytes.data() + 12, kEndianMarker);
    put16(bytes.data() + 16, kRecordBytes);
    put16(bytes.data() + 18, complete ? 1 : 0);
    put64(bytes.data() + 24, totalReads);
    put64(bytes.data() + 32, recordsWritten);
    put32(bytes.data() + 40, header_.featureCount);
    put32(bytes.data() + 44, static_cast<std::uint32_t>(header_.strand));
    bytes[48] = header_.crMultimapRescue ? 1 : 0;
    bytes[49] = header_.crIntronicFallback ? 1 : 0;
    if (!hexToBytes(header_.inputManifestSha256, bytes.data() + 64, SHA256_DIGEST_LENGTH)
        || !hexToBytes(header_.geneDictionarySha256, bytes.data() + 96, SHA256_DIGEST_LENGTH)) {
        error = "internal SHA-256 encoding failure";
        return false;
    }
    putText(bytes.data(), 128, 64, header_.starSuiteVersion);
    putText(bytes.data(), 192, 384, header_.sourceRevision);
    putText(bytes.data(), 576, 64, header_.featureType);
    putText(bytes.data(), 640, 256, header_.policy);
    return writeAt(binaryFd_, bytes.data(), bytes.size(), 0, error);
}

bool Writer::write(std::uint64_t ordinal, std::uint32_t laneIndex,
                   const std::string &normalizedReadName, const Record &record,
                   std::string &error)
{
    if (!isOpen() || writeFailed_.load()) {
        error = writeError_.empty() ? "spatial feature sidecar is not writable" : writeError_;
        return false;
    }
    if (normalizedReadName.empty()) return fail("empty normalized read name", error);
    if (!validateRecord(record, header_.featureCount, error)) return fail(error, error);
    const std::vector<unsigned char> bytes = encodeRecord(record);
    const std::uint64_t recordOffset = static_cast<std::uint64_t>(kHeaderBytes)
        + ordinal * static_cast<std::uint64_t>(kRecordBytes);
    if (recordOffset > static_cast<std::uint64_t>(std::numeric_limits<off_t>::max())) {
        return fail("spatial feature sidecar offset overflow", error);
    }
    if (!writeAt(binaryFd_, bytes.data(), bytes.size(), static_cast<off_t>(recordOffset), error)) {
        return fail(error, error);
    }

    unsigned char nameBytes[16] = {};
    put64(nameBytes, fnv1a64(normalizedReadName));
    put32(nameBytes + 8, laneIndex);
    put32(nameBytes + 12, kNamePresent);
    const std::uint64_t nameOffset = ordinal * sizeof(nameBytes);
    if (nameOffset > static_cast<std::uint64_t>(std::numeric_limits<off_t>::max())
        || !writeAt(namesFd_, nameBytes, sizeof(nameBytes), static_cast<off_t>(nameOffset), error)) {
        return fail(error.empty() ? "read-name scratch offset overflow" : error, error);
    }

    recordsWritten_.fetch_add(1);
    std::uint64_t first = firstOrdinal_.load();
    while (ordinal < first && !firstOrdinal_.compare_exchange_weak(first, ordinal)) {}
    std::uint64_t last = lastOrdinal_.load();
    while (ordinal > last && !lastOrdinal_.compare_exchange_weak(last, ordinal)) {}
    return true;
}

bool Writer::finalize(std::uint64_t totalReads, std::string &error)
{
    if (!isOpen() || namesFd_ < 0) {
        error = "spatial feature sidecar is not open";
        return false;
    }
    if (writeFailed_.load()) {
        error = writeError_.empty() ? "a spatial feature record write failed" : writeError_;
        return false;
    }
    const std::uint64_t written = recordsWritten_.load();
    if (totalReads == 0 || written != totalReads) {
        std::ostringstream message;
        message << "spatial feature record count mismatch: wrote " << written
                << ", expected " << totalReads;
        error = message.str();
        return false;
    }

    struct LaneBlock {
        std::uint32_t lane = 0;
        std::uint64_t first = 0;
        std::uint64_t count = 0;
        std::string digest;
        LaneBlock(std::uint32_t laneIn, std::uint64_t firstIn,
                  std::uint64_t countIn, const std::string &digestIn)
            : lane(laneIn), first(firstIn), count(countIn), digest(digestIn) {}
    };
    struct LaneTotal {
        std::uint32_t lane = 0;
        std::uint64_t first = 0;
        std::uint64_t count = 0;
        LaneTotal(std::uint32_t laneIn, std::uint64_t firstIn, std::uint64_t countIn)
            : lane(laneIn), first(firstIn), count(countIn) {}
    };
    std::vector<LaneBlock> blocks;
    std::vector<LaneTotal> lanes;
    std::uint32_t activeLane = UINT32_MAX;
    std::uint64_t activeBlockCount = 0;
    std::uint64_t activeBlockFirst = 0;
    SHA256_CTX blockContext;
    bool blockActive = false;
    auto finishBlock = [&]() {
        if (!blockActive) return;
        unsigned char digest[SHA256_DIGEST_LENGTH];
        SHA256_Final(digest, &blockContext);
        blocks.push_back(LaneBlock(activeLane, activeBlockFirst, activeBlockCount,
                                   hexDigest(digest)));
        blockActive = false;
        activeBlockCount = 0;
    };

    flagCounts_.clear();
    for (std::uint64_t ordinal = 0; ordinal < totalReads; ++ordinal) {
        unsigned char recordBytes[kRecordBytes];
        unsigned char nameBytes[16];
        if (!readAt(binaryFd_, recordBytes, sizeof(recordBytes),
                    static_cast<off_t>(kHeaderBytes + ordinal * kRecordBytes), error)
            || !readAt(namesFd_, nameBytes, sizeof(nameBytes),
                       static_cast<off_t>(ordinal * sizeof(nameBytes)), error)) {
            return false;
        }
        Record record;
        if (!decodeRecord(recordBytes, sizeof(recordBytes), record, error)
            || !validateRecord(record, header_.featureCount, error)) {
            error = "invalid record at ordinal " + std::to_string(ordinal) + ": " + error;
            return false;
        }
        if (get32(nameBytes + 12) != kNamePresent) {
            error = "missing read-name digest at ordinal " + std::to_string(ordinal);
            return false;
        }
        for (unsigned bit = 0; bit < 16; ++bit) {
            const std::uint16_t mask = static_cast<std::uint16_t>(1u << bit);
            if (record.statusFlags & mask) ++flagCounts_[mask];
        }
        const std::uint32_t lane = get32(nameBytes + 8);
        if (lane != activeLane) {
            finishBlock();
            if (!lanes.empty() && lane <= lanes.back().lane) {
                error = "input lane indices are not strictly increasing at ordinal "
                    + std::to_string(ordinal);
                return false;
            }
            activeLane = lane;
            lanes.push_back(LaneTotal(lane, ordinal, 0));
        }
        ++lanes.back().count;
        if (!blockActive) {
            SHA256_Init(&blockContext);
            activeBlockFirst = ordinal;
            activeBlockCount = 0;
            blockActive = true;
        }
        SHA256_Update(&blockContext, nameBytes, 8);
        ++activeBlockCount;
        if (activeBlockCount == kNameDigestBlockReads) finishBlock();
    }
    finishBlock();

    std::ostringstream digestTable;
    digestTable << "lane_index\tfirst_ordinal\tread_count\tsha256\n";
    for (const LaneBlock &block : blocks) {
        digestTable << block.lane << '\t' << block.first << '\t' << block.count
                    << '\t' << block.digest << '\n';
    }
    if (!writeAtomicText(digestPath_, digestTable.str(), error)) return false;

    if (!writeHeader(true, totalReads, written, error)) return false;
    const off_t finalSize = static_cast<off_t>(kHeaderBytes + totalReads * kRecordBytes);
    if (::ftruncate(binaryFd_, finalSize) != 0 || ::fsync(binaryFd_) != 0) {
        error = "could not flush finalized spatial feature sidecar: "
            + std::string(std::strerror(errno));
        return false;
    }
    if (::close(binaryFd_) != 0) {
        binaryFd_ = -1;
        error = "could not close finalized spatial feature sidecar";
        return false;
    }
    binaryFd_ = -1;
    ::close(namesFd_);
    namesFd_ = -1;
    ::unlink(namesTempPath_.c_str());
    const std::string binaryPath = config_.prefix + ".bin";
    if (::rename(binaryTempPath_.c_str(), binaryPath.c_str()) != 0) {
        error = "could not atomically commit spatial feature sidecar: "
            + std::string(std::strerror(errno));
        return false;
    }

    const std::string binarySha = sha256File(binaryPath, error);
    if (binarySha.empty()) return false;
    const std::string digestSha = sha256File(digestPath_, error);
    if (digestSha.empty()) return false;
    std::ostringstream summary;
    summary << "{\n"
            << "  \"schema\": \"star_suite.spatial_feature_sidecar.v1\",\n"
            << "  \"status\": \"complete\",\n"
            << "  \"total_reads\": " << totalReads << ",\n"
            << "  \"feature_count\": " << header_.featureCount << ",\n"
            << "  \"first_ordinal\": " << firstOrdinal_.load() << ",\n"
            << "  \"last_ordinal\": " << lastOrdinal_.load() << ",\n"
            << "  \"record_bytes\": " << kRecordBytes << ",\n"
            << "  \"name_digest_block_reads\": " << kNameDigestBlockReads << ",\n"
            << "  \"input_manifest_sha256\": \"" << header_.inputManifestSha256 << "\",\n"
            << "  \"gene_dictionary_sha256\": \"" << header_.geneDictionarySha256 << "\",\n"
            << "  \"binary_sha256\": \"" << binarySha << "\",\n"
            << "  \"read_name_digests_sha256\": \"" << digestSha << "\",\n"
            << "  \"star_suite_version\": \"" << jsonEscape(header_.starSuiteVersion) << "\",\n"
            << "  \"source_revision\": \"" << jsonEscape(header_.sourceRevision) << "\",\n"
            << "  \"feature_type\": \"" << jsonEscape(header_.featureType) << "\",\n"
            << "  \"policy\": \"" << jsonEscape(header_.policy) << "\",\n"
            << "  \"lane_boundaries\": [\n";
    for (std::size_t i = 0; i < lanes.size(); ++i) {
        summary << "    {\"lane_index\": " << lanes[i].lane
                << ", \"first_ordinal\": " << lanes[i].first
                << ", \"read_count\": " << lanes[i].count << "}"
                << (i + 1 == lanes.size() ? "\n" : ",\n");
    }
    summary << "  ],\n  \"status_flag_counts\": {";
    bool firstFlag = true;
    for (const auto &entry : flagCounts_) {
        if (!firstFlag) summary << ", ";
        summary << "\"" << entry.first << "\": " << entry.second;
        firstFlag = false;
    }
    summary << "}\n}\n";
    return writeAtomicText(summaryPath_, summary.str(), error);
}

Reader::Reader() : fd_(-1) {}

Reader::~Reader()
{
    if (fd_ >= 0) ::close(fd_);
}

bool Reader::open(const std::string &path, std::string &error)
{
    if (fd_ >= 0) {
        error = "spatial feature sidecar reader is already open";
        return false;
    }
    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        error = "cannot open spatial feature sidecar: " + path;
        return false;
    }
    unsigned char bytes[kHeaderBytes];
    if (!readAt(fd_, bytes, sizeof(bytes), 0, error)
        || !decodeHeader(bytes, sizeof(bytes), header_, error)) {
        return false;
    }
    if (!header_.complete || header_.totalReads == 0
        || header_.recordsWritten != header_.totalReads) {
        error = "spatial feature sidecar is incomplete";
        return false;
    }
    struct stat info;
    const std::uint64_t expected = static_cast<std::uint64_t>(header_.headerBytes)
        + header_.totalReads * header_.recordBytes;
    if (::fstat(fd_, &info) != 0 || static_cast<std::uint64_t>(info.st_size) != expected) {
        error = "spatial feature sidecar file size does not match header";
        return false;
    }
    return true;
}

bool Reader::read(std::uint64_t ordinal, Record &record, std::string &error) const
{
    if (fd_ < 0 || ordinal >= header_.totalReads) {
        error = "spatial feature ordinal is out of range";
        return false;
    }
    unsigned char bytes[kRecordBytes];
    if (!readAt(fd_, bytes, sizeof(bytes),
                static_cast<off_t>(header_.headerBytes + ordinal * header_.recordBytes), error)
        || !decodeRecord(bytes, sizeof(bytes), record, error)) {
        return false;
    }
    return validateRecord(record, header_.featureCount, error);
}

bool Reader::validateAll(std::string &error) const
{
    Record record;
    for (std::uint64_t ordinal = 0; ordinal < header_.totalReads; ++ordinal) {
        if (!read(ordinal, record, error)) {
            error = "invalid spatial feature record at ordinal " + std::to_string(ordinal)
                + ": " + error;
            return false;
        }
    }
    return true;
}

} // namespace spatial_feature_sidecar
