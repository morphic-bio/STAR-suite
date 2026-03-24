#ifndef H_SoloBinarySpool
#define H_SoloBinarySpool

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <limits>
#include <vector>
#include "IncludeDefine.h"
#include "SoloFeatureTypes.h"

namespace SoloBinarySpool {

static constexpr uint32_t kMagic = 0x4C4F5342u; // "BSOL"
static constexpr uint16_t kVersion = 1;
static constexpr uint8_t kFlagReadIndex = 1u;
static constexpr size_t kDefaultChunkBytes = 4u * 1024u * 1024u;

struct FileHeader {
    uint32_t magic;
    uint16_t version;
    int16_t featureType;
    uint8_t flags;
    uint8_t reserved[7];
};

struct MemoryChunk {
    std::unique_ptr<char[]> data;
    size_t used;

    MemoryChunk() : data(), used(0) {}
    explicit MemoryChunk(size_t capacity) : data(new char[capacity]), used(0) {}
};

struct MemoryBuffer {
    size_t chunkBytes;
    std::vector<MemoryChunk> chunks;
    size_t totalBytes;

    explicit MemoryBuffer(size_t chunkBytesIn = kDefaultChunkBytes)
        : chunkBytes(chunkBytesIn), chunks(), totalBytes(0) {}

    void clear() {
        chunks.clear();
        totalBytes = 0;
    }
};

inline MemoryChunk &ensureChunk(MemoryBuffer &buffer) {
    if (buffer.chunks.empty() || buffer.chunks.back().used == buffer.chunkBytes) {
        buffer.chunks.emplace_back(buffer.chunkBytes);
    }
    return buffer.chunks.back();
}

inline void writeBytes(MemoryBuffer &buffer, const void *srcVoid, size_t nBytes) {
    const char *src = reinterpret_cast<const char *>(srcVoid);
    size_t remaining = nBytes;
    while (remaining > 0) {
        MemoryChunk &chunk = ensureChunk(buffer);
        const size_t avail = buffer.chunkBytes - chunk.used;
        const size_t step = remaining < avail ? remaining : avail;
        std::memcpy(chunk.data.get() + chunk.used, src, step);
        chunk.used += step;
        buffer.totalBytes += step;
        src += step;
        remaining -= step;
    }
}

inline bool readBytes(const MemoryBuffer &buffer, size_t &offset, void *dstVoid, size_t nBytes) {
    if (offset + nBytes > buffer.totalBytes) {
        return false;
    }

    char *dst = reinterpret_cast<char *>(dstVoid);
    size_t remaining = nBytes;
    size_t localOffset = offset;
    while (remaining > 0) {
        const size_t chunkIndex = localOffset / buffer.chunkBytes;
        const size_t chunkOffset = localOffset % buffer.chunkBytes;
        const MemoryChunk &chunk = buffer.chunks[chunkIndex];
        const size_t avail = chunk.used - chunkOffset;
        const size_t step = remaining < avail ? remaining : avail;
        std::memcpy(dst, chunk.data.get() + chunkOffset, step);
        dst += step;
        localOffset += step;
        remaining -= step;
    }
    offset = localOffset;
    return true;
}

template <typename T>
inline void writePod(std::fstream &stream, const T &value) {
    stream.write(reinterpret_cast<const char *>(&value), sizeof(T));
}

template <typename T>
inline void writePod(MemoryBuffer &buffer, const T &value) {
    writeBytes(buffer, &value, sizeof(T));
}

template <typename T>
inline bool readPod(std::fstream &stream, T &value) {
    stream.read(reinterpret_cast<char *>(&value), sizeof(T));
    return stream.good();
}

template <typename T>
inline bool readPod(const MemoryBuffer &buffer, size_t &offset, T &value) {
    return readBytes(buffer, offset, &value, sizeof(T));
}

inline bool envEnabled() {
    const char *value = std::getenv("STAR_SOLO_BINARY_SPOOL");
    if (value == nullptr || value[0] == '\0') {
        value = std::getenv("STAR_SOLO_BINARY_SPOOL_IN_MEMORY");
    }
    if (value == nullptr || value[0] == '\0') {
        return false;
    }
    if (value[0] == '0' && value[1] == '\0') {
        return false;
    }
    if ((value[0] == 'n' || value[0] == 'N') && (value[1] == 'o' || value[1] == 'O') && value[2] == '\0') {
        return false;
    }
    return true;
}

inline bool envInMemoryEnabled() {
    const char *value = std::getenv("STAR_SOLO_BINARY_SPOOL_IN_MEMORY");
    if (value == nullptr || value[0] == '\0') {
        return false;
    }
    if (value[0] == '0' && value[1] == '\0') {
        return false;
    }
    if ((value[0] == 'n' || value[0] == 'N') && (value[1] == 'o' || value[1] == 'O') && value[2] == '\0') {
        return false;
    }
    return true;
}

inline uint64_t envInMemoryLimitBytes() {
    const char *value = std::getenv("STAR_SOLO_BINARY_SPOOL_IN_MEMORY_LIMIT_MB");
    if (value == nullptr || value[0] == '\0') {
        return 0;
    }
    char *end = nullptr;
    const unsigned long long limitMb = std::strtoull(value, &end, 10);
    if (end == value || *end != '\0') {
        return 0;
    }
    const unsigned long long bytes = limitMb * 1024ull * 1024ull;
    if (bytes > std::numeric_limits<uint64_t>::max()) {
        return std::numeric_limits<uint64_t>::max();
    }
    return static_cast<uint64_t>(bytes);
}

inline bool supportsFeature(const int32 featureType) {
    switch (featureType) {
        case SoloFeatureTypes::Gene:
        case SoloFeatureTypes::GeneFull:
        case SoloFeatureTypes::GeneFull_Ex50pAS:
        case SoloFeatureTypes::GeneFull_ExonOverIntron:
            return true;
        default:
            return false;
    }
}

inline void writeFileHeader(std::fstream &stream, const int32 featureType, const bool readIndexYes) {
    FileHeader header{};
    header.magic = kMagic;
    header.version = kVersion;
    header.featureType = static_cast<int16_t>(featureType);
    header.flags = readIndexYes ? kFlagReadIndex : 0u;
    stream.seekp(0, std::ios::beg);
    writePod(stream, header);
}

inline void writeFileHeader(MemoryBuffer &buffer, const int32 featureType, const bool readIndexYes) {
    buffer.clear();
    FileHeader header{};
    header.magic = kMagic;
    header.version = kVersion;
    header.featureType = static_cast<int16_t>(featureType);
    header.flags = readIndexYes ? kFlagReadIndex : 0u;
    writePod(buffer, header);
}

inline void flushToStream(const MemoryBuffer &buffer, std::fstream &stream) {
    for (const auto &chunk : buffer.chunks) {
        if (chunk.used > 0) {
            stream.write(chunk.data.get(), static_cast<std::streamsize>(chunk.used));
        }
    }
}

inline bool seekToRecords(std::fstream &stream, const int32 expectedFeatureType, const bool expectedReadIndexYes) {
    stream.clear();
    stream.seekg(0, std::ios::beg);
    FileHeader header{};
    if (!readPod(stream, header)) {
        return false;
    }
    return header.magic == kMagic &&
           header.version == kVersion &&
           header.featureType == static_cast<int16_t>(expectedFeatureType) &&
           ((header.flags & kFlagReadIndex) != 0u) == expectedReadIndexYes;
}

inline bool seekToRecords(const MemoryBuffer &buffer, size_t &offset, const int32 expectedFeatureType, const bool expectedReadIndexYes) {
    offset = 0;
    FileHeader header{};
    if (!readPod(buffer, offset, header)) {
        return false;
    }
    return header.magic == kMagic &&
           header.version == kVersion &&
           header.featureType == static_cast<int16_t>(expectedFeatureType) &&
           ((header.flags & kFlagReadIndex) != 0u) == expectedReadIndexYes;
}

inline void writeRecordHeader(std::fstream &stream, const bool readIndexYes, const uint64 umi,
                              const uint64 iRead, const uint32 readFlag, const uint32 feature,
                              const int32 cbmatch) {
    writePod(stream, umi);
    if (readIndexYes) {
        writePod(stream, iRead);
        writePod(stream, readFlag);
    }
    writePod(stream, feature);
    writePod(stream, cbmatch);
}

inline void writeRecordHeader(MemoryBuffer &buffer, const bool readIndexYes, const uint64 umi,
                              const uint64 iRead, const uint32 readFlag, const uint32 feature,
                              const int32 cbmatch) {
    writePod(buffer, umi);
    if (readIndexYes) {
        writePod(buffer, iRead);
        writePod(buffer, readFlag);
    }
    writePod(buffer, feature);
    writePod(buffer, cbmatch);
}

inline bool readRecordHeader(std::fstream &stream, const bool readIndexYes, uint64 &umi,
                             uint64 &iRead, uint32 &readFlag, uint32 &feature, int32 &cbmatch) {
    if (!readPod(stream, umi)) {
        return false;
    }
    if (readIndexYes) {
        if (!readPod(stream, iRead) || !readPod(stream, readFlag)) {
            return false;
        }
    } else {
        iRead = static_cast<uint64_t>(-1);
        readFlag = 0;
    }
    return readPod(stream, feature) && readPod(stream, cbmatch);
}

inline bool readRecordHeader(const MemoryBuffer &buffer, size_t &offset, const bool readIndexYes, uint64 &umi,
                             uint64 &iRead, uint32 &readFlag, uint32 &feature, int32 &cbmatch) {
    if (!readPod(buffer, offset, umi)) {
        return false;
    }
    if (readIndexYes) {
        if (!readPod(buffer, offset, iRead) || !readPod(buffer, offset, readFlag)) {
            return false;
        }
    } else {
        iRead = static_cast<uint64_t>(-1);
        readFlag = 0;
    }
    return readPod(buffer, offset, feature) && readPod(buffer, offset, cbmatch);
}

inline void writeResolvedCb(std::fstream &stream, const uint64 cb) {
    writePod(stream, cb);
}

inline void writeResolvedCb(MemoryBuffer &buffer, const uint64 cb) {
    writePod(buffer, cb);
}

inline bool readResolvedCb(std::fstream &stream, int64 &cb) {
    uint64_t raw = 0;
    if (!readPod(stream, raw)) {
        return false;
    }
    cb = static_cast<int64_t>(raw);
    return true;
}

inline bool readResolvedCb(const MemoryBuffer &buffer, size_t &offset, int64 &cb) {
    uint64_t raw = 0;
    if (!readPod(buffer, offset, raw)) {
        return false;
    }
    cb = static_cast<int64_t>(raw);
    return true;
}

inline void writeAmbiguousCandidate(std::fstream &stream, const uint32_t cbIdx, const uint8_t qual) {
    writePod(stream, cbIdx);
    writePod(stream, qual);
}

inline void writeAmbiguousCandidate(MemoryBuffer &buffer, const uint32_t cbIdx, const uint8_t qual) {
    writePod(buffer, cbIdx);
    writePod(buffer, qual);
}

inline bool readAmbiguousCandidate(std::fstream &stream, uint32_t &cbIdx, char &qual) {
    uint8_t qualByte = 0;
    if (!readPod(stream, cbIdx) || !readPod(stream, qualByte)) {
        return false;
    }
    qual = static_cast<char>(qualByte);
    return true;
}

inline bool readAmbiguousCandidate(const MemoryBuffer &buffer, size_t &offset, uint32_t &cbIdx, char &qual) {
    uint8_t qualByte = 0;
    if (!readPod(buffer, offset, cbIdx) || !readPod(buffer, offset, qualByte)) {
        return false;
    }
    qual = static_cast<char>(qualByte);
    return true;
}

inline bool skipCbPayload(std::fstream &stream, const int32_t cbmatch) {
    if (cbmatch <= 1) {
        uint64_t cb = 0;
        return readPod(stream, cb);
    }
    for (int32_t ii = 0; ii < cbmatch; ++ii) {
        uint32_t cbIdx = 0;
        uint8_t qual = 0;
        if (!readPod(stream, cbIdx) || !readPod(stream, qual)) {
            return false;
        }
    }
    return true;
}

inline bool skipCbPayload(const MemoryBuffer &buffer, size_t &offset, const int32_t cbmatch) {
    if (cbmatch <= 1) {
        uint64_t cb = 0;
        return readPod(buffer, offset, cb);
    }
    for (int32_t ii = 0; ii < cbmatch; ++ii) {
        uint32_t cbIdx = 0;
        uint8_t qual = 0;
        if (!readPod(buffer, offset, cbIdx) || !readPod(buffer, offset, qual)) {
            return false;
        }
    }
    return true;
}

} // namespace SoloBinarySpool

#endif
