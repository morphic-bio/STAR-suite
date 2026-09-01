#include "input/BgzfIndex.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <fcntl.h>
#include <limits>
#include <mutex>
#include <sstream>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <zlib.h>

namespace star {
namespace input {
namespace {

const unsigned char kBgzfEof[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};
const char kSidecarMagic[8] = {'S', 'T', 'A', 'R', 'B', 'G', 'Z', 'I'};
const uint32_t kSidecarVersion = 1;
const uint32_t kNoRecordOffset = std::numeric_limits<uint32_t>::max();

bool set_error(std::string* error, const std::string& message) {
    if (error != nullptr) {
        *error = message;
    }
    return false;
}

uint16_t load_u16(const unsigned char* data) {
    return static_cast<uint16_t>(data[0]) |
           static_cast<uint16_t>(static_cast<uint16_t>(data[1]) << 8);
}

uint32_t load_u32(const unsigned char* data) {
    return static_cast<uint32_t>(data[0]) |
           (static_cast<uint32_t>(data[1]) << 8) |
           (static_cast<uint32_t>(data[2]) << 16) |
           (static_cast<uint32_t>(data[3]) << 24);
}

bool parse_bgzf_header18(const unsigned char* header, uint32_t* compressed_size) {
    if (header[0] != 0x1f || header[1] != 0x8b || header[2] != 8 ||
        (header[3] & 0x04U) == 0) {
        return false;
    }
    const uint16_t xlen = load_u16(header + 10);
    const size_t available_extra = std::min<size_t>(xlen, 6);
    size_t pos = 12;
    const size_t end = 12 + available_extra;
    while (pos + 4 <= end) {
        const uint16_t slen = load_u16(header + pos + 2);
        pos += 4;
        if (pos + slen > end) {
            return false;
        }
        if (header[pos - 4] == 'B' && header[pos - 3] == 'C' && slen == 2) {
            if (compressed_size != nullptr) {
                *compressed_size = static_cast<uint32_t>(load_u16(header + pos)) + 1U;
            }
            return true;
        }
        pos += slen;
    }
    return false;
}

bool read_exact(std::ifstream& input,
                uint64_t offset,
                void* destination,
                size_t size) {
    if (offset > static_cast<uint64_t>(std::numeric_limits<std::streamoff>::max())) {
        return false;
    }
    input.clear();
    input.seekg(static_cast<std::streamoff>(offset), std::ios::beg);
    if (!input.good()) {
        return false;
    }
    input.read(static_cast<char*>(destination), static_cast<std::streamsize>(size));
    return input.good() || static_cast<size_t>(input.gcount()) == size;
}

bool read_exact_fd(int input_fd,
                   uint64_t offset,
                   void* destination,
                   size_t size,
                   std::string* error) {
    unsigned char* output = static_cast<unsigned char*>(destination);
    size_t completed = 0;
    while (completed < size) {
        if (offset + completed > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
            return set_error(error, "BGZF member offset exceeds platform file limits");
        }
        const ssize_t count = ::pread(input_fd,
                                      output + completed,
                                      size - completed,
                                      static_cast<off_t>(offset + completed));
        if (count < 0) {
            if (errno == EINTR) {
                continue;
            }
            return set_error(error, std::string("could not read BGZF member: ") +
                                    std::strerror(errno));
        }
        if (count == 0) {
            return set_error(error, "unexpected end of BGZF member");
        }
        completed += static_cast<size_t>(count);
    }
    return true;
}

struct FileIdentity {
    uint64_t size = 0;
    int64_t mtimeSeconds = 0;
    int64_t mtimeNanoseconds = 0;
};

bool file_identity(const std::string& path,
                   FileIdentity* identity,
                   std::string* error) {
    struct stat info;
    if (::stat(path.c_str(), &info) != 0) {
        return set_error(error, "could not stat BGZF input " + path + ": " + std::strerror(errno));
    }
    if (info.st_size < 0) {
        return set_error(error, "BGZF input has a negative file size: " + path);
    }
    identity->size = static_cast<uint64_t>(info.st_size);
    identity->mtimeSeconds = static_cast<int64_t>(info.st_mtime);
#if defined(__APPLE__)
    identity->mtimeNanoseconds = static_cast<int64_t>(info.st_mtimespec.tv_nsec);
#else
    identity->mtimeNanoseconds = static_cast<int64_t>(info.st_mtim.tv_nsec);
#endif
    return true;
}

bool detect_with_identity(const std::string& path,
                          const FileIdentity& identity,
                          BgzfDetection* detection,
                          std::string* error) {
    detection->isBgzf = false;
    detection->hasEofMarker = false;
    if (identity.size == 0) {
        return true;
    }
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input.good()) {
        return set_error(error, "could not open BGZF candidate " + path);
    }
    if (identity.size < 18) {
        return true;
    }
    unsigned char header[18];
    if (!read_exact(input, 0, header, sizeof(header))) {
        return set_error(error, "could not read gzip header from " + path);
    }
    detection->isBgzf = parse_bgzf_header18(header, nullptr);
    if (!detection->isBgzf || identity.size < sizeof(kBgzfEof)) {
        return true;
    }
    unsigned char eof[sizeof(kBgzfEof)];
    if (!read_exact(input, identity.size - sizeof(kBgzfEof), eof, sizeof(eof))) {
        return set_error(error, "could not read BGZF EOF marker from " + path);
    }
    detection->hasEofMarker = std::memcmp(eof, kBgzfEof, sizeof(eof)) == 0;
    return true;
}

bool scan_blocks(const std::string& path,
                 uint64_t file_size,
                 std::vector<BgzfBlock>* blocks,
                 std::string* error) {
    blocks->clear();
    std::ifstream input(path.c_str(), std::ios::binary);
    if (!input.good()) {
        return set_error(error, "could not open BGZF input " + path);
    }
    uint64_t offset = 0;
    while (offset < file_size) {
        if (file_size - offset < 18) {
            std::ostringstream message;
            message << "truncated BGZF header at block offset " << offset;
            return set_error(error, message.str());
        }
        unsigned char header[18];
        if (!read_exact(input, offset, header, sizeof(header))) {
            std::ostringstream message;
            message << "could not read BGZF header at block offset " << offset;
            return set_error(error, message.str());
        }
        uint32_t compressed_size = 0;
        if (!parse_bgzf_header18(header, &compressed_size)) {
            std::ostringstream message;
            message << "invalid BGZF header at block offset " << offset;
            return set_error(error, message.str());
        }
        if (compressed_size < 26 || compressed_size > 65536) {
            std::ostringstream message;
            message << "invalid BGZF member size " << compressed_size
                    << " at block offset " << offset;
            return set_error(error, message.str());
        }
        if (compressed_size > file_size - offset) {
            std::ostringstream message;
            message << "truncated BGZF member at block offset " << offset
                    << " (declared size " << compressed_size << ")";
            return set_error(error, message.str());
        }
        if (compressed_size == sizeof(kBgzfEof)) {
            unsigned char candidate[sizeof(kBgzfEof)];
            if (!read_exact(input, offset, candidate, sizeof(candidate))) {
                return set_error(error, "could not read BGZF EOF member");
            }
            if (std::memcmp(candidate, kBgzfEof, sizeof(candidate)) == 0) {
                if (offset + compressed_size != file_size) {
                    std::ostringstream message;
                    message << "BGZF EOF marker is not final at block offset " << offset;
                    return set_error(error, message.str());
                }
                return true;
            }
        }
        unsigned char trailer_isize[4];
        if (!read_exact(input, offset + compressed_size - 4,
                        trailer_isize, sizeof(trailer_isize))) {
            std::ostringstream message;
            message << "could not read BGZF trailer at block offset " << offset;
            return set_error(error, message.str());
        }
        BgzfBlock block;
        block.compressedOffset = offset;
        block.compressedSize = compressed_size;
        block.isize = load_u32(trailer_isize);
        blocks->push_back(block);
        offset += compressed_size;
    }
    return set_error(error, "detected BGZF input is missing the required EOF marker");
}

struct CountResult {
    uint64_t newlineCount = 0;
    std::array<uint32_t, 4> firstNewlineOffsets;
    uint32_t firstNewlineCount = 0;
    bool nonempty = false;
    bool endsWithNewline = false;
};

bool count_block(int input_fd,
                 const BgzfBlock& block,
                 CountResult* result,
                 std::string* error) {
    std::vector<unsigned char> inflated;
    if (!inflate_bgzf_block_fd(input_fd, block, true, &inflated, error)) {
        return false;
    }
    result->nonempty = !inflated.empty();
    result->endsWithNewline = !inflated.empty() && inflated.back() == '\n';
    for (size_t offset = 0; offset < inflated.size(); ++offset) {
        if (inflated[offset] != '\n') {
            continue;
        }
        ++result->newlineCount;
        if (result->firstNewlineCount < result->firstNewlineOffsets.size()) {
            result->firstNewlineOffsets[result->firstNewlineCount++] =
                static_cast<uint32_t>(offset);
        }
    }
    return true;
}

bool count_records(const std::string& path,
                   uint32_t reader_threads,
                   std::vector<BgzfBlock>* blocks,
                   uint64_t* record_count,
                   std::string* error) {
    std::vector<CountResult> counts(blocks->size());
    if (!blocks->empty()) {
        const uint32_t nthreads = std::max<uint32_t>(
            1, std::min<uint32_t>(reader_threads == 0 ? 1 : reader_threads,
                                  static_cast<uint32_t>(blocks->size())));
        std::atomic<size_t> next(0);
        std::atomic<bool> failed(false);
        std::mutex error_mutex;
        std::string worker_error;
        std::vector<std::thread> workers;
        workers.reserve(nthreads);
        for (uint32_t thread = 0; thread < nthreads; ++thread) {
            workers.emplace_back([&]() {
                const int input_fd = ::open(path.c_str(), O_RDONLY);
                if (input_fd < 0) {
                    failed.store(true);
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (worker_error.empty()) {
                        worker_error = "could not open BGZF input " + path + ": " +
                                       std::strerror(errno);
                    }
                    return;
                }
                while (!failed.load()) {
                    const size_t index = next.fetch_add(1);
                    if (index >= blocks->size()) {
                        break;
                    }
                    std::string local_error;
                    if (!count_block(input_fd, (*blocks)[index], &counts[index], &local_error)) {
                        failed.store(true);
                        std::lock_guard<std::mutex> lock(error_mutex);
                        if (worker_error.empty()) {
                            worker_error = local_error;
                        }
                        break;
                    }
                }
                ::close(input_fd);
            });
        }
        for (size_t thread = 0; thread < workers.size(); ++thread) {
            workers[thread].join();
        }
        if (failed.load()) {
            return set_error(error, worker_error.empty()
                ? "BGZF record-count worker failed" : worker_error);
        }
    }

    uint64_t cumulative_newlines = 0;
    uint64_t cumulative_records = 0;
    bool at_line_start = true;
    bool saw_payload = false;
    for (size_t index = 0; index < blocks->size(); ++index) {
        BgzfBlock& block = (*blocks)[index];
        const CountResult& count = counts[index];
        block.firstRecordOrdinal = cumulative_records;
        block.recordsStarting = 0;
        block.firstRecordOffset = 0;
        block.hasFirstRecordOffset = false;
        if (!count.nonempty) {
            continue;
        }
        saw_payload = true;
        const bool starts_at_zero = at_line_start && cumulative_newlines % 4 == 0;
        if (starts_at_zero) {
            block.firstRecordOffset = 0;
            block.hasFirstRecordOffset = true;
            ++block.recordsStarting;
        }

        const uint64_t phase = cumulative_newlines % 4;
        const uint64_t first_matching_newline = phase == 0 ? 4 : 4 - phase;
        uint64_t matching_newlines = 0;
        if (count.newlineCount >= first_matching_newline) {
            matching_newlines = 1 + (count.newlineCount - first_matching_newline) / 4;
        }
        const bool final_newline_matches = count.endsWithNewline &&
            count.newlineCount >= first_matching_newline &&
            (count.newlineCount - first_matching_newline) % 4 == 0;
        if (final_newline_matches) {
            --matching_newlines;
        }
        if (matching_newlines != 0) {
            if (!block.hasFirstRecordOffset) {
                if (first_matching_newline == 0 ||
                    first_matching_newline > count.firstNewlineCount) {
                    return set_error(error, "internal BGZF record-offset accounting error");
                }
                block.firstRecordOffset =
                    count.firstNewlineOffsets[static_cast<size_t>(first_matching_newline - 1)] + 1;
                block.hasFirstRecordOffset = true;
            }
            if (matching_newlines > std::numeric_limits<uint32_t>::max() - block.recordsStarting) {
                return set_error(error, "BGZF block record count exceeds uint32_t");
            }
            block.recordsStarting += static_cast<uint32_t>(matching_newlines);
        }
        cumulative_records += block.recordsStarting;
        cumulative_newlines += count.newlineCount;
        at_line_start = count.endsWithNewline;
    }

    const uint64_t logical_lines = cumulative_newlines +
        ((saw_payload && !at_line_start) ? 1U : 0U);
    if (logical_lines % 4 != 0) {
        std::ostringstream message;
        message << "BGZF FASTQ has " << logical_lines
                << " logical lines; expected a multiple of four";
        return set_error(error, message.str());
    }
    *record_count = cumulative_records;
    return true;
}

void write_u8(std::ostream& output, uint8_t value) {
    output.put(static_cast<char>(value));
}

void write_u32(std::ostream& output, uint32_t value) {
    unsigned char bytes[4];
    for (unsigned shift = 0; shift < 4; ++shift) {
        bytes[shift] = static_cast<unsigned char>((value >> (shift * 8)) & 0xffU);
    }
    output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

void write_u64(std::ostream& output, uint64_t value) {
    unsigned char bytes[8];
    for (unsigned shift = 0; shift < 8; ++shift) {
        bytes[shift] = static_cast<unsigned char>((value >> (shift * 8)) & 0xffU);
    }
    output.write(reinterpret_cast<const char*>(bytes), sizeof(bytes));
}

bool read_u8(std::istream& input, uint8_t* value) {
    const int byte = input.get();
    if (byte == std::char_traits<char>::eof()) {
        return false;
    }
    *value = static_cast<uint8_t>(byte);
    return true;
}

bool read_u32(std::istream& input, uint32_t* value) {
    unsigned char bytes[4];
    input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
    if (input.gcount() != static_cast<std::streamsize>(sizeof(bytes))) {
        return false;
    }
    *value = load_u32(bytes);
    return true;
}

bool read_u64(std::istream& input, uint64_t* value) {
    unsigned char bytes[8];
    input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
    if (input.gcount() != static_cast<std::streamsize>(sizeof(bytes))) {
        return false;
    }
    uint64_t result = 0;
    for (unsigned shift = 0; shift < 8; ++shift) {
        result |= static_cast<uint64_t>(bytes[shift]) << (shift * 8);
    }
    *value = result;
    return true;
}

bool validate_cached_blocks(const std::vector<BgzfBlock>& blocks,
                            uint64_t file_size,
                            uint64_t record_count) {
    uint64_t offset = 0;
    uint64_t records = 0;
    for (size_t index = 0; index < blocks.size(); ++index) {
        const BgzfBlock& block = blocks[index];
        if (offset > file_size || block.compressedOffset != offset || block.compressedSize < 26 ||
            block.compressedSize > 65536 ||
            block.compressedSize > file_size - offset ||
            block.firstRecordOrdinal != records ||
            (block.recordsStarting != 0 && !block.hasFirstRecordOffset) ||
            (block.hasFirstRecordOffset && block.firstRecordOffset >= block.isize)) {
            return false;
        }
        offset += block.compressedSize;
        records += block.recordsStarting;
    }
    return records == record_count && offset + sizeof(kBgzfEof) == file_size;
}

bool load_sidecar(const std::string& sidecar,
                  const FileIdentity& identity,
                  std::vector<BgzfBlock>* blocks,
                  uint64_t* record_count) {
    std::ifstream input(sidecar.c_str(), std::ios::binary);
    if (!input.good()) {
        return false;
    }
    char magic[sizeof(kSidecarMagic)];
    input.read(magic, sizeof(magic));
    if (input.gcount() != static_cast<std::streamsize>(sizeof(magic)) ||
        std::memcmp(magic, kSidecarMagic, sizeof(magic)) != 0) {
        return false;
    }
    uint32_t version = 0;
    uint64_t file_size = 0;
    uint64_t mtime_seconds = 0;
    uint64_t mtime_nanoseconds = 0;
    uint64_t block_count = 0;
    if (!read_u32(input, &version) || !read_u64(input, &file_size) ||
        !read_u64(input, &mtime_seconds) || !read_u64(input, &mtime_nanoseconds) ||
        !read_u64(input, record_count) || !read_u64(input, &block_count)) {
        return false;
    }
    if (version != kSidecarVersion || file_size != identity.size ||
        static_cast<int64_t>(mtime_seconds) != identity.mtimeSeconds ||
        static_cast<int64_t>(mtime_nanoseconds) != identity.mtimeNanoseconds ||
        block_count > identity.size / 26U + 1U ||
        block_count > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        return false;
    }
    std::vector<BgzfBlock> loaded(static_cast<size_t>(block_count));
    for (size_t index = 0; index < loaded.size(); ++index) {
        BgzfBlock& block = loaded[index];
        uint8_t has_offset = 0;
        if (!read_u64(input, &block.compressedOffset) ||
            !read_u32(input, &block.compressedSize) ||
            !read_u32(input, &block.isize) ||
            !read_u64(input, &block.firstRecordOrdinal) ||
            !read_u32(input, &block.recordsStarting) ||
            !read_u32(input, &block.firstRecordOffset) ||
            !read_u8(input, &has_offset)) {
            return false;
        }
        block.hasFirstRecordOffset = has_offset == 1;
        if (has_offset > 1) {
            return false;
        }
    }
    if (input.peek() != std::char_traits<char>::eof() ||
        !validate_cached_blocks(loaded, identity.size, *record_count)) {
        return false;
    }
    blocks->swap(loaded);
    return true;
}

void write_sidecar(const std::string& sidecar,
                   const FileIdentity& identity,
                   const std::vector<BgzfBlock>& blocks,
                   uint64_t record_count) {
    std::ostringstream temporary_name;
    temporary_name << sidecar << ".tmp." << static_cast<long long>(::getpid());
    const std::string temporary = temporary_name.str();
    std::ofstream output(temporary.c_str(), std::ios::binary | std::ios::trunc);
    if (!output.good()) {
        return;
    }
    output.write(kSidecarMagic, sizeof(kSidecarMagic));
    write_u32(output, kSidecarVersion);
    write_u64(output, identity.size);
    write_u64(output, static_cast<uint64_t>(identity.mtimeSeconds));
    write_u64(output, static_cast<uint64_t>(identity.mtimeNanoseconds));
    write_u64(output, record_count);
    write_u64(output, blocks.size());
    for (size_t index = 0; index < blocks.size(); ++index) {
        const BgzfBlock& block = blocks[index];
        write_u64(output, block.compressedOffset);
        write_u32(output, block.compressedSize);
        write_u32(output, block.isize);
        write_u64(output, block.firstRecordOrdinal);
        write_u32(output, block.recordsStarting);
        write_u32(output, block.hasFirstRecordOffset
            ? block.firstRecordOffset : kNoRecordOffset);
        write_u8(output, block.hasFirstRecordOffset ? 1 : 0);
    }
    output.close();
    if (!output.good() || std::rename(temporary.c_str(), sidecar.c_str()) != 0) {
        std::remove(temporary.c_str());
    }
}

} // namespace

const char* bgzf_cache_status_name(BgzfCacheStatus status) {
    return status == BgzfCacheStatus::Loaded ? "loaded" : "rebuilt";
}

bool inflate_bgzf_block(const std::string& path,
                        const BgzfBlock& block,
                        bool check_crc,
                        std::vector<unsigned char>* output,
                        std::string* error) {
    if (output == nullptr) {
        return set_error(error, "BGZF inflate output is null");
    }
    const int input_fd = ::open(path.c_str(), O_RDONLY);
    if (input_fd < 0) {
        return set_error(error, "could not open BGZF input " + path + ": " +
                                std::strerror(errno));
    }
    const bool result = inflate_bgzf_block_fd(input_fd, block, check_crc, output, error);
    ::close(input_fd);
    return result;
}

bool inflate_bgzf_block_fd(int input_fd,
                           const BgzfBlock& block,
                           bool check_crc,
                           std::vector<unsigned char>* output,
                           std::string* error) {
    if (output == nullptr) {
        return set_error(error, "BGZF inflate output is null");
    }
    if (input_fd < 0) {
        return set_error(error, "BGZF inflate input descriptor is invalid");
    }
    std::vector<unsigned char> member(block.compressedSize);
    std::string read_error;
    if (!read_exact_fd(input_fd, block.compressedOffset, member.data(), member.size(),
                       &read_error)) {
        std::ostringstream message;
        message << "could not read BGZF member at block offset " << block.compressedOffset
                << ": " << read_error;
        return set_error(error, message.str());
    }
    if (member.size() < 26) {
        std::ostringstream message;
        message << "invalid BGZF member at block offset " << block.compressedOffset;
        return set_error(error, message.str());
    }
    const size_t deflate_size = member.size() - 18 - 8;
    const uint32_t expected_crc = load_u32(member.data() + member.size() - 8);
    const uint32_t expected_isize = load_u32(member.data() + member.size() - 4);
    if (expected_isize != block.isize) {
        std::ostringstream message;
        message << "BGZF ISIZE changed at block offset " << block.compressedOffset;
        return set_error(error, message.str());
    }
    std::vector<unsigned char> inflated(std::max<size_t>(block.isize, 1));
    z_stream stream;
    std::memset(&stream, 0, sizeof(stream));
    stream.next_in = member.data() + 18;
    stream.avail_in = static_cast<uInt>(deflate_size);
    stream.next_out = inflated.data();
    stream.avail_out = static_cast<uInt>(inflated.size());
    const int init_status = inflateInit2(&stream, -15);
    if (init_status != Z_OK) {
        return set_error(error, "zlib could not initialize raw BGZF inflate");
    }
    const int inflate_status = ::inflate(&stream, Z_FINISH);
    inflateEnd(&stream);
    if (inflate_status != Z_STREAM_END || stream.total_out != block.isize ||
        stream.total_in != deflate_size) {
        std::ostringstream message;
        message << "raw inflate failed at block offset " << block.compressedOffset
                << " (zlib status " << inflate_status << ")";
        return set_error(error, message.str());
    }
    inflated.resize(block.isize);
    if (check_crc) {
        const Bytef* data = inflated.empty()
            ? reinterpret_cast<const Bytef*>("") : inflated.data();
        const uint32_t observed_crc = static_cast<uint32_t>(
            crc32(0L, data, static_cast<uInt>(inflated.size())));
        if (observed_crc != expected_crc) {
            std::ostringstream message;
            message << "BGZF CRC mismatch at block offset " << block.compressedOffset;
            return set_error(error, message.str());
        }
    }
    output->swap(inflated);
    return true;
}

bool BgzfIndex::detect(const std::string& path,
                       BgzfDetection* detection,
                       std::string* error) {
    if (detection == nullptr) {
        return set_error(error, "BGZF detection output is null");
    }
    FileIdentity identity;
    if (!file_identity(path, &identity, error)) {
        return false;
    }
    return detect_with_identity(path, identity, detection, error);
}

bool BgzfIndex::open(const std::string& path,
                     uint32_t reader_threads,
                     std::string* error) {
    path_.clear();
    blocks_.clear();
    recordCount_ = 0;
    FileIdentity identity;
    if (!file_identity(path, &identity, error)) {
        return false;
    }
    BgzfDetection detection;
    if (!detect_with_identity(path, identity, &detection, error)) {
        return false;
    }
    if (!detection.isBgzf) {
        return set_error(error, "input is not BGZF: " + path);
    }
    path_ = path;
    fileSize_ = identity.size;
    mtimeSeconds_ = identity.mtimeSeconds;
    mtimeNanoseconds_ = identity.mtimeNanoseconds;
    const std::string sidecar = path + ".bgzi";
    if (detection.hasEofMarker &&
        load_sidecar(sidecar, identity, &blocks_, &recordCount_)) {
        cacheStatus_ = BgzfCacheStatus::Loaded;
        return true;
    }
    if (!scan_blocks(path, identity.size, &blocks_, error) ||
        !count_records(path, reader_threads, &blocks_, &recordCount_, error)) {
        blocks_.clear();
        recordCount_ = 0;
        return false;
    }
    cacheStatus_ = BgzfCacheStatus::Rebuilt;
    write_sidecar(sidecar, identity, blocks_, recordCount_);
    return true;
}

const std::string& BgzfIndex::path() const {
    return path_;
}

const std::vector<BgzfBlock>& BgzfIndex::blocks() const {
    return blocks_;
}

uint64_t BgzfIndex::record_count() const {
    return recordCount_;
}

uint64_t BgzfIndex::file_size() const {
    return fileSize_;
}

BgzfCacheStatus BgzfIndex::cache_status() const {
    return cacheStatus_;
}

bool BgzfIndex::locate_record(uint64_t record_ordinal,
                              size_t* block_index,
                              uint32_t* decompressed_offset,
                              std::string* error) const {
    if (block_index == nullptr || decompressed_offset == nullptr) {
        return set_error(error, "BGZF record locator output is null");
    }
    if (record_ordinal >= recordCount_) {
        std::ostringstream message;
        message << "BGZF record ordinal " << record_ordinal
                << " is outside record count " << recordCount_;
        return set_error(error, message.str());
    }
    size_t low = 0;
    size_t high = blocks_.size();
    while (low < high) {
        const size_t middle = low + (high - low) / 2;
        const BgzfBlock& block = blocks_[middle];
        if (block.firstRecordOrdinal + block.recordsStarting <= record_ordinal) {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    if (low >= blocks_.size() || blocks_[low].recordsStarting == 0 ||
        !blocks_[low].hasFirstRecordOffset ||
        record_ordinal < blocks_[low].firstRecordOrdinal) {
        return set_error(error, "BGZF record ordinal could not be located in the index");
    }
    *block_index = low;
    *decompressed_offset = blocks_[low].firstRecordOffset;
    return true;
}

} // namespace input
} // namespace star
