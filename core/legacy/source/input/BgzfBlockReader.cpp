#include "input/BgzfBlockReader.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <sstream>
#include <sys/stat.h>
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
    size_t position = 12;
    const size_t end = 12 + available_extra;
    while (position + 4 <= end) {
        const uint16_t subfield_length = load_u16(header + position + 2);
        position += 4;
        if (position + subfield_length > end) {
            return false;
        }
        if (header[position - 4] == 'B' && header[position - 3] == 'C' &&
            subfield_length == 2) {
            if (compressed_size != nullptr) {
                *compressed_size = static_cast<uint32_t>(load_u16(header + position)) + 1U;
            }
            return true;
        }
        position += subfield_length;
    }
    return false;
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
        const ssize_t count = ::pread(input_fd, output + completed, size - completed,
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

bool file_size(const std::string& path, uint64_t* size, std::string* error) {
    struct stat info;
    if (::stat(path.c_str(), &info) != 0) {
        return set_error(error, "could not stat BGZF input " + path + ": " +
                                std::strerror(errno));
    }
    if (info.st_size < 0) {
        return set_error(error, "BGZF input has a negative file size: " + path);
    }
    *size = static_cast<uint64_t>(info.st_size);
    return true;
}

bool gzip_payload_offset(const unsigned char* member,
                         size_t member_size,
                         size_t* payload_offset,
                         std::string* error) {
    if (member == nullptr || member_size < 18 || member[0] != 0x1f || member[1] != 0x8b ||
        member[2] != 8) {
        return set_error(error, "invalid gzip header in BGZF member");
    }
    const unsigned char flags = member[3];
    size_t position = 10;
    if ((flags & 0x04U) != 0) {
        if (position + 2 > member_size) {
            return set_error(error, "truncated BGZF extra-field length");
        }
        const uint16_t xlen = load_u16(member + position);
        position += 2;
        if (xlen > member_size - position) {
            return set_error(error, "truncated BGZF extra field");
        }
        position += xlen;
    }
    const unsigned char zero_terminated_flags[2] = {0x08U, 0x10U};
    for (size_t flag_index = 0; flag_index < 2; ++flag_index) {
        if ((flags & zero_terminated_flags[flag_index]) == 0) {
            continue;
        }
        while (position < member_size && member[position] != 0) {
            ++position;
        }
        if (position == member_size) {
            return set_error(error, "unterminated optional BGZF header field");
        }
        ++position;
    }
    if ((flags & 0x02U) != 0) {
        if (position + 2 > member_size) {
            return set_error(error, "truncated BGZF header CRC");
        }
        position += 2;
    }
    if (position > member_size - 8) {
        return set_error(error, "BGZF header overlaps its trailer");
    }
    *payload_offset = position;
    return true;
}

} // namespace

bool detect_bgzf(const std::string& path,
                 BgzfDetection* detection,
                 std::string* error) {
    if (detection == nullptr) {
        return set_error(error, "BGZF detection output is null");
    }
    detection->isBgzf = false;
    detection->hasEofMarker = false;
    uint64_t size = 0;
    if (!file_size(path, &size, error)) {
        return false;
    }
    if (size == 0) {
        return true;
    }
    const int input_fd = ::open(path.c_str(), O_RDONLY);
    if (input_fd < 0) {
        return set_error(error, "could not open BGZF candidate " + path + ": " +
                                std::strerror(errno));
    }
    bool ok = true;
    if (size >= 18) {
        unsigned char header[18];
        ok = read_exact_fd(input_fd, 0, header, sizeof(header), error);
        if (ok) {
            detection->isBgzf = parse_bgzf_header18(header, nullptr);
        }
    }
    if (ok && detection->isBgzf && size >= sizeof(kBgzfEof)) {
        unsigned char eof[sizeof(kBgzfEof)];
        ok = read_exact_fd(input_fd, size - sizeof(kBgzfEof), eof, sizeof(eof), error);
        if (ok) {
            detection->hasEofMarker = std::memcmp(eof, kBgzfEof, sizeof(eof)) == 0;
        }
    }
    ::close(input_fd);
    return ok;
}

bool read_bgzf_block_header_fd(int input_fd,
                               uint64_t compressed_offset,
                               uint64_t range_end,
                               BgzfBlock* block,
                               bool* is_eof_marker,
                               std::string* error) {
    if (input_fd < 0 || block == nullptr || is_eof_marker == nullptr) {
        return set_error(error, "invalid BGZF block-header reader arguments");
    }
    *is_eof_marker = false;
    if (compressed_offset > range_end || range_end - compressed_offset < 18) {
        std::ostringstream message;
        message << "truncated BGZF header at block offset " << compressed_offset;
        return set_error(error, message.str());
    }
    unsigned char header[18];
    std::string read_error;
    if (!read_exact_fd(input_fd, compressed_offset, header, sizeof(header), &read_error)) {
        std::ostringstream message;
        message << "could not read BGZF header at block offset " << compressed_offset
                << ": " << read_error;
        return set_error(error, message.str());
    }
    uint32_t compressed_size = 0;
    if (!parse_bgzf_header18(header, &compressed_size)) {
        std::ostringstream message;
        message << "invalid BGZF header at block offset " << compressed_offset;
        return set_error(error, message.str());
    }
    if (compressed_size < 26 || compressed_size > 65536) {
        std::ostringstream message;
        message << "invalid BGZF member size " << compressed_size
                << " at block offset " << compressed_offset;
        return set_error(error, message.str());
    }
    if (compressed_size > range_end - compressed_offset) {
        std::ostringstream message;
        message << "truncated BGZF member at block offset " << compressed_offset
                << " (declared size " << compressed_size << ')';
        return set_error(error, message.str());
    }

    if (compressed_size == sizeof(kBgzfEof)) {
        unsigned char candidate[sizeof(kBgzfEof)];
        if (!read_exact_fd(input_fd, compressed_offset, candidate, sizeof(candidate),
                           &read_error)) {
            return set_error(error, read_error);
        }
        if (std::memcmp(candidate, kBgzfEof, sizeof(candidate)) == 0) {
            if (compressed_offset + compressed_size != range_end) {
                std::ostringstream message;
                message << "BGZF EOF marker is not final at block offset "
                        << compressed_offset;
                return set_error(error, message.str());
            }
            *is_eof_marker = true;
            block->compressedOffset = compressed_offset;
            block->compressedSize = compressed_size;
            block->isize = 0;
            return true;
        }
    }

    unsigned char trailer_isize[4];
    if (!read_exact_fd(input_fd, compressed_offset + compressed_size - 4,
                       trailer_isize, sizeof(trailer_isize), &read_error)) {
        std::ostringstream message;
        message << "could not read BGZF trailer at block offset " << compressed_offset
                << ": " << read_error;
        return set_error(error, message.str());
    }
    block->compressedOffset = compressed_offset;
    block->compressedSize = compressed_size;
    block->isize = load_u32(trailer_isize);
    return true;
}

bool read_bgzf_work_fd(int input_fd,
                       uint64_t compressed_offset,
                       uint64_t range_end,
                       uint64_t target_compressed_bytes,
                       size_t max_blocks,
                       std::vector<BgzfBlock>* blocks,
                       std::vector<unsigned char>* compressed_bytes,
                       bool* reached_end,
                       std::string* error) {
    if (input_fd < 0 || blocks == nullptr || compressed_bytes == nullptr ||
        reached_end == nullptr || max_blocks == 0) {
        return set_error(error, "invalid BGZF work-reader arguments");
    }
    blocks->clear();
    compressed_bytes->clear();
    *reached_end = false;
    if (compressed_offset > range_end) {
        return set_error(error, "BGZF work offset is beyond its range end");
    }
    if (compressed_offset == range_end) {
        *reached_end = true;
        return true;
    }

    const uint64_t remaining = range_end - compressed_offset;
    const uint64_t target = std::max<uint64_t>(target_compressed_bytes, 1);
    const uint64_t lookahead = 65536U;
    const uint64_t wanted = target > std::numeric_limits<uint64_t>::max() - lookahead
        ? remaining : target + lookahead;
    const uint64_t request64 = std::min<uint64_t>(remaining, wanted);
    if (request64 > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        return set_error(error, "BGZF work claim exceeds platform memory limits");
    }
    compressed_bytes->resize(static_cast<size_t>(request64));
    std::string read_error;
    if (!read_exact_fd(input_fd, compressed_offset, compressed_bytes->data(),
                       compressed_bytes->size(), &read_error)) {
        std::ostringstream message;
        message << "could not read BGZF work at block offset " << compressed_offset
                << ": " << read_error;
        return set_error(error, message.str());
    }

    size_t position = 0;
    while (blocks->size() < max_blocks &&
           (blocks->empty() || position < target)) {
        const uint64_t block_offset = compressed_offset + position;
        if (block_offset > range_end || range_end - block_offset < 18 ||
            compressed_bytes->size() - position < 18) {
            std::ostringstream message;
            message << "truncated BGZF header at block offset " << block_offset;
            return set_error(error, message.str());
        }
        const unsigned char* member = compressed_bytes->data() + position;
        uint32_t compressed_size = 0;
        if (!parse_bgzf_header18(member, &compressed_size)) {
            std::ostringstream message;
            message << "invalid BGZF header at block offset " << block_offset;
            return set_error(error, message.str());
        }
        if (compressed_size < 26 || compressed_size > 65536) {
            std::ostringstream message;
            message << "invalid BGZF member size " << compressed_size
                    << " at block offset " << block_offset;
            return set_error(error, message.str());
        }
        if (compressed_size > range_end - block_offset) {
            std::ostringstream message;
            message << "truncated BGZF member at block offset " << block_offset
                    << " (declared size " << compressed_size << ')';
            return set_error(error, message.str());
        }
        if (compressed_size > compressed_bytes->size() - position) {
            return set_error(error, "BGZF work lookahead did not cover a full member");
        }
        if (compressed_size == sizeof(kBgzfEof) &&
            std::memcmp(member, kBgzfEof, sizeof(kBgzfEof)) == 0) {
            if (block_offset + compressed_size != range_end) {
                std::ostringstream message;
                message << "BGZF EOF marker is not final at block offset "
                        << block_offset;
                return set_error(error, message.str());
            }
            *reached_end = true;
            compressed_bytes->resize(position);
            return true;
        }

        BgzfBlock block;
        block.compressedOffset = block_offset;
        block.compressedSize = compressed_size;
        block.isize = load_u32(member + compressed_size - 4);
        blocks->push_back(block);
        position += compressed_size;
        if (compressed_offset + position == range_end) {
            *reached_end = true;
            break;
        }
    }
    compressed_bytes->resize(position);
    return true;
}

bool inflate_bgzf_block_buffer(const unsigned char* member,
                               size_t member_size,
                               uint64_t compressed_offset,
                               bool check_crc,
                               unsigned char* output,
                               size_t output_size,
                               std::string* error) {
    if (member == nullptr || (output == nullptr && output_size != 0)) {
        return set_error(error, "invalid buffered BGZF inflate arguments");
    }
    if (member_size < 26) {
        std::ostringstream message;
        message << "invalid BGZF member at block offset " << compressed_offset;
        return set_error(error, message.str());
    }
    size_t payload_offset = 0;
    std::string header_error;
    if (!gzip_payload_offset(member, member_size, &payload_offset, &header_error)) {
        std::ostringstream message;
        message << header_error << " at block offset " << compressed_offset;
        return set_error(error, message.str());
    }
    const uint32_t expected_crc = load_u32(member + member_size - 8);
    const uint32_t expected_isize = load_u32(member + member_size - 4);
    if (output_size != expected_isize) {
        std::ostringstream message;
        message << "BGZF ISIZE mismatch at block offset " << compressed_offset;
        return set_error(error, message.str());
    }
    const size_t deflate_size = member_size - payload_offset - 8;

    unsigned char empty_output = 0;
    z_stream stream;
    std::memset(&stream, 0, sizeof(stream));
    stream.next_in = const_cast<Bytef*>(member + payload_offset);
    stream.avail_in = static_cast<uInt>(deflate_size);
    stream.next_out = output_size == 0 ? &empty_output : output;
    stream.avail_out = static_cast<uInt>(output_size == 0 ? 1 : output_size);
    const int init_status = inflateInit2(&stream, -15);
    if (init_status != Z_OK) {
        return set_error(error, "zlib could not initialize raw BGZF inflate");
    }
    const int inflate_status = ::inflate(&stream, Z_FINISH);
    inflateEnd(&stream);
    if (inflate_status != Z_STREAM_END || stream.total_out != expected_isize ||
        stream.total_in != deflate_size) {
        std::ostringstream message;
        message << "raw inflate failed at block offset " << compressed_offset
                << " (zlib status " << inflate_status << ')';
        return set_error(error, message.str());
    }
    if (check_crc) {
        const Bytef* data = output_size == 0
            ? reinterpret_cast<const Bytef*>("") : output;
        const uint32_t observed_crc = static_cast<uint32_t>(
            crc32(0L, data, static_cast<uInt>(output_size)));
        if (observed_crc != expected_crc) {
            std::ostringstream message;
            message << "BGZF CRC mismatch at block offset " << compressed_offset;
            return set_error(error, message.str());
        }
    }
    return true;
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
    const uint32_t expected_isize = load_u32(member.data() + member.size() - 4);
    if (expected_isize != block.isize) {
        std::ostringstream message;
        message << "BGZF ISIZE changed at block offset " << block.compressedOffset;
        return set_error(error, message.str());
    }

    output->resize(block.isize);
    return inflate_bgzf_block_buffer(member.data(), member.size(),
                                     block.compressedOffset, check_crc,
                                     output->empty() ? nullptr : output->data(),
                                     output->size(), error);
}

} // namespace input
} // namespace star
