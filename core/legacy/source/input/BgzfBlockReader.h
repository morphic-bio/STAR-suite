#ifndef CODE_input_BgzfBlockReader
#define CODE_input_BgzfBlockReader

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace star {
namespace input {

struct BgzfDetection {
    bool isBgzf = false;
    bool hasEofMarker = false;
};

struct BgzfBlock {
    uint64_t compressedOffset = 0;
    uint32_t compressedSize = 0;
    uint32_t isize = 0;
};

bool detect_bgzf(const std::string& path,
                 BgzfDetection* detection,
                 std::string* error);

// Reads the BC/BSIZE metadata for exactly one member at compressed_offset.
// range_end is the exclusive physical byte boundary owned by the stream.
bool read_bgzf_block_header_fd(int input_fd,
                               uint64_t compressed_offset,
                               uint64_t range_end,
                               BgzfBlock* block,
                               bool* is_eof_marker,
                               std::string* error);

// Reads one bounded contiguous claim at compressed_offset. Every member
// boundary is discovered from its inline BC/BSIZE field in the returned byte
// window; no index or pre-scan is involved. A final canonical EOF member is
// validated but excluded from blocks and compressed_bytes.
bool read_bgzf_work_fd(int input_fd,
                       uint64_t compressed_offset,
                       uint64_t range_end,
                       uint64_t target_compressed_bytes,
                       size_t max_blocks,
                       std::vector<BgzfBlock>* blocks,
                       std::vector<unsigned char>* compressed_bytes,
                       bool* reached_end,
                       std::string* error);

// Inflates a member already resident in a claimed compressed-byte window
// directly into caller-owned storage.
bool inflate_bgzf_block_buffer(const unsigned char* member,
                               size_t member_size,
                               uint64_t compressed_offset,
                               bool check_crc,
                               unsigned char* output,
                               size_t output_size,
                               std::string* error);

bool inflate_bgzf_block_fd(int input_fd,
                           const BgzfBlock& block,
                           bool check_crc,
                           std::vector<unsigned char>* output,
                           std::string* error);

} // namespace input
} // namespace star

#endif
