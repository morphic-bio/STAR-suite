#ifndef CODE_input_BgzfBlockReader
#define CODE_input_BgzfBlockReader

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

bool inflate_bgzf_block_fd(int input_fd,
                           const BgzfBlock& block,
                           bool check_crc,
                           std::vector<unsigned char>* output,
                           std::string* error);

} // namespace input
} // namespace star

#endif
