#ifndef CODE_input_BgzfIndex
#define CODE_input_BgzfIndex

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
    uint64_t firstRecordOrdinal = 0;
    uint32_t recordsStarting = 0;
    uint32_t firstRecordOffset = 0;
    bool hasFirstRecordOffset = false;
};

enum class BgzfCacheStatus {
    Rebuilt,
    Loaded
};

const char* bgzf_cache_status_name(BgzfCacheStatus status);

bool inflate_bgzf_block(const std::string& path,
                        const BgzfBlock& block,
                        bool check_crc,
                        std::vector<unsigned char>* output,
                        std::string* error);

bool inflate_bgzf_block_fd(int input_fd,
                           const BgzfBlock& block,
                           bool check_crc,
                           std::vector<unsigned char>* output,
                           std::string* error);

class BgzfIndex {
public:
    static bool detect(const std::string& path,
                       BgzfDetection* detection,
                       std::string* error);

    bool open(const std::string& path,
              uint32_t reader_threads,
              std::string* error);

    const std::string& path() const;
    const std::vector<BgzfBlock>& blocks() const;
    uint64_t record_count() const;
    uint64_t file_size() const;
    BgzfCacheStatus cache_status() const;

    bool locate_record(uint64_t record_ordinal,
                       size_t* block_index,
                       uint32_t* decompressed_offset,
                       std::string* error) const;

private:
    std::string path_;
    std::vector<BgzfBlock> blocks_;
    uint64_t recordCount_ = 0;
    uint64_t fileSize_ = 0;
    int64_t mtimeSeconds_ = 0;
    int64_t mtimeNanoseconds_ = 0;
    BgzfCacheStatus cacheStatus_ = BgzfCacheStatus::Rebuilt;
};

} // namespace input
} // namespace star

#endif
