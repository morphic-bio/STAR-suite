#ifndef CODE_input_CbqWriter
#define CODE_input_CbqWriter

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace star {
namespace input {

struct CbqWriterSegment {
    std::string header_payload;
    std::string sequence;
    std::string quality;
};

class CbqWriter {
public:
    CbqWriter();
    CbqWriter(const CbqWriter&) = delete;
    CbqWriter& operator=(const CbqWriter&) = delete;
    ~CbqWriter();

    bool open(const std::string& path,
              uint32_t mate_count,
              int compression_level,
              uint64_t block_size,
              std::string* error);
    bool add_record(const std::vector<CbqWriterSegment>& segments, std::string* error);
    bool finish(std::string* error);
    bool is_open() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace input
} // namespace star

#endif
