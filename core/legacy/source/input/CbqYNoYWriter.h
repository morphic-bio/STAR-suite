#ifndef CODE_input_CbqYNoYWriter
#define CODE_input_CbqYNoYWriter

#include "input/CbqWriter.h"

#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

class Parameters;

namespace star {
namespace input {

struct CbqYNoYRecord {
    bool is_y;
    std::vector<CbqWriterSegment> segments;
};

class CbqYNoYOrderedWriter {
public:
    CbqYNoYOrderedWriter();
    CbqYNoYOrderedWriter(const CbqYNoYOrderedWriter&) = delete;
    CbqYNoYOrderedWriter& operator=(const CbqYNoYOrderedWriter&) = delete;

    bool open(const std::string& y_path,
              const std::string& noy_path,
              uint32_t mate_count,
              int compression_level,
              uint64_t block_size,
              std::string* error);
    bool submit_chunk(uint64_t chunk_id,
                      const std::vector<CbqYNoYRecord>& records,
                      std::string* error);
    bool finish(std::string* error);

private:
    std::mutex mutex_;
    std::condition_variable cv_;
    uint64_t next_chunk_id_;
    bool failed_;
    std::string failure_;
    bool opened_;
    CbqWriter y_writer_;
    CbqWriter noy_writer_;
};

bool openGlobalCbqYNoYWriter(const Parameters& P, std::string* error);
bool finishGlobalCbqYNoYWriter(std::string* error);

} // namespace input
} // namespace star

#endif
