#ifndef STAR_SPATIAL_R1_FASTQ_TAP_H
#define STAR_SPATIAL_R1_FASTQ_TAP_H

#include <cstddef>
#include <cstdint>
#include <string>

namespace spatial_r1_fastq_tap {

class Writer {
  public:
    Writer();
    ~Writer();
    Writer(const Writer &) = delete;
    Writer &operator=(const Writer &) = delete;

    bool open(const std::string &path, bool requireFifo, std::string &error);
    bool write(const std::string &readName, const std::string &sequence,
               const std::string &quality, std::string &error);
    bool close(std::string &error);

    bool isOpen() const { return fd_ >= 0; }
    std::uint64_t recordsWritten() const { return recordsWritten_; }

  private:
    bool writeAll(const char *data, std::size_t size, std::string &error);

    int fd_;
    std::string path_;
    std::uint64_t recordsWritten_;
};

} // namespace spatial_r1_fastq_tap

#endif
