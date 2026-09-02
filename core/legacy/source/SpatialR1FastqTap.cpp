#include "SpatialR1FastqTap.h"

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

namespace spatial_r1_fastq_tap {

Writer::Writer() : fd_(-1), recordsWritten_(0) {}

Writer::~Writer()
{
    if (fd_ >= 0) ::close(fd_);
}

bool Writer::open(const std::string &path, bool requireFifo, std::string &error)
{
    if (isOpen() || path.empty() || path == "-" || path == "None") {
        error = "invalid or already-open raw-R1 FASTQ tap path";
        return false;
    }
    struct stat info;
    if (::lstat(path.c_str(), &info) != 0) {
        error = "cannot stat raw-R1 FASTQ tap " + path + ": " + std::strerror(errno);
        return false;
    }
    if (requireFifo && !S_ISFIFO(info.st_mode)) {
        error = "raw-R1 FASTQ tap is not a FIFO: " + path;
        return false;
    }
    const int flags = O_WRONLY
#ifdef O_CLOEXEC
        | O_CLOEXEC
#endif
        ;
    fd_ = ::open(path.c_str(), flags);
    if (fd_ < 0) {
        error = "cannot open raw-R1 FASTQ tap " + path + ": " + std::strerror(errno);
        return false;
    }
    path_ = path;
    recordsWritten_ = 0;
    return true;
}

bool Writer::writeAll(const char *data, std::size_t size, std::string &error)
{
    std::size_t offset = 0;
    while (offset < size) {
        const ssize_t count = ::write(fd_, data + offset, size - offset);
        if (count > 0) {
            offset += static_cast<std::size_t>(count);
            continue;
        }
        if (count < 0 && errno == EINTR) continue;
        error = "failed writing raw-R1 FASTQ tap " + path_ + ": "
            + (count == 0 ? std::string("short write") : std::strerror(errno));
        return false;
    }
    return true;
}

bool Writer::write(const std::string &readName, const std::string &sequence,
                   const std::string &quality, std::string &error)
{
    if (!isOpen()) {
        error = "raw-R1 FASTQ tap is not open";
        return false;
    }
    if (readName.empty() || readName.find_first_of("\t\r\n ") != std::string::npos) {
        error = "raw-R1 FASTQ tap received an invalid read name";
        return false;
    }
    if (sequence.empty() || sequence.size() != quality.size()
        || sequence.find_first_of("\r\n") != std::string::npos
        || quality.find_first_of("\r\n") != std::string::npos) {
        error = "raw-R1 FASTQ tap received invalid sequence/quality lengths";
        return false;
    }
    std::string record;
    record.reserve(readName.size() + sequence.size() + quality.size() + 8);
    record.push_back('@');
    record += readName;
    record.push_back('\n');
    record += sequence;
    record += "\n+\n";
    record += quality;
    record.push_back('\n');
    if (!writeAll(record.data(), record.size(), error)) return false;
    ++recordsWritten_;
    return true;
}

bool Writer::close(std::string &error)
{
    if (fd_ < 0) return true;
    const int descriptor = fd_;
    fd_ = -1;
    if (::close(descriptor) != 0) {
        error = "cannot close raw-R1 FASTQ tap " + path_ + ": " + std::strerror(errno);
        return false;
    }
    return true;
}

} // namespace spatial_r1_fastq_tap
