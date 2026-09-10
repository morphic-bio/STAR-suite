#ifndef STAR_OCM_COUNT_STORE_H
#define STAR_OCM_COUNT_STORE_H

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

namespace ocm {

// Process-local scratch records, not a persistent interchange format. Keeping
// source order also preserves the order of entries within each caller column.
struct CountRecord { uint32_t gene, cell, count; };
static_assert(sizeof(CountRecord) == 12, "OCM scratch record layout");

struct CountBufferBudget {
    explicit CountBufferBudget(uint64_t maximum) : limit(maximum) {}
    uint64_t limit, used = 0, peak = 0;
    bool reserve(uint64_t bytes) {
        if (bytes > limit - used) return false;
        used += bytes;
        peak = std::max(peak, used);
        return true;
    }
};

// Routing is serial. After finish(), a single sample task owns each store.
// Retained blocks stay charged until all sample tasks have joined; their budget
// is subtracted from matrix admission. Spilled stores use one fixed I/O buffer
// each and bulk I/O. Unlinked temporary files cannot become stale artifacts.
class CountStore {
public:
    enum : size_t { blockRecords = 4096, blockBytes = blockRecords * sizeof(CountRecord) };

    CountStore(CountBufferBudget& budget, const std::string& tempDir)
        : budget_(budget), tempDir_(tempDir) {}
    CountStore(const CountStore&) = delete;
    CountStore& operator=(const CountStore&) = delete;
    ~CountStore() {
        if (file_) std::fclose(file_);
        budget_.used -= blocks_.size() * blockBytes;
    }

    void append(const CountRecord& record) {
        if (finished_) throw std::logic_error("Append to finished OCM count store");
        if (size_ == std::numeric_limits<uint64_t>::max())
            throw std::overflow_error("OCM scratch record count overflow");
        if (!file_ && size_ % blockRecords == 0) {
            if (budget_.reserve(blockBytes)) {
                try {
                    std::unique_ptr<CountRecord[]> block(new CountRecord[blockRecords]);
                    blocks_.push_back(std::move(block));
                } catch (...) {
                    budget_.used -= blockBytes;
                    throw;
                }
            } else {
                spill();
            }
        }
        if (file_) {
            io_[pending_++] = record;
            if (pending_ == blockRecords) flush();
        } else {
            blocks_.back()[size_ % blockRecords] = record;
        }
        ++size_;
    }

    void finish() {
        if (finished_) return;
        if (file_) {
            flush();
            if (std::fflush(file_) != 0) throw std::runtime_error("Flush OCM count spill failed");
        }
        finished_ = true;
    }

    uint64_t size() const { return size_; }
    bool spilled() const { return file_ != nullptr; }

    template<class Function> void forEach(const Function& function) {
        if (!finished_) throw std::logic_error("Read unfinished OCM count store");
        if (file_) {
            if (fseeko(file_, 0, SEEK_SET) != 0)
                throw std::runtime_error("Rewind OCM count spill failed");
            for (uint64_t begin = 0; begin < size_;) {
                const size_t n = static_cast<size_t>(std::min<uint64_t>(blockRecords, size_ - begin));
                if (std::fread(io_.data(), sizeof(CountRecord), n, file_) != n)
                    throw std::runtime_error("Short read from OCM count spill");
                for (size_t i = 0; i < n; ++i) function(io_[i]);
                begin += n;
            }
        } else {
            uint64_t begin = 0;
            for (const auto& block : blocks_) {
                const size_t n = static_cast<size_t>(std::min<uint64_t>(blockRecords, size_ - begin));
                for (size_t i = 0; i < n; ++i) function(block[i]);
                begin += n;
            }
        }
    }

private:
    void write(const CountRecord* records, size_t n) {
        if (std::fwrite(records, sizeof(CountRecord), n, file_) != n)
            throw std::runtime_error("Write OCM count spill failed");
    }
    void flush() {
        write(io_.data(), pending_);
        pending_ = 0;
    }
    void spill() {
        std::string path = tempDir_ + "/ocm-counts-XXXXXX";
        std::vector<char> name(path.begin(), path.end());
        name.push_back('\0');
        const int fd = mkstemp(name.data());
        if (fd < 0) throw std::runtime_error("Create OCM count spill failed in " + tempDir_);
        if (unlink(name.data()) != 0) {
            close(fd);
            std::remove(name.data());
            throw std::runtime_error("Unlink OCM count spill failed");
        }
        file_ = fdopen(fd, "w+b");
        if (!file_) {
            close(fd);
            throw std::runtime_error("Open OCM count spill stream failed");
        }
        uint64_t begin = 0;
        for (const auto& block : blocks_) {
            const size_t n = static_cast<size_t>(std::min<uint64_t>(blockRecords, size_ - begin));
            write(block.get(), n);
            begin += n;
        }
        budget_.used -= blocks_.size() * blockBytes;
        blocks_.clear();
    }

    CountBufferBudget& budget_;
    std::string tempDir_;
    std::vector<std::unique_ptr<CountRecord[]>> blocks_;
    std::array<CountRecord, blockRecords> io_;
    FILE* file_ = nullptr;
    uint64_t size_ = 0;
    size_t pending_ = 0;
    bool finished_ = false;
};

} // namespace ocm
#endif
