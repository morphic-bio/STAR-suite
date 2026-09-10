#ifndef STAR_BGZF_PIPE_GROUP_H
#define STAR_BGZF_PIPE_GROUP_H
#include "input/BgzfRangeReader.h"
#include <atomic>
#include <cerrno>
#include <csignal>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <sstream>
#include <chrono>
#include <unistd.h>

namespace star { namespace input {
// Raw, ordered BGZF bytes feed the established STAR chunk parser. This keeps
// FASTA/FASTQ headers, lanes, reopening, limits, taps and batch markers intact.
class BgzfPipeGroup {
    std::vector<std::thread> producers_;
    std::mutex mutex_;
    std::string error_;
    std::atomic<bool> permitsEnabled_{false};
    BgzfWorkPermitHooks external_;
    std::atomic<uint64_t> decodeWork_{0}, decodeBytes_{0}, decodeNs_{0}, permitWaitNs_{0};
    static uint64_t acquire(void *p) {
        auto *self = static_cast<BgzfPipeGroup*>(p);
        return self->permitsEnabled_.load() && self->external_.enabled()
            ? self->external_.acquire(self->external_.context) : UINT64_MAX;
    }
    static void release(void *p, uint64_t wait, uint64_t units, uint64_t bytes, uint64_t ns) {
        auto *self = static_cast<BgzfPipeGroup*>(p);
        self->decodeWork_.fetch_add(units, std::memory_order_relaxed);
        self->decodeBytes_.fetch_add(bytes, std::memory_order_relaxed);
        self->decodeNs_.fetch_add(ns, std::memory_order_relaxed);
        if (wait != UINT64_MAX) self->permitWaitNs_.fetch_add(wait, std::memory_order_relaxed);
        if (wait != UINT64_MAX) self->external_.release(self->external_.context, wait, units, bytes, ns);
    }
    static bool writeAll(int fd, const char *p, size_t n) {
        while (n) {
            ssize_t k = ::write(fd, p, n);
            if (k < 0 && errno == EINTR) continue;
            if (k < 0 && errno == EPIPE) return false; // read limit / pass cancellation
            if (k <= 0) throw std::runtime_error(std::string("BGZF FIFO write failed: ") + std::strerror(errno));
            p += k; n -= size_t(k);
        }
        return true;
    }
public:
    explicit BgzfPipeGroup(const BgzfWorkPermitHooks& hooks) : external_(hooks) {}
    ~BgzfPipeGroup() { join(); }
    void enablePermits() { permitsEnabled_.store(true); }
    void join() { for (auto& thread : producers_) if (thread.joinable()) thread.join(); }
    std::string error() { std::lock_guard<std::mutex> lock(mutex_); return error_; }
    std::string summary() const {
        std::ostringstream out;
        out << "BGZF raw input totals: decode_work=" << decodeWork_.load()
            << " decode_bytes=" << decodeBytes_.load()
            << " inflate_seconds_sum=" << decodeNs_.load() / 1e9
            << " permit_wait_seconds_sum=" << permitWaitNs_.load() / 1e9 << "\n";
        return out.str();
    }
    void start(std::vector<std::string> paths, std::vector<bool> native,
               std::string fifo, unsigned threads, bool crc) {
        producers_.emplace_back([this, paths, native, fifo, threads, crc] {
            sigset_t blocked; sigemptyset(&blocked); sigaddset(&blocked, SIGPIPE);
            pthread_sigmask(SIG_BLOCK, &blocked, nullptr);
            int fd = -1;
            try {
                fd = ::open(fifo.c_str(), O_WRONLY);
                if (fd < 0) throw std::runtime_error("could not open BGZF FIFO " + fifo);
                bool live = true;
                for (size_t lane = 0; lane < paths.size() && live; ++lane) {
                    const std::string marker = "FILE " + std::to_string(lane) + "\n";
                    if (!writeAll(fd, marker.data(), marker.size())) break;
                    if (native[lane]) {
                        BgzfRangeReader reader;
                        std::string message;
                        BgzfWorkPermitHooks hooks; hooks.context = this; hooks.acquire = acquire; hooks.release = release;
                        if (!reader.open(paths[lane], 0, UINT64_MAX, threads, crc, &message, &hooks))
                            throw std::runtime_error(message);
                        const char *data; size_t size;
                        while (reader.next_bytes(&data, &size, &message)) {
                            if (!writeAll(fd, data, size)) { live = false; break; }
                        }
                        if (!message.empty()) throw std::runtime_error(message);
                    } else {
                        gzFile gz = gzopen(paths[lane].c_str(), "rb");
                        if (!gz) throw std::runtime_error("could not open mixed FASTQ lane " + paths[lane]);
                        std::vector<char> buffer(1<<20);
                        int n = 0;
                        try {
                            while ((n = gzread(gz, buffer.data(), buffer.size())) > 0)
                                if (!writeAll(fd, buffer.data(), n)) { live = false; break; }
                            if (n < 0) throw std::runtime_error("gzread failed for " + paths[lane]);
                        } catch (...) { gzclose(gz); throw; }
                        gzclose(gz);
                    }
                }
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lock(mutex_); if (error_.empty()) error_ = e.what();
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex_); if (error_.empty()) error_ = "unknown BGZF FIFO exception";
            }
            if (fd >= 0) ::close(fd);
        });
    }
};
}}
#endif
