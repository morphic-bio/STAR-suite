#ifndef STAR_POOLED_TRANSCRIPT_MAP_H
#define STAR_POOLED_TRANSCRIPT_MAP_H

#include "htslib/khash.h"
#include <cstdint>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

struct TranscriptRange { size_t offset; uint32_t count; };
KHASH_INIT(transcriptRanges, uint32_t, TranscriptRange, 1, __ac_Wang_hash, kh_int_hash_equal)

// Each UMI owns a slice in one per-cell transcript buffer. Intersection only
// shrinks a slice, so subsequent observations can update it without allocating.
template<class Transcript> class PooledTranscriptMap {
    khash_t(transcriptRanges)* table_ = nullptr;
    std::vector<Transcript> transcripts_;
    size_t liveSlots_ = 0;
    void initialize() {
        if (!table_) {
            table_ = kh_init(transcriptRanges);
            if (!table_) throw std::bad_alloc();
        }
    }
public:
    struct View {
        const Transcript* data;
        size_t count;
        bool empty() const { return count == 0; }
        size_t size() const { return count; }
        const Transcript& operator[](size_t i) const { return data[i]; }
        const Transcript* begin() const { return data; }
        const Transcript* end() const { return count ? data + count : data; }
    };
    class const_iterator {
        const PooledTranscriptMap* owner_;
        khint_t k_;
        void skip() {
            const auto* h = owner_->table_;
            while (h && k_ != kh_end(h) && !kh_exist(h, k_)) ++k_;
        }
    public:
        const_iterator(const PooledTranscriptMap* owner, khint_t k) : owner_(owner), k_(k) { skip(); }
        std::pair<uint32_t, View> operator*() const {
            const auto& range = kh_val(owner_->table_, k_);
            const Transcript* data = range.count ? owner_->transcripts_.data() + range.offset : nullptr;
            return {kh_key(owner_->table_, k_), {data, range.count}};
        }
        const_iterator& operator++() { ++k_; skip(); return *this; }
        bool operator!=(const const_iterator& other) const { return owner_ != other.owner_ || k_ != other.k_; }
    };
    PooledTranscriptMap() = default;
    PooledTranscriptMap(const PooledTranscriptMap&) = delete;
    PooledTranscriptMap& operator=(const PooledTranscriptMap&) = delete;
    PooledTranscriptMap(PooledTranscriptMap&& other) noexcept { swap(other); }
    PooledTranscriptMap& operator=(PooledTranscriptMap&& other) noexcept {
        if (this != &other) { clearAndFree(); swap(other); }
        return *this;
    }
    ~PooledTranscriptMap() { kh_destroy(transcriptRanges, table_); }
    void swap(PooledTranscriptMap& other) noexcept {
        std::swap(table_, other.table_);
        transcripts_.swap(other.transcripts_);
        std::swap(liveSlots_, other.liveSlots_);
    }
    void clearAndFree() {
        kh_destroy(transcriptRanges, table_); table_ = nullptr;
        std::vector<Transcript>().swap(transcripts_);
        liveSlots_ = 0;
    }
    size_t size() const { return table_ ? kh_size(table_) : 0; }
    size_t allocatedSlots() const { return transcripts_.size(); }
    size_t liveSlots() const { return liveSlots_; }
    bool rejected(uint32_t umi) const {
        if (!table_) return false;
        const khint_t k = kh_get(transcriptRanges, table_, umi);
        return k != kh_end(table_) && kh_val(table_, k).count == 0;
    }
    void reserve(size_t entries) {
        if (entries == 0) return;
        const size_t maxBuckets = size_t{1} << (sizeof(khint_t) * 8 - 1);
        if (entries > static_cast<size_t>(maxBuckets * .77) - 1)
            throw std::length_error("Transcript map capacity exceeds khash indices");
        initialize();
        const size_t buckets = static_cast<size_t>(entries / .77) + 1;
        if (buckets > kh_n_buckets(table_) && kh_resize(transcriptRanges, table_, buckets) < 0)
            throw std::bad_alloc();
    }
    void merge(uint32_t umi, const std::vector<Transcript>& incoming) {
        initialize();
        int absent;
        const khint_t k = kh_put(transcriptRanges, table_, umi, &absent);
        if (absent < 0) throw std::bad_alloc();
        auto& range = kh_val(table_, k);
        if (absent) {
            if (incoming.size() > UINT32_MAX || incoming.size() > transcripts_.max_size() - transcripts_.size())
                throw std::length_error("Transcript slice size overflow");
            range.offset = transcripts_.size();
            range.count = static_cast<uint32_t>(incoming.size());
            transcripts_.insert(transcripts_.end(), incoming.begin(), incoming.end());
            liveSlots_ += range.count;
            return;
        }
        if (range.count == 0) return; // Empty intersections stay rejected.
        size_t next = 0;
        uint32_t written = 0;
        for (uint32_t old = 0; old < range.count; ++old) {
            const Transcript previous = transcripts_[range.offset + old];
            while (next < incoming.size() && previous.tr > incoming[next].tr) ++next;
            if (next == incoming.size()) break;
            if (previous.tr == incoming[next].tr)
                transcripts_[range.offset + written++] = {incoming[next].tr,
                    static_cast<uint8_t>(previous.type | incoming[next].type)};
        }
        liveSlots_ -= range.count - written;
        range.count = written;
    }
    const_iterator begin() const { return {this, 0}; }
    const_iterator end() const { return {this, table_ ? kh_end(table_) : 0}; }
};
#endif
