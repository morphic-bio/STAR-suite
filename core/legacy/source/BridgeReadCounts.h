#ifndef STAR_BRIDGE_READ_COUNTS_H
#define STAR_BRIDGE_READ_COUNTS_H

#include "htslib/khash.h"
#include <cstdint>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

KHASH_INIT(bridgeReadCounts, uint32_t, uint64_t, 1, __ac_Wang_hash, kh_int_hash_equal)

// Integer read accounting, owned by one mapper or the gather thread. Values
// retain the existing low32=unique/high32=multi packing and arithmetic.
class BridgeReadCounts {
    khash_t(bridgeReadCounts)* table_ = nullptr;
    void initialize() {
        if (!table_) {
            table_ = kh_init(bridgeReadCounts);
            if (!table_) throw std::bad_alloc();
        }
    }
public:
    class const_iterator {
        const khash_t(bridgeReadCounts)* table_;
        khint_t index_;
        void skip() {
            while (table_ && index_ != kh_end(table_) && !kh_exist(table_, index_)) ++index_;
        }
    public:
        const_iterator(const khash_t(bridgeReadCounts)* table, khint_t index)
            : table_(table), index_(index) { skip(); }
        std::pair<uint32_t, uint64_t> operator*() const {
            return {kh_key(table_, index_), kh_val(table_, index_)};
        }
        const_iterator& operator++() { ++index_; skip(); return *this; }
        bool operator!=(const const_iterator& other) const {
            return table_ != other.table_ || index_ != other.index_;
        }
    };
    BridgeReadCounts() = default;
    ~BridgeReadCounts() { kh_destroy(bridgeReadCounts, table_); }
    BridgeReadCounts(const BridgeReadCounts& other) {
        try {
            reserve(other.size());
            for (const auto& entry : other) (*this)[entry.first] = entry.second;
        } catch (...) {
            kh_destroy(bridgeReadCounts, table_);
            throw;
        }
    }
    BridgeReadCounts(BridgeReadCounts&& other) noexcept : table_(other.table_) { other.table_ = nullptr; }
    BridgeReadCounts& operator=(BridgeReadCounts other) { swap(other); return *this; }
    void swap(BridgeReadCounts& other) noexcept { std::swap(table_, other.table_); }
    size_t size() const { return table_ ? kh_size(table_) : 0; }
    bool empty() const { return size() == 0; }
    void clear() { if (table_) kh_clear(bridgeReadCounts, table_); }
    void reserve(size_t entries) {
        if (entries == 0) return;
        const size_t maxBuckets = size_t{1} << (sizeof(khint_t) * 8 - 1);
        if (entries > static_cast<size_t>(maxBuckets * 0.77) - 1)
            throw std::length_error("Bridge read counter capacity exceeds khash indices");
        initialize();
        const size_t buckets = static_cast<size_t>(entries / 0.77) + 1;
        if (buckets > kh_n_buckets(table_) && kh_resize(bridgeReadCounts, table_, buckets) < 0)
            throw std::bad_alloc();
    }
    uint64_t& operator[](uint32_t key) {
        initialize();
        int absent;
        const khint_t k = kh_put(bridgeReadCounts, table_, key, &absent);
        if (absent < 0) throw std::bad_alloc();
        if (absent) kh_val(table_, k) = 0;
        return kh_val(table_, k);
    }
    const_iterator begin() const { return {table_, 0}; }
    const_iterator end() const { return {table_, table_ ? kh_end(table_) : 0}; }
};
#endif
