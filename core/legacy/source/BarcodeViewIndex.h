#ifndef STAR_BARCODE_VIEW_INDEX_H
#define STAR_BARCODE_VIEW_INDEX_H

#include "htslib/khash.h"
#include <cstdint>
#include <cstring>
#include <new>
#include <stdexcept>
#include <string>

struct BarcodeView {
    const char* data;
    size_t size;
    explicit BarcodeView(const std::string& s) : data(s.data()), size(s.size()) {}
    BarcodeView(const char* p, size_t n) : data(p), size(n) {}
};
inline khint_t barcodeViewHash(BarcodeView s) {
    uint64_t hash = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < s.size; ++i) {
        hash ^= static_cast<unsigned char>(s.data[i]);
        hash *= UINT64_C(1099511628211);
    }
    return static_cast<khint_t>(hash ^ (hash >> 32));
}
inline bool barcodeViewEqual(BarcodeView a, BarcodeView b) {
    return a.size == b.size && (a.size == 0 || std::memcmp(a.data, b.data, a.size) == 0);
}
KHASH_INIT(barcodeViewIndex, BarcodeView, uint32_t, 1, barcodeViewHash, barcodeViewEqual)

// Owns hash storage, borrows key bytes. Inserted strings must stay alive and
// unmoved/unmodified until clear/destruction; query strings need only live for
// find(). Length-aware keys retain general string semantics, including NULs.
class BarcodeViewIndex {
    khash_t(barcodeViewIndex)* table_ = nullptr;
    void initialize() {
        if (!table_) {
            table_ = kh_init(barcodeViewIndex);
            if (!table_) throw std::bad_alloc();
        }
    }
public:
    BarcodeViewIndex() = default;
    BarcodeViewIndex(const BarcodeViewIndex&) = delete;
    BarcodeViewIndex& operator=(const BarcodeViewIndex&) = delete;
    ~BarcodeViewIndex() { kh_destroy(barcodeViewIndex, table_); }
    size_t size() const { return table_ ? kh_size(table_) : 0; }
    void clear() { if (table_) kh_clear(barcodeViewIndex, table_); }
    void reserve(size_t entries) {
        if (entries == 0) return;
        const size_t maxBuckets = size_t{1} << (sizeof(khint_t) * 8 - 1);
        if (entries > static_cast<size_t>(maxBuckets * .77) - 1)
            throw std::length_error("Barcode index capacity exceeds khash indices");
        initialize();
        const size_t buckets = static_cast<size_t>(entries / .77) + 1;
        if (buckets > kh_n_buckets(table_) && kh_resize(barcodeViewIndex, table_, buckets) < 0)
            throw std::bad_alloc();
    }
    void insert(BarcodeView key, uint32_t value) {
        initialize();
        int absent;
        const khint_t k = kh_put(barcodeViewIndex, table_, key, &absent);
        if (absent < 0) throw std::bad_alloc();
        kh_val(table_, k) = value;
    }
    void insert(const std::string& key, uint32_t value) { insert(BarcodeView(key), value); }
    void insert(std::string&&, uint32_t) = delete;
    const uint32_t* find(BarcodeView key) const {
        if (!table_) return nullptr;
        const khint_t k = kh_get(barcodeViewIndex, table_, key);
        return k == kh_end(table_) ? nullptr : &kh_val(table_, k);
    }
    const uint32_t* find(const std::string& key) const { return find(BarcodeView(key)); }
};
#endif
