#ifndef STAR_BORROWED_BARCODE_INDEX_H
#define STAR_BORROWED_BARCODE_INDEX_H

#include "klib/khash.h"
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>

KHASH_MAP_INIT_STR(borrowed_barcode_cell, uint32_t)

// The caller owns the strings and must keep their character storage unchanged
// until this index is destroyed. After construction, concurrent lookup is safe.
class BorrowedBarcodeIndex {
public:
    explicit BorrowedBarcodeIndex(size_t expected) {
        if (expected > std::numeric_limits<khint_t>::max() / 2)
            throw std::length_error("Barcode index exceeds khash capacity");
        if (kh_resize(borrowed_barcode_cell, &table_, expected + expected / 2 + 1) < 0) {
            release();
            throw std::bad_alloc();
        }
    }
    ~BorrowedBarcodeIndex() { release(); }
    BorrowedBarcodeIndex(const BorrowedBarcodeIndex&) = delete;
    BorrowedBarcodeIndex& operator=(const BorrowedBarcodeIndex&) = delete;
    void insert(const std::string& barcode, uint32_t cell) {
        int inserted = 0;
        const auto at = kh_put(borrowed_barcode_cell, &table_, barcode.c_str(), &inserted);
        if (inserted < 0) throw std::bad_alloc();
        kh_val(&table_, at) = cell; // preserve last-value behavior for duplicate keys
    }
    uint32_t find(const std::string& barcode) const {
        const auto at = kh_get(borrowed_barcode_cell, &table_, barcode.c_str());
        return at == kh_end(&table_) ? UINT32_MAX : kh_val(&table_, at);
    }
private:
    khash_t(borrowed_barcode_cell) table_{};
    void release() {
        std::free(table_.flags);
        std::free(table_.keys);
        std::free(table_.vals);
    }
};
#endif
