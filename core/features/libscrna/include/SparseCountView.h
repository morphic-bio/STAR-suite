#ifndef STAR_SCRNA_SPARSE_COUNT_VIEW_H
#define STAR_SCRNA_SPARSE_COUNT_VIEW_H
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

// Borrowed, immutable words. Cell offsets and stride are in uint32_t words,
// relative to each gene/count base. The owner must outlive the synchronous call.
// Either 32-bit or 64-bit offsets can be supplied; no narrowing conversions.
struct SparseCountView {
    const uint32_t* genes = nullptr;
    const uint32_t* counts = nullptr;
    size_t geneWords = 0, countWords = 0, stride = 1, cells = 0;
    const uint32_t* offsets = nullptr;
    const uint64_t* offsets64 = nullptr;
    const uint32_t* entries = nullptr;
    uint64_t start(size_t cell) const { return offsets64 ? offsets64[cell] : offsets[cell]; }
    size_t position(size_t cell, size_t entry) const { return size_t(start(cell)) + entry * stride; }
    uint32_t gene(size_t cell, size_t entry) const { return genes[position(cell, entry)]; }
    uint32_t count(size_t cell, size_t entry) const { return counts[position(cell, entry)]; }
    void validate(size_t expectedCells) const {
        if (cells != expectedCells || !stride || (!offsets && !offsets64) || !entries)
            throw std::invalid_argument("Invalid sparse count view axes, offsets or stride");
        for (size_t c = 0; c < cells; ++c) {
            if (!entries[c]) continue;
            if (!genes || !counts || start(c) > std::numeric_limits<size_t>::max())
                throw std::invalid_argument("Missing sparse count data or offset overflow");
            const size_t begin = size_t(start(c)), n = entries[c] - 1;
            if (n > (std::numeric_limits<size_t>::max() - begin) / stride)
                throw std::overflow_error("Sparse count view offset overflow");
            const size_t last = begin + n * stride;
            if (last >= geneWords || last >= countWords)
                throw std::out_of_range("Sparse count view exceeds borrowed storage");
        }
    }
};
#endif
