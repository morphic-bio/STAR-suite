#ifndef STAR_BARCODE_EDIT_DISTANCE_H
#define STAR_BARCODE_EDIT_DISTANCE_H

#include <cstdint>

namespace barcode_edit {
// Unit-cost global Levenshtein distance for patterns of at most 64 bases.
struct Myers64Pattern {
    uint64_t peq[256] = {};
    int length;
    uint64_t mask;
    uint64_t last;

    Myers64Pattern(const char* pattern, int m)
        : length(m), mask(m == 64 ? ~0ULL : ((1ULL << m) - 1ULL)),
          last(m == 0 ? 0 : (1ULL << (m - 1))) {
        for (int i = 0; i < m; ++i) {
            peq[static_cast<unsigned char>(pattern[i])] |= (1ULL << i);
        }
    }

    int distance(const char* text, int n) const {
        if (length == 0) return n;
        if (n == 0) return length;
        uint64_t pv = ~0ULL;
        uint64_t mv = 0;
        int score = length;
        for (int j = 0; j < n; ++j) {
            const uint64_t eq = peq[static_cast<unsigned char>(text[j])];
            const uint64_t xv = eq | mv;
            const uint64_t xh = (((eq & pv) + pv) ^ pv) | eq;
            uint64_t ph = mv | ~(xh | pv);
            uint64_t mh = pv & xh;
            if (ph & last) {
                ++score;
            } else if (mh & last) {
                --score;
            }
            ph = (ph << 1) | 1ULL;
            mh <<= 1;
            pv = (mh | ~(xv | ph)) & mask;
            mv = (ph & xv) & mask;
        }
        return score;
    }
};
} // namespace barcode_edit

#endif
