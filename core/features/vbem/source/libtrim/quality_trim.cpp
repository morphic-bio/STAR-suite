#include "trim.h"
#include <cstring>
#include <algorithm>

// 3' quality trimming using cutadapt's modified Mott algorithm
// Returns the position to trim to (new length)
//
// Scoring follows cutadapt: score += (cutoff - q_phred), where low quality
// bases contribute positive values. We track the maximum score and the
// corresponding position; scanning stops when the cumulative score drops
// below zero (high-quality stretch reached).
uint32_t quality_trim_3p(const char* qual, uint32_t len, uint8_t cutoff) {
    if (len == 0) return 0;

    int64_t score = 0;
    int64_t max_score = 0;
    uint32_t cut_pos = len;  // Default: no trimming

    for (int32_t i = (int32_t)len - 1; i >= 0; i--) {
        int q = (int)(qual[i] - 33);
        score += (int)cutoff - q;  // positive when quality < cutoff

        if (score < 0) {
            // High-quality region reached; stop scanning
            break;
        }
        if (score > max_score) {
            max_score = score;
            cut_pos = (uint32_t)i;  // trim from here to end
        }
    }

    return cut_pos;
}

// 5' quality trimming (disabled by default, hook for future use)
// Mirrors 3' logic but scans 5' to 3'
uint32_t quality_trim_5p(const char* qual, uint32_t len, uint8_t cutoff) {
    if (len == 0) return 0;

    int64_t score = 0;
    int64_t max_score = 0;
    uint32_t cut_pos = 0;  // Default: no trimming

    for (uint32_t i = 0; i < len; i++) {
        int q = (int)(qual[i] - 33);
        score += (int)cutoff - q;  // positive when quality < cutoff

        if (score < 0) {
            // High-quality region reached; stop scanning
            break;
        }
        if (score > max_score) {
            max_score = score;
            cut_pos = i + 1;  // Trim up to (not including) this position
        }
    }
    return cut_pos;
}
