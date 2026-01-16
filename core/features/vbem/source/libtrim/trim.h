#ifndef LIBTRIM_TRIM_H
#define LIBTRIM_TRIM_H

#include <stdint.h>
#include <stdbool.h>

// Default TruSeq adapter sequences (Trim Galore defaults)
#define TRUSEQ_ADAPTER_R1 "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
#define TRUSEQ_ADAPTER_R2 "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

// Compatibility mode enum for adapter matching
enum TrimCompatMode {
    TRIM_COMPAT_OFF = 0,      // Default: cutadapt 5.1 parity
    TRIM_COMPAT_CUTADAPT3 = 1 // cutadapt 3.x compatibility mode
};

struct TrimParams {
    const char* adapter_r1;    // Default: TRUSEQ_ADAPTER_R1
    const char* adapter_r2;    // Default: TRUSEQ_ADAPTER_R2
    uint8_t quality_cutoff;    // Default: 20
    uint32_t min_length;       // Default: 20
    uint32_t min_overlap;      // Default: 1 (Trim Galore default --stringency 1)
    double max_error_rate;     // Default: 0.1 (mismatches = floor(overlap * 0.1))
    bool trim_5p_quality;      // Default: false (match Trim Galore)
    enum TrimCompatMode compat_mode; // Default: TRIM_COMPAT_OFF (cutadapt 5.1)
    // Future hooks (currently unused):
    // bool trim_5p_adapter;   // Reserved for 5' adapter trimming
    // const char** extra_adapters; // Reserved for multiple adapter support
};

struct TrimResult {
    uint32_t new_length;
    uint32_t qual_trimmed_3p;  // bases trimmed by quality (3')
    uint32_t qual_trimmed_5p;  // bases trimmed by quality (5', if enabled)
    uint32_t adapter_trimmed;  // bases trimmed by adapter
    bool dropped;              // true if new_length < min_length
};

struct TrimStats {
    uint64_t reads_processed;
    uint64_t reads_trimmed;     // reads with any trimming
    uint64_t reads_too_short;   // reads dropped due to min_length
    uint64_t bases_quality_trimmed;
    uint64_t bases_adapter_trimmed;
};

// Initialize params with Trim Galore defaults
void trim_params_init(struct TrimParams* params);

// Single read trimming (in-place modification of seq/qual)
// Returns result with counts of trimmed bases
// Pure function: no globals, thread-safe
struct TrimResult trim_read(char* seq, char* qual, uint32_t len, 
                            const char* adapter, const struct TrimParams* params);

// Paired-end orchestrator
// Applies trim_read to both mates using adapter_r1/r2
// Sets result.dropped=true for BOTH mates if either falls below min_length
void trim_pair(char* seq1, char* qual1, uint32_t* len1,
               char* seq2, char* qual2, uint32_t* len2,
               const struct TrimParams* params,
               struct TrimResult* result1, struct TrimResult* result2);

// Aggregate stats helper (thread-safe if caller synchronizes)
void trim_stats_add(struct TrimStats* total, const struct TrimResult* result);

#endif // LIBTRIM_TRIM_H
