#include "trim.h"
#include <cstring>
#include <algorithm>
#include <climits>
#include <cstdlib>
// #region agent log
#include <fstream>
#include <string>
static const char* DEBUG_LOG_PATH = "/mnt/pikachu/STAR-Flex/.cursor/debug.log";
static bool debug_log_enabled() {
    static int enabled = -1;
    if (enabled != -1) {
        return enabled == 1;
    }
    const char* v = getenv("STAR_TRIM_DEBUG_LOG");
    if (!v || v[0] == '\0' || v[0] == '0' || v[0] == 'N' || v[0] == 'n') {
        enabled = 0;
    } else {
        enabled = 1;
    }
    return enabled == 1;
}
static void debug_log(const char* hyp, const char* loc, const char* msg, const std::string& data) {
    if (!debug_log_enabled()) {
        return;
    }
    std::ofstream f(DEBUG_LOG_PATH, std::ios::app);
    if (f) f << "{\"hypothesisId\":\"" << hyp << "\",\"location\":\"" << loc << "\",\"message\":\"" << msg << "\",\"data\":" << data << "}\n";
}
// #endregion

// Forward declarations for internal functions (defined in quality_trim.cpp and adapter_trim.cpp)
uint32_t quality_trim_3p(const char* qual, uint32_t len, uint8_t cutoff);
uint32_t quality_trim_5p(const char* qual, uint32_t len, uint8_t cutoff);
uint32_t find_adapter_3p(const char* seq, uint32_t seq_len,
                        const char* adapter, uint32_t adapter_len,
                        uint32_t min_overlap, double max_error_rate);
uint32_t find_adapter_3p_compat3(const char* seq, uint32_t seq_len,
                                 const char* adapter, uint32_t adapter_len,
                                 uint32_t min_overlap, double max_error_rate);

// Initialize params with Trim Galore defaults
void trim_params_init(struct TrimParams* params) {
    params->adapter_r1 = TRUSEQ_ADAPTER_R1;
    params->adapter_r2 = TRUSEQ_ADAPTER_R2;
    params->quality_cutoff = 20;
    params->min_length = 20;
    params->min_overlap = 1;  // Trim Galore default (--stringency 1)
    params->max_error_rate = 0.1;
    params->trim_5p_quality = false;
    params->compat_mode = TRIM_COMPAT_OFF;  // Default: cutadapt 5.1 parity
}

// Single read trimming (in-place modification of seq/qual buffers)
// seq and qual are fixed-size buffers, len is updated via pointer
struct TrimResult trim_read(char* seq, char* qual, uint32_t len, 
                            const char* adapter, const struct TrimParams* params) {
    struct TrimResult result = {0, 0, 0, 0, false};
    result.new_length = len;
    uint32_t orig_len = len;
    
    if (len == 0) {
        result.dropped = true;
        return result;
    }
    
    // Step 1: Quality trim 3' (Mott algorithm)
    uint32_t qual_trim_pos_3p = quality_trim_3p(qual, len, params->quality_cutoff);
    if (qual_trim_pos_3p < len) {
        result.qual_trimmed_3p = len - qual_trim_pos_3p;
        len = qual_trim_pos_3p;
        // Buffers are truncated by updating len - caller will use new length
    }
    uint32_t len_after_qual = len;
    
    // Step 2: Quality trim 5' (if enabled - disabled by default)
    uint32_t qual_trim_pos_5p = 0;
    if (params->trim_5p_quality) {
        qual_trim_pos_5p = quality_trim_5p(qual, len, params->quality_cutoff);
        if (qual_trim_pos_5p > 0) {
            result.qual_trimmed_5p = qual_trim_pos_5p;
            uint32_t bases_to_remove = qual_trim_pos_5p;
            // Shift seq and qual left
            memmove(seq, seq + bases_to_remove, len - bases_to_remove);
            memmove(qual, qual + bases_to_remove, len - bases_to_remove);
            len -= bases_to_remove;
        }
    }
    
    // Step 3: Adapter search on quality-trimmed sequence
    uint32_t adapter_pos = len;
    if (adapter != nullptr && strlen(adapter) > 0) {
        uint32_t adapter_len = strlen(adapter);
        // Switch adapter matching algorithm based on compatibility mode
        if (params->compat_mode == TRIM_COMPAT_CUTADAPT3) {
            adapter_pos = find_adapter_3p_compat3(seq, len, adapter, adapter_len,
                                                 params->min_overlap, params->max_error_rate);
        } else {
            adapter_pos = find_adapter_3p(seq, len, adapter, adapter_len,
                                          params->min_overlap, params->max_error_rate);
        }
        if (adapter_pos < len) {
            result.adapter_trimmed = len - adapter_pos;
            len = adapter_pos;
            // Buffers are truncated by updating len
        }
    }
    
    // Step 4: Check min_length -> set dropped flag
    result.new_length = len;
    if (len < params->min_length) {
        result.dropped = true;
    }
    
    // #region agent log
    // Log when any trimming occurred (H1: quality, H2/H3/H4: adapter)
    if (debug_log_enabled() && (result.qual_trimmed_3p > 0 || result.adapter_trimmed > 0)) {
        std::string seq_prefix(seq, std::min(len, 30u));
        std::string last10 = (orig_len > 10) ? std::string(seq + orig_len - 10, 10) : std::string(seq, orig_len);
        debug_log("H1-H4", "trim.cpp:trim_read", "trim_applied", 
            "{\"orig_len\":" + std::to_string(orig_len) + 
            ",\"len_after_qual\":" + std::to_string(len_after_qual) +
            ",\"qual_trimmed\":" + std::to_string(result.qual_trimmed_3p) +
            ",\"adapter_pos\":" + std::to_string(adapter_pos) +
            ",\"adapter_trimmed\":" + std::to_string(result.adapter_trimmed) +
            ",\"final_len\":" + std::to_string(len) +
            ",\"seq_prefix\":\"" + seq_prefix + "\"" +
            ",\"last10_orig\":\"" + last10 + "\"}");
    }
    // #endregion
    
    return result;
}

// Paired-end orchestrator
// Modifies seq1/qual1/len1 and seq2/qual2/len2 in place
void trim_pair(char* seq1, char* qual1, uint32_t* len1,
               char* seq2, char* qual2, uint32_t* len2,
               const struct TrimParams* params,
               struct TrimResult* result1, struct TrimResult* result2) {
    // Trim each mate independently
    *result1 = trim_read(seq1, qual1, *len1, params->adapter_r1, params);
    *result2 = trim_read(seq2, qual2, *len2, params->adapter_r2, params);
    
    // Update lengths (trim_read modifies buffers in place, but we update len pointers)
    *len1 = result1->new_length;
    *len2 = result2->new_length;
    
    // If either mate is dropped, mark both as dropped (Trim Galore behavior)
    if (result1->dropped || result2->dropped) {
        result1->dropped = true;
        result2->dropped = true;
    }
}

// Aggregate stats helper - centralizes stats accumulation logic
// Used by STAR integration and available for standalone tools
void trim_stats_add(struct TrimStats* total, const struct TrimResult* result) {
    total->reads_processed++;
    if (result->dropped) {
        total->reads_too_short++;
    }
    if (result->qual_trimmed_3p > 0 || result->qual_trimmed_5p > 0 || result->adapter_trimmed > 0) {
        total->reads_trimmed++;
    }
    total->bases_quality_trimmed += result->qual_trimmed_3p + result->qual_trimmed_5p;
    total->bases_adapter_trimmed += result->adapter_trimmed;
}
