#ifndef ALIGNMENT_FILTER_H
#define ALIGNMENT_FILTER_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <limits>
#include <unordered_map>

// Mate status enum matching Salmon's pufferfish::util::MateStatus
enum class MateStatus : uint8_t {
    SINGLE_END = 0,
    PAIRED_END_LEFT = 1,
    PAIRED_END_RIGHT = 2,
    PAIRED_END_PAIRED = 3
};

// Raw alignment from aligner (mirrors Salmon's QuasiAlignment)
struct RawAlignment {
    uint32_t transcript_id;
    int32_t pos;
    int32_t score;              // Alignment score (SW or coverage)
    double est_aln_prob;        // Computed: exp(-scoreExp * (bestScore - score))
    double log_frag_prob;       // Log fragment length probability
    double log_compat_prob;     // Log alignment compatibility probability
    double err_like;            // Pre-computed error model likelihood (CIGAR-based)
    bool has_err_like;          // True if err_like was computed from CIGAR
    MateStatus mate_status;
    bool is_decoy;
    bool is_forward;            // True if read aligns to forward strand (not reverse complemented)
    
    // For paired-end
    int32_t mate_score;
    int32_t fragment_len;
    bool mate_is_forward;       // For paired-end pairs: orientation of mate (needed for observed format)
    int32_t mate_pos;           // For paired-end pairs: position of mate (needed for observed format)
    bool mate_fields_set;       // true when mate_pos/mate_is_forward are populated
    bool is_primary;            // true if alignment is primary
    
    RawAlignment()
        : transcript_id(0), pos(0), score(0), est_aln_prob(0.0),
          log_frag_prob(0.0), log_compat_prob(0.0),
          err_like(0.0), has_err_like(false),
          mate_status(MateStatus::SINGLE_END), is_decoy(false),
          is_forward(true), mate_score(0), fragment_len(0),
          mate_is_forward(false), mate_pos(-1), mate_fields_set(false),
          is_primary(true) {}
};

// Score tracking per read (mirrors Salmon's MappingScoreInfo)
struct MappingScoreInfo {
    int32_t best_score;
    int32_t second_best_score;
    int32_t best_decoy_score;
    double decoy_thresh;
    
    // Per-transcript best score tracking: (score, index)
    std::unordered_map<uint32_t, std::pair<int32_t, size_t>> best_score_per_transcript;
    
    MappingScoreInfo() 
        : best_score(std::numeric_limits<int32_t>::min()),
          second_best_score(std::numeric_limits<int32_t>::min()),
          best_decoy_score(std::numeric_limits<int32_t>::min()),
          decoy_thresh(1.0) {}
    
    MappingScoreInfo(double decoy_thresh_in)
        : best_score(std::numeric_limits<int32_t>::min()),
          second_best_score(std::numeric_limits<int32_t>::min()),
          best_decoy_score(std::numeric_limits<int32_t>::min()),
          decoy_thresh(decoy_thresh_in) {}
    
    void clear(size_t num_hits) {
        best_score = std::numeric_limits<int32_t>::min();
        second_best_score = std::numeric_limits<int32_t>::min();
        best_decoy_score = std::numeric_limits<int32_t>::min();
        best_score_per_transcript.clear();
    }
};

// Filter configuration matching Salmon's defaults
struct FilterParams {
    // Salmon defaults from SalmonDefaults.hpp
    double score_exp = 1.0;           // scoreExp: controls probability decay
    double min_aln_prob = 1e-5;      // minAlnProb: skip below this threshold
    double decoy_threshold = 1.0;    // decoyThreshold
    bool hard_filter = false;         // hardFilter: if true, keep only best score
    uint32_t max_read_occs = 200;    // maxReadOccs: discard if > this many hits
    double min_score_fraction = 0.0;  // minScoreFraction: 0 = disabled (default for STAR inline mode)
                                      // Set to 0.65 for Salmon-parity mode (drops alignments < 65% of best score)
};

// Update per-transcript best scores (called per alignment)
// Matches Salmon's updateRefMappings logic
void updateRefMappings(
    uint32_t tid, 
    int32_t score, 
    bool is_compat, 
    size_t idx,
    const std::vector<bool>& is_decoy,
    MappingScoreInfo& msi,
    std::vector<int32_t>& scores
);

// Filter and collect alignments matching Salmon's filterAndCollectAlignments
// Matches SalmonMappingUtils.hpp lines 271-390
std::vector<RawAlignment> filterAndCollectAlignments(
    std::vector<RawAlignment>& alignments,
    const FilterParams& params,
    MappingScoreInfo& msi
);

#endif // ALIGNMENT_FILTER_H
