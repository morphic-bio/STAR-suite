// Native Visium HD R1 anchored+indel decoder.
//
// Production-speed C++ port of scripts/probe_hd_r1_anchored_full_decode.py::decode_record.
// Reads R1 FASTQ(.gz), scans the chemistry full-barcode start offset window, anchors BC1/BC2
// against the full slide oligos (variable length 15/16 and 14/15) with Hamming + deletion-signature
// (indel) candidate generation, scores candidates, and emits decode_reads.tsv (read_id ->
// universe_unique_assigned, row2, col2) for the existing R1-rewrite preprocessor.
//
// Parity reference: the Python decoder is the oracle. This tool must reproduce its per-read
// (universe_unique_assigned, row2, col2) assignment on frozen fixtures. col2 == BC1/X, row2 == BC2/Y.
//
// This file vendors the validated companion decoder implementation at
// 5e9af58f9b86e9d3b95f8612a4a77e40f9a3ed86. The STAR-owned library wrapper
// is below; define STAR_SPATIAL_R1_DECODER_STANDALONE to retain the original
// command-line tool for parity work.

#include "SpatialR1Decoder.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <sys/stat.h>
#include <htslib/khash.h>
#include <zlib.h>

namespace {

// ---- distance helpers (match Python edit_distance exactly for accepted candidates) ----
constexpr int kMaxLen = 64;
int edit_distance_myers64(const char* pattern, int m, const char* text, int n) {
    if (m == 0) return n;
    if (n == 0) return m;
    uint64_t peq[256] = {0};
    for (int i = 0; i < m; ++i) {
        peq[static_cast<unsigned char>(pattern[i])] |= (1ULL << i);
    }

    const uint64_t mask = (m == 64) ? ~0ULL : ((1ULL << m) - 1ULL);
    const uint64_t last = 1ULL << (m - 1);
    uint64_t pv = ~0ULL;
    uint64_t mv = 0;
    int score = m;

    for (int j = 0; j < n; ++j) {
        uint64_t eq = peq[static_cast<unsigned char>(text[j])];
        uint64_t xv = eq | mv;
        uint64_t xh = (((eq & pv) + pv) ^ pv) | eq;
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

int edit_distance_bounded(const char* a, int la, const char* b, int lb, int maxk) {
    if (la > kMaxLen) la = kMaxLen;
    if (lb > kMaxLen) lb = kMaxLen;
    if (la - lb > maxk || lb - la > maxk) return maxk + 1;
    int d = edit_distance_myers64(b, lb, a, la);
    return d <= maxk ? d : maxk + 1;
}

struct TerminalNormalizedEdit {
    int edit_distance;
    int terminal_delta;

    TerminalNormalizedEdit(
        int editDistance = std::numeric_limits<int>::max() / 4,
        int terminalDelta = 0)
        : edit_distance(editDistance), terminal_delta(terminalDelta) {}
};

TerminalNormalizedEdit terminal_normalized_edit_distance(const char* observed, int observed_len,
                                                         const char* target, int target_len,
                                                         int maxk) {
    const int raw = edit_distance_bounded(observed, observed_len, target, target_len, maxk);
    TerminalNormalizedEdit best{raw, std::abs(observed_len - target_len)};
    if (best.edit_distance == 0) return best;

    for (int obs_left = 0; obs_left <= maxk; ++obs_left) {
        for (int obs_right = 0; obs_left + obs_right <= maxk; ++obs_right) {
            const int after_obs_trim = observed_len - obs_left - obs_right;
            if (after_obs_trim <= 0) continue;
            for (int target_left = 0; obs_left + obs_right + target_left <= maxk;
                 ++target_left) {
                for (int target_right = 0;
                     obs_left + obs_right + target_left + target_right <= maxk;
                     ++target_right) {
                    const int terminal_delta =
                        obs_left + obs_right + target_left + target_right;
                    if (terminal_delta == 0) continue;
                    if ((obs_left + obs_right) > 0 &&
                        (target_left + target_right) > 0) {
                        continue;
                    }
                    const int after_target_trim = target_len - target_left - target_right;
                    if (after_target_trim <= 0) continue;
                    if (after_obs_trim != after_target_trim) continue;
                    if (std::memcmp(observed + obs_left, target + target_left,
                                    static_cast<size_t>(after_obs_trim)) != 0) {
                        continue;
                    }
                    if (std::make_pair(0, terminal_delta) <
                            std::make_pair(best.edit_distance, best.terminal_delta)) {
                        best.edit_distance = 0;
                        best.terminal_delta = terminal_delta;
                        if (best.terminal_delta == 1) return best;
                    }
                }
            }
        }
    }
    return best;
}

const char BASES[4] = {'A', 'C', 'G', 'T'};

int base_bits(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

bool pack_sequence(const char* seq, int len, uint64_t& packed) {
    if (len < 0 || len > 28) return false;  // 56 payload bits; high 8 bits store length.
    uint64_t v = 0;
    for (int i = 0; i < len; ++i) {
        int b = base_bits(seq[i]);
        if (b < 0) return false;
        v = (v << 2) | static_cast<uint64_t>(b);
    }
    packed = v;
    return true;
}

uint64_t make_packed_key(int len, uint64_t packed) {
    return (static_cast<uint64_t>(len) << 56) | packed;
}

int packed_key_len(uint64_t key) {
    return static_cast<int>(key >> 56);
}

uint64_t packed_key_payload(uint64_t key) {
    return key & 0x00FFFFFFFFFFFFFFULL;
}

int packed_hamming(uint64_t a, uint64_t b, int len) {
    uint64_t x = a ^ b;
    int d = 0;
    for (int i = 0; i < len; ++i) {
        d += ((x & 0x3ULL) != 0);
        x >>= 2;
    }
    return d;
}

void unpack_payload(uint64_t packed, int len, char* out) {
    for (int i = len - 1; i >= 0; --i) {
        out[i] = BASES[packed & 0x3ULL];
        packed >>= 2;
    }
}

// All sequences within Hamming distance <= max_d of `seq` (includes seq).
void hamming_variants(const std::string& seq, int max_d, std::vector<std::string>& out) {
    out.clear();
    out.push_back(seq);
    if (max_d >= 1) {
        for (size_t i = 0; i < seq.size(); ++i) {
            std::string v = seq;
            for (char b : BASES) {
                if (b == seq[i]) continue;
                v[i] = b;
                out.push_back(v);
            }
        }
    }
    if (max_d >= 2) {
        for (size_t i = 0; i < seq.size(); ++i) {
            for (size_t j = i + 1; j < seq.size(); ++j) {
                std::string v = seq;
                for (char bi : BASES) {
                    if (bi == seq[i]) continue;
                    v[i] = bi;
                    for (char bj : BASES) {
                        if (bj == seq[j]) continue;
                        v[j] = bj;
                        out.push_back(v);
                    }
                }
                // restore handled by reassignment v=seq next iter
            }
        }
    }
}

// All sequences within edit distance <= 1 of `seq` (includes seq), deduped.
void edit1_variants(const std::string& seq, std::vector<std::string>& out) {
    out.clear();
    out.push_back(seq);
    for (size_t i = 0; i < seq.size(); ++i) {
        std::string v = seq;
        for (char b : BASES) {
            if (b == seq[i]) continue;
            v[i] = b;
            out.push_back(v);
        }
    }
    for (size_t i = 0; i < seq.size(); ++i) {
        std::string v;
        v.reserve(seq.size() - 1);
        v.append(seq.data(), i);
        v.append(seq.data() + i + 1, seq.size() - i - 1);
        out.push_back(v);
    }
    for (size_t i = 0; i <= seq.size(); ++i) {
        for (char b : BASES) {
            std::string v;
            v.reserve(seq.size() + 1);
            v.append(seq.data(), i);
            v.push_back(b);
            v.append(seq.data() + i, seq.size() - i);
            out.push_back(v);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

// All deletion signatures (original + up to max_del positions removed), deduped.
void deletion_signatures(const std::string& seq, int max_del, std::vector<std::string>& out) {
    out.clear();
    out.push_back(seq);
    int md = std::min<int>(max_del, static_cast<int>(seq.size()));
    // combinations of positions to delete
    for (int k = 1; k <= md; ++k) {
        std::vector<int> idx(k);
        for (int i = 0; i < k; ++i) idx[i] = i;
        while (true) {
            std::string s;
            s.reserve(seq.size() - k);
            int p = 0;
            for (size_t i = 0; i < seq.size(); ++i) {
                if (p < k && idx[p] == static_cast<int>(i)) { ++p; continue; }
                s.push_back(seq[i]);
            }
            out.push_back(s);
            // next combination
            int t = k - 1;
            while (t >= 0 && idx[t] == static_cast<int>(seq.size()) - k + t) --t;
            if (t < 0) break;
            ++idx[t];
            for (int i = t + 1; i < k; ++i) idx[i] = idx[i - 1] + 1;
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

bool packed_deletion_signature_keys(const char* seq, int len, int max_del, std::vector<uint64_t>& out) {
    out.clear();
    uint64_t packed = 0;
    if (!pack_sequence(seq, len, packed)) return false;
    out.push_back(make_packed_key(len, packed));

    int md = std::min(max_del, len);
    for (int k = 1; k <= md; ++k) {
        std::vector<int> idx(k);
        for (int i = 0; i < k; ++i) idx[i] = i;
        while (true) {
            uint64_t sig = 0;
            int p = 0;
            for (int i = 0; i < len; ++i) {
                if (p < k && idx[p] == i) {
                    ++p;
                    continue;
                }
                int b = base_bits(seq[i]);
                if (b < 0) return false;
                sig = (sig << 2) | static_cast<uint64_t>(b);
            }
            out.push_back(make_packed_key(len - k, sig));

            int t = k - 1;
            while (t >= 0 && idx[t] == len - k + t) --t;
            if (t < 0) break;
            ++idx[t];
            for (int i = t + 1; i < k; ++i) idx[i] = idx[i - 1] + 1;
        }
    }

    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return true;
}

struct PackedHit {
    uint16_t idx;
    uint8_t distance;
    uint8_t reserved;

    PackedHit(uint16_t index = 0, uint8_t editDistance = 0,
              uint8_t reservedValue = 0)
        : idx(index), distance(editDistance), reserved(reservedValue) {}
};

struct HitSpan {
    const PackedHit* data;
    uint16_t count;
    bool found;

    HitSpan(const PackedHit* hitData = nullptr, uint16_t hitCount = 0,
            bool wasFound = false)
        : data(hitData), count(hitCount), found(wasFound) {}
};

struct PackedEntry {
    uint64_t key;
    uint32_t offset;
    uint16_t count;
    uint16_t reserved;

    PackedEntry(uint64_t entryKey = 0, uint32_t entryOffset = 0,
                uint16_t entryCount = 0, uint16_t reservedValue = 0)
        : key(entryKey), offset(entryOffset), count(entryCount),
          reserved(reservedValue) {}
};

struct MemoValue {
    uint32_t offset;
    uint16_t count;

    MemoValue(uint32_t valueOffset = 0, uint16_t valueCount = 0)
        : offset(valueOffset), count(valueCount) {}
};

KHASH_MAP_INIT_INT64(memo64, MemoValue)

struct PackedCandidate {
    uint64_t key;
    uint16_t idx;

    PackedCandidate(uint64_t candidateKey = 0, uint16_t index = 0)
        : key(candidateKey), idx(index) {}
};

struct PackedScoredCandidate {
    uint64_t key;
    uint16_t idx;
    uint8_t distance;
    uint8_t raw_distance;

    PackedScoredCandidate(uint64_t candidateKey = 0, uint16_t index = 0,
                          uint8_t editDistance = 0, uint8_t rawDistance = 0)
        : key(candidateKey), idx(index), distance(editDistance),
          raw_distance(rawDistance) {}
};

struct PackedSizedCandidate {
    uint64_t key;
    uint16_t idx;
    uint8_t distance;
    uint8_t observed_len;

    PackedSizedCandidate(uint64_t candidateKey = 0, uint16_t index = 0,
                         uint8_t editDistance = 0, uint8_t observedLength = 0)
        : key(candidateKey), idx(index), distance(editDistance),
          observed_len(observedLength) {}
};

struct LocalEditVariant {
    uint64_t key;
    uint8_t distance;

    LocalEditVariant(uint64_t variantKey = 0, uint8_t editDistance = 0)
        : key(variantKey), distance(editDistance) {}
};

enum class TargetSliceMode { Full, Prefix, Suffix };

struct DecodeStats {
    uint64_t anchor_lookup_calls = 0;
    uint64_t anchor_hamming_keys_found = 0;
    uint64_t anchor_hamming_hits_returned = 0;
    uint64_t anchor_deletion_dp_calls = 0;
    uint64_t anchor_deletion_dp_accepts = 0;
    uint64_t full_cache_lookup_calls = 0;
    uint64_t full_cache_keys_found = 0;
    uint64_t full_cache_target_hits = 0;
    uint64_t full_same_length_target_miss_dp = 0;
    uint64_t full_length_shift_dp = 0;
    uint64_t full_dp_calls = 0;
    uint64_t full_dp_accepts = 0;

    void add(const DecodeStats& other) {
        anchor_lookup_calls += other.anchor_lookup_calls;
        anchor_hamming_keys_found += other.anchor_hamming_keys_found;
        anchor_hamming_hits_returned += other.anchor_hamming_hits_returned;
        anchor_deletion_dp_calls += other.anchor_deletion_dp_calls;
        anchor_deletion_dp_accepts += other.anchor_deletion_dp_accepts;
        full_cache_lookup_calls += other.full_cache_lookup_calls;
        full_cache_keys_found += other.full_cache_keys_found;
        full_cache_target_hits += other.full_cache_target_hits;
        full_same_length_target_miss_dp += other.full_same_length_target_miss_dp;
        full_length_shift_dp += other.full_length_shift_dp;
        full_dp_calls += other.full_dp_calls;
        full_dp_accepts += other.full_dp_accepts;
    }
};

struct ObservedCacheStats {
    uint64_t lookup_calls = 0;
    uint64_t memo_hits = 0;
    uint64_t memo_misses = 0;
    uint64_t e01_keys_found = 0;
    uint64_t e01_hits_returned = 0;
    uint64_t h2_keys_found = 0;
    uint64_t h2_hits_returned = 0;
    uint64_t deletion_signature_candidates = 0;
    uint64_t deletion_candidates_after_dedupe = 0;
    uint64_t skipped_already_resolved = 0;
    uint64_t dp_calls = 0;
    uint64_t dp_accepts = 0;
    uint64_t non_acgt_keys = 0;
    uint64_t non_acgt_dp_calls = 0;
    uint64_t answer_hits_stored = 0;

    void add(const ObservedCacheStats& other) {
        lookup_calls += other.lookup_calls;
        memo_hits += other.memo_hits;
        memo_misses += other.memo_misses;
        e01_keys_found += other.e01_keys_found;
        e01_hits_returned += other.e01_hits_returned;
        h2_keys_found += other.h2_keys_found;
        h2_hits_returned += other.h2_hits_returned;
        deletion_signature_candidates += other.deletion_signature_candidates;
        deletion_candidates_after_dedupe += other.deletion_candidates_after_dedupe;
        skipped_already_resolved += other.skipped_already_resolved;
        dp_calls += other.dp_calls;
        dp_accepts += other.dp_accepts;
        non_acgt_keys += other.non_acgt_keys;
        non_acgt_dp_calls += other.non_acgt_dp_calls;
        answer_hits_stored += other.answer_hits_stored;
    }
};

struct PackedCacheHeader {
    char magic[16];
    uint32_t version = 1;
    uint32_t score_mode = 0;
    uint32_t max_distance = 0;
    uint32_t target_count = 0;
    uint64_t target_checksum = 0;
    uint64_t entry_count = 0;
    uint64_t hit_count = 0;
};

static_assert(sizeof(PackedHit) == 4, "PackedHit must stay serialized as 4 bytes");
static_assert(sizeof(PackedEntry) == 16, "PackedEntry must stay serialized as 16 bytes");

uint64_t fnv1a_update(uint64_t h, const void* data, size_t n) {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < n; ++i) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

uint64_t checksum_targets(const std::vector<std::string>& targets) {
    uint64_t h = 1469598103934665603ULL;
    for (auto& target : targets) {
        uint64_t len = static_cast<uint64_t>(target.size());
        h = fnv1a_update(h, &len, sizeof(len));
        h = fnv1a_update(h, target.data(), target.size());
    }
    return h;
}

bool ensure_directory(const std::string& path) {
    if (path.empty()) return true;
    if (::mkdir(path.c_str(), 0775) == 0) return true;
    if (errno == EEXIST) return true;
    std::cerr << "ERROR: cannot create cache dir " << path << ": " << std::strerror(errno) << "\n";
    return false;
}

std::string cache_path_for(const std::string& cache_dir, const std::string& label) {
    if (cache_dir.empty()) return std::string();
    return cache_dir + "/" + label + ".hd_r1_cache.bin";
}

double pct(uint64_t num, uint64_t den) {
    return den ? (100.0 * static_cast<double>(num) / static_cast<double>(den)) : 0.0;
}

enum class CandidateScoreMode { LegacyEdit, AffineGap };
enum class DirectBestResolutionMode { StrictBoundary, OverlapSpanDp };

struct AffineScoringConfig {
    int substitution_cost = 5;
    int gap_open_cost = 2;
    int gap_extend_cost = 1;
    int n_cost = 1;
};

struct AlignmentMetrics {
    int cost = -1;
    int substitutions = 0;
    int insertions = 0;
    int deletions = 0;
    int n_edits = 0;
};

struct DpState {
    int cost = std::numeric_limits<int>::max() / 4;
    int substitutions = 0;
    int insertions = 0;
    int deletions = 0;
    int n_edits = 0;
};

bool dp_state_valid(const DpState& s) {
    return s.cost < std::numeric_limits<int>::max() / 8;
}

int dp_state_ops(const DpState& s) {
    return s.substitutions + s.insertions + s.deletions + s.n_edits;
}

bool dp_state_less(const DpState& a, const DpState& b) {
    if (!dp_state_valid(a)) return false;
    if (!dp_state_valid(b)) return true;
    if (a.cost != b.cost) return a.cost < b.cost;
    if (dp_state_ops(a) != dp_state_ops(b)) return dp_state_ops(a) < dp_state_ops(b);
    if (a.substitutions != b.substitutions) return a.substitutions < b.substitutions;
    const int a_gaps = a.insertions + a.deletions;
    const int b_gaps = b.insertions + b.deletions;
    if (a_gaps != b_gaps) return a_gaps < b_gaps;
    if (a.n_edits != b.n_edits) return a.n_edits < b.n_edits;
    if (a.insertions != b.insertions) return a.insertions < b.insertions;
    return a.deletions < b.deletions;
}

DpState best_dp_state(const DpState& a, const DpState& b) {
    return dp_state_less(b, a) ? b : a;
}

DpState add_pair_score(DpState s, char observed, char target,
                       const AffineScoringConfig& cfg) {
    if (!dp_state_valid(s)) return s;
    const bool observed_acgt = base_bits(observed) >= 0;
    const bool target_acgt = base_bits(target) >= 0;
    if (!observed_acgt || !target_acgt) {
        s.cost += cfg.n_cost;
        ++s.n_edits;
    } else if (observed != target) {
        s.cost += cfg.substitution_cost;
        ++s.substitutions;
    }
    return s;
}

DpState add_insertion_score(DpState s, int gap_cost) {
    if (!dp_state_valid(s)) return s;
    s.cost += gap_cost;
    ++s.insertions;
    return s;
}

DpState add_deletion_score(DpState s, int gap_cost) {
    if (!dp_state_valid(s)) return s;
    s.cost += gap_cost;
    ++s.deletions;
    return s;
}

AlignmentMetrics affine_alignment_metrics(const char* observed, int observed_len,
                                          const char* target, int target_len,
                                          const AffineScoringConfig& cfg) {
    AlignmentMetrics out;
    if (observed_len < 0 || target_len < 0 ||
        observed_len > kMaxLen || target_len > kMaxLen) {
        out.cost = std::numeric_limits<int>::max() / 4;
        return out;
    }

    DpState invalid;
    DpState prev_m[kMaxLen + 1], prev_i[kMaxLen + 1], prev_d[kMaxLen + 1];
    DpState cur_m[kMaxLen + 1], cur_i[kMaxLen + 1], cur_d[kMaxLen + 1];
    for (int j = 0; j <= target_len; ++j) {
        prev_m[j] = invalid;
        prev_i[j] = invalid;
        prev_d[j] = invalid;
    }
    prev_m[0].cost = 0;
    for (int j = 1; j <= target_len; ++j) {
        prev_d[j] = add_deletion_score(prev_d[j - 1],
                                       j == 1 ? cfg.gap_open_cost : cfg.gap_extend_cost);
        if (j == 1) prev_d[j] = add_deletion_score(prev_m[j - 1], cfg.gap_open_cost);
    }

    for (int i = 1; i <= observed_len; ++i) {
        for (int j = 0; j <= target_len; ++j) {
            cur_m[j] = invalid;
            cur_i[j] = invalid;
            cur_d[j] = invalid;
        }

        cur_i[0] = add_insertion_score(prev_i[0],
                                       i == 1 ? cfg.gap_open_cost : cfg.gap_extend_cost);
        if (i == 1) cur_i[0] = add_insertion_score(prev_m[0], cfg.gap_open_cost);

        for (int j = 1; j <= target_len; ++j) {
            DpState from_diag = best_dp_state(prev_m[j - 1], prev_i[j - 1]);
            from_diag = best_dp_state(from_diag, prev_d[j - 1]);
            cur_m[j] = add_pair_score(from_diag, observed[i - 1], target[j - 1], cfg);

            DpState ins_from_m = add_insertion_score(prev_m[j], cfg.gap_open_cost);
            DpState ins_from_i = add_insertion_score(prev_i[j], cfg.gap_extend_cost);
            DpState ins_from_d = add_insertion_score(prev_d[j], cfg.gap_open_cost);
            cur_i[j] = best_dp_state(best_dp_state(ins_from_m, ins_from_i), ins_from_d);

            DpState del_from_m = add_deletion_score(cur_m[j - 1], cfg.gap_open_cost);
            DpState del_from_d = add_deletion_score(cur_d[j - 1], cfg.gap_extend_cost);
            DpState del_from_i = add_deletion_score(cur_i[j - 1], cfg.gap_open_cost);
            cur_d[j] = best_dp_state(best_dp_state(del_from_m, del_from_d), del_from_i);
        }

        for (int j = 0; j <= target_len; ++j) {
            prev_m[j] = cur_m[j];
            prev_i[j] = cur_i[j];
            prev_d[j] = cur_d[j];
        }
    }

    DpState best = best_dp_state(best_dp_state(prev_m[target_len], prev_i[target_len]),
                                 prev_d[target_len]);
    out.cost = best.cost;
    out.substitutions = best.substitutions;
    out.insertions = best.insertions;
    out.deletions = best.deletions;
    out.n_edits = best.n_edits;
    return out;
}

const char* score_mode_name(CandidateScoreMode mode) {
    return mode == CandidateScoreMode::AffineGap ? "affine" : "legacy";
}

const char* direct_best_resolution_name(DirectBestResolutionMode mode) {
    return mode == DirectBestResolutionMode::OverlapSpanDp
        ? "overlap_span_dp"
        : "strict_boundary";
}

void print_decode_stats(const DecodeStats& stats) {
    std::cerr << "decode_stats:"
              << " anchor_lookup_calls=" << stats.anchor_lookup_calls
              << " anchor_hamming_keys_found=" << stats.anchor_hamming_keys_found
              << " (" << pct(stats.anchor_hamming_keys_found, stats.anchor_lookup_calls) << "%)"
              << " anchor_hamming_hits_returned=" << stats.anchor_hamming_hits_returned
              << " anchor_deletion_dp_calls=" << stats.anchor_deletion_dp_calls
              << " anchor_deletion_dp_accepts=" << stats.anchor_deletion_dp_accepts
              << " (" << pct(stats.anchor_deletion_dp_accepts, stats.anchor_deletion_dp_calls) << "%)"
              << "\n";
    std::cerr << "decode_stats:"
              << " full_cache_lookup_calls=" << stats.full_cache_lookup_calls
              << " full_cache_keys_found=" << stats.full_cache_keys_found
              << " (" << pct(stats.full_cache_keys_found, stats.full_cache_lookup_calls) << "%)"
              << " full_cache_target_hits=" << stats.full_cache_target_hits
              << " full_same_length_target_miss_dp=" << stats.full_same_length_target_miss_dp
              << " full_length_shift_dp=" << stats.full_length_shift_dp
              << " full_dp_calls=" << stats.full_dp_calls
              << " full_dp_accepts=" << stats.full_dp_accepts
              << " (" << pct(stats.full_dp_accepts, stats.full_dp_calls) << "%)"
              << "\n";
}

void maybe_add_edit_variant(const std::string& observed, int distance,
                            const bool* allowed_lengths,
                            std::vector<LocalEditVariant>& out) {
    const int len = static_cast<int>(observed.size());
    if (len < 0 || len > kMaxLen || !allowed_lengths[len]) return;
    uint64_t packed = 0;
    if (!pack_sequence(observed.data(), len, packed)) return;
    out.push_back({make_packed_key(len, packed), static_cast<uint8_t>(distance)});
}

void enumerate_edit_variants_rec(const std::string& target, int pos, int edits,
                                 int max_dist, const bool* allowed_lengths,
                                 int min_allowed_len, int max_allowed_len,
                                 std::string& observed,
                                 std::vector<LocalEditVariant>& out) {
    const int observed_len = static_cast<int>(observed.size());
    const int remaining = static_cast<int>(target.size()) - pos;
    const int edits_left = max_dist - edits;
    const int min_final_len = observed_len + std::max(0, remaining - edits_left);
    const int max_final_len = observed_len + remaining + edits_left;
    if (observed_len > max_allowed_len ||
        min_final_len > max_allowed_len ||
        max_final_len < min_allowed_len) {
        return;
    }

    if (pos == static_cast<int>(target.size())) {
        maybe_add_edit_variant(observed, edits, allowed_lengths, out);
        if (edits >= max_dist || observed_len >= max_allowed_len) return;
        for (char base : BASES) {
            observed.push_back(base);
            enumerate_edit_variants_rec(target, pos, edits + 1, max_dist,
                                        allowed_lengths, min_allowed_len,
                                        max_allowed_len, observed, out);
            observed.pop_back();
        }
        return;
    }

    const char target_base = target[pos];
    observed.push_back(target_base);
    enumerate_edit_variants_rec(target, pos + 1, edits, max_dist,
                                allowed_lengths, min_allowed_len,
                                max_allowed_len, observed, out);
    observed.pop_back();

    if (edits >= max_dist) return;

    for (char base : BASES) {
        if (base == target_base) continue;
        observed.push_back(base);
        enumerate_edit_variants_rec(target, pos + 1, edits + 1, max_dist,
                                    allowed_lengths, min_allowed_len,
                                    max_allowed_len, observed, out);
        observed.pop_back();
    }

    enumerate_edit_variants_rec(target, pos + 1, edits + 1, max_dist,
                                allowed_lengths, min_allowed_len,
                                max_allowed_len, observed, out);

    if (observed_len < max_allowed_len) {
        for (char base : BASES) {
            observed.push_back(base);
            enumerate_edit_variants_rec(target, pos, edits + 1, max_dist,
                                        allowed_lengths, min_allowed_len,
                                        max_allowed_len, observed, out);
            observed.pop_back();
        }
    }
}

class PackedAnswerCache {
  public:
    enum class ScoreMode { Hamming, Edit };

    PackedAnswerCache() : index_(kh_init(memo64)) {
        if (!index_) {
            std::cerr << "ERROR: failed to allocate packed answer khash index\n";
            std::exit(2);
        }
    }

    ~PackedAnswerCache() {
        kh_destroy(memo64, index_);
    }

    PackedAnswerCache(const PackedAnswerCache&) = delete;
    PackedAnswerCache& operator=(const PackedAnswerCache&) = delete;

    void build_hamming_universe(const std::vector<std::string>& targets, int max_dist,
                                ScoreMode score_mode, const char* label,
                                const std::string& cache_path = std::string()) {
        entries_.clear();
        hits_.clear();
        clear_index();
        max_distance_ = max_dist;
        target_count_ = targets.size();
        score_mode_ = score_mode;
        enabled_ = false;

        if (targets.empty()) return;
        if (targets.size() > static_cast<size_t>(std::numeric_limits<uint16_t>::max())) {
            std::cerr << "ERROR: packed HD cache target count exceeds uint16_t index range: "
                      << targets.size() << "\n";
            std::exit(2);
        }
        if (max_dist < 0 || max_dist > 2) {
            std::cerr << "ERROR: packed HD cache currently supports H0/H1/H2 only; requested H"
                      << max_dist << " for " << label << "\n";
            std::exit(2);
        }

        const uint64_t target_checksum = checksum_targets(targets);
        if (!cache_path.empty() &&
            load_cache(cache_path, targets.size(), target_checksum, max_dist, score_mode, label)) {
            enabled_ = true;
            return;
        }

        std::vector<uint64_t> target_packed(targets.size(), 0);
        for (size_t i = 0; i < targets.size(); ++i) {
            if (!pack_sequence(targets[i].data(), static_cast<int>(targets[i].size()), target_packed[i])) {
                std::cerr << "ERROR: non-ACGT target in packed HD cache " << label
                          << " at index " << i << "\n";
                std::exit(2);
            }
        }

        size_t max_len = 0;
        for (auto& t : targets) max_len = std::max(max_len, t.size());
        size_t variants_per_target = 1;
        if (max_dist >= 1) variants_per_target += 3 * max_len;
        if (max_dist >= 2) variants_per_target += 9 * max_len * (max_len - 1) / 2;

        std::vector<PackedCandidate> candidates;
        candidates.reserve(targets.size() * variants_per_target);
        std::vector<std::string> variants;

        for (size_t idx = 0; idx < targets.size(); ++idx) {
            hamming_variants(targets[idx], max_dist, variants);
            for (auto& observed : variants) {
                uint64_t packed = 0;
                if (!pack_sequence(observed.data(), static_cast<int>(observed.size()), packed)) continue;
                candidates.push_back({make_packed_key(static_cast<int>(observed.size()), packed),
                                      static_cast<uint16_t>(idx)});
            }
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const PackedCandidate& a, const PackedCandidate& b) {
                      if (a.key != b.key) return a.key < b.key;
                      return a.idx < b.idx;
                  });

        char observed_buf[kMaxLen + 1];
        size_t i = 0;
        while (i < candidates.size()) {
            const uint64_t key = candidates[i].key;
            const int len = packed_key_len(key);
            const uint64_t payload = packed_key_payload(key);
            if (score_mode == ScoreMode::Edit) unpack_payload(payload, len, observed_buf);

            const uint32_t offset = static_cast<uint32_t>(hits_.size());
            uint16_t count = 0;
            uint16_t last_idx = std::numeric_limits<uint16_t>::max();

            while (i < candidates.size() && candidates[i].key == key) {
                uint16_t idx = candidates[i].idx;
                ++i;
                if (idx == last_idx) continue;
                last_idx = idx;

                int d = max_dist + 1;
                if (score_mode == ScoreMode::Hamming) {
                    d = packed_hamming(payload, target_packed[idx], len);
                } else {
                    d = edit_distance_bounded(observed_buf, len, targets[idx].data(),
                                              static_cast<int>(targets[idx].size()), max_dist);
                }
                if (d <= max_dist) {
                    if (count == std::numeric_limits<uint16_t>::max()) {
                        std::cerr << "ERROR: too many packed hits for one HD cache key in "
                                  << label << "\n";
                        std::exit(2);
                    }
                    hits_.push_back({idx, static_cast<uint8_t>(d), 0});
                    ++count;
                }
            }

            if (count > 0) entries_.push_back({key, offset, count, 0});
            else hits_.resize(offset);
        }

        candidates.clear();
        candidates.shrink_to_fit();
        rebuild_index_from_entries(label);
        enabled_ = true;
        std::cerr << "built packed " << label << " cache: keys=" << entries_.size()
                  << " hits=" << hits_.size() << " max_h=" << max_dist << "\n";
        if (!cache_path.empty()) save_cache(cache_path, target_checksum, label);
    }

    void build_edit1_universe(const std::vector<std::string>& targets, const char* label) {
        entries_.clear();
        hits_.clear();
        clear_index();
        max_distance_ = 1;
        target_count_ = targets.size();
        score_mode_ = ScoreMode::Edit;
        enabled_ = false;

        if (targets.empty()) return;
        if (targets.size() > static_cast<size_t>(std::numeric_limits<uint16_t>::max())) {
            std::cerr << "ERROR: packed HD E0/E1 cache target count exceeds uint16_t index range: "
                      << targets.size() << "\n";
            std::exit(2);
        }

        size_t max_variants = 0;
        for (auto& target : targets) {
            size_t n = target.size();
            max_variants = std::max(max_variants, 1 + 3 * n + n + 4 * (n + 1));
        }
        std::vector<PackedCandidate> candidates;
        candidates.reserve(targets.size() * max_variants);
        std::vector<std::string> variants;
        for (size_t idx = 0; idx < targets.size(); ++idx) {
            edit1_variants(targets[idx], variants);
            for (auto& observed : variants) {
                uint64_t packed = 0;
                if (!pack_sequence(observed.data(), static_cast<int>(observed.size()), packed)) continue;
                candidates.push_back({make_packed_key(static_cast<int>(observed.size()), packed),
                                      static_cast<uint16_t>(idx)});
            }
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const PackedCandidate& a, const PackedCandidate& b) {
                      if (a.key != b.key) return a.key < b.key;
                      return a.idx < b.idx;
                  });

        char observed_buf[kMaxLen + 1];
        size_t i = 0;
        while (i < candidates.size()) {
            const uint64_t key = candidates[i].key;
            const int len = packed_key_len(key);
            const uint64_t payload = packed_key_payload(key);
            unpack_payload(payload, len, observed_buf);

            const uint32_t offset = static_cast<uint32_t>(hits_.size());
            uint16_t count = 0;
            uint16_t last_idx = std::numeric_limits<uint16_t>::max();
            while (i < candidates.size() && candidates[i].key == key) {
                uint16_t idx = candidates[i].idx;
                ++i;
                if (idx == last_idx) continue;
                last_idx = idx;
                int d = edit_distance_bounded(observed_buf, len, targets[idx].data(),
                                              static_cast<int>(targets[idx].size()), 1);
                if (d <= 1) {
                    hits_.push_back({idx, static_cast<uint8_t>(d), 0});
                    ++count;
                }
            }
            if (count > 0) entries_.push_back({key, offset, count, 0});
            else hits_.resize(offset);
        }

        rebuild_index_from_entries(label);
        enabled_ = true;
        std::cerr << "built packed " << label << " E0/E1 cache: keys=" << entries_.size()
                  << " hits=" << hits_.size() << "\n";
    }

    void build_suffix_wildcard_edit_universe_for_fixed_width(
        const std::vector<std::string>& targets,
        int fixed_width,
        int max_dist,
        const char* label) {
        entries_.clear();
        hits_.clear();
        clear_index();
        max_distance_ = max_dist;
        target_count_ = targets.size();
        score_mode_ = ScoreMode::Edit;
        enabled_ = false;

        if (targets.empty()) return;
        if (targets.size() > static_cast<size_t>(std::numeric_limits<uint16_t>::max())) {
            std::cerr << "ERROR: packed fixed-width cache target count exceeds uint16_t index range: "
                      << targets.size() << "\n";
            std::exit(2);
        }
        if (fixed_width <= 0 || fixed_width > 28) {
            std::cerr << "ERROR: invalid fixed-width packed key length " << fixed_width
                      << " for " << label << "\n";
            std::exit(2);
        }
        if (max_dist < 0 || max_dist > 2) {
            std::cerr << "ERROR: packed fixed-width cache currently supports E0/E1/E2 only; requested E"
                      << max_dist << " for " << label << "\n";
            std::exit(2);
        }

        std::vector<PackedSizedCandidate> candidates;
        std::vector<LocalEditVariant> local_variants;
        std::string observed;

        for (size_t idx = 0; idx < targets.size(); ++idx) {
            bool allowed_lengths[kMaxLen + 1] = {false};
            int min_allowed_len = kMaxLen;
            int max_allowed_len = 0;
            const int target_len = static_cast<int>(targets[idx].size());
            for (int delta = -max_dist; delta <= max_dist; ++delta) {
                const int observed_len = target_len + delta;
                if (observed_len <= 0 || observed_len > fixed_width ||
                    observed_len > kMaxLen) {
                    continue;
                }
                allowed_lengths[observed_len] = true;
                min_allowed_len = std::min(min_allowed_len, observed_len);
                max_allowed_len = std::max(max_allowed_len, observed_len);
            }
            if (max_allowed_len == 0) continue;

            local_variants.clear();
            observed.clear();
            enumerate_edit_variants_rec(targets[idx], 0, 0, max_dist,
                                        allowed_lengths, min_allowed_len,
                                        max_allowed_len, observed, local_variants);
            std::sort(local_variants.begin(), local_variants.end(),
                      [](const LocalEditVariant& a, const LocalEditVariant& b) {
                          if (a.key != b.key) return a.key < b.key;
                          return a.distance < b.distance;
                      });
            size_t i = 0;
            while (i < local_variants.size()) {
                const uint64_t observed_key = local_variants[i].key;
                uint8_t best_distance = local_variants[i].distance;
                ++i;
                while (i < local_variants.size() && local_variants[i].key == observed_key) {
                    best_distance = std::min(best_distance, local_variants[i].distance);
                    ++i;
                }

                const int observed_len = packed_key_len(observed_key);
                const int suffix_len = fixed_width - observed_len;
                if (suffix_len < 0) continue;
                const uint64_t observed_payload = packed_key_payload(observed_key);
                const uint64_t suffix_count = 1ULL << (2 * suffix_len);
                const uint64_t fixed_prefix = observed_payload << (2 * suffix_len);
                for (uint64_t suffix = 0; suffix < suffix_count; ++suffix) {
                    candidates.push_back({
                        make_packed_key(fixed_width, fixed_prefix | suffix),
                        static_cast<uint16_t>(idx),
                        best_distance,
                        static_cast<uint8_t>(observed_len),
                    });
                }
            }
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const PackedSizedCandidate& a, const PackedSizedCandidate& b) {
                      if (a.key != b.key) return a.key < b.key;
                      if (a.idx != b.idx) return a.idx < b.idx;
                      if (a.observed_len != b.observed_len) {
                          return a.observed_len < b.observed_len;
                      }
                      return a.distance < b.distance;
                  });

        size_t i = 0;
        while (i < candidates.size()) {
            const uint64_t key = candidates[i].key;
            const uint32_t offset = static_cast<uint32_t>(hits_.size());
            uint16_t count = 0;
            while (i < candidates.size() && candidates[i].key == key) {
                const uint16_t idx = candidates[i].idx;
                const uint8_t observed_len = candidates[i].observed_len;
                uint8_t best_distance = candidates[i].distance;
                ++i;
                while (i < candidates.size() && candidates[i].key == key &&
                       candidates[i].idx == idx &&
                       candidates[i].observed_len == observed_len) {
                    best_distance = std::min(best_distance, candidates[i].distance);
                    ++i;
                }
                if (count == std::numeric_limits<uint16_t>::max()) {
                    std::cerr << "ERROR: too many packed hits for one fixed-width key in "
                              << label << "\n";
                    std::exit(2);
                }
                hits_.push_back({idx, best_distance, observed_len});
                ++count;
            }
            if (count > 0) entries_.push_back({key, offset, count, 0});
            else hits_.resize(offset);
        }

        candidates.clear();
        candidates.shrink_to_fit();
        rebuild_index_from_entries(label);
        enabled_ = true;
        std::cerr << "built khash " << label
                  << " fixed-width suffix-wildcard search: keys=" << kh_size(index_)
                  << " hits=" << hits_.size()
                  << " max_e=" << max_dist
                  << " fixed_width=" << fixed_width
                  << "\n";
    }

    void build_best_edit_universe_for_query_lengths(const std::vector<std::string>& targets,
                                                    const std::vector<int>& query_lengths,
                                                    int max_dist,
                                                    const char* label,
                                                    TargetSliceMode slice_mode = TargetSliceMode::Full,
                                                    bool normalize_terminal_indels = false,
                                                    bool keep_only_best = true) {
        entries_.clear();
        hits_.clear();
        clear_index();
        max_distance_ = max_dist;
        target_count_ = targets.size();
        score_mode_ = ScoreMode::Edit;
        enabled_ = false;

        if (targets.empty() || query_lengths.empty()) return;
        if (targets.size() > static_cast<size_t>(std::numeric_limits<uint16_t>::max())) {
            std::cerr << "ERROR: packed HD best-edit cache target count exceeds uint16_t index range: "
                      << targets.size() << "\n";
            std::exit(2);
        }
        if (max_dist < 0 || max_dist > 2) {
            std::cerr << "ERROR: packed HD best-edit cache currently supports E0/E1/E2 only; requested E"
                      << max_dist << " for " << label << "\n";
            std::exit(2);
        }

        bool allowed_lengths[kMaxLen + 1] = {false};
        int min_allowed_len = kMaxLen;
        int max_allowed_len = 0;
        for (int len : query_lengths) {
            if (len <= 0 || len > kMaxLen) {
                std::cerr << "ERROR: invalid query length " << len
                          << " for packed HD best-edit cache " << label << "\n";
                std::exit(2);
            }
            allowed_lengths[len] = true;
            min_allowed_len = std::min(min_allowed_len, len);
            max_allowed_len = std::max(max_allowed_len, len);
        }

        std::vector<PackedScoredCandidate> candidates;
        std::vector<LocalEditVariant> local_variants;
        std::string observed;
        std::string target_slice;
        char observed_buf[kMaxLen + 1];
        for (size_t idx = 0; idx < targets.size(); ++idx) {
            local_variants.clear();
            if (slice_mode == TargetSliceMode::Full) {
                observed.clear();
                enumerate_edit_variants_rec(targets[idx], 0, 0, max_dist,
                                            allowed_lengths, min_allowed_len,
                                            max_allowed_len, observed, local_variants);
            } else {
                for (int query_len : query_lengths) {
                    if (static_cast<int>(targets[idx].size()) < query_len) continue;
                    if (slice_mode == TargetSliceMode::Prefix) {
                        target_slice.assign(targets[idx].data(), query_len);
                    } else {
                        target_slice.assign(targets[idx].data() + targets[idx].size() - query_len,
                                            query_len);
                    }
                    bool single_allowed_lengths[kMaxLen + 1] = {false};
                    single_allowed_lengths[query_len] = true;
                    observed.clear();
                    enumerate_edit_variants_rec(target_slice, 0, 0, max_dist,
                                                single_allowed_lengths, query_len,
                                                query_len, observed, local_variants);
                }
            }
            std::sort(local_variants.begin(), local_variants.end(),
                      [](const LocalEditVariant& a, const LocalEditVariant& b) {
                          if (a.key != b.key) return a.key < b.key;
                          return a.distance < b.distance;
                      });
            size_t i = 0;
            while (i < local_variants.size()) {
                const uint64_t key = local_variants[i].key;
                uint8_t best_raw_distance = local_variants[i].distance;
                ++i;
                while (i < local_variants.size() && local_variants[i].key == key) {
                    best_raw_distance = std::min(best_raw_distance, local_variants[i].distance);
                    ++i;
                }
                uint8_t best_distance = best_raw_distance;
                if (normalize_terminal_indels && best_raw_distance > 0) {
                    const int len = packed_key_len(key);
                    unpack_payload(packed_key_payload(key), len, observed_buf);
                    TerminalNormalizedEdit normalized = terminal_normalized_edit_distance(
                        observed_buf, len, targets[idx].data(),
                        static_cast<int>(targets[idx].size()), max_dist);
                    if (normalized.edit_distance < best_raw_distance) {
                        best_distance = static_cast<uint8_t>(normalized.edit_distance);
                    }
                }
                candidates.push_back({key, static_cast<uint16_t>(idx), best_distance,
                                      best_raw_distance});
            }
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const PackedScoredCandidate& a, const PackedScoredCandidate& b) {
                      if (a.key != b.key) return a.key < b.key;
                      if (a.idx != b.idx) return a.idx < b.idx;
                      return a.distance < b.distance;
                  });

        size_t i = 0;
        while (i < candidates.size()) {
            const uint64_t key = candidates[i].key;
            const uint32_t offset = static_cast<uint32_t>(hits_.size());
            const size_t key_start = i;
            uint8_t best_key_distance = std::numeric_limits<uint8_t>::max();
            uint8_t best_key_raw_distance = std::numeric_limits<uint8_t>::max();

            while (i < candidates.size() && candidates[i].key == key) {
                const uint16_t idx = candidates[i].idx;
                uint8_t best_idx_distance = candidates[i].distance;
                uint8_t best_idx_raw_distance = candidates[i].raw_distance;
                ++i;
                while (i < candidates.size() && candidates[i].key == key &&
                       candidates[i].idx == idx) {
                    best_idx_distance = std::min(best_idx_distance, candidates[i].distance);
                    best_idx_raw_distance =
                        std::min(best_idx_raw_distance, candidates[i].raw_distance);
                    ++i;
                }
                best_key_distance = std::min(best_key_distance, best_idx_distance);
                best_key_raw_distance = std::min(best_key_raw_distance, best_idx_raw_distance);
            }

            uint16_t count = 0;
            size_t j = key_start;
            while (j < i) {
                const uint16_t idx = candidates[j].idx;
                uint8_t best_idx_distance = candidates[j].distance;
                uint8_t best_idx_raw_distance = candidates[j].raw_distance;
                ++j;
                while (j < i && candidates[j].key == key && candidates[j].idx == idx) {
                    best_idx_distance = std::min(best_idx_distance, candidates[j].distance);
                    best_idx_raw_distance =
                        std::min(best_idx_raw_distance, candidates[j].raw_distance);
                    ++j;
                }
                if (!keep_only_best ||
                    best_idx_distance == best_key_distance ||
                    (normalize_terminal_indels &&
                     best_idx_raw_distance == best_key_raw_distance)) {
                    if (count == std::numeric_limits<uint16_t>::max()) {
                        std::cerr << "ERROR: too many packed hits for one HD best-edit key in "
                                  << label << "\n";
                        std::exit(2);
                    }
                    hits_.push_back({idx, best_idx_distance, 0});
                    ++count;
                }
            }

            if (count > 0) entries_.push_back({key, offset, count, 0});
            else hits_.resize(offset);
        }

        candidates.clear();
        candidates.shrink_to_fit();
        rebuild_index_from_entries(label);
        enabled_ = true;
        std::cerr << "built khash " << label << " best-edit search: keys=" << kh_size(index_)
                  << " hits=" << hits_.size() << " max_e=" << max_dist
                  << " hit_policy=" << (keep_only_best ? "best" : "all")
                  << " query_lengths=";
        for (size_t j = 0; j < query_lengths.size(); ++j) {
            if (j) std::cerr << ",";
            std::cerr << query_lengths[j];
        }
        std::cerr << "\n";
    }

    bool lookup(const char* seq, int len, std::vector<PackedHit>& out) const {
        out.clear();
        HitSpan hits = lookup_span(seq, len);
        if (!hits.found) return false;
        if (hits.count > 0) out.insert(out.end(), hits.data, hits.data + hits.count);
        return true;
    }

    HitSpan lookup_span(const char* seq, int len) const {
        if (!enabled_) return {};
        uint64_t packed = 0;
        if (!pack_sequence(seq, len, packed)) return {};
        const uint64_t key = make_packed_key(len, packed);
        khiter_t it = kh_get(memo64, index_, static_cast<khint64_t>(key));
        if (it == kh_end(index_)) return {};
        MemoValue mv = kh_value(index_, it);
        return {hits_.data() + mv.offset, mv.count, true};
    }

    size_t entry_count() const { return kh_size(index_); }
    size_t hit_count() const { return hits_.size(); }

  private:
    void clear_index() {
        kh_clear(memo64, index_);
    }

    void reserve_index(size_t n, const char* label) {
        if (n == 0) return;
        if (kh_resize(memo64, index_, static_cast<khint_t>(n)) < 0) {
            std::cerr << "ERROR: failed to reserve packed answer khash index for "
                      << label << "\n";
            std::exit(2);
        }
    }

    void add_index_entry(uint64_t key, uint32_t offset, uint16_t count, const char* label) {
        int ret = 0;
        khiter_t it = kh_put(memo64, index_, static_cast<khint64_t>(key), &ret);
        if (ret < 0) {
            std::cerr << "ERROR: failed to grow packed answer khash index for "
                      << label << "\n";
            std::exit(2);
        }
        if (ret == 0) {
            std::cerr << "ERROR: duplicate packed answer khash key while building "
                      << label << "\n";
            std::exit(2);
        }
        kh_value(index_, it) = MemoValue{offset, count};
    }

    void rebuild_index_from_entries(const char* label) {
        clear_index();
        reserve_index(entries_.size(), label);
        for (const PackedEntry& entry : entries_) {
            add_index_entry(entry.key, entry.offset, entry.count, label);
        }
    }

    bool load_cache(const std::string& path, size_t target_count, uint64_t target_checksum,
                    int max_dist, ScoreMode score_mode, const char* label) {
        std::ifstream in(path, std::ios::binary);
        if (!in) return false;

        PackedCacheHeader header;
        in.read(reinterpret_cast<char*>(&header), sizeof(header));
        const char expected_magic[16] = {'H', 'D', 'R', '1', 'C', 'A', 'C', 'H', 'E', 'v', '1', 0};
        if (!in || std::memcmp(header.magic, expected_magic, sizeof(expected_magic)) != 0 ||
            header.version != 1 ||
            header.score_mode != static_cast<uint32_t>(score_mode) ||
            header.max_distance != static_cast<uint32_t>(max_dist) ||
            header.target_count != static_cast<uint32_t>(target_count) ||
            header.target_checksum != target_checksum) {
            return false;
        }

        entries_.resize(static_cast<size_t>(header.entry_count));
        hits_.resize(static_cast<size_t>(header.hit_count));
        in.read(reinterpret_cast<char*>(entries_.data()),
                static_cast<std::streamsize>(entries_.size() * sizeof(PackedEntry)));
        in.read(reinterpret_cast<char*>(hits_.data()),
                static_cast<std::streamsize>(hits_.size() * sizeof(PackedHit)));
        if (!in) {
            entries_.clear();
            hits_.clear();
            clear_index();
            return false;
        }
        rebuild_index_from_entries(label);
        std::cerr << "loaded packed " << label << " cache: keys=" << entries_.size()
                  << " hits=" << hits_.size() << " khash_keys=" << kh_size(index_)
                  << " max_h=" << max_dist << "\n";
        return true;
    }

    void save_cache(const std::string& path, uint64_t target_checksum, const char* label) const {
        const std::string tmp_path = path + ".tmp";
        std::ofstream out(tmp_path, std::ios::binary);
        if (!out) {
            std::cerr << "WARNING: cannot write packed " << label << " cache to "
                      << tmp_path << "\n";
            return;
        }
        PackedCacheHeader header;
        const char magic[16] = {'H', 'D', 'R', '1', 'C', 'A', 'C', 'H', 'E', 'v', '1', 0};
        std::memcpy(header.magic, magic, sizeof(magic));
        header.score_mode = static_cast<uint32_t>(score_mode_);
        header.max_distance = static_cast<uint32_t>(max_distance_);
        header.target_count = static_cast<uint32_t>(target_count_);
        header.target_checksum = target_checksum;
        header.entry_count = static_cast<uint64_t>(entries_.size());
        header.hit_count = static_cast<uint64_t>(hits_.size());
        out.write(reinterpret_cast<const char*>(&header), sizeof(header));
        out.write(reinterpret_cast<const char*>(entries_.data()),
                  static_cast<std::streamsize>(entries_.size() * sizeof(PackedEntry)));
        out.write(reinterpret_cast<const char*>(hits_.data()),
                  static_cast<std::streamsize>(hits_.size() * sizeof(PackedHit)));
        out.close();
        if (!out) {
            std::cerr << "WARNING: failed while writing packed " << label << " cache to "
                      << tmp_path << "\n";
            std::remove(tmp_path.c_str());
            return;
        }
        if (std::rename(tmp_path.c_str(), path.c_str()) != 0) {
            std::cerr << "WARNING: cannot install packed " << label << " cache at "
                      << path << ": " << std::strerror(errno) << "\n";
            std::remove(tmp_path.c_str());
        }
    }

    bool enabled_ = false;
    int max_distance_ = 0;
    size_t target_count_ = 0;
    ScoreMode score_mode_ = ScoreMode::Hamming;
    std::vector<PackedEntry> entries_;
    std::vector<PackedHit> hits_;
    khash_t(memo64)* index_ = nullptr;
};

int hit_distance_for_idx(const std::vector<PackedHit>& hits, int target_idx, int missing) {
    auto it = std::lower_bound(hits.begin(), hits.end(), target_idx,
                               [](const PackedHit& h, int idx) { return h.idx < idx; });
    if (it != hits.end() && it->idx == target_idx) return it->distance;
    return missing;
}

int hit_distance_for_idx(HitSpan hits, int target_idx, int missing) {
    if (hits.count == 0) return missing;
    const PackedHit* begin = hits.data;
    const PackedHit* end = hits.data + hits.count;
    auto it = std::lower_bound(begin, end, target_idx,
                               [](const PackedHit& h, int idx) { return h.idx < idx; });
    if (it != end && it->idx == target_idx) return it->distance;
    return missing;
}

struct ObservedLookupScratch {
    std::vector<PackedHit> tmp;
    std::vector<PackedHit> computed;
    std::vector<int> best;
    std::vector<int> touched;
    std::vector<int> seen_stamp;
    std::vector<uint16_t> seen_touched;
    int seen_generation = 1;
    std::vector<uint64_t> sigbuf;

    void ensure_target_count(size_t n) {
        if (best.size() == n && seen_stamp.size() == n) return;
        best.assign(n, std::numeric_limits<int>::max());
        seen_stamp.assign(n, 0);
        touched.clear();
        seen_touched.clear();
        seen_generation = 1;
    }

    void clear_best() {
        for (int idx : touched) best[idx] = std::numeric_limits<int>::max();
        touched.clear();
    }

    void add_best(int idx, int distance) {
        if (best[idx] == std::numeric_limits<int>::max()) touched.push_back(idx);
        if (distance < best[idx]) best[idx] = distance;
    }

    void start_seen() {
        seen_touched.clear();
        ++seen_generation;
        if (seen_generation == std::numeric_limits<int>::max()) {
            std::fill(seen_stamp.begin(), seen_stamp.end(), 0);
            seen_generation = 1;
        }
    }
};

struct DpOnlyMemoStats {
    uint64_t lookup_calls = 0;
    uint64_t memo_hits = 0;
    uint64_t memo_misses = 0;
    uint64_t dp_calls = 0;
    uint64_t dp_accepts = 0;
    uint64_t answer_hits_stored = 0;

    void add(const DpOnlyMemoStats& other) {
        lookup_calls += other.lookup_calls;
        memo_hits += other.memo_hits;
        memo_misses += other.memo_misses;
        dp_calls += other.dp_calls;
        dp_accepts += other.dp_accepts;
        answer_hits_stored += other.answer_hits_stored;
    }
};

struct DpOnlyMemo {
    khash_t(memo64)* memo = nullptr;
    std::unordered_map<std::string, MemoValue> string_memo;
    std::vector<PackedHit> hit_pool;
    DpOnlyMemoStats stats;

    DpOnlyMemo() : memo(kh_init(memo64)) {
        if (!memo) {
            std::cerr << "ERROR: failed to allocate DP miss memo\n";
            std::exit(2);
        }
    }

    ~DpOnlyMemo() {
        kh_destroy(memo64, memo);
    }

    DpOnlyMemo(const DpOnlyMemo&) = delete;
    DpOnlyMemo& operator=(const DpOnlyMemo&) = delete;
};

bool lookup_dp_local_packed(DpOnlyMemo& memo, uint64_t key, MemoValue& value) {
    khiter_t it = kh_get(memo64, memo.memo, static_cast<khint64_t>(key));
    if (it == kh_end(memo.memo)) return false;
    value = kh_value(memo.memo, it);
    return true;
}

void remember_dp_local_packed(DpOnlyMemo& memo, uint64_t key,
                              const std::vector<PackedHit>& hits) {
    int ret = 0;
    khiter_t it = kh_put(memo64, memo.memo, static_cast<khint64_t>(key), &ret);
    if (ret == 0) return;
    if (ret < 0) {
        std::cerr << "ERROR: failed to grow DP miss memo\n";
        std::exit(2);
    }
    MemoValue mv{static_cast<uint32_t>(memo.hit_pool.size()), static_cast<uint16_t>(hits.size())};
    if (!hits.empty()) memo.hit_pool.insert(memo.hit_pool.end(), hits.begin(), hits.end());
    kh_value(memo.memo, it) = mv;
}

void remember_dp_local_string(DpOnlyMemo& memo, const std::string& key,
                              const std::vector<PackedHit>& hits) {
    if (memo.string_memo.find(key) != memo.string_memo.end()) return;
    MemoValue mv{static_cast<uint32_t>(memo.hit_pool.size()), static_cast<uint16_t>(hits.size())};
    if (!hits.empty()) memo.hit_pool.insert(memo.hit_pool.end(), hits.begin(), hits.end());
    memo.string_memo[key] = mv;
}

struct DeletionCandidateIndex;

void lookup_dp_only(const char* query, int len, const std::vector<std::string>& targets,
                    const DeletionCandidateIndex* deletion_index,
                    int max_dist, ObservedLookupScratch& scratch,
                    DpOnlyMemo& memo, std::vector<PackedHit>& out);

class PrebuiltObservedAnswerCache {
  public:
    PrebuiltObservedAnswerCache() : memo_(kh_init(memo64)) {
        if (!memo_) {
            std::cerr << "ERROR: failed to allocate prebuilt observed cache\n";
            std::exit(2);
        }
    }

    ~PrebuiltObservedAnswerCache() {
        kh_destroy(memo64, memo_);
    }

    PrebuiltObservedAnswerCache(const PrebuiltObservedAnswerCache&) = delete;
    PrebuiltObservedAnswerCache& operator=(const PrebuiltObservedAnswerCache&) = delete;

    void build_from(const std::unordered_map<uint64_t, MemoValue>& memo,
                    const std::unordered_map<std::string, MemoValue>& string_memo,
                    const std::vector<PackedHit>& hits,
                    const char* label,
                    const ObservedCacheStats& stats) {
        label_ = label;
        stats_ = stats;
        hits_ = hits;
        string_memo_ = string_memo;
        kh_clear(memo64, memo_);
        if (memo.size() > 0 && kh_resize(memo64, memo_, static_cast<khint_t>(memo.size())) < 0) {
            std::cerr << "ERROR: failed to reserve prebuilt observed cache " << label_ << "\n";
            std::exit(2);
        }
        for (auto& kv : memo) {
            int ret = 0;
            khiter_t it = kh_put(memo64, memo_, static_cast<khint64_t>(kv.first), &ret);
            if (ret < 0) {
                std::cerr << "ERROR: failed to build prebuilt observed cache " << label_ << "\n";
                std::exit(2);
            }
            kh_value(memo_, it) = kv.second;
        }
        std::cerr << "prebuilt observed " << label_ << " cache: keys=" << kh_size(memo_)
                  << " string_keys=" << string_memo_.size()
                  << " hits=" << hits_.size() << "\n";
    }

    HitSpan lookup_span(const char* query, int len) const {
        uint64_t packed = 0;
        if (!pack_sequence(query, len, packed)) {
            auto it = string_memo_.find(std::string(query, len));
            if (it == string_memo_.end()) return {};
            if (it->second.count == 0) return {nullptr, 0, true};
            return {hits_.data() + it->second.offset, it->second.count, true};
        }
        uint64_t key = make_packed_key(len, packed);
        khiter_t it = kh_get(memo64, memo_, static_cast<khint64_t>(key));
        if (it == kh_end(memo_)) return {};
        MemoValue mv = kh_value(memo_, it);
        if (mv.count == 0) return {nullptr, 0, true};
        return {hits_.data() + mv.offset, mv.count, true};
    }

    size_t entry_count() const { return kh_size(memo_); }
    size_t string_entry_count() const { return string_memo_.size(); }
    size_t hit_count() const { return hits_.size(); }
    const char* label() const { return label_; }
    const ObservedCacheStats& stats() const { return stats_; }

  private:
    const char* label_ = "";
    ObservedCacheStats stats_;
    khash_t(memo64)* memo_ = nullptr;
    std::unordered_map<std::string, MemoValue> string_memo_;
    std::vector<PackedHit> hits_;
};

struct DeletionCandidateIndex {
    const std::vector<std::string>* targets = nullptr;
    int max_deletions = 2;
    std::unordered_map<uint64_t, std::vector<uint16_t>> index;
    size_t indexed_hits = 0;

    void build(const std::vector<std::string>& t, int max_del, const char* label) {
        targets = &t;
        max_deletions = max_del;
        index.clear();
        indexed_hits = 0;
        std::vector<uint64_t> sigs;
        for (int idx = 0; idx < static_cast<int>(t.size()); ++idx) {
            if (!packed_deletion_signature_keys(t[idx].data(), static_cast<int>(t[idx].size()),
                                                max_deletions, sigs)) {
                std::cerr << "ERROR: non-ACGT target in packed deletion-signature index "
                          << label << " at index " << idx << "\n";
                std::exit(2);
            }
            for (uint64_t sig : sigs) {
                index[sig].push_back(static_cast<uint16_t>(idx));
                ++indexed_hits;
            }
        }
        for (auto& kv : index) {
            auto& v = kv.second;
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        }
        std::cerr << "built " << label << " deletion-signature index: keys=" << index.size()
                  << " hits=" << indexed_hits << " max_del=" << max_deletions << "\n";
    }
};

void lookup_dp_only(const char* query, int len, const std::vector<std::string>& targets,
                    const DeletionCandidateIndex* deletion_index,
                    int max_dist, ObservedLookupScratch& scratch,
                    DpOnlyMemo& memo, std::vector<PackedHit>& out) {
    out.clear();
    ++memo.stats.lookup_calls;
    uint64_t packed = 0;
    if (pack_sequence(query, len, packed)) {
        uint64_t key = make_packed_key(len, packed);
        MemoValue mv;
        if (lookup_dp_local_packed(memo, key, mv)) {
            ++memo.stats.memo_hits;
            if (mv.count > 0) {
                out.insert(out.end(), memo.hit_pool.data() + mv.offset,
                           memo.hit_pool.data() + mv.offset + mv.count);
            }
            return;
        }
        ++memo.stats.memo_misses;
        if (deletion_index && packed_deletion_signature_keys(query, len,
                                                             deletion_index->max_deletions,
                                                             scratch.sigbuf)) {
            scratch.ensure_target_count(targets.size());
            scratch.start_seen();
            for (uint64_t sig : scratch.sigbuf) {
                auto sig_it = deletion_index->index.find(sig);
                if (sig_it == deletion_index->index.end()) continue;
                for (uint16_t idx : sig_it->second) {
                    if (scratch.seen_stamp[idx] == scratch.seen_generation) continue;
                    scratch.seen_stamp[idx] = scratch.seen_generation;
                    scratch.seen_touched.push_back(idx);
                }
            }
            for (uint16_t idx : scratch.seen_touched) {
                const std::string& target = targets[idx];
                ++memo.stats.dp_calls;
                int d = edit_distance_bounded(query, len, target.data(),
                                              static_cast<int>(target.size()), max_dist);
                if (d <= max_dist) {
                    ++memo.stats.dp_accepts;
                    out.push_back({idx, static_cast<uint8_t>(d), 0});
                }
            }
            std::sort(out.begin(), out.end(),
                      [](const PackedHit& a, const PackedHit& b) { return a.idx < b.idx; });
        } else {
            for (int idx = 0; idx < static_cast<int>(targets.size()); ++idx) {
                const std::string& target = targets[idx];
                if (std::abs(len - static_cast<int>(target.size())) > max_dist) continue;
                ++memo.stats.dp_calls;
                int d = edit_distance_bounded(query, len, target.data(),
                                              static_cast<int>(target.size()), max_dist);
                if (d <= max_dist) {
                    ++memo.stats.dp_accepts;
                    out.push_back({static_cast<uint16_t>(idx), static_cast<uint8_t>(d), 0});
                }
            }
        }
        remember_dp_local_packed(memo, key, out);
        memo.stats.answer_hits_stored += out.size();
        return;
    }

    std::string key(query, len);
    auto it = memo.string_memo.find(key);
    if (it != memo.string_memo.end()) {
        ++memo.stats.memo_hits;
        MemoValue mv = it->second;
        if (mv.count > 0) {
            out.insert(out.end(), memo.hit_pool.data() + mv.offset,
                       memo.hit_pool.data() + mv.offset + mv.count);
        }
        return;
    }
    ++memo.stats.memo_misses;
    for (int idx = 0; idx < static_cast<int>(targets.size()); ++idx) {
        const std::string& target = targets[idx];
        if (std::abs(len - static_cast<int>(target.size())) > max_dist) continue;
        ++memo.stats.dp_calls;
        int d = edit_distance_bounded(query, len, target.data(),
                                      static_cast<int>(target.size()), max_dist);
        if (d <= max_dist) {
            ++memo.stats.dp_accepts;
            out.push_back({static_cast<uint16_t>(idx), static_cast<uint8_t>(d), 0});
        }
    }
    remember_dp_local_string(memo, key, out);
    memo.stats.answer_hits_stored += out.size();
}

class ObservedAnswerCache {
  public:
    struct LocalCache {
        std::unordered_map<uint64_t, MemoValue> memo;
        std::unordered_map<std::string, MemoValue> string_memo;
        std::vector<PackedHit> hit_pool;
        std::vector<uint64_t> pending_keys;
        std::vector<std::string> pending_string_keys;
        ObservedCacheStats stats;
        uint64_t lookups_since_flush = 0;
        uint64_t lookups_since_sync = 0;
        size_t packed_sync_cursor = 0;
        size_t string_sync_cursor = 0;

        size_t pending_count() const {
            return pending_keys.size() + pending_string_keys.size();
        }
    };

    void init(const std::vector<std::string>& targets, const PackedAnswerCache& e01,
              const PackedAnswerCache& h2, const DeletionCandidateIndex& deletion_index,
              int max_dist, const char* label, bool use_e01) {
        targets_ = &targets;
        e01_ = &e01;
        h2_ = &h2;
        deletion_index_ = &deletion_index;
        max_distance_ = max_dist;
        label_ = label;
        use_e01_ = use_e01;
    }

    void lookup(const char* query, int len, LocalCache& local,
                ObservedLookupScratch& scratch, std::vector<PackedHit>& out) {
        out.clear();
        ++local.stats.lookup_calls;
        ++local.lookups_since_flush;
        ++local.lookups_since_sync;
        uint64_t packed = 0;
        if (!pack_sequence(query, len, packed)) {
            lookup_non_acgt(query, len, local, out);
            return;
        }

        uint64_t key = make_packed_key(len, packed);
        auto local_it = local.memo.find(key);
        if (local_it != local.memo.end()) {
            ++local.stats.memo_hits;
            append_local_hits(local, local_it->second, out);
            return;
        }

        if (should_sync(local)) {
            sync_from_global(local);
            local_it = local.memo.find(key);
            if (local_it != local.memo.end()) {
                ++local.stats.memo_hits;
                append_local_hits(local, local_it->second, out);
                return;
            }
        }

        ++local.stats.memo_misses;
        compute_acgt(query, len, scratch, local.stats);
        remember_local_packed(local, key, scratch.computed, true);
        out = scratch.computed;
        maybe_flush(local);
    }

    void flush(LocalCache& local) {
        if (local.pending_count() == 0 &&
            local.lookups_since_flush == 0 &&
            stats_empty(local.stats)) {
            return;
        }

        std::lock_guard<std::mutex> lock(mutex_);
        for (uint64_t key : local.pending_keys) {
            if (memo_.find(key) != memo_.end()) continue;
            auto local_it = local.memo.find(key);
            if (local_it == local.memo.end()) continue;
            MemoValue src = local_it->second;
            MemoValue dst{static_cast<uint32_t>(hit_pool_.size()), src.count};
            if (src.count > 0) {
                hit_pool_.insert(hit_pool_.end(), local.hit_pool.data() + src.offset,
                                 local.hit_pool.data() + src.offset + src.count);
            }
            memo_[key] = dst;
            memo_log_.push_back(key);
        }
        for (const std::string& key : local.pending_string_keys) {
            if (string_memo_.find(key) != string_memo_.end()) continue;
            auto local_it = local.string_memo.find(key);
            if (local_it == local.string_memo.end()) continue;
            MemoValue src = local_it->second;
            MemoValue dst{static_cast<uint32_t>(hit_pool_.size()), src.count};
            if (src.count > 0) {
                hit_pool_.insert(hit_pool_.end(), local.hit_pool.data() + src.offset,
                                 local.hit_pool.data() + src.offset + src.count);
            }
            string_memo_[key] = dst;
            string_memo_log_.push_back(key);
        }
        stats.add(local.stats);
        local.stats = ObservedCacheStats();
        local.pending_keys.clear();
        local.pending_string_keys.clear();
        local.lookups_since_flush = 0;
        local.lookups_since_sync = 0;
    }

    size_t memo_size() const { return memo_.size(); }
    size_t string_memo_size() const { return string_memo_.size(); }
    size_t hit_pool_size() const { return hit_pool_.size(); }
    const char* label() const { return label_; }

    void freeze_into(PrebuiltObservedAnswerCache& out) const {
        out.build_from(memo_, string_memo_, hit_pool_, label_, stats);
    }

    ObservedCacheStats stats;

  private:
    static constexpr size_t kFlushPendingThreshold = 4096;
    static constexpr uint64_t kFlushLookupInterval = 65536;
    static constexpr uint64_t kSyncLookupInterval = 65536;
    mutable std::mutex mutex_;

    static bool stats_empty(const ObservedCacheStats& s) {
        return s.lookup_calls == 0 &&
               s.memo_hits == 0 &&
               s.memo_misses == 0 &&
               s.e01_keys_found == 0 &&
               s.e01_hits_returned == 0 &&
               s.h2_keys_found == 0 &&
               s.h2_hits_returned == 0 &&
               s.deletion_signature_candidates == 0 &&
               s.deletion_candidates_after_dedupe == 0 &&
               s.skipped_already_resolved == 0 &&
               s.dp_calls == 0 &&
               s.dp_accepts == 0 &&
               s.non_acgt_keys == 0 &&
               s.non_acgt_dp_calls == 0 &&
               s.answer_hits_stored == 0;
    }

    static void append_local_hits(const LocalCache& local, MemoValue mv,
                                  std::vector<PackedHit>& out) {
        if (mv.count == 0) return;
        out.insert(out.end(), local.hit_pool.data() + mv.offset,
                   local.hit_pool.data() + mv.offset + mv.count);
    }

    void remember_local_packed(LocalCache& local, uint64_t key,
                               const std::vector<PackedHit>& hits, bool pending) const {
        if (local.memo.find(key) != local.memo.end()) return;
        MemoValue mv{static_cast<uint32_t>(local.hit_pool.size()), static_cast<uint16_t>(hits.size())};
        local.hit_pool.insert(local.hit_pool.end(), hits.begin(), hits.end());
        local.memo[key] = mv;
        if (pending) {
            local.pending_keys.push_back(key);
            local.stats.answer_hits_stored += hits.size();
        }
    }

    void remember_local_string(LocalCache& local, const std::string& key,
                               const std::vector<PackedHit>& hits, bool pending) const {
        if (local.string_memo.find(key) != local.string_memo.end()) return;
        MemoValue mv{static_cast<uint32_t>(local.hit_pool.size()), static_cast<uint16_t>(hits.size())};
        local.hit_pool.insert(local.hit_pool.end(), hits.begin(), hits.end());
        local.string_memo[key] = mv;
        if (pending) {
            local.pending_string_keys.push_back(key);
            local.stats.answer_hits_stored += hits.size();
        }
    }

    static bool should_sync(const LocalCache& local) {
        return local.lookups_since_sync >= kSyncLookupInterval;
    }

    void sync_from_global(LocalCache& local) const {
        std::lock_guard<std::mutex> lock(mutex_);
        for (size_t i = local.packed_sync_cursor; i < memo_log_.size(); ++i) {
            uint64_t key = memo_log_[i];
            if (local.memo.find(key) != local.memo.end()) continue;
            auto memo_it = memo_.find(key);
            if (memo_it == memo_.end()) continue;
            MemoValue src = memo_it->second;
            MemoValue dst{static_cast<uint32_t>(local.hit_pool.size()), src.count};
            if (src.count > 0) {
                local.hit_pool.insert(local.hit_pool.end(), hit_pool_.data() + src.offset,
                                      hit_pool_.data() + src.offset + src.count);
            }
            local.memo[key] = dst;
        }
        for (size_t i = local.string_sync_cursor; i < string_memo_log_.size(); ++i) {
            const std::string& key = string_memo_log_[i];
            if (local.string_memo.find(key) != local.string_memo.end()) continue;
            auto memo_it = string_memo_.find(key);
            if (memo_it == string_memo_.end()) continue;
            MemoValue src = memo_it->second;
            MemoValue dst{static_cast<uint32_t>(local.hit_pool.size()), src.count};
            if (src.count > 0) {
                local.hit_pool.insert(local.hit_pool.end(), hit_pool_.data() + src.offset,
                                      hit_pool_.data() + src.offset + src.count);
            }
            local.string_memo[key] = dst;
        }
        local.packed_sync_cursor = memo_log_.size();
        local.string_sync_cursor = string_memo_log_.size();
        local.lookups_since_sync = 0;
    }

    void maybe_flush(LocalCache& local) {
        if (local.pending_count() == 0) return;
        if (local.pending_count() >= kFlushPendingThreshold ||
            local.lookups_since_flush >= kFlushLookupInterval) {
            flush(local);
        }
    }

    void compute_acgt(const char* query, int len, ObservedLookupScratch& scratch,
                      ObservedCacheStats& local_stats) const {
        scratch.ensure_target_count(targets_->size());
        scratch.clear_best();
        if (use_e01_ && e01_->lookup(query, len, scratch.tmp)) {
            ++local_stats.e01_keys_found;
            local_stats.e01_hits_returned += scratch.tmp.size();
            for (auto& hit : scratch.tmp) scratch.add_best(hit.idx, hit.distance);
        }
        if (h2_->lookup(query, len, scratch.tmp)) {
            ++local_stats.h2_keys_found;
            local_stats.h2_hits_returned += scratch.tmp.size();
            for (auto& hit : scratch.tmp) scratch.add_best(hit.idx, hit.distance);
        }

        packed_deletion_signature_keys(query, len, deletion_index_->max_deletions, scratch.sigbuf);
        scratch.start_seen();
        for (uint64_t sig : scratch.sigbuf) {
            auto it = deletion_index_->index.find(sig);
            if (it == deletion_index_->index.end()) continue;
            local_stats.deletion_signature_candidates += it->second.size();
            for (uint16_t idx : it->second) {
                if (scratch.seen_stamp[idx] == scratch.seen_generation) continue;
                scratch.seen_stamp[idx] = scratch.seen_generation;
                scratch.seen_touched.push_back(idx);
            }
        }

        local_stats.deletion_candidates_after_dedupe += scratch.seen_touched.size();
        for (uint16_t idx : scratch.seen_touched) {
            if (scratch.best[idx] <= max_distance_) {
                ++local_stats.skipped_already_resolved;
                continue;
            }
            ++local_stats.dp_calls;
            const std::string& target = (*targets_)[idx];
            int d = edit_distance_bounded(query, len, target.data(),
                                          static_cast<int>(target.size()), max_distance_);
            if (d <= max_distance_) {
                ++local_stats.dp_accepts;
                scratch.add_best(idx, d);
            }
        }

        scratch.computed.clear();
        scratch.computed.reserve(scratch.touched.size());
        for (int idx : scratch.touched) {
            scratch.computed.push_back({static_cast<uint16_t>(idx),
                                        static_cast<uint8_t>(scratch.best[idx]), 0});
        }
        std::sort(scratch.computed.begin(), scratch.computed.end(),
                  [](const PackedHit& a, const PackedHit& b) { return a.idx < b.idx; });
        if (scratch.computed.size() > std::numeric_limits<uint16_t>::max()) {
            std::cerr << "ERROR: too many observed-cache hits for " << label_ << "\n";
            std::exit(2);
        }
    }

    void lookup_non_acgt(const char* query, int len, LocalCache& local,
                         std::vector<PackedHit>& out) {
        std::string key(query, len);
        auto local_it = local.string_memo.find(key);
        if (local_it != local.string_memo.end()) {
            ++local.stats.memo_hits;
            append_local_hits(local, local_it->second, out);
            return;
        }
        if (should_sync(local)) {
            sync_from_global(local);
            local_it = local.string_memo.find(key);
            if (local_it != local.string_memo.end()) {
                ++local.stats.memo_hits;
                append_local_hits(local, local_it->second, out);
                return;
            }
        }

        ++local.stats.memo_misses;
        ++local.stats.non_acgt_keys;
        std::vector<PackedHit> computed;
        for (int idx = 0; idx < static_cast<int>(targets_->size()); ++idx) {
            const std::string& target = (*targets_)[idx];
            if (std::abs(len - static_cast<int>(target.size())) > max_distance_) continue;
            ++local.stats.dp_calls;
            ++local.stats.non_acgt_dp_calls;
            int d = edit_distance_bounded(query, len, target.data(),
                                          static_cast<int>(target.size()), max_distance_);
            if (d <= max_distance_) {
                ++local.stats.dp_accepts;
                computed.push_back({static_cast<uint16_t>(idx), static_cast<uint8_t>(d), 0});
            }
        }
        if (computed.size() > std::numeric_limits<uint16_t>::max()) {
            std::cerr << "ERROR: too many non-ACGT observed-cache hits for " << label_ << "\n";
            std::exit(2);
        }
        remember_local_string(local, key, computed, true);
        out.swap(computed);
        maybe_flush(local);
    }

    const std::vector<std::string>* targets_ = nullptr;
    const PackedAnswerCache* e01_ = nullptr;
    const PackedAnswerCache* h2_ = nullptr;
    const DeletionCandidateIndex* deletion_index_ = nullptr;
    int max_distance_ = 2;
    const char* label_ = "";
    bool use_e01_ = true;
    std::unordered_map<uint64_t, MemoValue> memo_;
    std::unordered_map<std::string, MemoValue> string_memo_;
    std::vector<PackedHit> hit_pool_;
    std::vector<uint64_t> memo_log_;
    std::vector<std::string> string_memo_log_;
};

void print_observed_cache_stats(const ObservedAnswerCache& cache) {
    const auto& s = cache.stats;
    std::cerr << "observed_cache_stats " << cache.label()
              << ": memo_entries=" << cache.memo_size()
              << " string_memo_entries=" << cache.string_memo_size()
              << " memo_hits_pool=" << cache.hit_pool_size()
              << " lookup_calls=" << s.lookup_calls
              << " memo_hits=" << s.memo_hits
              << " memo_misses=" << s.memo_misses
              << " e01_keys_found=" << s.e01_keys_found
              << " e01_hits_returned=" << s.e01_hits_returned
              << " h2_keys_found=" << s.h2_keys_found
              << " h2_hits_returned=" << s.h2_hits_returned
              << " deletion_signature_candidates=" << s.deletion_signature_candidates
              << " deletion_candidates_after_dedupe=" << s.deletion_candidates_after_dedupe
              << " skipped_already_resolved=" << s.skipped_already_resolved
              << " dp_calls=" << s.dp_calls
              << " dp_accepts=" << s.dp_accepts
              << " (" << pct(s.dp_accepts, s.dp_calls) << "%)"
              << " non_acgt_keys=" << s.non_acgt_keys
              << " non_acgt_dp_calls=" << s.non_acgt_dp_calls
              << " answer_hits_stored=" << s.answer_hits_stored
              << "\n";
}

void print_prebuilt_observed_cache_stats(const PrebuiltObservedAnswerCache& cache) {
    const auto& s = cache.stats();
    std::cerr << "prebuilt_observed_cache_stats " << cache.label()
              << ": memo_entries=" << cache.entry_count()
              << " string_memo_entries=" << cache.string_entry_count()
              << " memo_hits_pool=" << cache.hit_count()
              << " build_lookup_calls=" << s.lookup_calls
              << " build_memo_hits=" << s.memo_hits
              << " build_memo_misses=" << s.memo_misses
              << " h2_keys_found=" << s.h2_keys_found
              << " h2_hits_returned=" << s.h2_hits_returned
              << " deletion_signature_candidates=" << s.deletion_signature_candidates
              << " deletion_candidates_after_dedupe=" << s.deletion_candidates_after_dedupe
              << " skipped_already_resolved=" << s.skipped_already_resolved
              << " dp_calls=" << s.dp_calls
              << " dp_accepts=" << s.dp_accepts
              << " (" << pct(s.dp_accepts, s.dp_calls) << "%)"
              << " non_acgt_keys=" << s.non_acgt_keys
              << " non_acgt_dp_calls=" << s.non_acgt_dp_calls
              << " answer_hits_stored=" << s.answer_hits_stored
              << "\n";
}

struct LookupScratch {
    std::vector<int> best;
    std::vector<int> touched;
    std::vector<int> seen_stamp;
    std::vector<int> seen_touched;
    int seen_generation = 1;
    std::vector<std::string> sigbuf;

    void ensure(size_t n) {
        if (best.size() < n) best.assign(n, std::numeric_limits<int>::max());
        if (seen_stamp.size() < n) seen_stamp.assign(n, 0);
    }

    void clear_best() {
        for (int idx : touched) best[idx] = std::numeric_limits<int>::max();
        touched.clear();
    }

    void add_hit(int idx, int distance) {
        if (best[idx] == std::numeric_limits<int>::max()) touched.push_back(idx);
        if (distance < best[idx]) best[idx] = distance;
    }

    void start_seen() {
        seen_touched.clear();
        ++seen_generation;
        if (seen_generation == std::numeric_limits<int>::max()) {
            std::fill(seen_stamp.begin(), seen_stamp.end(), 0);
            seen_generation = 1;
        }
    }

    bool mark_seen(int idx) {
        if (seen_stamp[idx] == seen_generation) return false;
        seen_stamp[idx] = seen_generation;
        seen_touched.push_back(idx);
        return true;
    }
};

// ---- anchor lookup (hamming index and/or deletion-signature index) ----
struct AnchorLookup {
    const std::vector<std::string>* anchors = nullptr;
    int max_distance = 0;
    bool use_hamming = false;
    bool use_deletion = false;
    int signature_deletions = 0;
    PackedAnswerCache hamming_cache;
    std::unordered_map<std::string, std::vector<int>> deletion_index;

    void build(const std::vector<std::string>& anc, int max_dist, const std::string& mode, int sig_del,
               const char* label, const std::string& cache_path = std::string()) {
        anchors = &anc;
        max_distance = max_dist;
        signature_deletions = sig_del;
        use_hamming = (mode == "hamming" || mode == "hamming-deletion-hash");
        use_deletion = (mode == "deletion-hash" || mode == "hamming-deletion-hash");
        std::vector<std::string> tmp;
        if (use_hamming) {
            hamming_cache.build_hamming_universe(anc, max_distance,
                                                 PackedAnswerCache::ScoreMode::Hamming, label,
                                                 cache_path);
        }
        if (use_deletion) {
            for (int idx = 0; idx < static_cast<int>(anc.size()); ++idx) {
                deletion_signatures(anc[idx], signature_deletions, tmp);
                for (auto& s : tmp) deletion_index[s].push_back(idx);
            }
            for (auto& kv : deletion_index) {
                auto& v = kv.second;
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
            }
        }
    }

    // returns (idx, distance) with idx ascending; distance is best for that idx.
    void lookup(const char* query, int len, std::vector<PackedHit>& out,
                LookupScratch& scratch, std::vector<PackedHit>& cache_hits,
                DecodeStats* stats) const {
        out.clear();
        if (stats) ++stats->anchor_lookup_calls;
        scratch.ensure(anchors->size());
        scratch.clear_best();
        if (use_hamming) {
            if (hamming_cache.lookup(query, len, cache_hits)) {
                if (stats) {
                    ++stats->anchor_hamming_keys_found;
                    stats->anchor_hamming_hits_returned += cache_hits.size();
                }
                for (auto& hit : cache_hits) scratch.add_hit(hit.idx, hit.distance);
            }
        }
        if (use_deletion) {
            std::string query_string(query, len);
            deletion_signatures(query_string, signature_deletions, scratch.sigbuf);
            scratch.start_seen();
            for (auto& sig : scratch.sigbuf) {
                auto it = deletion_index.find(sig);
                if (it == deletion_index.end()) continue;
                for (int idx : it->second) scratch.mark_seen(idx);
            }
            for (int idx : scratch.seen_touched) {
                if (stats) ++stats->anchor_deletion_dp_calls;
                int d = edit_distance_bounded(query, len, (*anchors)[idx].data(),
                                              static_cast<int>((*anchors)[idx].size()),
                                              max_distance);
                if (d <= max_distance) {
                    if (stats) ++stats->anchor_deletion_dp_accepts;
                    scratch.add_hit(idx, d);
                }
            }
        }
        out.reserve(scratch.touched.size());
        for (int idx : scratch.touched) out.push_back({static_cast<uint16_t>(idx),
                                                       static_cast<uint8_t>(scratch.best[idx]), 0});
        std::sort(out.begin(), out.end(),
                  [](const PackedHit& a, const PackedHit& b) { return a.idx < b.idx; });
    }
};

int cached_edit_for_target_from_hits(bool cache_found, const std::vector<PackedHit>& cache_hits,
                                     const std::string& seq, int start, int length,
                                     const std::vector<std::string>& targets, int target_idx, int maxk,
                                     DecodeStats* stats) {
    const std::string& target = targets[target_idx];
    if (cache_found && length == static_cast<int>(target.size())) {
        auto it = std::lower_bound(cache_hits.begin(), cache_hits.end(), target_idx,
                                   [](const PackedHit& h, int idx) { return h.idx < idx; });
        if (it != cache_hits.end() && it->idx == target_idx) {
            if (stats) ++stats->full_cache_target_hits;
            return it->distance;
        }
    }
    if (stats) {
        ++stats->full_dp_calls;
        if (length == static_cast<int>(target.size())) {
            ++stats->full_same_length_target_miss_dp;
        } else {
            ++stats->full_length_shift_dp;
        }
    }
    int d = edit_distance_bounded(seq.data() + start, length, target.data(),
                                  static_cast<int>(target.size()), maxk);
    if (stats && d <= maxk) ++stats->full_dp_accepts;
    return d;
}

int cached_edit_for_target(const PackedAnswerCache& cache, const std::string& seq, int start, int length,
                           const std::vector<std::string>& targets, int target_idx, int maxk,
                           std::vector<PackedHit>& cache_hits, DecodeStats* stats) {
    if (start < 0 || length <= 0 || start + length > static_cast<int>(seq.size())) return maxk + 1;
    bool cache_found = false;
    if (length == static_cast<int>(targets[target_idx].size())) {
        if (stats) ++stats->full_cache_lookup_calls;
        cache_found = cache.lookup(seq.data() + start, length, cache_hits);
        if (stats && cache_found) ++stats->full_cache_keys_found;
    }
    return cached_edit_for_target_from_hits(cache_found, cache_hits, seq, start, length,
                                            targets, target_idx, maxk, stats);
}

// best_observed_edit_distance: min (distance, abs(length_delta)) over length deltas in [-window, window].
int best_observed_edit_distance(const std::string& seq, int start, int target_idx,
                                const std::vector<std::string>& targets,
                                const PackedAnswerCache& cache,
                                int window, int maxk, std::vector<PackedHit>& cache_hits,
                                DecodeStats* stats) {
    const std::string& target = targets[target_idx];
    int best_distance = -1;
    int best_abs_len_delta = 0;
    for (int ld = -window; ld <= window; ++ld) {
        int length = static_cast<int>(target.size()) + ld;
        if (length <= 0) continue;
        int end = start + length;
        if (start < 0 || end > static_cast<int>(seq.size())) continue;
        int d = cached_edit_for_target(cache, seq, start, length, targets, target_idx,
                                       maxk, cache_hits, stats);
        if (best_distance < 0 ||
            std::make_pair(d, std::abs(ld)) <
                std::make_pair(best_distance, best_abs_len_delta)) {
            best_distance = d;
            best_abs_len_delta = std::abs(ld);
        }
    }
    if (best_distance < 0) return std::max(1, static_cast<int>(target.size()));
    return best_distance;
}

int best_observed_edit_distance_observed(const std::string& seq, int start, int target_idx,
                                         const std::vector<std::string>& targets,
                                         ObservedAnswerCache& cache,
                                         ObservedAnswerCache::LocalCache& local_cache,
                                         ObservedLookupScratch& scratch,
                                         int window, int maxk, std::vector<PackedHit>& observed_hits) {
    const std::string& target = targets[target_idx];
    int best_distance = -1;
    int best_abs_len_delta = 0;
    for (int ld = -window; ld <= window; ++ld) {
        int length = static_cast<int>(target.size()) + ld;
        if (length <= 0) continue;
        int end = start + length;
        if (start < 0 || end > static_cast<int>(seq.size())) continue;
        cache.lookup(seq.data() + start, length, local_cache, scratch, observed_hits);
        int d = hit_distance_for_idx(observed_hits, target_idx, maxk + 1);
        if (best_distance < 0 ||
            std::make_pair(d, std::abs(ld)) <
                std::make_pair(best_distance, best_abs_len_delta)) {
            best_distance = d;
            best_abs_len_delta = std::abs(ld);
        }
    }
    if (best_distance < 0) return std::max(1, static_cast<int>(target.size()));
    return best_distance;
}

int best_observed_edit_distance_prebuilt(const std::string& seq, int start, int target_idx,
                                         const std::vector<std::string>& targets,
                                         const PrebuiltObservedAnswerCache& cache,
                                         const DeletionCandidateIndex& deletion_index,
                                         ObservedLookupScratch& dp_scratch,
                                         DpOnlyMemo& dp_memo,
                                         std::vector<PackedHit>& dp_hits,
                                         int window, int maxk, bool fallback_on_cache_miss,
                                         bool* cache_miss) {
    const std::string& target = targets[target_idx];
    int best_distance = -1;
    int best_abs_len_delta = 0;
    for (int ld = -window; ld <= window; ++ld) {
        int length = static_cast<int>(target.size()) + ld;
        if (length <= 0) continue;
        int end = start + length;
        if (start < 0 || end > static_cast<int>(seq.size())) continue;
        HitSpan hits = cache.lookup_span(seq.data() + start, length);
        if (!hits.found) {
            if (fallback_on_cache_miss) {
                if (cache_miss) *cache_miss = true;
                return maxk + 1;
            }
            lookup_dp_only(seq.data() + start, length, targets, &deletion_index,
                           maxk, dp_scratch, dp_memo, dp_hits);
            hits = {dp_hits.data(), static_cast<uint16_t>(dp_hits.size()), true};
        }
        int d = hit_distance_for_idx(hits, target_idx, maxk + 1);
        if (best_distance < 0 ||
            std::make_pair(d, std::abs(ld)) <
                std::make_pair(best_distance, best_abs_len_delta)) {
            best_distance = d;
            best_abs_len_delta = std::abs(ld);
        }
    }
    if (best_distance < 0) return std::max(1, static_cast<int>(target.size()));
    return best_distance;
}

struct ObservedAlignmentChoice {
    int edit_distance = std::numeric_limits<int>::max() / 4;
    int obs_len = 0;
    int abs_len_delta = 0;
    AlignmentMetrics affine;
};

ObservedAlignmentChoice best_observed_alignment_choice(const std::string& seq, int start,
                                                       const std::string& target,
                                                       int window, int maxk,
                                                       const AffineScoringConfig& affine_cfg,
                                                       bool rank_by_affine,
                                                       bool compute_affine,
                                                       bool normalize_terminal_indels = false) {
    compute_affine = compute_affine || rank_by_affine;
    ObservedAlignmentChoice best;
    bool have_best = false;
    for (int ld = -window; ld <= window; ++ld) {
        const int length = static_cast<int>(target.size()) + ld;
        if (length <= 0) continue;
        const int end = start + length;
        if (start < 0 || end > static_cast<int>(seq.size())) continue;
        TerminalNormalizedEdit normalized;
        if (normalize_terminal_indels) {
            normalized = terminal_normalized_edit_distance(
                seq.data() + start, length, target.data(),
                static_cast<int>(target.size()), maxk);
        } else {
            normalized.edit_distance = edit_distance_bounded(
                seq.data() + start, length, target.data(),
                static_cast<int>(target.size()), maxk);
            normalized.terminal_delta = std::abs(ld);
        }
        const int edit = normalized.edit_distance;
        if (edit > maxk) continue;

        ObservedAlignmentChoice candidate;
        candidate.edit_distance = edit;
        candidate.obs_len = length;
        candidate.abs_len_delta = normalize_terminal_indels
            ? normalized.terminal_delta
            : std::abs(ld);
        if (compute_affine) {
            candidate.affine = affine_alignment_metrics(seq.data() + start, length,
                                                        target.data(),
                                                        static_cast<int>(target.size()),
                                                        affine_cfg);
        }

        bool better = false;
        if (!have_best) {
            better = true;
        } else if (rank_by_affine) {
            const long key[3] = {candidate.affine.cost, candidate.edit_distance,
                                 candidate.abs_len_delta};
            const long best_key[3] = {best.affine.cost, best.edit_distance,
                                      best.abs_len_delta};
            for (int i = 0; i < 3; ++i) {
                if (key[i] != best_key[i]) {
                    better = key[i] < best_key[i];
                    break;
                }
            }
        } else if (std::make_pair(candidate.edit_distance, candidate.abs_len_delta) <
                   std::make_pair(best.edit_distance, best.abs_len_delta)) {
            better = true;
        }
        if (better) {
            have_best = true;
            best = candidate;
        }
    }
    if (!have_best) {
        best.edit_distance = std::max(1, static_cast<int>(target.size()));
        best.obs_len = static_cast<int>(target.size());
        best.abs_len_delta = 0;
    }
    return best;
}

struct Config {
    std::vector<std::string> bc1_oligos, bc2_oligos;
    std::vector<std::string> bc1_anchors, bc2_anchors;
    std::vector<int> bc1_len, bc2_len;
    int anchor_len = 14;      // bc1/bc2 segment length
    int bc2_seg_len = 14;
    int full_start_min = 8, full_start_max = 12;
    int bc1_anchor_hamming = 2, bc2_anchor_hamming = 2;
    int bc1_full_hamming = 2, bc2_full_hamming = 2;
    int bc1_anchor_offset_window = 2, bc2_offset_window = 2;
    int grid_rows = 3350, grid_cols = 3350;
    AnchorLookup bc1_lookup, bc2_lookup;
    PackedAnswerCache bc1_full_lookup, bc2_full_lookup;
    PackedAnswerCache bc1_direct_best_edit_lookup, bc2_direct_best_edit_lookup;
    PackedAnswerCache bc1_tiered_exact_lookup, bc2_tiered_exact_lookup;
    PackedAnswerCache bc1_tiered_e1_lookup, bc2_tiered_e1_lookup;
    PackedAnswerCache bc1_tiered_fixed_lookup, bc2_tiered_fixed_lookup;
    PackedAnswerCache bc1_tiered_h2_lookup, bc2_tiered_h2_lookup;
    bool observed_cache_experiment = false;
    bool observed_cache_prebuild = false;
    bool observed_cache_prebuild_fallback = false;
    bool observed_cache_use_e01 = true;
    long observed_cache_prebuild_sample_reads = 0;
    long decode_skip_reads = 0;
    long decode_read_limit = 0;
    PackedAnswerCache bc1_anchor_e01, bc2_anchor_e01, bc1_full_e01, bc2_full_e01;
    DeletionCandidateIndex bc1_anchor_deletion2, bc2_anchor_deletion2, bc1_full_deletion2, bc2_full_deletion2;
    ObservedAnswerCache bc1_anchor_observed, bc2_anchor_observed, bc1_full_observed, bc2_full_observed;
    PrebuiltObservedAnswerCache bc1_anchor_prebuilt, bc2_anchor_prebuilt, bc1_full_prebuilt, bc2_full_prebuilt;
    bool direct_full_decode = false;
    bool direct_full_exact_lengths = false;
    bool direct_best_edit_decode = false;
    bool direct_tiered_decode = false;
    bool direct_tiered_fixed_decode = false;
    bool direct_tiered_h2_decode = false;
    bool frozen_support_consensus = false;
    double frozen_support_min_odds = 1.0;
    std::vector<uint64_t> bc1_frozen_h0;
    std::vector<uint64_t> bc2_frozen_h0;
    std::vector<uint64_t> bc1_frozen_exposure;
    std::vector<uint64_t> bc2_frozen_exposure;
    bool direct_best_non_acgt_dp_fallback = false;
    DirectBestResolutionMode direct_best_resolution = DirectBestResolutionMode::StrictBoundary;
    bool overlap_dp_length_guard = true;
    int overlap_dp_bc1_min_obs_len = 13;
    int overlap_dp_bc2_min_obs_len = 13;
    int overlap_dp_bc1_max_len_delta = 2;
    int overlap_dp_bc2_max_len_delta = 2;
    bool overlap_dp_anchor_gate_only = false;
    bool overlap_dp_offset_include_terminal_delta = false;
    bool overlap_dp_use_edit_hash = true;
    std::string non_acgt_reads_out;
    long progress_log_every_batches = 0;
    CandidateScoreMode score_mode = CandidateScoreMode::LegacyEdit;
    AffineScoringConfig affine;
    bool emit_decode_annotations = false;
    std::string sr_oracle_tsv;
    std::string candidate_audit_out;
    std::string candidate_preserving_out;
    // Optional research sidecar containing every assignment-eligible H0--H4
    // source candidate found by the direct tiered H2 decoder.  This is kept
    // separate from best_candidate_traces so enabling it cannot change the
    // frozen compatibility decision.
    std::string latent_candidate_out;
    bool suppress_bc1_overlong_prefix_e1 = false;
    bool suppress_bc2_short_terminal_e2 = true;
    int audit_full_start_padding = 2;
    int audit_bc2_start_padding = 2;
    int audit_edit_max = 4;
};

struct DirectHalfCandidate {
    uint16_t idx;
    uint8_t distance;
    uint8_t start;
    uint8_t full_start;
    uint8_t boundary;
    uint8_t offset;

    DirectHalfCandidate(uint16_t index = 0, uint8_t editDistance = 0,
                        uint8_t observedStart = 0, uint8_t fullStart = 0,
                        uint8_t observedBoundary = 0, uint8_t observedOffset = 0)
        : idx(index), distance(editDistance), start(observedStart),
          full_start(fullStart), boundary(observedBoundary), offset(observedOffset) {}
};

struct HalfSpanChoice {
    bool found = false;
    int edit_distance = std::numeric_limits<int>::max() / 4;
    int raw_edit_distance = std::numeric_limits<int>::max() / 4;
    int terminal_delta = 0;
    int start = -1;
    int obs_len = -1;
    AlignmentMetrics affine;
};

struct DirectDpHalfCandidate {
    uint16_t idx;
    HalfSpanChoice choice;

    DirectDpHalfCandidate(uint16_t index = 0,
                          const HalfSpanChoice& spanChoice = HalfSpanChoice())
        : idx(index), choice(spanChoice) {}
};

constexpr int kDirectBestBc1QueryLengths[] = {13, 14, 15, 16, 17, 18};
constexpr int kDirectBestBc2QueryLengths[] = {12, 13, 14, 15, 16, 17};
constexpr int kDirectBestBc1MinQueryLength = 13;
constexpr int kDirectBestBc1MaxQueryLength = 18;
constexpr int kDirectBestBc2MinQueryLength = 12;
constexpr int kDirectBestBc2MaxQueryLength = 17;
constexpr int kDirectBestMaxBc2QueryLength = 17;

void dedupe_direct_half_candidates(std::vector<DirectHalfCandidate>& candidates) {
    std::sort(candidates.begin(), candidates.end(),
              [](const DirectHalfCandidate& a, const DirectHalfCandidate& b) {
                  return std::make_tuple(a.idx, a.start, a.full_start, a.boundary,
                                         a.offset, a.distance) <
                         std::make_tuple(b.idx, b.start, b.full_start, b.boundary,
                                         b.offset, b.distance);
              });
    size_t out = 0;
    for (size_t i = 0; i < candidates.size();) {
        DirectHalfCandidate best = candidates[i];
        ++i;
        while (i < candidates.size() &&
               std::make_tuple(candidates[i].idx, candidates[i].start,
                               candidates[i].full_start, candidates[i].boundary,
                               candidates[i].offset) ==
                   std::make_tuple(best.idx, best.start, best.full_start,
                                   best.boundary, best.offset)) {
            if (candidates[i].distance < best.distance) best = candidates[i];
            ++i;
        }
        candidates[out++] = best;
    }
    candidates.resize(out);
}

struct ObservedDecodeCacheLocals {
    ObservedAnswerCache::LocalCache bc1_anchor;
    ObservedAnswerCache::LocalCache bc2_anchor;
    ObservedAnswerCache::LocalCache bc1_full;
    ObservedAnswerCache::LocalCache bc2_full;
    ObservedLookupScratch scratch;
    DpOnlyMemo bc1_anchor_dp_miss;
    DpOnlyMemo bc2_anchor_dp_miss;
    DpOnlyMemo bc1_full_dp_miss;
    DpOnlyMemo bc2_full_dp_miss;
    std::vector<PackedHit> bc1_anchor_dp_hits;
    std::vector<PackedHit> bc2_anchor_dp_hits;
    std::vector<PackedHit> bc1_full_dp_hits;
    std::vector<PackedHit> bc2_full_dp_hits;
    std::vector<DirectHalfCandidate> direct_bc1;
    std::vector<DirectHalfCandidate> direct_bc2;
    std::vector<DirectDpHalfCandidate> overlap_bc1;
    std::vector<DirectDpHalfCandidate> overlap_bc2_by_boundary[kMaxLen];

    void flush(Config& cfg) {
        cfg.bc1_anchor_observed.flush(bc1_anchor);
        cfg.bc2_anchor_observed.flush(bc2_anchor);
        cfg.bc1_full_observed.flush(bc1_full);
        cfg.bc2_full_observed.flush(bc2_full);
    }

    DpOnlyMemoStats dp_miss_stats() const {
        DpOnlyMemoStats out;
        out.add(bc1_anchor_dp_miss.stats);
        out.add(bc2_anchor_dp_miss.stats);
        out.add(bc1_full_dp_miss.stats);
        out.add(bc2_full_dp_miss.stats);
        return out;
    }
};

bool is_acgt(const std::string& s, int start, int len) {
    for (int i = start; i < start + len; ++i) {
        char c = s[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    }
    return true;
}

bool has_non_acgt(const std::string& s, int start, int len) {
    if (start < 0) start = 0;
    const int end = std::min<int>(static_cast<int>(s.size()), start + len);
    for (int i = start; i < end; ++i) {
        char c = s[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return true;
    }
    return false;
}

bool has_non_acgt(const char* s, int len) {
    for (int i = 0; i < len; ++i) {
        char c = s[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return true;
    }
    return false;
}

bool use_full_bc1_span_anchor(const Config& cfg, int bc1_obs_len, int col2) {
    if (col2 < 0 || col2 >= static_cast<int>(cfg.bc1_len.size())) return false;
    return bc1_obs_len > 0 &&
           bc1_obs_len <= cfg.anchor_len &&
           bc1_obs_len < cfg.bc1_len[col2];
}

bool use_full_bc2_span_anchor(const Config& cfg, int bc2_obs_len, int row2) {
    if (row2 < 0 || row2 >= static_cast<int>(cfg.bc2_len.size())) return false;
    return bc2_obs_len > 0 && bc2_obs_len != cfg.bc2_len[row2];
}

bool should_suppress_bc2_short_terminal_e2(const Config& cfg, int row2,
                                           int bc2_edit, int bc2_obs_len) {
    if (!cfg.suppress_bc2_short_terminal_e2) return false;
    if (row2 < 0 || row2 >= static_cast<int>(cfg.bc2_len.size())) return false;
    if (cfg.bc2_full_hamming < 2 || bc2_edit < cfg.bc2_full_hamming) return false;
    const int target_len = cfg.bc2_len[row2];
    return target_len == cfg.bc2_seg_len &&
           bc2_obs_len <= target_len - cfg.bc2_full_hamming;
}

struct DecodeWindowClassification {
    uint64_t n_mask;
    uint32_t n_count;
    bool unsupported;

    DecodeWindowClassification()
        : n_mask(0), n_count(0), unsupported(false) {}
};

DecodeWindowClassification classify_direct_best_decode_window(
    const std::string& seq, const Config& cfg) {
    DecodeWindowClassification classification;
    const int start = std::max(0, cfg.full_start_min);
    const int max_bc2_start =
        cfg.full_start_max + cfg.anchor_len + 2 + cfg.bc1_full_hamming;
    const int end = std::min<int>(
        static_cast<int>(seq.size()), max_bc2_start + kDirectBestMaxBc2QueryLength);
    for (int index = start; index < end; ++index) {
        const char base = static_cast<char>(
            std::toupper(static_cast<unsigned char>(seq[static_cast<size_t>(index)])));
        if (base == 'N') {
            ++classification.n_count;
            if (index < 64) classification.n_mask |= uint64_t(1) << index;
        } else if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            classification.unsupported = true;
        }
    }
    return classification;
}

bool direct_best_decode_window_has_non_acgt(const std::string& seq, const Config& cfg) {
    const DecodeWindowClassification classification =
        classify_direct_best_decode_window(seq, cfg);
    return classification.n_count != 0 || classification.unsupported;
}

bool has_overlong_exact_prefix_overlap(const char* observed, int observed_len,
                                       const std::string& target,
                                       int min_overlap, int edit_distance) {
    if (edit_distance <= 0 || observed_len <= min_overlap) return false;
    const int target_len = static_cast<int>(target.size());
    const int max_overlap = std::min(observed_len - 1, target_len);
    for (int overlap = max_overlap; overlap >= min_overlap; --overlap) {
        if (std::memcmp(observed, target.data(), static_cast<size_t>(overlap)) == 0) {
            return true;
        }
    }
    return false;
}

bool half_span_choice_less(const HalfSpanChoice& a, const HalfSpanChoice& b,
                           int expected_start, int expected_len) {
    if (!a.found) return false;
    if (!b.found) return true;
    return std::make_tuple(a.edit_distance, a.terminal_delta,
                           std::abs(a.start - expected_start),
                           std::abs(a.obs_len - expected_len),
                           a.raw_edit_distance, a.start, a.obs_len) <
           std::make_tuple(b.edit_distance, b.terminal_delta,
                           std::abs(b.start - expected_start),
                           std::abs(b.obs_len - expected_len),
                           b.raw_edit_distance, b.start, b.obs_len);
}

struct WindowEditCell {
    int cost = std::numeric_limits<int>::max() / 4;
    int start = 0;
};

bool window_edit_cell_less(const WindowEditCell& a, const WindowEditCell& b) {
    if (a.cost != b.cost) return a.cost < b.cost;
    return a.start < b.start;
}

WindowEditCell best_window_edit_cell(WindowEditCell a, WindowEditCell b) {
    return window_edit_cell_less(b, a) ? b : a;
}

HalfSpanChoice best_half_window_dp_choice(const std::string& seq,
                                          int window_start,
                                          int window_end,
                                          const std::string& target,
                                          int max_edit,
                                          int expected_start,
                                          int expected_len,
                                          int min_obs_len,
                                          int max_len_delta,
                                          bool free_observed_prefix,
                                          const AffineScoringConfig& affine_cfg,
                                          bool compute_affine) {
    HalfSpanChoice best;
    const int slen = static_cast<int>(seq.size());
    window_start = std::max(0, window_start);
    window_end = std::min(slen, window_end);
    if (window_start >= window_end) return best;
    const int observed_len = window_end - window_start;
    const int target_len = static_cast<int>(target.size());
    if (observed_len > kMaxLen || target_len > kMaxLen) return best;

    WindowEditCell prev[kMaxLen + 1], cur[kMaxLen + 1];
    prev[0].cost = 0;
    prev[0].start = 0;
    for (int j = 1; j <= target_len; ++j) {
        prev[j].cost = j;
        prev[j].start = 0;
    }

    auto maybe_store_best = [&](const WindowEditCell& cell, int end_offset) {
        if (cell.cost > max_edit) return;
        const int start_offset = cell.start;
        const int obs_len = end_offset - start_offset;
        if (obs_len <= 0) return;
        if (min_obs_len >= 0 && obs_len < min_obs_len) return;
        if (max_len_delta >= 0 && std::abs(obs_len - expected_len) > max_len_delta) {
            return;
        }
        const int start = window_start + start_offset;
        HalfSpanChoice candidate;
        candidate.found = true;
        candidate.edit_distance = cell.cost;
        candidate.raw_edit_distance = cell.cost;
        candidate.terminal_delta = std::abs(obs_len - target_len);
        candidate.start = start;
        candidate.obs_len = obs_len;
        if (compute_affine) {
            candidate.affine = affine_alignment_metrics(
                seq.data() + start, obs_len, target.data(), target_len, affine_cfg);
        }
        if (half_span_choice_less(candidate, best, expected_start, expected_len)) {
            best = candidate;
        }
    };

    maybe_store_best(prev[target_len], 0);
    for (int i = 1; i <= observed_len; ++i) {
        if (free_observed_prefix) {
            cur[0].cost = 0;
            cur[0].start = i;
        } else {
            cur[0].cost = i;
            cur[0].start = 0;
        }
        const char observed_base = seq[window_start + i - 1];
        for (int j = 1; j <= target_len; ++j) {
            WindowEditCell subst = prev[j - 1];
            subst.cost += (observed_base == target[j - 1]) ? 0 : 1;
            WindowEditCell ins = prev[j];
            ++ins.cost;
            WindowEditCell del = cur[j - 1];
            ++del.cost;
            cur[j] = best_window_edit_cell(best_window_edit_cell(subst, ins), del);
        }
        maybe_store_best(cur[target_len], i);
        for (int j = 0; j <= target_len; ++j) prev[j] = cur[j];
    }
    return best;
}

struct OracleTarget {
    bool resolved = false;
    std::string cb;
    int row2 = -1;
    int col2 = -1;
    std::string status;
};

struct OracleTargetAudit {
    bool possible = false;
    int min_total_edit = -1;
    int bc1_edit = -1;
    int bc2_edit = -1;
    int anchor_score = -1;
    int offset_sum = -1;
    int full_start = -1;
    int bc2_start = -1;
    int bc1_obs_len = -1;
    int bc2_obs_len = -1;
    bool observed_non_acgt = false;
    std::string anchor_failure_class;
    int anchor_bc2_edit = -1;
    int anchor_current_bc1_edit = -1;
    int anchor_current_d1 = 0;
    int anchor_current_d2 = 0;
    int anchor_any_bc1_edit = -1;
    int anchor_any_d1 = 0;
    int anchor_any_d2 = 0;
};

struct AnchorLimitProbe {
    bool evaluated = false;
    int bc2_anchor_edit = -1;
    int best_current_bc1_edit = -1;
    int best_current_d1 = 0;
    int best_current_d2 = 0;
    int best_any_bc1_edit = -1;
    int best_any_d1 = 0;
    int best_any_d2 = 0;
    bool current_offsets_available = false;
    bool any_offsets_available = false;
    std::string failure_class;
};

struct CandidateTrace {
    int row2;
    int col2;
    int full_start;
    int bc1_edit;
    int bc2_edit;
    int anchor_score;
    int offset_sum;
    int bc1_obs_len;
    int bc2_obs_len;

    CandidateTrace(int row = -1, int column = -1, int fullStart = -1,
                   int bc1Edit = -1, int bc2Edit = -1,
                   int anchorScore = -1, int offsetSum = -1,
                   int bc1ObservedLength = -1, int bc2ObservedLength = -1)
        : row2(row), col2(column), full_start(fullStart), bc1_edit(bc1Edit),
          bc2_edit(bc2Edit), anchor_score(anchorScore), offset_sum(offsetSum),
          bc1_obs_len(bc1ObservedLength), bc2_obs_len(bc2ObservedLength) {}
};

bool same_candidate_coord(const CandidateTrace& trace,
                          int row2, int col2, int full_start) {
    return trace.row2 == row2 && trace.col2 == col2 &&
           trace.full_start == full_start;
}

void append_candidate_trace(std::vector<CandidateTrace>& traces,
                            int row2, int col2, int full_start,
                            int bc1_edit, int bc2_edit,
                            int anchor_score, int offset_sum,
                            int bc1_obs_len, int bc2_obs_len) {
    // Candidate-preserving mode uses this same trace vector as its transparent
    // sequence-only ledger.  Do not cap it: silent truncation would turn an
    // ambiguous read into an artificially confident candidate set.
    for (const CandidateTrace& trace : traces) {
        if (same_candidate_coord(trace, row2, col2, full_start)) return;
    }
    traces.push_back({row2, col2, full_start, bc1_edit, bc2_edit,
                      anchor_score, offset_sum, bc1_obs_len, bc2_obs_len});
}

// Decode one read -> (assigned, row2, col2, full_start). Mirrors decode_record candidate gen
// + choose_unique_best (min score, unique coord).
struct Decoded {
    bool assigned = false;
    bool sequestered_non_acgt = false;
    bool non_acgt_dp_checked = false;
    uint32_t non_acgt_null_queries = 0;
    int row2 = 0;
    int col2 = 0;
    int full_start = -1;
    bool have_best = false;
    bool best_unique = true;
    uint64_t candidate_pairs = 0;
    uint32_t best_tie_candidates = 0;
    long score_key[4] = {-1, -1, -1, -1};
    std::vector<CandidateTrace> best_candidate_traces;
    std::vector<CandidateTrace> latent_candidate_traces;
    bool frozen_support_attempted = false;
    bool frozen_support_resolved = false;
    std::string frozen_support_status;
    double frozen_h0_odds = 0.0;
    double frozen_exposure_odds = 0.0;
    bool have_suppressed_best = false;
    long suppressed_score_key[4] = {-1, -1, -1, -1};
    int bc1_edit = -1;
    int bc2_edit = -1;
    int anchor_score = -1;
    int offset_sum = -1;
    int bc1_obs_len = -1;
    int bc2_obs_len = -1;
    AlignmentMetrics bc1_affine;
    AlignmentMetrics bc2_affine;
    bool oracle_candidate_present = false;
    bool oracle_bc1_seed_seen = false;
    bool oracle_bc2_seed_seen = false;
    bool oracle_pair_seen = false;
    bool oracle_pair_bc1_pass = false;
    bool oracle_pair_bc2_pass = false;
    bool oracle_pair_anchor_pass = false;
    long oracle_score_key[4] = {-1, -1, -1, -1};
    int oracle_full_start = -1;
    int oracle_bc1_edit = -1;
    int oracle_bc2_edit = -1;
    int oracle_anchor_score = -1;
    int oracle_offset_sum = -1;
    int oracle_bc1_obs_len = -1;
    int oracle_bc2_obs_len = -1;
    AlignmentMetrics oracle_bc1_affine;
    AlignmentMetrics oracle_bc2_affine;
};

struct FullCandidate {
    uint16_t idx;
    uint8_t anchor_distance;
    uint8_t full_distance;
    uint8_t obs_len;
    AlignmentMetrics affine;

    FullCandidate(uint16_t index = 0, uint8_t anchorDistance = 0,
                  uint8_t fullDistance = 0, uint8_t observedLength = 0,
                  const AlignmentMetrics& affineMetrics = AlignmentMetrics())
        : idx(index), anchor_distance(anchorDistance), full_distance(fullDistance),
          obs_len(observedLength), affine(affineMetrics) {}
};

void make_candidate_score_key(const Config& cfg,
                              int bc1_edit, int bc2_edit,
                              const AlignmentMetrics& bc1_affine,
                              const AlignmentMetrics& bc2_affine,
                              int anchor_score, int offset_sum,
                              int full_start, long key[4]) {
    key[0] = (cfg.score_mode == CandidateScoreMode::AffineGap)
        ? static_cast<long>(bc1_affine.cost) + static_cast<long>(bc2_affine.cost)
        : static_cast<long>(bc1_edit) + static_cast<long>(bc2_edit);
    key[1] = (cfg.direct_best_resolution == DirectBestResolutionMode::OverlapSpanDp &&
              cfg.overlap_dp_anchor_gate_only)
        ? 0
        : anchor_score;
    key[2] = offset_sum;
    key[3] = full_start;
}

int compare_score_keys(const long a[4], const long b[4]) {
    for (int i = 0; i < 4; ++i) {
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    }
    return 0;
}

void store_candidate_trace(Decoded& out, const long key[4],
                           int row2, int col2, int full_start,
                           int bc1_edit, int bc2_edit,
                           const AlignmentMetrics& bc1_affine,
                           const AlignmentMetrics& bc2_affine,
                           int anchor_score, int offset_sum,
                           int bc1_obs_len, int bc2_obs_len) {
    std::memcpy(out.score_key, key, sizeof(out.score_key));
    out.row2 = row2;
    out.col2 = col2;
    out.full_start = full_start;
    out.bc1_edit = bc1_edit;
    out.bc2_edit = bc2_edit;
    out.bc1_affine = bc1_affine;
    out.bc2_affine = bc2_affine;
    out.anchor_score = anchor_score;
    out.offset_sum = offset_sum;
    out.bc1_obs_len = bc1_obs_len;
    out.bc2_obs_len = bc2_obs_len;
}

void consider_candidate(Decoded& out, const Config& cfg,
                        int bc1_edit, int bc2_edit,
                        const AlignmentMetrics& bc1_affine,
                        const AlignmentMetrics& bc2_affine,
                        int anchor_score, int offset_sum,
                        int full_start, int row2, int col2,
                        int bc1_obs_len, int bc2_obs_len,
                        const OracleTarget* oracle = nullptr,
                        bool assignment_eligible = true) {
    ++out.candidate_pairs;
    long key[4] = {0, 0, 0, 0};
    make_candidate_score_key(cfg, bc1_edit, bc2_edit, bc1_affine, bc2_affine,
                             anchor_score, offset_sum, full_start, key);
    if (oracle && oracle->resolved && row2 == oracle->row2 && col2 == oracle->col2) {
        if (!out.oracle_candidate_present ||
            compare_score_keys(key, out.oracle_score_key) < 0) {
            out.oracle_candidate_present = true;
            std::memcpy(out.oracle_score_key, key, sizeof(out.oracle_score_key));
            out.oracle_full_start = full_start;
            out.oracle_bc1_edit = bc1_edit;
            out.oracle_bc2_edit = bc2_edit;
            out.oracle_bc1_affine = bc1_affine;
            out.oracle_bc2_affine = bc2_affine;
            out.oracle_anchor_score = anchor_score;
            out.oracle_offset_sum = offset_sum;
            out.oracle_bc1_obs_len = bc1_obs_len;
            out.oracle_bc2_obs_len = bc2_obs_len;
        }
    }
    if (!assignment_eligible) {
        if (!out.have_suppressed_best ||
            compare_score_keys(key, out.suppressed_score_key) < 0) {
            out.have_suppressed_best = true;
            std::memcpy(out.suppressed_score_key, key, sizeof(out.suppressed_score_key));
        }
        return;
    }
    if (!out.have_best) {
        out.have_best = true;
        out.best_unique = true;
        out.best_tie_candidates = 0;
        store_candidate_trace(out, key, row2, col2, full_start, bc1_edit, bc2_edit,
                              bc1_affine, bc2_affine, anchor_score, offset_sum,
                              bc1_obs_len, bc2_obs_len);
        out.best_candidate_traces.clear();
        append_candidate_trace(out.best_candidate_traces, row2, col2, full_start,
                               bc1_edit, bc2_edit, anchor_score, offset_sum,
                               bc1_obs_len, bc2_obs_len);
        return;
    }

    const int cmp = compare_score_keys(key, out.score_key);
    if (cmp < 0) {
        out.best_unique = true;
        out.best_tie_candidates = 0;
        store_candidate_trace(out, key, row2, col2, full_start, bc1_edit, bc2_edit,
                              bc1_affine, bc2_affine, anchor_score, offset_sum,
                              bc1_obs_len, bc2_obs_len);
        out.best_candidate_traces.clear();
        append_candidate_trace(out.best_candidate_traces, row2, col2, full_start,
                               bc1_edit, bc2_edit, anchor_score, offset_sum,
                               bc1_obs_len, bc2_obs_len);
    } else if (cmp == 0) {
        ++out.best_tie_candidates;
        if (row2 != out.row2 || col2 != out.col2) out.best_unique = false;
        append_candidate_trace(out.best_candidate_traces, row2, col2, full_start,
                               bc1_edit, bc2_edit, anchor_score, offset_sum,
                               bc1_obs_len, bc2_obs_len);
    }
}

void finalize_decoded_assignment(Decoded& out) {
    if (out.have_best && out.best_unique &&
        (!out.have_suppressed_best ||
         compare_score_keys(out.score_key, out.suppressed_score_key) < 0)) {
        out.assigned = true;
    }
}

bool direct_best_anchor_score_for_pair(const std::string& seq, const Config& cfg,
                                       int full_start, int bc2_boundary,
                                       int bc2_align_start,
                                       int bc2_obs_len, int row2, int col2,
                                       int& anchor_score, int& offset_sum,
                                       bool allow_non_acgt_anchor = false,
                                       bool include_terminal_delta_in_offset = true) {
    const int slen = static_cast<int>(seq.size());
    if (row2 < 0 || row2 >= static_cast<int>(cfg.bc2_anchors.size()) ||
        col2 < 0 || col2 >= static_cast<int>(cfg.bc1_anchors.size())) {
        return false;
    }
    if (bc2_boundary < 0 || bc2_align_start < 0) return false;
    const bool use_full_bc2 = use_full_bc2_span_anchor(cfg, bc2_obs_len, row2);
    int bc2_anchor_d = cfg.bc2_anchor_hamming + 1;
    const int bc2_start_offset = std::abs(bc2_align_start - bc2_boundary);
    int bc2_span_offset = bc2_start_offset;
    if (use_full_bc2) {
        if (bc2_align_start + bc2_obs_len > slen) return false;
        if (!allow_non_acgt_anchor &&
            !is_acgt(seq, bc2_align_start, bc2_obs_len)) {
            return false;
        }
        bc2_anchor_d = edit_distance_bounded(
            seq.data() + bc2_align_start, bc2_obs_len,
            cfg.bc2_oligos[row2].data(), static_cast<int>(cfg.bc2_oligos[row2].size()),
            cfg.bc2_anchor_hamming);
        if (include_terminal_delta_in_offset) {
            bc2_span_offset += std::abs(bc2_obs_len - cfg.bc2_len[row2]);
        }
    } else {
        if (bc2_align_start + cfg.bc2_seg_len > slen) return false;
        if (!allow_non_acgt_anchor &&
            !is_acgt(seq, bc2_align_start, cfg.bc2_seg_len)) {
            return false;
        }
        bc2_anchor_d = edit_distance_bounded(
            seq.data() + bc2_align_start, cfg.bc2_seg_len,
            cfg.bc2_anchors[row2].data(), static_cast<int>(cfg.bc2_anchors[row2].size()),
            cfg.bc2_anchor_hamming);
    }
    if (bc2_anchor_d > cfg.bc2_anchor_hamming) return false;

    const int bc1_obs_len = bc2_boundary - full_start;
    if (use_full_bc1_span_anchor(cfg, bc1_obs_len, col2)) {
        if (full_start < 0 || full_start + bc1_obs_len > slen) return false;
        if (!allow_non_acgt_anchor && !is_acgt(seq, full_start, bc1_obs_len)) {
            return false;
        }
        const int bc1_anchor_d = edit_distance_bounded(
            seq.data() + full_start, bc1_obs_len,
            cfg.bc1_oligos[col2].data(), static_cast<int>(cfg.bc1_oligos[col2].size()),
            cfg.bc1_anchor_hamming);
        if (bc1_anchor_d > cfg.bc1_anchor_hamming) return false;
        const int boundary_delta = bc2_boundary - (full_start + cfg.bc1_len[col2]);
        if (std::abs(boundary_delta) > cfg.bc2_offset_window) return false;
        anchor_score = bc1_anchor_d + bc2_anchor_d;
        offset_sum = (include_terminal_delta_in_offset ? std::abs(boundary_delta) : 0) +
            bc2_span_offset;
        return true;
    }

    const int extra = cfg.bc1_len[col2] - cfg.anchor_len;
    const int bc1_anchor_expected_start = full_start + extra;
    bool have_best_anchor = false;
    int best_anchor_score = 0;
    int best_offset_sum = 0;
    const int min_bc1_anchor_obs_len = std::max(1, cfg.anchor_len - cfg.bc1_anchor_hamming);
    const int max_bc1_anchor_obs_len = cfg.anchor_len + cfg.bc1_anchor_hamming;
    for (int bc1_anchor_obs_len = min_bc1_anchor_obs_len;
         bc1_anchor_obs_len <= max_bc1_anchor_obs_len; ++bc1_anchor_obs_len) {
        for (int d1 = -cfg.bc1_anchor_offset_window; d1 <= cfg.bc1_anchor_offset_window; ++d1) {
            const int bc1_anchor_start = bc1_anchor_expected_start + d1;
            const int bc1_anchor_end = bc1_anchor_start + bc1_anchor_obs_len;
            if (bc1_anchor_start < full_start || bc1_anchor_end > slen) continue;
            if (!allow_non_acgt_anchor &&
                !is_acgt(seq, bc1_anchor_start, bc1_anchor_obs_len)) {
                continue;
            }
            const int d2 = bc2_boundary - bc1_anchor_end;
            if (d2 < -cfg.bc2_offset_window || d2 > cfg.bc2_offset_window) continue;
            const int bc1_anchor_d = edit_distance_bounded(
                seq.data() + bc1_anchor_start, bc1_anchor_obs_len,
                cfg.bc1_anchors[col2].data(),
                static_cast<int>(cfg.bc1_anchors[col2].size()),
                cfg.bc1_anchor_hamming);
            if (bc1_anchor_d > cfg.bc1_anchor_hamming) continue;

            const int candidate_anchor_score = bc1_anchor_d + bc2_anchor_d;
            const int candidate_offset_sum = std::abs(d1) + std::abs(d2) +
                (include_terminal_delta_in_offset
                     ? std::abs(bc1_anchor_obs_len - cfg.anchor_len)
                     : 0) +
                bc2_span_offset;
            if (!have_best_anchor ||
                std::make_pair(candidate_anchor_score, candidate_offset_sum) <
                    std::make_pair(best_anchor_score, best_offset_sum)) {
                have_best_anchor = true;
                best_anchor_score = candidate_anchor_score;
                best_offset_sum = candidate_offset_sum;
            }
        }
    }
    if (!have_best_anchor) return false;
    anchor_score = best_anchor_score;
    offset_sum = best_offset_sum;
    return true;
}

bool anchor_probe_tuple_less(int edit, int d1, int d2,
                             int best_edit, int best_d1, int best_d2) {
    return std::make_tuple(edit, std::abs(d1) + std::abs(d2), std::abs(d1),
                           std::abs(d2), d1, d2) <
           std::make_tuple(best_edit, std::abs(best_d1) + std::abs(best_d2),
                           std::abs(best_d1), std::abs(best_d2), best_d1, best_d2);
}

AnchorLimitProbe probe_anchor_limits_for_pair(const std::string& seq, const Config& cfg,
                                              int full_start, int bc2_boundary,
                                              int bc2_align_start,
                                              int bc2_obs_len, int row2, int col2) {
    AnchorLimitProbe probe;
    const int slen = static_cast<int>(seq.size());
    if (row2 < 0 || row2 >= static_cast<int>(cfg.bc2_anchors.size()) ||
        col2 < 0 || col2 >= static_cast<int>(cfg.bc1_anchors.size()) ||
        bc2_boundary < 0 || bc2_align_start < 0) {
        probe.failure_class = "anchor_probe_out_of_bounds";
        return probe;
    }
    probe.evaluated = true;
    if (use_full_bc2_span_anchor(cfg, bc2_obs_len, row2)) {
        if (bc2_align_start + bc2_obs_len > slen) {
            probe.failure_class = "anchor_probe_out_of_bounds";
            return probe;
        }
        probe.bc2_anchor_edit = edit_distance_bounded(
            seq.data() + bc2_align_start, bc2_obs_len,
            cfg.bc2_oligos[row2].data(), static_cast<int>(cfg.bc2_oligos[row2].size()),
            std::max(bc2_obs_len, cfg.bc2_len[row2]));
    } else {
        if (bc2_align_start + cfg.bc2_seg_len > slen) {
            probe.failure_class = "anchor_probe_out_of_bounds";
            return probe;
        }
        probe.bc2_anchor_edit = edit_distance_bounded(
            seq.data() + bc2_align_start, cfg.bc2_seg_len,
            cfg.bc2_anchors[row2].data(), static_cast<int>(cfg.bc2_anchors[row2].size()),
            cfg.bc2_seg_len);
    }

    const int bc1_obs_len = bc2_boundary - full_start;
    if (use_full_bc1_span_anchor(cfg, bc1_obs_len, col2)) {
        if (full_start < 0 || full_start + bc1_obs_len > slen) {
            probe.failure_class = "anchor_probe_out_of_bounds";
            return probe;
        }
        const int bc1_edit = edit_distance_bounded(
            seq.data() + full_start, bc1_obs_len,
            cfg.bc1_oligos[col2].data(), static_cast<int>(cfg.bc1_oligos[col2].size()),
            cfg.bc1_len[col2]);
        const int boundary_delta = bc2_boundary - (full_start + cfg.bc1_len[col2]);
        probe.any_offsets_available = true;
        probe.best_any_bc1_edit = bc1_edit;
        probe.best_any_d1 = 0;
        probe.best_any_d2 = boundary_delta;
        if (std::abs(boundary_delta) <= cfg.bc2_offset_window) {
            probe.current_offsets_available = true;
            probe.best_current_bc1_edit = bc1_edit;
            probe.best_current_d1 = 0;
            probe.best_current_d2 = boundary_delta;
        }
        if (probe.bc2_anchor_edit > cfg.bc2_anchor_hamming) {
            probe.failure_class = "bc2_anchor_edit";
        } else if (!probe.current_offsets_available) {
            probe.failure_class = bc1_edit <= cfg.bc1_anchor_hamming
                ? "bc1_offset_window"
                : "bc1_offset_window_and_bc1_anchor_edit";
        } else if (bc1_edit > cfg.bc1_anchor_hamming) {
            probe.failure_class = "bc1_anchor_edit";
        } else {
            probe.failure_class = "strict_anchor_pass";
        }
        return probe;
    }

    const int extra = cfg.bc1_len[col2] - cfg.anchor_len;
    const int bc1_anchor_expected_start = full_start + extra;
    const int min_bc1_anchor_obs_len = std::max(1, cfg.anchor_len - cfg.bc1_anchor_hamming);
    const int max_bc1_anchor_obs_len = cfg.anchor_len + cfg.bc1_anchor_hamming;
    for (int bc1_anchor_obs_len = min_bc1_anchor_obs_len;
         bc1_anchor_obs_len <= max_bc1_anchor_obs_len; ++bc1_anchor_obs_len) {
        const int min_anchor_start = full_start;
        const int max_anchor_start = bc2_boundary - bc1_anchor_obs_len;
        for (int bc1_anchor_start = min_anchor_start; bc1_anchor_start <= max_anchor_start;
             ++bc1_anchor_start) {
            const int bc1_anchor_end = bc1_anchor_start + bc1_anchor_obs_len;
            if (bc1_anchor_start < 0 || bc1_anchor_end > slen) continue;
            const int d1 = bc1_anchor_start - bc1_anchor_expected_start;
            const int d2 = bc2_boundary - bc1_anchor_end;
            const int bc1_edit = edit_distance_bounded(
                seq.data() + bc1_anchor_start, bc1_anchor_obs_len,
                cfg.bc1_anchors[col2].data(),
                static_cast<int>(cfg.bc1_anchors[col2].size()),
                cfg.anchor_len);

            if (!probe.any_offsets_available ||
                anchor_probe_tuple_less(bc1_edit, d1, d2, probe.best_any_bc1_edit,
                                        probe.best_any_d1, probe.best_any_d2)) {
                probe.any_offsets_available = true;
                probe.best_any_bc1_edit = bc1_edit;
                probe.best_any_d1 = d1;
                probe.best_any_d2 = d2;
            }

            if (std::abs(d1) <= cfg.bc1_anchor_offset_window &&
                std::abs(d2) <= cfg.bc2_offset_window) {
                if (!probe.current_offsets_available ||
                    anchor_probe_tuple_less(bc1_edit, d1, d2, probe.best_current_bc1_edit,
                                            probe.best_current_d1, probe.best_current_d2)) {
                    probe.current_offsets_available = true;
                    probe.best_current_bc1_edit = bc1_edit;
                    probe.best_current_d1 = d1;
                    probe.best_current_d2 = d2;
                }
            }
        }
    }

    if (probe.bc2_anchor_edit > cfg.bc2_anchor_hamming) {
        probe.failure_class = "bc2_anchor_edit";
    } else if (!probe.current_offsets_available) {
        if (probe.any_offsets_available && probe.best_any_bc1_edit <= cfg.bc1_anchor_hamming) {
            probe.failure_class = "bc1_offset_window";
        } else if (probe.any_offsets_available) {
            probe.failure_class = "bc1_offset_window_and_bc1_anchor_edit";
        } else {
            probe.failure_class = "bc1_anchor_segment_out_of_bounds";
        }
    } else if (probe.best_current_bc1_edit > cfg.bc1_anchor_hamming) {
        probe.failure_class = "bc1_anchor_edit";
    } else {
        probe.failure_class = "strict_anchor_pass";
    }
    return probe;
}

OracleTargetAudit audit_oracle_target_envelope(const std::string& seq, const Config& cfg,
                                               const OracleTarget& oracle,
                                               int full_start_min, int full_start_max,
                                               int bc2_start_padding,
                                               int edit_max,
                                               int bc2_length_window,
                                               bool require_anchor) {
    OracleTargetAudit best;
    if (!oracle.resolved) return best;
    if (oracle.row2 < 0 || oracle.row2 >= static_cast<int>(cfg.bc2_oligos.size()) ||
        oracle.col2 < 0 || oracle.col2 >= static_cast<int>(cfg.bc1_oligos.size())) {
        return best;
    }
    const int slen = static_cast<int>(seq.size());
    const std::string& bc1_target = cfg.bc1_oligos[oracle.col2];
    const std::string& bc2_target = cfg.bc2_oligos[oracle.row2];

    for (int full_start = full_start_min; full_start <= full_start_max; ++full_start) {
        if (full_start < 0 || full_start >= slen) continue;
        const int min_bc2_start =
            full_start + cfg.anchor_len - edit_max - bc2_start_padding;
        const int max_bc2_start =
            full_start + cfg.anchor_len + 2 + edit_max + bc2_start_padding;
        for (int bc2_start = min_bc2_start; bc2_start <= max_bc2_start; ++bc2_start) {
            const int bc1_obs_len = bc2_start - full_start;
            if (bc1_obs_len <= 0 || full_start + bc1_obs_len > slen) continue;
            const int bc1_edit = edit_distance_bounded(
                seq.data() + full_start, bc1_obs_len,
                bc1_target.data(), static_cast<int>(bc1_target.size()), edit_max);
            if (bc1_edit > edit_max) continue;

            int best_bc2_edit = edit_max + 1;
            int best_bc2_obs_len = -1;
            int best_bc2_align_start = bc2_start;
            int best_bc2_span_delta = 0;
            for (int ld = -bc2_length_window; ld <= bc2_length_window; ++ld) {
                const int bc2_obs_len = static_cast<int>(bc2_target.size()) + ld;
                if (bc2_obs_len <= 0 || bc2_start + bc2_obs_len > slen) continue;
                const int bc2_edit = edit_distance_bounded(
                    seq.data() + bc2_start, bc2_obs_len,
                    bc2_target.data(), static_cast<int>(bc2_target.size()), edit_max);
                const int span_delta = std::abs(ld);
                if (std::make_pair(bc2_edit, span_delta) <
                    std::make_pair(best_bc2_edit, best_bc2_span_delta)) {
                    best_bc2_edit = bc2_edit;
                    best_bc2_obs_len = bc2_obs_len;
                    best_bc2_span_delta = span_delta;
                }
            }
            if (best_bc2_edit > edit_max) continue;

            int anchor_score = -1;
            int offset_sum = -1;
            bool anchor_ok = direct_best_anchor_score_for_pair(
                seq, cfg, full_start, bc2_start, best_bc2_align_start,
                best_bc2_obs_len,
                oracle.row2, oracle.col2,
                anchor_score, offset_sum,
                cfg.direct_best_resolution == DirectBestResolutionMode::OverlapSpanDp,
                cfg.direct_best_resolution == DirectBestResolutionMode::OverlapSpanDp
                    ? cfg.overlap_dp_offset_include_terminal_delta
                    : true);
            AnchorLimitProbe anchor_probe = probe_anchor_limits_for_pair(
                seq, cfg, full_start, bc2_start, best_bc2_align_start,
                best_bc2_obs_len,
                oracle.row2, oracle.col2);
            if (require_anchor && !anchor_ok) continue;

            const int total_edit = bc1_edit + best_bc2_edit;
            const int probe_bc2 = anchor_probe.evaluated ? anchor_probe.bc2_anchor_edit : 9999;
            const int probe_bc1 = anchor_probe.any_offsets_available
                ? anchor_probe.best_any_bc1_edit
                : 9999;
            const int probe_offset_sum = anchor_probe.any_offsets_available
                ? std::abs(anchor_probe.best_any_d1) + std::abs(anchor_probe.best_any_d2)
                : 9999;
            if (!best.possible ||
                std::make_tuple(total_edit, anchor_ok ? anchor_score : 9999,
                                anchor_ok ? offset_sum : 9999, probe_bc2,
                                probe_bc1, probe_offset_sum, full_start, bc2_start) <
                    std::make_tuple(best.min_total_edit,
                                    best.anchor_score >= 0 ? best.anchor_score : 9999,
                                    best.offset_sum >= 0 ? best.offset_sum : 9999,
                                    best.anchor_bc2_edit >= 0 ? best.anchor_bc2_edit : 9999,
                                    best.anchor_any_bc1_edit >= 0 ? best.anchor_any_bc1_edit : 9999,
                                    best.anchor_any_bc1_edit >= 0
                                        ? std::abs(best.anchor_any_d1) + std::abs(best.anchor_any_d2)
                                        : 9999,
                                    best.full_start, best.bc2_start)) {
                best.possible = true;
                best.min_total_edit = total_edit;
                best.bc1_edit = bc1_edit;
                best.bc2_edit = best_bc2_edit;
                best.anchor_score = anchor_ok ? anchor_score : -1;
                best.offset_sum = anchor_ok ? offset_sum : -1;
                best.full_start = full_start;
                best.bc2_start = bc2_start;
                best.bc1_obs_len = bc1_obs_len;
                best.bc2_obs_len = best_bc2_obs_len;
                best.observed_non_acgt =
                    has_non_acgt(seq, full_start, bc1_obs_len) ||
                    has_non_acgt(seq, bc2_start, best_bc2_obs_len);
                best.anchor_failure_class = anchor_probe.failure_class;
                best.anchor_bc2_edit = anchor_probe.bc2_anchor_edit;
                best.anchor_current_bc1_edit = anchor_probe.best_current_bc1_edit;
                best.anchor_current_d1 = anchor_probe.best_current_d1;
                best.anchor_current_d2 = anchor_probe.best_current_d2;
                best.anchor_any_bc1_edit = anchor_probe.best_any_bc1_edit;
                best.anchor_any_d1 = anchor_probe.best_any_d1;
                best.anchor_any_d2 = anchor_probe.best_any_d2;
            }
        }
    }
    return best;
}

Decoded decode_record(const std::string& seq, Config& cfg,
                      std::vector<PackedHit>& lk1,
                      std::vector<PackedHit>& lk2,
                      LookupScratch& lookup_scratch,
                      std::vector<PackedHit>& cache_hits,
                      std::vector<PackedHit>& bc1_full_hits,
                      std::vector<PackedHit>& bc2_full_hits,
                      std::vector<FullCandidate>& bc1_valid,
                      std::vector<FullCandidate>& bc2_valid,
                      ObservedDecodeCacheLocals& observed_locals,
                      bool use_observed_cache,
                      DecodeStats* stats,
                      const OracleTarget* oracle = nullptr) {
    const int anchor_len = cfg.anchor_len;
    const int slen = static_cast<int>(seq.size());
    const bool need_affine =
        cfg.score_mode == CandidateScoreMode::AffineGap || cfg.emit_decode_annotations;
    Decoded out;

    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        if (full_start < 0) continue;
        for (int extra = 1; extra <= 2; ++extra) {
            int bc1_len = anchor_len + extra;
            int bc1_anchor_expected_start = full_start + extra;
            for (int d1 = -cfg.bc1_anchor_offset_window; d1 <= cfg.bc1_anchor_offset_window; ++d1) {
                int bc1_anchor_start = bc1_anchor_expected_start + d1;
                int bc1_anchor_end = bc1_anchor_start + anchor_len;
                if (bc1_anchor_start < full_start || slen < bc1_anchor_end) continue;
                if (!is_acgt(seq, bc1_anchor_start, anchor_len)) continue;
                if (use_observed_cache) {
                    cfg.bc1_anchor_observed.lookup(seq.data() + bc1_anchor_start, anchor_len,
                                                   observed_locals.bc1_anchor,
                                                   observed_locals.scratch, lk1);
                } else {
                    cfg.bc1_lookup.lookup(seq.data() + bc1_anchor_start, anchor_len,
                                          lk1, lookup_scratch, cache_hits, stats);
                }
                if (lk1.empty()) continue;
                int bc2_expected_start = bc1_anchor_end;
                for (int d2 = -cfg.bc2_offset_window; d2 <= cfg.bc2_offset_window; ++d2) {
                    int bc2_anchor_start = bc2_expected_start + d2;
                    int bc2_anchor_end = bc2_anchor_start + cfg.bc2_seg_len;
                    if (bc2_anchor_start < 0 || slen < bc2_anchor_end) continue;
                    int bc1_obs_len = bc2_anchor_start - full_start;
                    if (bc1_obs_len <= 0 || full_start + bc1_obs_len > slen) continue;
                    bool bc1_cache_found = false;
                    if (use_observed_cache) {
                        cfg.bc1_full_observed.lookup(seq.data() + full_start, bc1_obs_len,
                                                     observed_locals.bc1_full,
                                                     observed_locals.scratch, bc1_full_hits);
                    } else if (bc1_obs_len == bc1_len) {
                        if (stats) ++stats->full_cache_lookup_calls;
                        bc1_cache_found = cfg.bc1_full_lookup.lookup(seq.data() + full_start,
                                                                      bc1_obs_len, bc1_full_hits);
                        if (stats && bc1_cache_found) ++stats->full_cache_keys_found;
                    }

                    bc1_valid.clear();
                    for (auto& c1 : lk1) {
                        int bc1_idx = c1.idx, bc1_anchor_h = c1.distance;
                        if (cfg.bc1_len[bc1_idx] != bc1_len) continue;  // variable-length filter
                        int bc1_full_h = use_observed_cache
                            ? hit_distance_for_idx(bc1_full_hits, bc1_idx, cfg.bc1_full_hamming + 1)
                            : cached_edit_for_target_from_hits(
                                  bc1_cache_found, bc1_full_hits, seq, full_start, bc1_obs_len,
                                  cfg.bc1_oligos, bc1_idx, cfg.bc1_full_hamming, stats);
                        if (bc1_full_h > cfg.bc1_full_hamming) continue;
                        AlignmentMetrics bc1_affine;
                        if (need_affine) {
                            bc1_affine = affine_alignment_metrics(
                                seq.data() + full_start, bc1_obs_len,
                                cfg.bc1_oligos[bc1_idx].data(),
                                static_cast<int>(cfg.bc1_oligos[bc1_idx].size()),
                                cfg.affine);
                        }
                        bc1_valid.push_back({static_cast<uint16_t>(bc1_idx),
                                             static_cast<uint8_t>(bc1_anchor_h),
                                             static_cast<uint8_t>(bc1_full_h),
                                             static_cast<uint8_t>(bc1_obs_len),
                                             bc1_affine});
                    }
                    if (bc1_valid.empty()) continue;

                    if (!is_acgt(seq, bc2_anchor_start, cfg.bc2_seg_len)) continue;
                    if (use_observed_cache) {
                        cfg.bc2_anchor_observed.lookup(seq.data() + bc2_anchor_start,
                                                       cfg.bc2_seg_len,
                                                       observed_locals.bc2_anchor,
                                                       observed_locals.scratch, lk2);
                    } else {
                        cfg.bc2_lookup.lookup(seq.data() + bc2_anchor_start, cfg.bc2_seg_len,
                                              lk2, lookup_scratch, cache_hits, stats);
                    }
                    if (lk2.empty()) continue;

                    bc2_valid.clear();
                    for (auto& c2 : lk2) {
                        int bc2_idx = c2.idx;
                        int bc2_full_h = cfg.bc2_full_hamming + 1;
                        int bc2_obs_len = 0;
                        AlignmentMetrics bc2_affine;
                        if (need_affine) {
                            ObservedAlignmentChoice choice = best_observed_alignment_choice(
                                seq, bc2_anchor_start, cfg.bc2_oligos[bc2_idx],
                                cfg.bc2_offset_window, cfg.bc2_full_hamming, cfg.affine,
                                cfg.score_mode == CandidateScoreMode::AffineGap,
                                need_affine);
                            bc2_full_h = choice.edit_distance;
                            bc2_obs_len = choice.obs_len;
                            bc2_affine = choice.affine;
                        } else {
                            bc2_full_h = use_observed_cache
                                ? best_observed_edit_distance_observed(
                                      seq, bc2_anchor_start, bc2_idx, cfg.bc2_oligos,
                                      cfg.bc2_full_observed, observed_locals.bc2_full,
                                      observed_locals.scratch, cfg.bc2_offset_window,
                                      cfg.bc2_full_hamming, bc2_full_hits)
                                : best_observed_edit_distance(
                                      seq, bc2_anchor_start, bc2_idx, cfg.bc2_oligos,
                                      cfg.bc2_full_lookup, cfg.bc2_offset_window,
                                      cfg.bc2_full_hamming, bc2_full_hits, stats);
                            bc2_obs_len = static_cast<int>(cfg.bc2_oligos[bc2_idx].size());
                        }
                        if (bc2_full_h <= cfg.bc2_full_hamming) {
                            bc2_valid.push_back({static_cast<uint16_t>(bc2_idx), c2.distance,
                                                 static_cast<uint8_t>(bc2_full_h),
                                                 static_cast<uint8_t>(bc2_obs_len),
                                                 bc2_affine});
                        }
                    }
                    if (bc2_valid.empty()) continue;

                    for (auto& c1 : bc1_valid) {
                        int bc1_idx = c1.idx;
                        int bc1_anchor_h = c1.anchor_distance;
                        int bc1_full_h = c1.full_distance;
                        for (auto& c2 : bc2_valid) {
                            int bc2_idx = c2.idx, bc2_anchor_h = c2.anchor_distance;
                            int row2 = bc2_idx, col2 = bc1_idx;
                            if (row2 < 0 || row2 >= cfg.grid_rows || col2 < 0 || col2 >= cfg.grid_cols) continue;
                            int bc2_full_h = c2.full_distance;
                            int anchor_h = bc1_anchor_h + bc2_anchor_h;
                            int offset_sum = std::abs(d1) + std::abs(d2);
                            consider_candidate(out, cfg, bc1_full_h, bc2_full_h,
                                               c1.affine, c2.affine, anchor_h,
                                               offset_sum, full_start, row2, col2,
                                               c1.obs_len, c2.obs_len, oracle);
                        }
                    }
                }
            }
        }
    }
    finalize_decoded_assignment(out);
    return out;
}

HitSpan lookup_full_prebuilt_or_dp(const std::string& seq, int start, int length,
                                   const std::vector<std::string>& targets,
                                   const PrebuiltObservedAnswerCache& cache,
                                   const DeletionCandidateIndex& deletion_index,
                                   ObservedLookupScratch& scratch,
                                   DpOnlyMemo& dp_memo,
                                   std::vector<PackedHit>& dp_hits,
                                   int maxk,
                                   bool fallback_on_cache_miss,
                                   bool* cache_miss) {
    if (start < 0 || length <= 0 || start + length > static_cast<int>(seq.size())) return {};
    HitSpan hits = cache.lookup_span(seq.data() + start, length);
    if (hits.found) return hits;
    if (fallback_on_cache_miss) {
        if (cache_miss) *cache_miss = true;
        return {};
    }
    lookup_dp_only(seq.data() + start, length, targets, &deletion_index,
                   maxk, scratch, dp_memo, dp_hits);
    return {dp_hits.data(), static_cast<uint16_t>(dp_hits.size()), true};
}

void append_direct_candidates(HitSpan hits,
                              const std::vector<int>& target_lengths,
                              int expected_target_length,
                              int start,
                              int boundary,
                              int max_distance,
                              std::vector<DirectHalfCandidate>& out) {
    if (!hits.found || hits.count == 0) return;
    for (uint16_t i = 0; i < hits.count; ++i) {
        const PackedHit& hit = hits.data[i];
        if (hit.distance > max_distance) continue;
        if (expected_target_length > 0 && target_lengths[hit.idx] != expected_target_length) continue;
        out.push_back({hit.idx, hit.distance, static_cast<uint8_t>(start),
                       static_cast<uint8_t>(start), static_cast<uint8_t>(boundary), 0});
    }
}

Decoded decode_record_direct_prebuilt(const std::string& seq, Config& cfg,
                                      ObservedDecodeCacheLocals& observed_locals,
                                      bool fallback_on_cache_miss,
                                      bool* cache_miss) {
    const int slen = static_cast<int>(seq.size());
    bool have_best = false;
    long best_key[2] = {0, 0};
    int best_row = -1, best_col = -1;
    bool best_unique = true;

    auto consider = [&](int full_h, int full_start, int row2, int col2) {
        long key[2] = {full_h, full_start};
        if (!have_best) {
            have_best = true;
            std::memcpy(best_key, key, sizeof(best_key));
            best_row = row2;
            best_col = col2;
            best_unique = true;
            return;
        }
        int cmp = 0;
        for (int i = 0; i < 2; ++i) {
            if (key[i] != best_key[i]) {
                cmp = (key[i] < best_key[i]) ? -1 : 1;
                break;
            }
        }
        if (cmp < 0) {
            std::memcpy(best_key, key, sizeof(best_key));
            best_row = row2;
            best_col = col2;
            best_unique = true;
        } else if (cmp == 0) {
            if (row2 != best_row || col2 != best_col) best_unique = false;
        }
    };

    const int bc1_obs_len_min = cfg.direct_full_exact_lengths
        ? cfg.anchor_len + 1
        : std::max(1, cfg.anchor_len + 1 - cfg.bc1_full_hamming);
    const int bc1_obs_len_max = cfg.direct_full_exact_lengths
        ? cfg.anchor_len + 2
        : cfg.anchor_len + 2 + cfg.bc1_full_hamming;
    const int bc2_obs_len_min = cfg.direct_full_exact_lengths
        ? cfg.bc2_seg_len
        : std::max(1, cfg.bc2_seg_len - cfg.bc2_full_hamming);
    const int bc2_obs_len_max = cfg.direct_full_exact_lengths
        ? cfg.bc2_seg_len + 1
        : cfg.bc2_seg_len + 1 + cfg.bc2_full_hamming;
    observed_locals.direct_bc1.clear();
    observed_locals.direct_bc2.clear();

    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        if (full_start < 0) continue;
        for (int bc1_obs_len = bc1_obs_len_min; bc1_obs_len <= bc1_obs_len_max; ++bc1_obs_len) {
            const int boundary = full_start + bc1_obs_len;
            if (boundary > slen) continue;
            HitSpan hits = lookup_full_prebuilt_or_dp(
                seq, full_start, bc1_obs_len, cfg.bc1_oligos, cfg.bc1_full_prebuilt,
                cfg.bc1_full_deletion2, observed_locals.scratch,
                observed_locals.bc1_full_dp_miss, observed_locals.bc1_full_dp_hits,
                cfg.bc1_full_hamming, fallback_on_cache_miss, cache_miss);
            if (cache_miss && *cache_miss) return {};
            const int expected_bc1_len = cfg.direct_full_exact_lengths ? bc1_obs_len : 0;
            append_direct_candidates(hits, cfg.bc1_len, expected_bc1_len, full_start, boundary,
                                     cfg.bc1_full_hamming, observed_locals.direct_bc1);
        }
    }
    if (observed_locals.direct_bc1.empty()) return {};

    bool candidate_bc2_starts[kMaxLen] = {false};
    int min_bc2_start = slen;
    int max_bc2_start = -1;
    for (const DirectHalfCandidate& c1 : observed_locals.direct_bc1) {
        if (c1.boundary >= kMaxLen) continue;
        candidate_bc2_starts[c1.boundary] = true;
        min_bc2_start = std::min<int>(min_bc2_start, c1.boundary);
        max_bc2_start = std::max<int>(max_bc2_start, c1.boundary);
    }
    if (max_bc2_start < min_bc2_start) return {};

    for (int bc2_start = min_bc2_start; bc2_start <= max_bc2_start; ++bc2_start) {
        if (bc2_start < 0 || bc2_start >= kMaxLen || !candidate_bc2_starts[bc2_start]) continue;
        for (int bc2_obs_len = bc2_obs_len_min; bc2_obs_len <= bc2_obs_len_max; ++bc2_obs_len) {
            if (bc2_start + bc2_obs_len > slen) continue;
            HitSpan hits = lookup_full_prebuilt_or_dp(
                seq, bc2_start, bc2_obs_len, cfg.bc2_oligos, cfg.bc2_full_prebuilt,
                cfg.bc2_full_deletion2, observed_locals.scratch,
                observed_locals.bc2_full_dp_miss, observed_locals.bc2_full_dp_hits,
                cfg.bc2_full_hamming, fallback_on_cache_miss, cache_miss);
            if (cache_miss && *cache_miss) return {};
            const int expected_bc2_len = cfg.direct_full_exact_lengths ? bc2_obs_len : 0;
            append_direct_candidates(hits, cfg.bc2_len, expected_bc2_len, bc2_start, bc2_start,
                                     cfg.bc2_full_hamming, observed_locals.direct_bc2);
        }
    }
    if (observed_locals.direct_bc2.empty()) return {};

    for (const DirectHalfCandidate& c1 : observed_locals.direct_bc1) {
        for (const DirectHalfCandidate& c2 : observed_locals.direct_bc2) {
            if (c1.boundary != c2.start) continue;
            int row2 = c2.idx;
            int col2 = c1.idx;
            if (row2 < 0 || row2 >= cfg.grid_rows || col2 < 0 || col2 >= cfg.grid_cols) continue;
            consider(static_cast<int>(c1.distance) + static_cast<int>(c2.distance),
                     c1.start, row2, col2);
        }
    }

    Decoded out;
    if (have_best && best_unique) {
        out.assigned = true;
        out.row2 = best_row;
        out.col2 = best_col;
        out.full_start = static_cast<int>(best_key[1]);
    }
    return out;
}

HitSpan lookup_direct_best_span_or_non_acgt_null(const PackedAnswerCache& cache,
                                                 const char* query,
                                                 int len,
                                                 bool non_acgt_null_enabled,
                                                 uint32_t& non_acgt_null_queries) {
    uint64_t packed = 0;
    if (pack_sequence(query, len, packed)) {
        return cache.lookup_span(query, len);
    }
    if (!has_non_acgt(query, len)) return {};
    if (non_acgt_null_enabled) {
        ++non_acgt_null_queries;
        return {nullptr, 0, true};
    }
    return {};
}

void dedupe_overlap_dp_half_candidates(std::vector<DirectDpHalfCandidate>& candidates,
                                       const std::vector<int>& target_lengths,
                                       int expected_start) {
    std::sort(candidates.begin(), candidates.end(),
              [](const DirectDpHalfCandidate& a, const DirectDpHalfCandidate& b) {
                  if (a.idx != b.idx) return a.idx < b.idx;
                  return std::make_tuple(a.choice.edit_distance,
                                         a.choice.terminal_delta,
                                         a.choice.start,
                                         a.choice.obs_len) <
                         std::make_tuple(b.choice.edit_distance,
                                         b.choice.terminal_delta,
                                         b.choice.start,
                                         b.choice.obs_len);
              });
    size_t out = 0;
    for (size_t i = 0; i < candidates.size();) {
        DirectDpHalfCandidate best = candidates[i];
        ++i;
        while (i < candidates.size() && candidates[i].idx == best.idx) {
            const int idx = candidates[i].idx;
            const int expected_len =
                (idx >= 0 && idx < static_cast<int>(target_lengths.size()))
                    ? target_lengths[idx]
                    : 0;
            if (half_span_choice_less(candidates[i].choice, best.choice,
                                      expected_start, expected_len)) {
                best = candidates[i];
            }
            ++i;
        }
        candidates[out++] = best;
    }
    candidates.resize(out);
}

HalfSpanChoice overlap_choice_from_hash_hit(const std::string& seq,
                                            int start,
                                            int obs_len,
                                            const std::string& target,
                                            uint8_t distance,
                                            const Config& cfg,
                                            bool need_affine) {
    HalfSpanChoice choice;
    choice.found = true;
    choice.edit_distance = distance;
    choice.raw_edit_distance = distance;
    choice.terminal_delta = std::abs(obs_len - static_cast<int>(target.size()));
    choice.start = start;
    choice.obs_len = obs_len;
    if (need_affine) {
        choice.affine = affine_alignment_metrics(
            seq.data() + start, obs_len, target.data(),
            static_cast<int>(target.size()), cfg.affine);
    }
    return choice;
}

void append_overlap_hash_candidates(HitSpan hits,
                                    const std::string& seq,
                                    int start,
                                    int obs_len,
                                    const std::vector<std::string>& targets,
                                    const std::vector<int>& target_lengths,
                                    int min_obs_len,
                                    int max_len_delta,
                                    int max_distance,
                                    const Config& cfg,
                                    bool need_affine,
                                    const OracleTarget* oracle,
                                    bool oracle_is_bc1,
                                    Decoded& out,
                                    std::vector<DirectDpHalfCandidate>& candidates) {
    if (!hits.found || hits.count == 0) return;
    if (min_obs_len >= 0 && obs_len < min_obs_len) return;
    for (uint16_t i = 0; i < hits.count; ++i) {
        const PackedHit& hit = hits.data[i];
        if (hit.distance > max_distance) continue;
        if (hit.idx >= target_lengths.size() || hit.idx >= targets.size()) continue;
        if (max_len_delta >= 0 &&
            std::abs(obs_len - target_lengths[hit.idx]) > max_len_delta) {
            continue;
        }
        if (oracle && oracle->resolved) {
            if (oracle_is_bc1 && hit.idx == oracle->col2) {
                out.oracle_bc1_seed_seen = true;
            } else if (!oracle_is_bc1 && hit.idx == oracle->row2) {
                out.oracle_bc2_seed_seen = true;
            }
        }
        HalfSpanChoice choice = overlap_choice_from_hash_hit(
            seq, start, obs_len, targets[hit.idx], hit.distance, cfg, need_affine);
        candidates.push_back({hit.idx, choice});
    }
}

Decoded decode_record_overlap_span_dp_target_scan(const std::string& seq, Config& cfg,
                                                  ObservedDecodeCacheLocals& observed_locals,
                                                  const OracleTarget* oracle = nullptr) {
    const int slen = static_cast<int>(seq.size());
    Decoded out;
    const bool need_affine =
        cfg.score_mode == CandidateScoreMode::AffineGap || cfg.emit_decode_annotations;

    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        if (full_start < 0 || full_start >= slen) continue;
        observed_locals.overlap_bc1.clear();
        for (int col2 = 0; col2 < static_cast<int>(cfg.bc1_oligos.size()); ++col2) {
            const int bc1_max_len_delta = cfg.overlap_dp_length_guard
                ? cfg.overlap_dp_bc1_max_len_delta
                : -1;
            const int bc1_min_obs_len = cfg.overlap_dp_length_guard
                ? cfg.overlap_dp_bc1_min_obs_len
                : -1;
            HalfSpanChoice choice = best_half_window_dp_choice(
                seq, full_start, full_start + kDirectBestBc1MaxQueryLength,
                cfg.bc1_oligos[col2], cfg.bc1_full_hamming,
                full_start, cfg.bc1_len[col2], bc1_min_obs_len,
                bc1_max_len_delta, false, cfg.affine, need_affine);
            if (!choice.found) continue;
            if (oracle && oracle->resolved && col2 == oracle->col2) {
                out.oracle_bc1_seed_seen = true;
            }
            observed_locals.overlap_bc1.push_back(
                {static_cast<uint16_t>(col2), choice});
        }
        if (observed_locals.overlap_bc1.empty()) continue;

        bool boundary_needed[kMaxLen] = {false};
        std::vector<int> boundaries;
        for (const DirectDpHalfCandidate& c1 : observed_locals.overlap_bc1) {
            const int boundary = full_start + c1.choice.obs_len;
            if (boundary < 0 || boundary >= kMaxLen || boundary > slen) continue;
            if (!boundary_needed[boundary]) {
                boundary_needed[boundary] = true;
                boundaries.push_back(boundary);
            }
        }
        for (int boundary : boundaries) {
            observed_locals.overlap_bc2_by_boundary[boundary].clear();
            const int bc2_window_start =
                full_start + cfg.anchor_len - cfg.bc1_full_hamming;
            const int bc2_window_end =
                full_start + cfg.anchor_len + 2 + cfg.bc1_full_hamming +
                kDirectBestBc2MaxQueryLength;
            for (int row2 = 0; row2 < static_cast<int>(cfg.bc2_oligos.size()); ++row2) {
                const int bc2_max_len_delta = cfg.overlap_dp_length_guard
                    ? cfg.overlap_dp_bc2_max_len_delta
                    : -1;
                const int bc2_min_obs_len = cfg.overlap_dp_length_guard
                    ? cfg.overlap_dp_bc2_min_obs_len
                    : -1;
                HalfSpanChoice choice = best_half_window_dp_choice(
                    seq, bc2_window_start, bc2_window_end, cfg.bc2_oligos[row2],
                    cfg.bc2_full_hamming, boundary, cfg.bc2_len[row2],
                    bc2_min_obs_len, bc2_max_len_delta, true, cfg.affine, need_affine);
                if (!choice.found) continue;
                if (oracle && oracle->resolved && row2 == oracle->row2) {
                    out.oracle_bc2_seed_seen = true;
                }
                observed_locals.overlap_bc2_by_boundary[boundary].push_back(
                    {static_cast<uint16_t>(row2), choice});
            }
        }

        for (const DirectDpHalfCandidate& c1 : observed_locals.overlap_bc1) {
            const int col2 = c1.idx;
            const int bc2_boundary = full_start + c1.choice.obs_len;
            if (bc2_boundary < 0 || bc2_boundary >= kMaxLen) continue;
            const auto& bc2_candidates =
                observed_locals.overlap_bc2_by_boundary[bc2_boundary];
            if (bc2_candidates.empty()) continue;
            for (const DirectDpHalfCandidate& c2 : bc2_candidates) {
                const int row2 = c2.idx;
                if (row2 < 0 || row2 >= cfg.grid_rows ||
                    col2 < 0 || col2 >= cfg.grid_cols) {
                    continue;
                }
                const bool is_oracle_pair =
                    oracle && oracle->resolved &&
                    row2 == oracle->row2 && col2 == oracle->col2;
                if (is_oracle_pair) {
                    out.oracle_pair_seen = true;
                    out.oracle_pair_bc1_pass = true;
                    out.oracle_pair_bc2_pass = true;
                }

                int anchor_score = 0;
                int offset_sum = 0;
                if (!direct_best_anchor_score_for_pair(
                        seq, cfg, full_start, bc2_boundary, c2.choice.start,
                        c2.choice.obs_len, row2, col2, anchor_score, offset_sum,
                        true, cfg.overlap_dp_offset_include_terminal_delta)) {
                    continue;
                }
                if (is_oracle_pair) out.oracle_pair_anchor_pass = true;
                consider_candidate(out, cfg,
                                   c1.choice.edit_distance,
                                   c2.choice.edit_distance,
                                   c1.choice.affine, c2.choice.affine,
                                   anchor_score, offset_sum, full_start,
                                   row2, col2, c1.choice.obs_len,
                                   c2.choice.obs_len, oracle, true);
            }
        }
    }

    finalize_decoded_assignment(out);
    return out;
}

Decoded decode_record_overlap_span_dp_hash(const std::string& seq, Config& cfg,
                                           ObservedDecodeCacheLocals& observed_locals,
                                           const OracleTarget* oracle = nullptr) {
    if (direct_best_decode_window_has_non_acgt(seq, cfg)) {
        Decoded out =
            decode_record_overlap_span_dp_target_scan(seq, cfg, observed_locals, oracle);
        out.non_acgt_dp_checked = true;
        return out;
    }

    const int slen = static_cast<int>(seq.size());
    Decoded out;
    const bool need_affine =
        cfg.score_mode == CandidateScoreMode::AffineGap || cfg.emit_decode_annotations;
    const int bc1_max_len_delta = cfg.overlap_dp_length_guard
        ? cfg.overlap_dp_bc1_max_len_delta
        : -1;
    const int bc1_min_obs_len = cfg.overlap_dp_length_guard
        ? cfg.overlap_dp_bc1_min_obs_len
        : -1;
    const int bc2_max_len_delta = cfg.overlap_dp_length_guard
        ? cfg.overlap_dp_bc2_max_len_delta
        : -1;
    const int bc2_min_obs_len = cfg.overlap_dp_length_guard
        ? cfg.overlap_dp_bc2_min_obs_len
        : -1;

    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        if (full_start < 0 || full_start >= slen) continue;
        observed_locals.overlap_bc1.clear();
        for (int bc1_query_len : kDirectBestBc1QueryLengths) {
            if (full_start + bc1_query_len > slen) continue;
            HitSpan bc1_hits =
                cfg.bc1_direct_best_edit_lookup.lookup_span(seq.data() + full_start,
                                                            bc1_query_len);
            append_overlap_hash_candidates(
                bc1_hits, seq, full_start, bc1_query_len, cfg.bc1_oligos,
                cfg.bc1_len, bc1_min_obs_len, bc1_max_len_delta,
                cfg.bc1_full_hamming, cfg, need_affine, oracle, true, out,
                observed_locals.overlap_bc1);
        }
        dedupe_overlap_dp_half_candidates(observed_locals.overlap_bc1,
                                          cfg.bc1_len, full_start);
        if (observed_locals.overlap_bc1.empty()) continue;

        bool boundary_needed[kMaxLen] = {false};
        std::vector<int> boundaries;
        for (const DirectDpHalfCandidate& c1 : observed_locals.overlap_bc1) {
            const int boundary = full_start + c1.choice.obs_len;
            if (boundary < 0 || boundary >= kMaxLen || boundary > slen) continue;
            if (!boundary_needed[boundary]) {
                boundary_needed[boundary] = true;
                boundaries.push_back(boundary);
            }
        }

        for (int boundary : boundaries) {
            observed_locals.overlap_bc2_by_boundary[boundary].clear();
            const int bc2_window_start =
                std::max(0, full_start + cfg.anchor_len - cfg.bc1_full_hamming);
            const int bc2_window_end = std::min(
                slen,
                full_start + cfg.anchor_len + 2 + cfg.bc1_full_hamming +
                    kDirectBestBc2MaxQueryLength);
            for (int bc2_start = bc2_window_start; bc2_start < bc2_window_end; ++bc2_start) {
                for (int bc2_query_len : kDirectBestBc2QueryLengths) {
                    if (bc2_start + bc2_query_len > bc2_window_end) continue;
                    HitSpan bc2_hits =
                        cfg.bc2_direct_best_edit_lookup.lookup_span(seq.data() + bc2_start,
                                                                    bc2_query_len);
                    append_overlap_hash_candidates(
                        bc2_hits, seq, bc2_start, bc2_query_len, cfg.bc2_oligos,
                        cfg.bc2_len, bc2_min_obs_len, bc2_max_len_delta,
                        cfg.bc2_full_hamming, cfg, need_affine, oracle, false,
                        out, observed_locals.overlap_bc2_by_boundary[boundary]);
                }
            }
            dedupe_overlap_dp_half_candidates(
                observed_locals.overlap_bc2_by_boundary[boundary],
                cfg.bc2_len, boundary);
        }

        for (const DirectDpHalfCandidate& c1 : observed_locals.overlap_bc1) {
            const int col2 = c1.idx;
            const int bc2_boundary = full_start + c1.choice.obs_len;
            if (bc2_boundary < 0 || bc2_boundary >= kMaxLen) continue;
            const auto& bc2_candidates =
                observed_locals.overlap_bc2_by_boundary[bc2_boundary];
            if (bc2_candidates.empty()) continue;
            for (const DirectDpHalfCandidate& c2 : bc2_candidates) {
                const int row2 = c2.idx;
                if (row2 < 0 || row2 >= cfg.grid_rows ||
                    col2 < 0 || col2 >= cfg.grid_cols) {
                    continue;
                }
                const bool is_oracle_pair =
                    oracle && oracle->resolved &&
                    row2 == oracle->row2 && col2 == oracle->col2;
                if (is_oracle_pair) {
                    out.oracle_pair_seen = true;
                    out.oracle_pair_bc1_pass = true;
                    out.oracle_pair_bc2_pass = true;
                }

                int anchor_score = 0;
                int offset_sum = 0;
                if (!direct_best_anchor_score_for_pair(
                        seq, cfg, full_start, bc2_boundary, c2.choice.start,
                        c2.choice.obs_len, row2, col2, anchor_score, offset_sum,
                        true, cfg.overlap_dp_offset_include_terminal_delta)) {
                    continue;
                }
                if (is_oracle_pair) out.oracle_pair_anchor_pass = true;
                consider_candidate(out, cfg,
                                   c1.choice.edit_distance,
                                   c2.choice.edit_distance,
                                   c1.choice.affine, c2.choice.affine,
                                   anchor_score, offset_sum, full_start,
                                   row2, col2, c1.choice.obs_len,
                                   c2.choice.obs_len, oracle, true);
            }
        }
    }

    finalize_decoded_assignment(out);
    return out;
}

Decoded decode_record_direct_best_edit(const std::string& seq, Config& cfg,
                                       ObservedDecodeCacheLocals& observed_locals,
                                       const OracleTarget* oracle = nullptr) {
    const int slen = static_cast<int>(seq.size());
    Decoded out;
    if (cfg.direct_best_resolution == DirectBestResolutionMode::OverlapSpanDp) {
        if (cfg.overlap_dp_use_edit_hash) {
            return decode_record_overlap_span_dp_hash(seq, cfg, observed_locals, oracle);
        }
        return decode_record_overlap_span_dp_target_scan(seq, cfg, observed_locals, oracle);
    }
    const bool non_acgt_window = direct_best_decode_window_has_non_acgt(seq, cfg);
    if (non_acgt_window && !cfg.direct_best_non_acgt_dp_fallback) {
        out.sequestered_non_acgt = true;
        return out;
    }

    uint32_t non_acgt_null_queries = 0;

    auto anchor_score_for_pair = [&](int full_start, int bc2_boundary,
                                     int bc2_align_start, int bc2_obs_len,
                                     int row2, int col2, int& anchor_score,
                                     int& offset_sum) {
        return direct_best_anchor_score_for_pair(
            seq, cfg, full_start, bc2_boundary, bc2_align_start, bc2_obs_len,
            row2, col2, anchor_score, offset_sum);
    };

    observed_locals.direct_bc1.clear();
    observed_locals.direct_bc2.clear();

    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        if (full_start < 0) continue;
        for (int bc1_query_len : kDirectBestBc1QueryLengths) {
            if (full_start + bc1_query_len > slen) continue;
            HitSpan bc1_hits = lookup_direct_best_span_or_non_acgt_null(
                cfg.bc1_direct_best_edit_lookup, seq.data() + full_start, bc1_query_len,
                cfg.direct_best_non_acgt_dp_fallback, non_acgt_null_queries);
            if (!bc1_hits.found || bc1_hits.count == 0) continue;
            for (uint16_t i = 0; i < bc1_hits.count; ++i) {
                const PackedHit& hit = bc1_hits.data[i];
                if (oracle && oracle->resolved && hit.idx == oracle->col2) {
                    out.oracle_bc1_seed_seen = true;
                }
                observed_locals.direct_bc1.push_back(
                    {hit.idx, hit.distance, static_cast<uint8_t>(full_start),
                     static_cast<uint8_t>(full_start), 0, 0});
            }
        }
    }
    dedupe_direct_half_candidates(observed_locals.direct_bc1);
    if (observed_locals.direct_bc1.empty()) {
        out.non_acgt_null_queries = non_acgt_null_queries;
        return out;
    }

    bool full_start_has_bc1[kMaxLen] = {false};
    for (const DirectHalfCandidate& c1 : observed_locals.direct_bc1) {
        const int full_start = c1.full_start;
        if (full_start < 0 || full_start >= kMaxLen) continue;
        full_start_has_bc1[full_start] = true;
    }

    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        if (full_start < 0 || full_start >= kMaxLen || !full_start_has_bc1[full_start]) continue;
        const int min_bc2_start = full_start + cfg.anchor_len - cfg.bc1_full_hamming;
        const int max_bc2_start = full_start + cfg.anchor_len + 2 + cfg.bc1_full_hamming;
        for (int bc2_boundary = min_bc2_start; bc2_boundary <= max_bc2_start; ++bc2_boundary) {
            if (bc2_boundary < 0 || bc2_boundary > slen) continue;
            for (int bc2_query_len : kDirectBestBc2QueryLengths) {
                if (bc2_boundary + bc2_query_len > slen) continue;
                HitSpan bc2_hits = lookup_direct_best_span_or_non_acgt_null(
                    cfg.bc2_direct_best_edit_lookup, seq.data() + bc2_boundary,
                    bc2_query_len, cfg.direct_best_non_acgt_dp_fallback,
                    non_acgt_null_queries);
                if (!bc2_hits.found || bc2_hits.count == 0) continue;
                for (uint16_t i = 0; i < bc2_hits.count; ++i) {
                    const PackedHit& hit = bc2_hits.data[i];
                    if (oracle && oracle->resolved && hit.idx == oracle->row2) {
                        out.oracle_bc2_seed_seen = true;
                    }
                    observed_locals.direct_bc2.push_back(
                        {hit.idx, hit.distance, static_cast<uint8_t>(bc2_boundary),
                         static_cast<uint8_t>(full_start),
                         static_cast<uint8_t>(bc2_boundary), 0});
                }
            }
        }
    }
    dedupe_direct_half_candidates(observed_locals.direct_bc2);
    if (observed_locals.direct_bc2.empty()) {
        out.non_acgt_null_queries = non_acgt_null_queries;
        return out;
    }

    for (const DirectHalfCandidate& c1 : observed_locals.direct_bc1) {
        const int full_start = c1.full_start;
        for (const DirectHalfCandidate& c2 : observed_locals.direct_bc2) {
            if (c2.full_start != c1.full_start) continue;
            const int row2 = c2.idx;
            const int col2 = c1.idx;
            if (row2 < 0 || row2 >= cfg.grid_rows ||
                col2 < 0 || col2 >= cfg.grid_cols) {
                continue;
            }
            const bool is_oracle_pair =
                oracle && oracle->resolved && row2 == oracle->row2 && col2 == oracle->col2;
            if (is_oracle_pair) out.oracle_pair_seen = true;

            const std::string& bc1_target = cfg.bc1_oligos[col2];
            const std::string& bc2_target = cfg.bc2_oligos[row2];
            int bc2_boundary = static_cast<int>(c2.boundary);
            int bc2_align_start = static_cast<int>(c2.start);
            int bc1_obs_len = bc2_boundary - full_start;
            int bc2_obs_len = -1;
            int bc1_full_d = cfg.bc1_full_hamming + 1;
            int bc2_full_d = cfg.bc2_full_hamming + 1;
            AlignmentMetrics bc1_affine;
            AlignmentMetrics bc2_affine;

            if (bc1_obs_len <= 0 || full_start + bc1_obs_len > slen) continue;
            bc1_full_d = edit_distance_bounded(
                seq.data() + full_start, bc1_obs_len,
                bc1_target.data(), static_cast<int>(bc1_target.size()),
                cfg.bc1_full_hamming);
            if (bc1_full_d > cfg.bc1_full_hamming) continue;
            if (cfg.score_mode == CandidateScoreMode::AffineGap ||
                cfg.emit_decode_annotations) {
                bc1_affine = affine_alignment_metrics(
                    seq.data() + full_start, bc1_obs_len,
                    bc1_target.data(), static_cast<int>(bc1_target.size()),
                    cfg.affine);
            }

            ObservedAlignmentChoice bc2_choice = best_observed_alignment_choice(
                seq, c2.start, bc2_target, cfg.bc2_offset_window,
                cfg.bc2_full_hamming, cfg.affine,
                cfg.score_mode == CandidateScoreMode::AffineGap,
                cfg.score_mode == CandidateScoreMode::AffineGap ||
                    cfg.emit_decode_annotations);
            bc2_full_d = bc2_choice.edit_distance;
            if (bc2_full_d > cfg.bc2_full_hamming) continue;
            bc2_obs_len = bc2_choice.obs_len;
            bc2_affine = bc2_choice.affine;

            if (bc1_obs_len <= 0 || full_start + bc1_obs_len > slen ||
                bc2_obs_len <= 0 || bc2_align_start + bc2_obs_len > slen) {
                continue;
            }
            if (is_oracle_pair) {
                out.oracle_pair_bc1_pass = true;
                out.oracle_pair_bc2_pass = true;
            }

            int anchor_score = 0;
            int offset_sum = 0;
            if (!anchor_score_for_pair(full_start, bc2_boundary, bc2_align_start,
                                       bc2_obs_len, row2, col2,
                                       anchor_score, offset_sum)) {
                continue;
            }
            offset_sum += c2.offset;
            if (is_oracle_pair) out.oracle_pair_anchor_pass = true;

            const bool bc1_overlong_exact_prefix = has_overlong_exact_prefix_overlap(
                seq.data() + full_start, bc1_obs_len, bc1_target,
                std::max(1, static_cast<int>(bc1_target.size()) -
                             cfg.bc1_full_hamming),
                bc1_full_d);
            const bool suppress_overlong_prefix_rescue =
                cfg.suppress_bc1_overlong_prefix_e1 &&
                bc1_overlong_exact_prefix && (bc1_full_d + bc2_full_d <= 1);
            const bool suppress_bc2_short_terminal_rescue =
                should_suppress_bc2_short_terminal_e2(
                    cfg, row2, bc2_full_d, bc2_obs_len);
            consider_candidate(out, cfg, bc1_full_d, bc2_full_d,
                               bc1_affine, bc2_affine, anchor_score,
                               offset_sum, full_start, row2, col2,
                               bc1_obs_len, bc2_obs_len, oracle,
                               !(suppress_overlong_prefix_rescue ||
                                 suppress_bc2_short_terminal_rescue));
        }
    }

    out.non_acgt_null_queries = non_acgt_null_queries;
    finalize_decoded_assignment(out);
    return out;
}

struct TieredCandidate {
    int tier;
    int row2;
    int col2;
    int full_start;
    int bc1_edit;
    int bc2_edit;
    int bc1_obs_len;
    int bc2_obs_len;

    TieredCandidate(int candidateTier = 0, int row = -1, int column = -1,
                    int fullStart = -1, int bc1Edit = -1, int bc2Edit = -1,
                    int bc1ObservedLength = -1, int bc2ObservedLength = -1)
        : tier(candidateTier), row2(row), col2(column), full_start(fullStart),
          bc1_edit(bc1Edit), bc2_edit(bc2Edit),
          bc1_obs_len(bc1ObservedLength), bc2_obs_len(bc2ObservedLength) {}
};

std::vector<int> unique_sorted_lengths(const std::vector<int>& lengths) {
    std::vector<int> out;
    out.reserve(lengths.size());
    for (int length : lengths) {
        if (length > 0) out.push_back(length);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<int> one_edit_query_lengths(const std::vector<int>& target_lengths) {
    std::vector<int> out;
    out.reserve(target_lengths.size() * 3);
    for (int length : target_lengths) {
        for (int delta = -1; delta <= 1; ++delta) {
            const int observed_len = length + delta;
            if (observed_len > 0) out.push_back(observed_len);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<int> edit_query_lengths(const std::vector<int>& target_lengths,
                                    int max_dist) {
    std::vector<int> out;
    out.reserve(target_lengths.size() * static_cast<size_t>(max_dist * 2 + 1));
    for (int length : target_lengths) {
        for (int delta = -max_dist; delta <= max_dist; ++delta) {
            const int observed_len = length + delta;
            if (observed_len > 0) out.push_back(observed_len);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

bool tiered_candidate_coord_seen(const std::vector<TieredCandidate>& candidates,
                                 int tier, int row2, int col2) {
    for (const TieredCandidate& candidate : candidates) {
        if (candidate.tier == tier && candidate.row2 == row2 && candidate.col2 == col2) {
            return true;
        }
    }
    return false;
}

void set_tiered_score_key(int bc1_edit, int bc2_edit, int tier, long key[4]) {
    key[0] = static_cast<long>(bc1_edit) + static_cast<long>(bc2_edit);
    key[1] = tier;
    key[2] = 0;
    key[3] = 0;
}

// Encode (substitutions, deletions, insertions) for edit distances up to two.
// Equal-length edit-distance-two alignments are the only case with two possible
// operation-count signatures: SS and DI.
uint8_t edit_operation_signature(int substitutions, int deletions, int insertions) {
    return static_cast<uint8_t>(substitutions * 9 + deletions * 3 + insertions);
}

bool has_one_delete_one_insert_alignment(const std::string& observed,
                                         const std::string& target) {
    if (observed.size() != target.size() || observed.empty()) return false;
    for (size_t observed_skip = 0; observed_skip < observed.size(); ++observed_skip) {
        for (size_t target_skip = 0; target_skip < target.size(); ++target_skip) {
            size_t oi = 0;
            size_t ti = 0;
            bool equal = true;
            while (oi < observed.size() && ti < target.size()) {
                if (oi == observed_skip) ++oi;
                if (ti == target_skip) ++ti;
                if (oi == observed.size() || ti == target.size()) break;
                if (observed[oi] != target[ti]) {
                    equal = false;
                    break;
                }
                ++oi;
                ++ti;
            }
            if (equal) return true;
        }
    }
    return false;
}

std::vector<uint8_t> optimal_edit_operation_signatures(const std::string& observed,
                                                       const std::string& target,
                                                       int edit_distance) {
    std::vector<uint8_t> signatures;
    const int length_delta = static_cast<int>(observed.size()) -
                             static_cast<int>(target.size());
    if (edit_distance == 0 && length_delta == 0) {
        signatures.push_back(edit_operation_signature(0, 0, 0));
    } else if (edit_distance == 1) {
        if (length_delta == 0) signatures.push_back(edit_operation_signature(1, 0, 0));
        if (length_delta == -1) signatures.push_back(edit_operation_signature(0, 1, 0));
        if (length_delta == 1) signatures.push_back(edit_operation_signature(0, 0, 1));
    } else if (edit_distance == 2) {
        if (length_delta == -2) signatures.push_back(edit_operation_signature(0, 2, 0));
        if (length_delta == -1) signatures.push_back(edit_operation_signature(1, 1, 0));
        if (length_delta == 0) {
            int hamming = 0;
            for (size_t i = 0; i < observed.size(); ++i) {
                hamming += observed[i] != target[i];
            }
            if (hamming == 2) signatures.push_back(edit_operation_signature(2, 0, 0));
            if (has_one_delete_one_insert_alignment(observed, target)) {
                signatures.push_back(edit_operation_signature(0, 1, 1));
            }
        }
        if (length_delta == 1) signatures.push_back(edit_operation_signature(1, 0, 1));
        if (length_delta == 2) signatures.push_back(edit_operation_signature(0, 0, 2));
    }
    std::sort(signatures.begin(), signatures.end());
    signatures.erase(std::unique(signatures.begin(), signatures.end()), signatures.end());
    return signatures;
}

bool candidate_group_has_ambiguous_registration(
    const std::string& seq, const Config& cfg,
    const std::vector<TieredCandidate>& candidates) {
    std::vector<std::pair<int, int>> coords;
    for (const TieredCandidate& candidate : candidates) {
        const std::pair<int, int> coord(candidate.row2, candidate.col2);
        if (std::find(coords.begin(), coords.end(), coord) == coords.end()) {
            coords.push_back(coord);
        }
    }

    for (const auto& coord : coords) {
        std::vector<std::pair<uint8_t, uint8_t>> combined_signatures;
        for (const TieredCandidate& candidate : candidates) {
            if (candidate.row2 != coord.first || candidate.col2 != coord.second) continue;
            const int bc2_start = candidate.full_start + candidate.bc1_obs_len;
            if (candidate.full_start < 0 || candidate.bc1_obs_len <= 0 ||
                candidate.bc2_obs_len <= 0 || bc2_start < 0 ||
                bc2_start + candidate.bc2_obs_len > static_cast<int>(seq.size())) {
                return true;
            }
            const std::string observed_bc1 = seq.substr(
                static_cast<size_t>(candidate.full_start),
                static_cast<size_t>(candidate.bc1_obs_len));
            const std::string observed_bc2 = seq.substr(
                static_cast<size_t>(bc2_start),
                static_cast<size_t>(candidate.bc2_obs_len));
            const auto bc1_signatures = optimal_edit_operation_signatures(
                observed_bc1, cfg.bc1_oligos[candidate.col2], candidate.bc1_edit);
            const auto bc2_signatures = optimal_edit_operation_signatures(
                observed_bc2, cfg.bc2_oligos[candidate.row2], candidate.bc2_edit);
            if (bc1_signatures.empty() || bc2_signatures.empty()) return true;
            for (uint8_t bc1_signature : bc1_signatures) {
                for (uint8_t bc2_signature : bc2_signatures) {
                    const std::pair<uint8_t, uint8_t> signature(
                        bc1_signature, bc2_signature);
                    if (std::find(combined_signatures.begin(), combined_signatures.end(),
                                  signature) == combined_signatures.end()) {
                        combined_signatures.push_back(signature);
                        if (combined_signatures.size() > 1) return true;
                    }
                }
            }
        }
        if (combined_signatures.size() != 1) return true;
    }
    return false;
}

using SupportProduct = unsigned __int128;

struct FrozenSupportRanking {
    bool unique = false;
    int row2 = -1;
    int col2 = -1;
    double odds = 0.0;
};

FrozenSupportRanking rank_frozen_support(
    const std::vector<TieredCandidate>& unique_coords,
    const std::vector<uint64_t>& bc1_counts,
    const std::vector<uint64_t>& bc2_counts) {
    FrozenSupportRanking ranking;
    SupportProduct best = 0;
    SupportProduct second = 0;
    int best_count = 0;
    for (const TieredCandidate& candidate : unique_coords) {
        const SupportProduct score =
            static_cast<SupportProduct>(bc1_counts[candidate.col2] + 1) *
            static_cast<SupportProduct>(bc2_counts[candidate.row2] + 1);
        if (score > best) {
            second = best;
            best = score;
            ranking.row2 = candidate.row2;
            ranking.col2 = candidate.col2;
            best_count = 1;
        } else if (score == best) {
            ++best_count;
        } else if (score > second) {
            second = score;
        }
    }
    ranking.unique = best_count == 1;
    if (ranking.unique && second > 0) {
        ranking.odds = static_cast<double>(static_cast<long double>(best) /
                                           static_cast<long double>(second));
    }
    return ranking;
}

void apply_frozen_support_consensus(
    Decoded& out, const std::string& seq, const Config& cfg,
    const std::vector<TieredCandidate>& all_candidates,
    const std::vector<TieredCandidate>& unique_coords) {
    out.frozen_support_attempted = true;
    if (candidate_group_has_ambiguous_registration(seq, cfg, all_candidates)) {
        out.frozen_support_status = "registration_ambiguous";
        return;
    }
    const FrozenSupportRanking h0 = rank_frozen_support(
        unique_coords, cfg.bc1_frozen_h0, cfg.bc2_frozen_h0);
    if (!h0.unique) {
        out.frozen_support_status = "h0_tie";
        return;
    }
    const FrozenSupportRanking exposure = rank_frozen_support(
        unique_coords, cfg.bc1_frozen_exposure, cfg.bc2_frozen_exposure);
    if (!exposure.unique) {
        out.frozen_support_status = "exposure_tie";
        return;
    }
    out.frozen_h0_odds = h0.odds;
    out.frozen_exposure_odds = exposure.odds;
    if (h0.row2 != exposure.row2 || h0.col2 != exposure.col2) {
        out.frozen_support_status = "support_conflict";
        return;
    }
    if (!(h0.odds > cfg.frozen_support_min_odds) ||
        !(exposure.odds > cfg.frozen_support_min_odds)) {
        out.frozen_support_status = "below_min_odds";
        return;
    }
    const TieredCandidate* winner = nullptr;
    for (const TieredCandidate& candidate : all_candidates) {
        if (candidate.row2 == h0.row2 && candidate.col2 == h0.col2) {
            winner = &candidate;
            break;
        }
    }
    if (!winner) {
        out.frozen_support_status = "internal_missing_winner";
        return;
    }
    out.assigned = true;
    out.frozen_support_resolved = true;
    out.frozen_support_status = "resolved";
    out.row2 = winner->row2;
    out.col2 = winner->col2;
    out.full_start = winner->full_start;
    out.bc1_edit = winner->bc1_edit;
    out.bc2_edit = winner->bc2_edit;
    out.bc1_obs_len = winner->bc1_obs_len;
    out.bc2_obs_len = winner->bc2_obs_len;
}

void collect_suffix_tolerant_fixed_hits(const PackedAnswerCache& cache,
                                        const std::string& seq,
                                        int start,
                                        int fixed_width,
                                        std::vector<PackedHit>& out) {
    out.clear();
    const int slen = static_cast<int>(seq.size());
    if (start < 0 || start >= slen || fixed_width <= 0 || fixed_width > kMaxLen) return;

    const int available = std::min(fixed_width, slen - start);
    int prefix_len = 0;
    while (prefix_len < available && base_bits(seq[start + prefix_len]) >= 0) {
        ++prefix_len;
    }

    if (prefix_len == fixed_width) {
        HitSpan span = cache.lookup_span(seq.data() + start, fixed_width);
        if (span.found && span.count > 0) {
            out.insert(out.end(), span.data, span.data + span.count);
        }
        return;
    }

    const int suffix_len = fixed_width - prefix_len;
    if (suffix_len < 0 || suffix_len > 6) return;

    char query[kMaxLen + 1];
    if (prefix_len > 0) {
        std::memcpy(query, seq.data() + start, static_cast<size_t>(prefix_len));
    }
    const uint64_t suffix_count = 1ULL << (2 * suffix_len);
    for (uint64_t suffix = 0; suffix < suffix_count; ++suffix) {
        uint64_t suffix_bits = suffix;
        for (int pos = suffix_len - 1; pos >= 0; --pos) {
            query[prefix_len + pos] = BASES[suffix_bits & 0x3ULL];
            suffix_bits >>= 2;
        }
        HitSpan span = cache.lookup_span(query, fixed_width);
        if (!span.found || span.count == 0) continue;
        for (uint16_t i = 0; i < span.count; ++i) {
            if (static_cast<int>(span.data[i].reserved) <= prefix_len) {
                out.push_back(span.data[i]);
            }
        }
    }
    if (out.empty()) return;

    std::sort(out.begin(), out.end(),
              [](const PackedHit& a, const PackedHit& b) {
                  if (a.idx != b.idx) return a.idx < b.idx;
                  if (a.reserved != b.reserved) return a.reserved < b.reserved;
                  return a.distance < b.distance;
              });
    std::vector<PackedHit> deduped;
    deduped.reserve(out.size());
    size_t i = 0;
    while (i < out.size()) {
        PackedHit best = out[i];
        ++i;
        while (i < out.size() &&
               out[i].idx == best.idx &&
               out[i].reserved == best.reserved) {
            best.distance = std::min(best.distance, out[i].distance);
            ++i;
        }
        deduped.push_back(best);
    }
    out.swap(deduped);
}

Decoded decode_record_direct_tiered(const std::string& seq, Config& cfg,
                                    const std::vector<int>& bc1_exact_lengths,
                                    const std::vector<int>& bc2_exact_lengths,
                                    const std::vector<int>& bc1_e1_query_lengths,
                                    const std::vector<int>& bc2_e1_query_lengths,
                                    const OracleTarget* oracle = nullptr) {
    const int slen = static_cast<int>(seq.size());
    Decoded out;
    std::vector<TieredCandidate> candidates;

    auto add_candidate = [&](int tier, int row2, int col2, int full_start,
                             int bc1_edit, int bc2_edit,
                             int bc1_obs_len, int bc2_obs_len) {
        if (row2 < 0 || row2 >= cfg.grid_rows || col2 < 0 || col2 >= cfg.grid_cols) {
            return;
        }
        ++out.candidate_pairs;
        candidates.push_back({tier, row2, col2, full_start, bc1_edit, bc2_edit,
                              bc1_obs_len, bc2_obs_len});
        if (oracle && oracle->resolved && row2 == oracle->row2 && col2 == oracle->col2) {
            long key[4] = {0, 0, 0, 0};
            set_tiered_score_key(bc1_edit, bc2_edit, tier, key);
            if (!out.oracle_candidate_present ||
                compare_score_keys(key, out.oracle_score_key) < 0) {
                out.oracle_candidate_present = true;
                std::memcpy(out.oracle_score_key, key, sizeof(out.oracle_score_key));
                out.oracle_full_start = full_start;
                out.oracle_bc1_edit = bc1_edit;
                out.oracle_bc2_edit = bc2_edit;
                out.oracle_anchor_score = 0;
                out.oracle_offset_sum = 0;
                out.oracle_bc1_obs_len = bc1_obs_len;
                out.oracle_bc2_obs_len = bc2_obs_len;
            }
        }
    };

    auto add_exact_bc2_candidates = [&](int tier, int full_start, int bc1_obs_len,
                                        int col2, int bc1_edit, int bc2_start) {
        for (int bc2_obs_len : bc2_exact_lengths) {
            if (bc2_start < 0 || bc2_start + bc2_obs_len > slen) continue;
            HitSpan bc2_hits = cfg.bc2_tiered_exact_lookup.lookup_span(
                seq.data() + bc2_start, bc2_obs_len);
            if (!bc2_hits.found || bc2_hits.count == 0) continue;
            for (uint16_t i = 0; i < bc2_hits.count; ++i) {
                const PackedHit& c2 = bc2_hits.data[i];
                if (c2.distance != 0) continue;
                add_candidate(tier, c2.idx, col2, full_start, bc1_edit, 0,
                              bc1_obs_len, bc2_obs_len);
            }
        }
    };

    auto add_e1_bc2_candidates = [&](int tier, int full_start, int bc1_obs_len,
                                     int col2, int bc2_start) {
        for (int bc2_obs_len : bc2_e1_query_lengths) {
            if (bc2_start < 0 || bc2_start + bc2_obs_len > slen) continue;
            HitSpan bc2_hits = cfg.bc2_tiered_e1_lookup.lookup_span(
                seq.data() + bc2_start, bc2_obs_len);
            if (!bc2_hits.found || bc2_hits.count == 0) continue;
            for (uint16_t i = 0; i < bc2_hits.count; ++i) {
                const PackedHit& c2 = bc2_hits.data[i];
                if (c2.distance != 1) continue;
                add_candidate(tier, c2.idx, col2, full_start, 0, 1,
                              bc1_obs_len, bc2_obs_len);
            }
        }
    };

    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        if (full_start < 0 || full_start >= slen) continue;

        // Tier 0 and the 0+1 side of Tier 1: exact BC1 followed immediately by BC2.
        for (int bc1_obs_len : bc1_exact_lengths) {
            if (full_start + bc1_obs_len > slen) continue;
            HitSpan bc1_hits = cfg.bc1_tiered_exact_lookup.lookup_span(
                seq.data() + full_start, bc1_obs_len);
            if (!bc1_hits.found || bc1_hits.count == 0) continue;
            const int bc2_start = full_start + bc1_obs_len;
            for (uint16_t i = 0; i < bc1_hits.count; ++i) {
                const PackedHit& c1 = bc1_hits.data[i];
                if (c1.distance != 0) continue;
                add_exact_bc2_candidates(0, full_start, bc1_obs_len, c1.idx, 0, bc2_start);
                add_e1_bc2_candidates(1, full_start, bc1_obs_len, c1.idx, bc2_start);
            }
        }

        // The 1+0 side of Tier 1: one-edit BC1 followed immediately by exact BC2.
        for (int bc1_obs_len : bc1_e1_query_lengths) {
            if (full_start + bc1_obs_len > slen) continue;
            HitSpan bc1_hits = cfg.bc1_tiered_e1_lookup.lookup_span(
                seq.data() + full_start, bc1_obs_len);
            if (!bc1_hits.found || bc1_hits.count == 0) continue;
            const int bc2_start = full_start + bc1_obs_len;
            for (uint16_t i = 0; i < bc1_hits.count; ++i) {
                const PackedHit& c1 = bc1_hits.data[i];
                if (c1.distance != 1) continue;
                add_exact_bc2_candidates(1, full_start, bc1_obs_len, c1.idx, 1, bc2_start);
            }
        }
    }

    if (candidates.empty()) return out;

    const TieredCandidate* best = nullptr;
    int best_tier = std::numeric_limits<int>::max();
    for (const TieredCandidate& candidate : candidates) {
        if (candidate.tier < best_tier) {
            best_tier = candidate.tier;
            best = &candidate;
        }
    }
    if (!best) return out;

    out.have_best = true;
    set_tiered_score_key(best->bc1_edit, best->bc2_edit, best->tier, out.score_key);
    out.row2 = best->row2;
    out.col2 = best->col2;
    out.full_start = best->full_start;
    out.bc1_edit = best->bc1_edit;
    out.bc2_edit = best->bc2_edit;
    out.anchor_score = 0;
    out.offset_sum = 0;
    out.bc1_obs_len = best->bc1_obs_len;
    out.bc2_obs_len = best->bc2_obs_len;

    std::vector<TieredCandidate> unique_coords;
    for (const TieredCandidate& candidate : candidates) {
        if (candidate.tier != best_tier) continue;
        append_candidate_trace(out.best_candidate_traces, candidate.row2, candidate.col2,
                               candidate.full_start, candidate.bc1_edit, candidate.bc2_edit,
                               0, 0, candidate.bc1_obs_len, candidate.bc2_obs_len);
        if (!tiered_candidate_coord_seen(unique_coords, candidate.tier,
                                         candidate.row2, candidate.col2)) {
            unique_coords.push_back(candidate);
        }
    }

    out.best_unique = unique_coords.size() == 1;
    out.best_tie_candidates =
        unique_coords.empty() ? 0 : static_cast<uint32_t>(unique_coords.size() - 1);
    if (out.best_unique) {
        const TieredCandidate& unique = unique_coords.front();
        out.assigned = true;
        out.row2 = unique.row2;
        out.col2 = unique.col2;
        out.full_start = unique.full_start;
        out.bc1_edit = unique.bc1_edit;
        out.bc2_edit = unique.bc2_edit;
        out.bc1_obs_len = unique.bc1_obs_len;
        out.bc2_obs_len = unique.bc2_obs_len;
        set_tiered_score_key(unique.bc1_edit, unique.bc2_edit, unique.tier,
                             out.score_key);
    }
    return out;
}

Decoded decode_record_direct_tiered_fixed(const std::string& seq, Config& cfg,
                                          int bc1_fixed_width,
                                          int bc2_fixed_width,
                                          const OracleTarget* oracle = nullptr) {
    const int slen = static_cast<int>(seq.size());
    Decoded out;
    std::vector<TieredCandidate> candidates;

    auto add_candidate = [&](int tier, int row2, int col2, int full_start,
                             int bc1_edit, int bc2_edit,
                             int bc1_obs_len, int bc2_obs_len) {
        if (row2 < 0 || row2 >= cfg.grid_rows || col2 < 0 || col2 >= cfg.grid_cols) {
            return;
        }
        ++out.candidate_pairs;
        candidates.push_back({tier, row2, col2, full_start, bc1_edit, bc2_edit,
                              bc1_obs_len, bc2_obs_len});
        if (oracle && oracle->resolved && row2 == oracle->row2 && col2 == oracle->col2) {
            long key[4] = {0, 0, 0, 0};
            set_tiered_score_key(bc1_edit, bc2_edit, tier, key);
            if (!out.oracle_candidate_present ||
                compare_score_keys(key, out.oracle_score_key) < 0) {
                out.oracle_candidate_present = true;
                std::memcpy(out.oracle_score_key, key, sizeof(out.oracle_score_key));
                out.oracle_full_start = full_start;
                out.oracle_bc1_edit = bc1_edit;
                out.oracle_bc2_edit = bc2_edit;
                out.oracle_anchor_score = 0;
                out.oracle_offset_sum = 0;
                out.oracle_bc1_obs_len = bc1_obs_len;
                out.oracle_bc2_obs_len = bc2_obs_len;
            }
        }
    };

    std::vector<PackedHit> bc1_fixed_hits;
    std::vector<PackedHit> bc2_fixed_hits;
    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        collect_suffix_tolerant_fixed_hits(cfg.bc1_tiered_fixed_lookup, seq,
                                           full_start, bc1_fixed_width,
                                           bc1_fixed_hits);
        if (bc1_fixed_hits.empty()) continue;
        for (const PackedHit& c1 : bc1_fixed_hits) {
            const int bc1_edit = static_cast<int>(c1.distance);
            const int bc1_obs_len = static_cast<int>(c1.reserved);
            if (bc1_edit > 1 || bc1_obs_len <= 0) continue;
            const int bc2_start = full_start + bc1_obs_len;
            collect_suffix_tolerant_fixed_hits(cfg.bc2_tiered_fixed_lookup, seq,
                                               bc2_start, bc2_fixed_width,
                                               bc2_fixed_hits);
            if (bc2_fixed_hits.empty()) continue;
            for (const PackedHit& c2 : bc2_fixed_hits) {
                const int bc2_edit = static_cast<int>(c2.distance);
                const int bc2_obs_len = static_cast<int>(c2.reserved);
                if (bc2_edit > 1 || bc2_obs_len <= 0) continue;
                const int tier = bc1_edit + bc2_edit;
                if (tier > 1) continue;
                add_candidate(tier, c2.idx, c1.idx, full_start,
                              bc1_edit, bc2_edit, bc1_obs_len, bc2_obs_len);
            }
        }
    }

    if (candidates.empty()) return out;

    const TieredCandidate* best = nullptr;
    int best_tier = std::numeric_limits<int>::max();
    for (const TieredCandidate& candidate : candidates) {
        if (candidate.tier < best_tier) {
            best_tier = candidate.tier;
            best = &candidate;
        }
    }
    if (!best) return out;

    out.have_best = true;
    set_tiered_score_key(best->bc1_edit, best->bc2_edit, best->tier, out.score_key);
    out.row2 = best->row2;
    out.col2 = best->col2;
    out.full_start = best->full_start;
    out.bc1_edit = best->bc1_edit;
    out.bc2_edit = best->bc2_edit;
    out.anchor_score = 0;
    out.offset_sum = 0;
    out.bc1_obs_len = best->bc1_obs_len;
    out.bc2_obs_len = best->bc2_obs_len;

    std::vector<TieredCandidate> unique_coords;
    for (const TieredCandidate& candidate : candidates) {
        if (candidate.tier != best_tier) continue;
        append_candidate_trace(out.best_candidate_traces, candidate.row2, candidate.col2,
                               candidate.full_start, candidate.bc1_edit, candidate.bc2_edit,
                               0, 0, candidate.bc1_obs_len, candidate.bc2_obs_len);
        if (!tiered_candidate_coord_seen(unique_coords, candidate.tier,
                                         candidate.row2, candidate.col2)) {
            unique_coords.push_back(candidate);
        }
    }

    out.best_unique = unique_coords.size() == 1;
    out.best_tie_candidates =
        unique_coords.empty() ? 0 : static_cast<uint32_t>(unique_coords.size() - 1);
    if (out.best_unique) {
        const TieredCandidate& unique = unique_coords.front();
        out.assigned = true;
        out.row2 = unique.row2;
        out.col2 = unique.col2;
        out.full_start = unique.full_start;
        out.bc1_edit = unique.bc1_edit;
        out.bc2_edit = unique.bc2_edit;
        out.bc1_obs_len = unique.bc1_obs_len;
        out.bc2_obs_len = unique.bc2_obs_len;
        set_tiered_score_key(unique.bc1_edit, unique.bc2_edit, unique.tier,
                             out.score_key);
    }
    return out;
}

HitSpan tiered_target_scan_span(
    const std::string& seq, int start, int observed_len,
    const std::vector<std::string>& targets, int max_distance,
    std::map<std::pair<int, int>, std::vector<PackedHit> >& cache) {
    const std::pair<int, int> key(start, observed_len);
    std::map<std::pair<int, int>, std::vector<PackedHit> >::iterator found =
        cache.find(key);
    if (found == cache.end()) {
        std::vector<PackedHit> hits;
        if (start >= 0 && observed_len > 0
            && start + observed_len <= static_cast<int>(seq.size())) {
            for (size_t index = 0; index < targets.size(); ++index) {
                const std::string& target = targets[index];
                if (std::abs(observed_len - static_cast<int>(target.size()))
                    > max_distance) {
                    continue;
                }
                const int distance = edit_distance_bounded(
                    seq.data() + start, observed_len, target.data(),
                    static_cast<int>(target.size()), max_distance);
                if (distance <= max_distance) {
                    hits.push_back(PackedHit(static_cast<uint16_t>(index),
                                             static_cast<uint8_t>(distance), 0));
                }
            }
        }
        found = cache.insert(std::make_pair(key, hits)).first;
    }
    return HitSpan(found->second.empty() ? nullptr : found->second.data(),
                   static_cast<uint16_t>(found->second.size()), true);
}

Decoded decode_record_direct_tiered_h2(const std::string& seq, Config& cfg,
                                       const std::vector<int>& bc1_query_lengths,
                                       const std::vector<int>& bc2_query_lengths,
                                       const OracleTarget* oracle = nullptr,
                                       const DecodeWindowClassification*
                                           window_classification = nullptr) {
    const int slen = static_cast<int>(seq.size());
    Decoded out;
    const bool audit_oracle = oracle && oracle->resolved;
    const bool target_scan = window_classification != nullptr
        ? (window_classification->n_count != 0 ||
           window_classification->unsupported)
        : direct_best_decode_window_has_non_acgt(seq, cfg);
    std::map<std::pair<int, int>, std::vector<PackedHit> > bc1_target_cache;
    std::map<std::pair<int, int>, std::vector<PackedHit> > bc2_target_cache;
    if (target_scan) out.non_acgt_dp_checked = true;

    auto note_oracle_candidate = [&](int tier, int row2, int col2, int full_start,
                                     int bc1_edit, int bc2_edit,
                                     int bc1_obs_len, int bc2_obs_len) {
        if (!audit_oracle || row2 != oracle->row2 || col2 != oracle->col2) return;
        long key[4] = {0, 0, 0, 0};
        set_tiered_score_key(bc1_edit, bc2_edit, tier, key);
        if (!out.oracle_candidate_present ||
            compare_score_keys(key, out.oracle_score_key) < 0) {
            out.oracle_candidate_present = true;
            std::memcpy(out.oracle_score_key, key, sizeof(out.oracle_score_key));
            out.oracle_full_start = full_start;
            out.oracle_bc1_edit = bc1_edit;
            out.oracle_bc2_edit = bc2_edit;
            out.oracle_anchor_score = 0;
            out.oracle_offset_sum = 0;
            out.oracle_bc1_obs_len = bc1_obs_len;
            out.oracle_bc2_obs_len = bc2_obs_len;
        }
    };

    auto add_candidate = [&](std::vector<TieredCandidate>& candidates,
                             bool collect_for_assignment,
                             int tier, int row2, int col2, int full_start,
                             int bc1_edit, int bc2_edit,
                             int bc1_obs_len, int bc2_obs_len) {
        if (row2 < 0 || row2 >= cfg.grid_rows || col2 < 0 || col2 >= cfg.grid_cols) {
            return;
        }
        if (!cfg.latent_candidate_out.empty()) {
            append_candidate_trace(out.latent_candidate_traces, row2, col2, full_start,
                                   bc1_edit, bc2_edit, 0, 0,
                                   bc1_obs_len, bc2_obs_len);
        }
        if (collect_for_assignment) {
            ++out.candidate_pairs;
            candidates.push_back({tier, row2, col2, full_start, bc1_edit, bc2_edit,
                                  bc1_obs_len, bc2_obs_len});
        }
        note_oracle_candidate(tier, row2, col2, full_start, bc1_edit, bc2_edit,
                              bc1_obs_len, bc2_obs_len);
    };

    std::vector<TieredCandidate> best_tier_candidates;
    bool have_assignment_tier = false;
    for (int tier = 0; tier <= 4; ++tier) {
        std::vector<TieredCandidate> tier_candidates;
        const bool collect_for_assignment = !have_assignment_tier;
        for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max;
             ++full_start) {
            if (full_start < 0 || full_start >= slen) continue;
            for (int bc1_obs_len : bc1_query_lengths) {
                if (full_start + bc1_obs_len > slen) continue;
                HitSpan bc1_hits = target_scan
                    ? tiered_target_scan_span(
                          seq, full_start, bc1_obs_len, cfg.bc1_oligos, 2,
                          bc1_target_cache)
                    : cfg.bc1_tiered_h2_lookup.lookup_span(
                          seq.data() + full_start, bc1_obs_len);
                if (!bc1_hits.found || bc1_hits.count == 0) continue;
                const int bc2_start = full_start + bc1_obs_len;
                if (bc2_start < 0 || bc2_start >= slen) continue;
                for (uint16_t i = 0; i < bc1_hits.count; ++i) {
                    const PackedHit& c1 = bc1_hits.data[i];
                    const int bc1_edit = static_cast<int>(c1.distance);
                    if (bc1_edit < 0 || bc1_edit > 2 || bc1_edit > tier) continue;
                    const int needed_bc2_edit = tier - bc1_edit;
                    if (needed_bc2_edit < 0 || needed_bc2_edit > 2) continue;
                    for (int bc2_obs_len : bc2_query_lengths) {
                        if (bc2_start + bc2_obs_len > slen) continue;
                        HitSpan bc2_hits = target_scan
                            ? tiered_target_scan_span(
                                  seq, bc2_start, bc2_obs_len, cfg.bc2_oligos, 2,
                                  bc2_target_cache)
                            : cfg.bc2_tiered_h2_lookup.lookup_span(
                                  seq.data() + bc2_start, bc2_obs_len);
                        if (!bc2_hits.found || bc2_hits.count == 0) continue;
                        for (uint16_t j = 0; j < bc2_hits.count; ++j) {
                            const PackedHit& c2 = bc2_hits.data[j];
                            const int bc2_edit = static_cast<int>(c2.distance);
                            if (bc2_edit != needed_bc2_edit) continue;
                            add_candidate(tier_candidates, collect_for_assignment,
                                          tier, c2.idx, c1.idx,
                                          full_start, bc1_edit, bc2_edit,
                                          bc1_obs_len, bc2_obs_len);
                        }
                    }
                }
            }
        }
        if (!have_assignment_tier && !tier_candidates.empty()) {
            best_tier_candidates.swap(tier_candidates);
            have_assignment_tier = true;
            if (!audit_oracle && cfg.latent_candidate_out.empty()) break;
        }
    }

    if (best_tier_candidates.empty()) return out;

    const TieredCandidate& first = best_tier_candidates.front();
    out.have_best = true;
    set_tiered_score_key(first.bc1_edit, first.bc2_edit, first.tier, out.score_key);
    out.row2 = first.row2;
    out.col2 = first.col2;
    out.full_start = first.full_start;
    out.bc1_edit = first.bc1_edit;
    out.bc2_edit = first.bc2_edit;
    out.anchor_score = 0;
    out.offset_sum = 0;
    out.bc1_obs_len = first.bc1_obs_len;
    out.bc2_obs_len = first.bc2_obs_len;

    std::vector<TieredCandidate> unique_coords;
    for (const TieredCandidate& candidate : best_tier_candidates) {
        append_candidate_trace(out.best_candidate_traces, candidate.row2, candidate.col2,
                               candidate.full_start, candidate.bc1_edit, candidate.bc2_edit,
                               0, 0, candidate.bc1_obs_len, candidate.bc2_obs_len);
        if (!tiered_candidate_coord_seen(unique_coords, candidate.tier,
                                         candidate.row2, candidate.col2)) {
            unique_coords.push_back(candidate);
        }
    }

    out.best_unique = unique_coords.size() == 1;
    out.best_tie_candidates =
        unique_coords.empty() ? 0 : static_cast<uint32_t>(unique_coords.size() - 1);
    if (out.best_unique) {
        const TieredCandidate& unique = unique_coords.front();
        out.assigned = true;
        out.row2 = unique.row2;
        out.col2 = unique.col2;
        out.full_start = unique.full_start;
        out.bc1_edit = unique.bc1_edit;
        out.bc2_edit = unique.bc2_edit;
        out.bc1_obs_len = unique.bc1_obs_len;
        out.bc2_obs_len = unique.bc2_obs_len;
        set_tiered_score_key(unique.bc1_edit, unique.bc2_edit, unique.tier,
                             out.score_key);
    } else if (cfg.frozen_support_consensus) {
        apply_frozen_support_consensus(out, seq, cfg, best_tier_candidates,
                                       unique_coords);
    }
    return out;
}

Decoded decode_record_prebuilt(const std::string& seq, Config& cfg,
                               std::vector<FullCandidate>& bc1_valid,
                               std::vector<FullCandidate>& bc2_valid,
                               ObservedDecodeCacheLocals& observed_locals,
                               bool fallback_on_cache_miss,
                               bool* cache_miss) {
    const int anchor_len = cfg.anchor_len;
    const int slen = static_cast<int>(seq.size());
    bool have_best = false;
    long best_key[4] = {0, 0, 0, 0};
    int best_row = -1, best_col = -1;
    bool best_unique = true;

    auto consider = [&](int full_h, int anchor_h, int offset_sum, int full_start, int row2, int col2) {
        long key[4] = {full_h, anchor_h, offset_sum, full_start};
        if (!have_best) {
            have_best = true;
            std::memcpy(best_key, key, sizeof(best_key));
            best_row = row2;
            best_col = col2;
            best_unique = true;
            return;
        }
        int cmp = 0;
        for (int i = 0; i < 4; ++i) {
            if (key[i] != best_key[i]) {
                cmp = (key[i] < best_key[i]) ? -1 : 1;
                break;
            }
        }
        if (cmp < 0) {
            std::memcpy(best_key, key, sizeof(best_key));
            best_row = row2;
            best_col = col2;
            best_unique = true;
        } else if (cmp == 0) {
            if (row2 != best_row || col2 != best_col) best_unique = false;
        }
    };

    for (int full_start = cfg.full_start_min; full_start <= cfg.full_start_max; ++full_start) {
        if (full_start < 0) continue;
        for (int extra = 1; extra <= 2; ++extra) {
            int bc1_len = anchor_len + extra;
            int bc1_anchor_expected_start = full_start + extra;
            for (int d1 = -cfg.bc1_anchor_offset_window; d1 <= cfg.bc1_anchor_offset_window; ++d1) {
                int bc1_anchor_start = bc1_anchor_expected_start + d1;
                int bc1_anchor_end = bc1_anchor_start + anchor_len;
                if (bc1_anchor_start < full_start || slen < bc1_anchor_end) continue;
                if (!is_acgt(seq, bc1_anchor_start, anchor_len)) continue;
                HitSpan lk1 = cfg.bc1_anchor_prebuilt.lookup_span(seq.data() + bc1_anchor_start, anchor_len);
                if (!lk1.found) {
                    if (fallback_on_cache_miss) {
                        if (cache_miss) *cache_miss = true;
                        return {};
                    }
                    lookup_dp_only(seq.data() + bc1_anchor_start, anchor_len,
                                   cfg.bc1_anchors, &cfg.bc1_anchor_deletion2,
                                   cfg.bc1_anchor_hamming, observed_locals.scratch,
                                   observed_locals.bc1_anchor_dp_miss,
                                   observed_locals.bc1_anchor_dp_hits);
                    lk1 = {observed_locals.bc1_anchor_dp_hits.data(),
                           static_cast<uint16_t>(observed_locals.bc1_anchor_dp_hits.size()),
                           true};
                }
                if (lk1.count == 0) continue;

                int bc2_expected_start = bc1_anchor_end;
                for (int d2 = -cfg.bc2_offset_window; d2 <= cfg.bc2_offset_window; ++d2) {
                    int bc2_anchor_start = bc2_expected_start + d2;
                    int bc2_anchor_end = bc2_anchor_start + cfg.bc2_seg_len;
                    if (bc2_anchor_start < 0 || slen < bc2_anchor_end) continue;
                    int bc1_obs_len = bc2_anchor_start - full_start;
                    if (bc1_obs_len <= 0 || full_start + bc1_obs_len > slen) continue;

                    HitSpan bc1_full_hits = cfg.bc1_full_prebuilt.lookup_span(seq.data() + full_start,
                                                                              bc1_obs_len);
                    if (!bc1_full_hits.found) {
                        if (fallback_on_cache_miss) {
                            if (cache_miss) *cache_miss = true;
                            return {};
                        }
                        lookup_dp_only(seq.data() + full_start, bc1_obs_len,
                                       cfg.bc1_oligos, &cfg.bc1_full_deletion2,
                                       cfg.bc1_full_hamming, observed_locals.scratch,
                                       observed_locals.bc1_full_dp_miss,
                                       observed_locals.bc1_full_dp_hits);
                        bc1_full_hits = {observed_locals.bc1_full_dp_hits.data(),
                                         static_cast<uint16_t>(observed_locals.bc1_full_dp_hits.size()),
                                         true};
                    }
                    bc1_valid.clear();
                    for (uint16_t i = 0; i < lk1.count; ++i) {
                        const PackedHit& c1 = lk1.data[i];
                        int bc1_idx = c1.idx;
                        if (cfg.bc1_len[bc1_idx] != bc1_len) continue;
                        int bc1_full_h = hit_distance_for_idx(bc1_full_hits, bc1_idx,
                                                              cfg.bc1_full_hamming + 1);
                        if (bc1_full_h > cfg.bc1_full_hamming) continue;
                        bc1_valid.push_back({static_cast<uint16_t>(bc1_idx), c1.distance,
                                             static_cast<uint8_t>(bc1_full_h)});
                    }
                    if (bc1_valid.empty()) continue;

                    if (!is_acgt(seq, bc2_anchor_start, cfg.bc2_seg_len)) continue;
                    HitSpan lk2 = cfg.bc2_anchor_prebuilt.lookup_span(seq.data() + bc2_anchor_start,
                                                                      cfg.bc2_seg_len);
                    if (!lk2.found) {
                        if (fallback_on_cache_miss) {
                            if (cache_miss) *cache_miss = true;
                            return {};
                        }
                        lookup_dp_only(seq.data() + bc2_anchor_start, cfg.bc2_seg_len,
                                       cfg.bc2_anchors, &cfg.bc2_anchor_deletion2,
                                       cfg.bc2_anchor_hamming, observed_locals.scratch,
                                       observed_locals.bc2_anchor_dp_miss,
                                       observed_locals.bc2_anchor_dp_hits);
                        lk2 = {observed_locals.bc2_anchor_dp_hits.data(),
                               static_cast<uint16_t>(observed_locals.bc2_anchor_dp_hits.size()),
                               true};
                    }
                    if (lk2.count == 0) continue;

                    bc2_valid.clear();
                    for (uint16_t i = 0; i < lk2.count; ++i) {
                        const PackedHit& c2 = lk2.data[i];
                        int bc2_idx = c2.idx;
                        int bc2_full_h = best_observed_edit_distance_prebuilt(
                            seq, bc2_anchor_start, bc2_idx, cfg.bc2_oligos,
                            cfg.bc2_full_prebuilt, cfg.bc2_full_deletion2,
                            observed_locals.scratch, observed_locals.bc2_full_dp_miss,
                            observed_locals.bc2_full_dp_hits, cfg.bc2_offset_window,
                            cfg.bc2_full_hamming, fallback_on_cache_miss, cache_miss);
                        if (cache_miss && *cache_miss) return {};
                        if (bc2_full_h <= cfg.bc2_full_hamming) {
                            bc2_valid.push_back({static_cast<uint16_t>(bc2_idx), c2.distance,
                                                 static_cast<uint8_t>(bc2_full_h)});
                        }
                    }
                    if (bc2_valid.empty()) continue;

                    for (auto& c1 : bc1_valid) {
                        int bc1_idx = c1.idx;
                        int bc1_anchor_h = c1.anchor_distance;
                        int bc1_full_h = c1.full_distance;
                        for (auto& c2 : bc2_valid) {
                            int bc2_idx = c2.idx;
                            int row2 = bc2_idx, col2 = bc1_idx;
                            if (row2 < 0 || row2 >= cfg.grid_rows || col2 < 0 || col2 >= cfg.grid_cols) continue;
                            int full_h = bc1_full_h + c2.full_distance;
                            int anchor_h = bc1_anchor_h + c2.anchor_distance;
                            int offset_sum = std::abs(d1) + std::abs(d2);
                            consider(full_h, anchor_h, offset_sum, full_start, row2, col2);
                        }
                    }
                }
            }
        }
    }

    Decoded out;
    if (have_best && best_unique) {
        out.assigned = true;
        out.row2 = best_row;
        out.col2 = best_col;
        out.full_start = static_cast<int>(best_key[3]);
    }
    return out;
}

// ---- FASTQ reading (gzip), batched ----
struct Record {
    uint64_t global_ordinal = 0;
    std::string read_id;
    std::string seq;
    std::string qual;
};

std::string read_id_from_header(const char* line) {
    // line starts with '@'; id is up to first whitespace
    const char* p = line + 1;
    const char* q = p;
    while (*q && *q != ' ' && *q != '\t' && *q != '\n' && *q != '\r') ++q;
    return std::string(p, q - p);
}

void load_lines(const std::vector<std::string>& bc_paths, std::vector<std::string>& out) {
    for (auto& path : bc_paths) {
        std::ifstream f(path);
        std::string line;
        while (std::getline(f, line)) {
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back();
            if (!line.empty()) out.push_back(line);
        }
    }
}

std::vector<std::string> split_tab_line(const std::string& line) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= line.size()) {
        size_t pos = line.find('\t', start);
        if (pos == std::string::npos) {
            out.push_back(line.substr(start));
            break;
        }
        out.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return out;
}

bool load_frozen_support_counts(const std::string& path,
                                const std::string& count_column,
                                const std::vector<std::string>& bc1_oligos,
                                const std::vector<std::string>& bc2_oligos,
                                std::vector<uint64_t>& bc1_counts,
                                std::vector<uint64_t>& bc2_counts) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "ERROR: cannot open frozen support TSV " << path << "\n";
        return false;
    }
    std::string header;
    if (!std::getline(in, header)) {
        std::cerr << "ERROR: empty frozen support TSV " << path << "\n";
        return false;
    }
    const auto columns = split_tab_line(header);
    std::unordered_map<std::string, size_t> column_index;
    for (size_t i = 0; i < columns.size(); ++i) column_index[columns[i]] = i;
    for (const std::string& required :
         {std::string("barcode_half"), std::string("oligo_index"),
          std::string("oligo_sequence"), count_column}) {
        if (column_index.find(required) == column_index.end()) {
            std::cerr << "ERROR: frozen support TSV " << path
                      << " is missing column " << required << "\n";
            return false;
        }
    }

    bc1_counts.assign(bc1_oligos.size(), 0);
    bc2_counts.assign(bc2_oligos.size(), 0);
    std::vector<bool> bc1_seen(bc1_oligos.size(), false);
    std::vector<bool> bc2_seen(bc2_oligos.size(), false);
    std::string line;
    size_t line_number = 1;
    while (std::getline(in, line)) {
        ++line_number;
        if (line.empty()) continue;
        const auto fields = split_tab_line(line);
        if (fields.size() < columns.size()) {
            std::cerr << "ERROR: short row in frozen support TSV " << path
                      << " at line " << line_number << "\n";
            return false;
        }
        const std::string& half = fields[column_index["barcode_half"]];
        const std::vector<std::string>* oligos = nullptr;
        std::vector<uint64_t>* counts = nullptr;
        std::vector<bool>* seen = nullptr;
        if (half == "BC1") {
            oligos = &bc1_oligos;
            counts = &bc1_counts;
            seen = &bc1_seen;
        } else if (half == "BC2") {
            oligos = &bc2_oligos;
            counts = &bc2_counts;
            seen = &bc2_seen;
        } else {
            std::cerr << "ERROR: invalid barcode_half in frozen support TSV " << path
                      << " at line " << line_number << "\n";
            return false;
        }
        try {
            const size_t index = static_cast<size_t>(
                std::stoull(fields[column_index["oligo_index"]]));
            const uint64_t count = std::stoull(fields[column_index[count_column]]);
            if (index >= oligos->size() || (*seen)[index]) {
                std::cerr << "ERROR: out-of-range or duplicate oligo index in " << path
                          << " at line " << line_number << "\n";
                return false;
            }
            if (fields[column_index["oligo_sequence"]] != (*oligos)[index]) {
                std::cerr << "ERROR: oligo sequence/index mismatch in " << path
                          << " at line " << line_number << "\n";
                return false;
            }
            (*counts)[index] = count;
            (*seen)[index] = true;
        } catch (const std::exception&) {
            std::cerr << "ERROR: invalid integer in frozen support TSV " << path
                      << " at line " << line_number << "\n";
            return false;
        }
    }
    if (std::find(bc1_seen.begin(), bc1_seen.end(), false) != bc1_seen.end() ||
        std::find(bc2_seen.begin(), bc2_seen.end(), false) != bc2_seen.end()) {
        std::cerr << "ERROR: frozen support TSV does not cover every oligo: " << path
                  << "\n";
        return false;
    }
    return true;
}

std::string strip_10x_barcode_suffix(std::string value) {
    const size_t dash = value.find('-');
    if (dash != std::string::npos) value.resize(dash);
    return value;
}

bool parse_int_field(const std::string& value, int& out) {
    if (value.empty()) return false;
    char* end = nullptr;
    long parsed = std::strtol(value.c_str(), &end, 10);
    if (end == value.c_str() || *end != '\0' ||
        parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max()) {
        return false;
    }
    out = static_cast<int>(parsed);
    return true;
}

bool resolve_cb_to_oracle_target(const std::string& raw_cb, const Config& cfg,
                                 const std::unordered_map<std::string, int>& bc1_by_seq,
                                 const std::unordered_map<std::string, int>& bc2_by_seq,
                                 OracleTarget& target) {
    target.cb = raw_cb;
    const std::string cb = strip_10x_barcode_suffix(raw_cb);
    if (cb.empty()) {
        target.status = "missing_cb";
        return false;
    }
    const std::string spatial_prefix = "s_002um_";
    if (cb.rfind(spatial_prefix, 0) == 0) {
        const std::string rest = cb.substr(spatial_prefix.size());
        const size_t sep = rest.find('_');
        if (sep != std::string::npos) {
            int row2 = -1;
            int col2 = -1;
            if (parse_int_field(rest.substr(0, sep), row2) &&
                parse_int_field(rest.substr(sep + 1), col2) &&
                row2 >= 0 && row2 < static_cast<int>(cfg.bc2_oligos.size()) &&
                col2 >= 0 && col2 < static_cast<int>(cfg.bc1_oligos.size())) {
                target.resolved = true;
                target.row2 = row2;
                target.col2 = col2;
                target.status = "resolved_from_spatial_cb";
                return true;
            }
        }
        target.status = "invalid_spatial_cb";
        return false;
    }

    int matches = 0;
    int matched_row = -1;
    int matched_col = -1;
    bool tried_len[kMaxLen + 1] = {false};
    for (int bc1_len : cfg.bc1_len) {
        if (bc1_len <= 0 || bc1_len >= static_cast<int>(cb.size())) continue;
        if (bc1_len <= kMaxLen && tried_len[bc1_len]) continue;
        if (bc1_len <= kMaxLen) tried_len[bc1_len] = true;
        const std::string bc1 = cb.substr(0, static_cast<size_t>(bc1_len));
        const std::string bc2 = cb.substr(static_cast<size_t>(bc1_len));
        auto bc1_it = bc1_by_seq.find(bc1);
        if (bc1_it == bc1_by_seq.end()) continue;
        auto bc2_it = bc2_by_seq.find(bc2);
        if (bc2_it == bc2_by_seq.end()) continue;
        ++matches;
        matched_col = bc1_it->second;
        matched_row = bc2_it->second;
    }
    if (matches == 1) {
        target.resolved = true;
        target.row2 = matched_row;
        target.col2 = matched_col;
        target.status = "resolved_from_cb";
        return true;
    }
    target.status = matches == 0 ? "cb_not_in_oligo_universe" : "cb_ambiguous_split";
    return false;
}

std::unordered_map<std::string, OracleTarget> load_sr_oracle_tsv(const std::string& path,
                                                                 const Config& cfg) {
    std::unordered_map<std::string, OracleTarget> out;
    if (path.empty()) return out;
    std::ifstream in(path);
    if (!in) {
        std::cerr << "ERROR: cannot open SR oracle TSV " << path << "\n";
        std::exit(1);
    }

    std::unordered_map<std::string, int> bc1_by_seq;
    std::unordered_map<std::string, int> bc2_by_seq;
    for (int i = 0; i < static_cast<int>(cfg.bc1_oligos.size()); ++i) {
        bc1_by_seq.emplace(cfg.bc1_oligos[i], i);
    }
    for (int i = 0; i < static_cast<int>(cfg.bc2_oligos.size()); ++i) {
        bc2_by_seq.emplace(cfg.bc2_oligos[i], i);
    }

    std::string header;
    if (!std::getline(in, header)) return out;
    auto columns = split_tab_line(header);
    std::unordered_map<std::string, size_t> column_idx;
    for (size_t i = 0; i < columns.size(); ++i) column_idx[columns[i]] = i;
    auto has_col = [&](const std::string& name) { return column_idx.find(name) != column_idx.end(); };
    if (!has_col("read_id")) {
        std::cerr << "ERROR: SR oracle TSV is missing read_id column: " << path << "\n";
        std::exit(1);
    }

    const bool has_cb = has_col("cb") || has_col("sr_cb");
    const bool has_row_col = has_col("row2") && has_col("col2");
    const std::string cb_col = has_col("cb") ? "cb" : "sr_cb";
    if (!has_cb && !has_row_col) {
        std::cerr << "ERROR: SR oracle TSV needs cb/sr_cb or row2+col2 columns: "
                  << path << "\n";
        std::exit(1);
    }

    uint64_t rows = 0;
    uint64_t resolved = 0;
    uint64_t duplicates = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        auto fields = split_tab_line(line);
        if (fields.size() < columns.size()) fields.resize(columns.size());
        const std::string& read_id = fields[column_idx["read_id"]];
        if (read_id.empty()) continue;
        ++rows;

        OracleTarget target;
        bool ok = false;
        if (has_row_col) {
            int row2 = -1;
            int col2 = -1;
            if (parse_int_field(fields[column_idx["row2"]], row2) &&
                parse_int_field(fields[column_idx["col2"]], col2) &&
                row2 >= 0 && row2 < static_cast<int>(cfg.bc2_oligos.size()) &&
                col2 >= 0 && col2 < static_cast<int>(cfg.bc1_oligos.size())) {
                target.resolved = true;
                target.row2 = row2;
                target.col2 = col2;
                target.cb = cfg.bc1_oligos[col2] + cfg.bc2_oligos[row2];
                target.status = "resolved_from_row_col";
                ok = true;
            }
        }
        if (!ok && has_cb) {
            resolve_cb_to_oracle_target(fields[column_idx[cb_col]], cfg,
                                        bc1_by_seq, bc2_by_seq, target);
        }
        if (target.resolved) ++resolved;
        auto existing = out.find(read_id);
        if (existing == out.end() || (!existing->second.resolved && target.resolved)) {
            out[read_id] = std::move(target);
        } else {
            ++duplicates;
        }
    }
    std::cerr << "loaded SR oracle TSV: rows=" << rows
              << " read_ids=" << out.size()
              << " resolved_rows=" << resolved
              << " duplicates_ignored=" << duplicates
              << " path=" << path << "\n";
    return out;
}

void write_decode_annotation_header(std::ofstream& out, const Config& cfg) {
    out << "\tscore_mode"
        << "\tdirect_best_resolution"
        << "\tcandidate_pairs"
        << "\tbest_unique"
        << "\tbest_tie_candidates"
        << "\tscore0"
        << "\tscore1"
        << "\tscore2"
        << "\tscore3"
        << "\tbc1_edit"
        << "\tbc2_edit"
        << "\tbc1_affine_cost"
        << "\tbc2_affine_cost"
        << "\tbc1_substitutions"
        << "\tbc1_insertions"
        << "\tbc1_deletions"
        << "\tbc1_n_edits"
        << "\tbc2_substitutions"
        << "\tbc2_insertions"
        << "\tbc2_deletions"
        << "\tbc2_n_edits"
        << "\tanchor_score"
        << "\toffset_sum"
        << "\tbc1_obs_len"
        << "\tbc2_obs_len";
    if (cfg.frozen_support_consensus) {
        out << "\tfrozen_support_attempted"
            << "\tfrozen_support_resolved"
            << "\tfrozen_support_status"
            << "\tfrozen_h0_odds_top_over_second"
            << "\tfrozen_exposure_odds_top_over_second";
    }
}

void append_long_field(std::string& line, long value, bool present) {
    line += '\t';
    if (present) line += std::to_string(value);
}

void append_int_field(std::string& line, int value) {
    line += '\t';
    if (value >= 0) line += std::to_string(value);
}

void append_signed_int_field(std::string& line, int value, bool present) {
    line += '\t';
    if (present) line += std::to_string(value);
}

void append_alignment_fields(std::string& line, const AlignmentMetrics& metrics) {
    append_int_field(line, metrics.cost);
}

void append_alignment_op_field(std::string& line, const AlignmentMetrics& metrics, int value) {
    line += '\t';
    if (metrics.cost >= 0) line += std::to_string(value);
}

void append_bool_value(std::string& line, bool value);

void append_decode_annotations(std::string& line, const Decoded& d,
                               const Config& cfg) {
    line += '\t';
    line += score_mode_name(cfg.score_mode);
    line += '\t';
    line += direct_best_resolution_name(cfg.direct_best_resolution);
    line += '\t';
    line += std::to_string(d.candidate_pairs);
    line += '\t';
    if (d.have_best) line += d.best_unique ? "True" : "False";
    line += '\t';
    if (d.have_best) line += std::to_string(d.best_tie_candidates);
    for (int i = 0; i < 4; ++i) append_long_field(line, d.score_key[i], d.have_best);
    append_int_field(line, d.bc1_edit);
    append_int_field(line, d.bc2_edit);
    append_alignment_fields(line, d.bc1_affine);
    append_alignment_fields(line, d.bc2_affine);
    append_alignment_op_field(line, d.bc1_affine, d.bc1_affine.substitutions);
    append_alignment_op_field(line, d.bc1_affine, d.bc1_affine.insertions);
    append_alignment_op_field(line, d.bc1_affine, d.bc1_affine.deletions);
    append_alignment_op_field(line, d.bc1_affine, d.bc1_affine.n_edits);
    append_alignment_op_field(line, d.bc2_affine, d.bc2_affine.substitutions);
    append_alignment_op_field(line, d.bc2_affine, d.bc2_affine.insertions);
    append_alignment_op_field(line, d.bc2_affine, d.bc2_affine.deletions);
    append_alignment_op_field(line, d.bc2_affine, d.bc2_affine.n_edits);
    append_int_field(line, d.anchor_score);
    append_int_field(line, d.offset_sum);
    append_int_field(line, d.bc1_obs_len);
    append_int_field(line, d.bc2_obs_len);
    if (cfg.frozen_support_consensus) {
        append_bool_value(line, d.frozen_support_attempted);
        append_bool_value(line, d.frozen_support_resolved);
        line += '\t';
        line += d.frozen_support_status;
        line += '\t';
        if (d.frozen_h0_odds > 0.0) line += std::to_string(d.frozen_h0_odds);
        line += '\t';
        if (d.frozen_exposure_odds > 0.0) {
            line += std::to_string(d.frozen_exposure_odds);
        }
    }
}

void append_bool_value(std::string& line, bool value) {
    line += '\t';
    line += value ? "True" : "False";
}

void write_candidate_audit_header(std::ofstream& out) {
    out << "read_id"
        << "\tsr_oracle_status"
        << "\tsr_cb"
        << "\tsr_row2"
        << "\tsr_col2"
        << "\taudit_bucket"
        << "\tdecoder_assigned"
        << "\tdecoder_row2"
        << "\tdecoder_col2"
        << "\tdecoder_full_start"
        << "\tcandidate_pairs"
        << "\tbest_unique"
        << "\tbest_tie_candidates"
        << "\tbest_score0"
        << "\tbest_score1"
        << "\tbest_score2"
        << "\tbest_score3"
        << "\tsr_candidate_present"
        << "\tsr_score0"
        << "\tsr_score1"
        << "\tsr_score2"
        << "\tsr_score3"
        << "\tsr_bc1_seed_seen"
        << "\tsr_bc2_seed_seen"
        << "\tsr_pair_seen"
        << "\tsr_pair_bc1_pass"
        << "\tsr_pair_bc2_pass"
        << "\tsr_pair_anchor_pass"
        << "\tsr_bc1_edit"
        << "\tsr_bc2_edit"
        << "\tsr_anchor_score"
        << "\tsr_offset_sum"
        << "\tsr_full_start"
        << "\tsr_bc1_obs_len"
        << "\tsr_bc2_obs_len"
        << "\tcurrent_target_possible"
        << "\tcurrent_min_total_edit"
        << "\tcurrent_bc1_edit"
        << "\tcurrent_bc2_edit"
        << "\tcurrent_anchor_score"
        << "\tcurrent_offset_sum"
        << "\tcurrent_full_start"
        << "\tcurrent_bc2_start"
        << "\tcurrent_bc1_obs_len"
        << "\tcurrent_bc2_obs_len"
        << "\tcurrent_no_anchor_possible"
        << "\tcurrent_no_anchor_failure_class"
        << "\tcurrent_no_anchor_bc2_anchor_edit"
        << "\tcurrent_no_anchor_current_bc1_anchor_edit"
        << "\tcurrent_no_anchor_current_d1"
        << "\tcurrent_no_anchor_current_d2"
        << "\tcurrent_no_anchor_any_bc1_anchor_edit"
        << "\tcurrent_no_anchor_any_d1"
        << "\tcurrent_no_anchor_any_d2"
        << "\tbroad_target_possible"
        << "\tbroad_min_total_edit"
        << "\tbroad_bc1_edit"
        << "\tbroad_bc2_edit"
        << "\tbroad_anchor_score"
        << "\tbroad_offset_sum"
        << "\tbroad_full_start"
        << "\tbroad_bc2_start"
        << "\tbroad_bc1_obs_len"
        << "\tbroad_bc2_obs_len"
        << "\tbroad_no_anchor_possible"
        << "\tbroad_no_anchor_failure_class"
        << "\tbroad_no_anchor_bc2_anchor_edit"
        << "\tbroad_no_anchor_current_bc1_anchor_edit"
        << "\tbroad_no_anchor_current_d1"
        << "\tbroad_no_anchor_current_d2"
        << "\tbroad_no_anchor_any_bc1_anchor_edit"
        << "\tbroad_no_anchor_any_d1"
        << "\tbroad_no_anchor_any_d2"
        << "\tnon_acgt_sequestered"
        << "\tdirect_best_resolution"
        << "\tbest_candidate_count"
        << "\tbest_candidate_traces";
}

std::string candidate_audit_bucket(const Decoded& d, const OracleTarget& oracle,
                                   const OracleTargetAudit& current,
                                   const OracleTargetAudit& current_no_anchor,
                                   const OracleTargetAudit& broad,
                                   const OracleTargetAudit& broad_no_anchor,
                                   const Config& cfg) {
    if (!oracle.resolved) return "sr_oracle_unresolved";
    if (d.oracle_candidate_present) {
        if (d.assigned && d.row2 == oracle.row2 && d.col2 == oracle.col2) {
            return "sr_candidate_present_selected";
        }
        if (d.have_best && compare_score_keys(d.oracle_score_key, d.score_key) == 0) {
            return "sr_candidate_present_in_best_tie";
        }
        return "sr_candidate_present_ranked_lower";
    }
    if (d.sequestered_non_acgt) return "non_acgt_or_quality_edge";
    if (current.possible && current.observed_non_acgt) return "non_acgt_or_quality_edge";
    if (current.possible) return "sr_candidate_absent_within_current_target_possible";
    if (current_no_anchor.possible && !current.possible) return "sr_candidate_absent_anchor_gated";
    if (d.oracle_pair_seen && !d.oracle_pair_bc1_pass) return "sr_candidate_absent_bc1_full_gated";
    if (d.oracle_pair_bc1_pass && !d.oracle_pair_bc2_pass) return "sr_candidate_absent_bc2_full_gated";
    if (d.oracle_pair_bc2_pass && !d.oracle_pair_anchor_pass) return "sr_candidate_absent_anchor_gated";
    if (broad.possible &&
        broad.bc1_edit <= cfg.bc1_full_hamming &&
        broad.bc2_edit <= cfg.bc2_full_hamming) {
        return "sr_candidate_absent_windowed_out";
    }
    if (!broad.possible && broad_no_anchor.possible) return "sr_candidate_absent_broad_anchor_gated";
    if (broad.possible) return "sr_candidate_absent_needs_edit_3_or_4";
    return "sr_candidate_absent_no_broad_match";
}

void append_audit_envelope(std::string& line, const OracleTargetAudit& audit) {
    append_bool_value(line, audit.possible);
    append_int_field(line, audit.min_total_edit);
    append_int_field(line, audit.bc1_edit);
    append_int_field(line, audit.bc2_edit);
    append_int_field(line, audit.anchor_score);
    append_int_field(line, audit.offset_sum);
    append_int_field(line, audit.full_start);
    append_int_field(line, audit.bc2_start);
    append_int_field(line, audit.bc1_obs_len);
    append_int_field(line, audit.bc2_obs_len);
}

void append_anchor_limit_fields(std::string& line, const OracleTargetAudit& audit) {
    line += '\t';
    line += audit.anchor_failure_class;
    append_int_field(line, audit.anchor_bc2_edit);
    append_int_field(line, audit.anchor_current_bc1_edit);
    append_signed_int_field(line, audit.anchor_current_d1, audit.anchor_current_bc1_edit >= 0);
    append_signed_int_field(line, audit.anchor_current_d2, audit.anchor_current_bc1_edit >= 0);
    append_int_field(line, audit.anchor_any_bc1_edit);
    append_signed_int_field(line, audit.anchor_any_d1, audit.anchor_any_bc1_edit >= 0);
    append_signed_int_field(line, audit.anchor_any_d2, audit.anchor_any_bc1_edit >= 0);
}

void append_best_candidate_traces(std::string& line, const Decoded& d) {
    line += '\t';
    line += std::to_string(d.best_candidate_traces.size());
    line += '\t';
    bool first = true;
    for (const CandidateTrace& trace : d.best_candidate_traces) {
        if (!first) line += ';';
        first = false;
        line += std::to_string(trace.row2);
        line += ':';
        line += std::to_string(trace.col2);
        line += ':';
        line += std::to_string(trace.full_start);
        line += ':';
        line += std::to_string(trace.bc1_edit);
        line += ':';
        line += std::to_string(trace.bc2_edit);
        line += ':';
        line += std::to_string(trace.anchor_score);
        line += ':';
        line += std::to_string(trace.offset_sum);
        line += ':';
        line += std::to_string(trace.bc1_obs_len);
        line += ':';
        line += std::to_string(trace.bc2_obs_len);
    }
}

void append_candidate_trace_list(std::string& line,
                                 const std::vector<CandidateTrace>& traces) {
    bool first = true;
    for (const CandidateTrace& trace : traces) {
        if (!first) line += ';';
        first = false;
        line += std::to_string(trace.row2);
        line += ':';
        line += std::to_string(trace.col2);
        line += ':';
        line += std::to_string(trace.full_start);
        line += ':';
        line += std::to_string(trace.bc1_edit);
        line += ':';
        line += std::to_string(trace.bc2_edit);
        line += ':';
        line += std::to_string(trace.anchor_score);
        line += ':';
        line += std::to_string(trace.offset_sum);
        line += ':';
        line += std::to_string(trace.bc1_obs_len);
        line += ':';
        line += std::to_string(trace.bc2_obs_len);
    }
}

std::string make_latent_candidate_line(const std::string& read_id,
                                       const Decoded& decoded) {
    std::string line;
    line.reserve(read_id.size() + decoded.latent_candidate_traces.size() * 36 + 32);
    line += read_id;
    line += '\t';
    if (decoded.have_best) line += std::to_string(decoded.score_key[0]);
    line += '\t';
    line += std::to_string(decoded.latent_candidate_traces.size());
    line += '\t';
    append_candidate_trace_list(line, decoded.latent_candidate_traces);
    return line;
}

std::string make_candidate_audit_line(const std::string& read_id,
                                      const OracleTarget& oracle,
                                      const Decoded& d,
                                      const OracleTargetAudit& current,
                                      const OracleTargetAudit& current_no_anchor,
                                      const OracleTargetAudit& broad,
                                      const OracleTargetAudit& broad_no_anchor,
                                      const Config& cfg) {
    std::string line;
    line.reserve(read_id.size() + oracle.cb.size() + 512);
    const std::string bucket =
        candidate_audit_bucket(d, oracle, current, current_no_anchor, broad,
                               broad_no_anchor, cfg);
    line += read_id;
    line += '\t';
    line += oracle.status;
    line += '\t';
    line += oracle.cb;
    append_int_field(line, oracle.row2);
    append_int_field(line, oracle.col2);
    line += '\t';
    line += bucket;
    append_bool_value(line, d.assigned);
    append_int_field(line, d.assigned ? d.row2 : -1);
    append_int_field(line, d.assigned ? d.col2 : -1);
    append_int_field(line, d.assigned ? d.full_start : -1);
    line += '\t';
    line += std::to_string(d.candidate_pairs);
    append_bool_value(line, d.have_best && d.best_unique);
    append_int_field(line, d.have_best ? static_cast<int>(d.best_tie_candidates) : -1);
    for (int i = 0; i < 4; ++i) append_long_field(line, d.score_key[i], d.have_best);
    append_bool_value(line, d.oracle_candidate_present);
    for (int i = 0; i < 4; ++i) {
        append_long_field(line, d.oracle_score_key[i], d.oracle_candidate_present);
    }
    append_bool_value(line, d.oracle_bc1_seed_seen);
    append_bool_value(line, d.oracle_bc2_seed_seen);
    append_bool_value(line, d.oracle_pair_seen);
    append_bool_value(line, d.oracle_pair_bc1_pass);
    append_bool_value(line, d.oracle_pair_bc2_pass);
    append_bool_value(line, d.oracle_pair_anchor_pass);
    append_int_field(line, d.oracle_bc1_edit);
    append_int_field(line, d.oracle_bc2_edit);
    append_int_field(line, d.oracle_anchor_score);
    append_int_field(line, d.oracle_offset_sum);
    append_int_field(line, d.oracle_full_start);
    append_int_field(line, d.oracle_bc1_obs_len);
    append_int_field(line, d.oracle_bc2_obs_len);
    append_audit_envelope(line, current);
    append_bool_value(line, current_no_anchor.possible);
    append_anchor_limit_fields(line, current_no_anchor);
    append_audit_envelope(line, broad);
    append_bool_value(line, broad_no_anchor.possible);
    append_anchor_limit_fields(line, broad_no_anchor);
    append_bool_value(line, d.sequestered_non_acgt);
    line += '\t';
    line += direct_best_resolution_name(cfg.direct_best_resolution);
    append_best_candidate_traces(line, d);
    return line;
}

double raw_candidate_error_probability(int phred) {
    const double probability =
        std::pow(10.0, -static_cast<double>(std::max(0, phred)) / 10.0);
    return std::max(1.0e-10, std::min(0.75, probability));
}

double raw_candidate_pair_log_likelihood(char observed, char candidate, int phred) {
    const char obs = static_cast<char>(std::toupper(static_cast<unsigned char>(observed)));
    const char truth = static_cast<char>(std::toupper(static_cast<unsigned char>(candidate)));
    if ((obs != 'A' && obs != 'C' && obs != 'G' && obs != 'T') ||
        (truth != 'A' && truth != 'C' && truth != 'G' && truth != 'T')) {
        return std::log(0.25);
    }
    const double error = raw_candidate_error_probability(phred);
    return obs == truth ? std::log1p(-error) : std::log(error / 3.0);
}

double raw_candidate_alignment_log_likelihood(
    const std::string& observed, const std::string& qualities,
    const std::string& candidate, int gap_q = 30, int missing_q = 30) {
    std::vector<double> previous(candidate.size() + 1, 0.0);
    std::vector<double> current(candidate.size() + 1, 0.0);
    const double deletion = std::log(raw_candidate_error_probability(gap_q));
    for (size_t j = 1; j <= candidate.size(); ++j) previous[j] = previous[j - 1] + deletion;
    for (size_t i = 1; i <= observed.size(); ++i) {
        const int q = i - 1 < qualities.size()
            ? std::max(0, static_cast<int>(static_cast<unsigned char>(qualities[i - 1])) - 33)
            : missing_q;
        const double insertion = std::log(raw_candidate_error_probability(q) / 4.0);
        current[0] = previous[0] + insertion;
        for (size_t j = 1; j <= candidate.size(); ++j) {
            current[j] = std::max({
                previous[j - 1] + raw_candidate_pair_log_likelihood(
                    observed[i - 1], candidate[j - 1], q),
                previous[j] + insertion,
                current[j - 1] + deletion,
            });
        }
        previous.swap(current);
    }
    return previous[candidate.size()];
}

std::string make_candidate_preserving_lines(
    const Record& record, const Decoded& decoded, const Config& cfg) {
    struct Choice {
        CandidateTrace trace;
        double log_likelihood;

        Choice(const CandidateTrace& candidateTrace = CandidateTrace(),
               double likelihood = -std::numeric_limits<double>::infinity())
            : trace(candidateTrace), log_likelihood(likelihood) {}
    };
    std::map<std::pair<int, int>, Choice> choices;
    for (const CandidateTrace& trace : decoded.best_candidate_traces) {
        const int observed_length = trace.bc1_obs_len + trace.bc2_obs_len;
        if (trace.row2 < 0 || trace.col2 < 0 ||
            trace.row2 >= static_cast<int>(cfg.bc2_oligos.size()) ||
            trace.col2 >= static_cast<int>(cfg.bc1_oligos.size()) ||
            trace.full_start < 0 || observed_length <= 0 ||
            trace.full_start + observed_length > static_cast<int>(record.seq.size())) {
            continue;
        }
        const std::string observed = record.seq.substr(
            static_cast<size_t>(trace.full_start), static_cast<size_t>(observed_length));
        const std::string qualities = trace.full_start < static_cast<int>(record.qual.size())
            ? record.qual.substr(static_cast<size_t>(trace.full_start),
                                 static_cast<size_t>(observed_length))
            : std::string();
        const std::string candidate = cfg.bc1_oligos[static_cast<size_t>(trace.col2)] +
                                      cfg.bc2_oligos[static_cast<size_t>(trace.row2)];
        const double likelihood = raw_candidate_alignment_log_likelihood(
            observed, qualities, candidate);
        const std::pair<int, int> coordinate(trace.row2, trace.col2);
        auto found = choices.find(coordinate);
        if (found == choices.end() || likelihood > found->second.log_likelihood ||
            (likelihood == found->second.log_likelihood &&
             std::tie(trace.full_start, trace.bc1_obs_len, trace.bc2_obs_len) <
             std::tie(found->second.trace.full_start, found->second.trace.bc1_obs_len,
                      found->second.trace.bc2_obs_len))) {
            choices[coordinate] = {trace, likelihood};
        }
    }
    std::ostringstream output;
    const std::string raw_umi = record.seq.substr(0, std::min<size_t>(9, record.seq.size()));
    const int candidate_count = static_cast<int>(choices.size());
    for (const auto& item : choices) {
        const Choice& choice = item.second;
        const CandidateTrace& trace = choice.trace;
        output << record.global_ordinal << '\t' << record.read_id << '\t' << raw_umi << '\t'
               << candidate_count << '\t'
               << trace.row2 << '\t' << trace.col2 << '\t'
               << (trace.bc1_edit + trace.bc2_edit) << '\t'
               << trace.bc1_edit << '\t' << trace.bc2_edit << '\t'
               << trace.bc1_obs_len << '\t' << trace.bc2_obs_len << '\t'
               << trace.full_start << '\t' << std::setprecision(17)
               << choice.log_likelihood << '\n';
    }
    return output.str();
}

using InternalDecoderConfig = Config;
using InternalDecoded = Decoded;
using InternalCandidateTrace = CandidateTrace;

}  // namespace

#ifdef STAR_SPATIAL_R1_DECODER_STANDALONE
int main(int argc, char** argv) {
    Config cfg;
    std::vector<std::string> r1_files;
    std::string bc1_path, bc2_path, out_path, mode = "hamming-deletion-hash";
    std::string cache_dir;
    std::string oligo_mutation_stats_out;
    std::string frozen_h0_prior_path;
    std::string frozen_oligo_stats_path;
    int sig_del = 1;
    bool emit_stats = false;
    bool suppress_read_output = false;
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads < 1) threads = 4;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--r1-fastq") r1_files.push_back(next());
        else if (a == "--bc1-oligos") bc1_path = next();
        else if (a == "--bc2-oligos") bc2_path = next();
        else if (a == "--out") out_path = next();
        else if (a == "--cache-dir") cache_dir = next();
        else if (a == "--stats") emit_stats = true;
        else if (a == "--oligo-mutation-stats-out") oligo_mutation_stats_out = next();
        else if (a == "--suppress-read-output") suppress_read_output = true;
        else if (a == "--observed-cache-experiment") cfg.observed_cache_experiment = true;
        else if (a == "--observed-cache-prebuild") {
            cfg.observed_cache_experiment = true;
            cfg.observed_cache_prebuild = true;
        }
        else if (a == "--observed-cache-prebuild-fallback") {
            cfg.observed_cache_prebuild_fallback = true;
        }
        else if (a == "--observed-cache-prebuild-sample-reads" ||
                 a == "--observed-cache-prebuild-read-limit") {
            cfg.observed_cache_prebuild_sample_reads = std::stol(next());
        }
        else if (a == "--decode-skip-reads") cfg.decode_skip_reads = std::stol(next());
        else if (a == "--decode-read-limit") cfg.decode_read_limit = std::stol(next());
        else if (a == "--observed-cache-no-e01") cfg.observed_cache_use_e01 = false;
        else if (a == "--direct-full-decode") cfg.direct_full_decode = true;
        else if (a == "--direct-full-exact-lengths") cfg.direct_full_exact_lengths = true;
        else if (a == "--direct-best-edit-decode") cfg.direct_best_edit_decode = true;
        else if (a == "--direct-tiered-decode") cfg.direct_tiered_decode = true;
        else if (a == "--direct-tiered-fixed-decode") cfg.direct_tiered_fixed_decode = true;
        else if (a == "--direct-tiered-h2-decode") cfg.direct_tiered_h2_decode = true;
        else if (a == "--tie-resolver") {
            const std::string value = next();
            if (value != "frozen-support-consensus") {
                std::cerr << "ERROR: --tie-resolver must be frozen-support-consensus\n";
                return 2;
            }
            cfg.frozen_support_consensus = true;
        }
        else if (a == "--frozen-h0-prior-tsv") frozen_h0_prior_path = next();
        else if (a == "--frozen-oligo-stats-tsv") frozen_oligo_stats_path = next();
        else if (a == "--frozen-support-min-odds") {
            cfg.frozen_support_min_odds = std::stod(next());
        }
        else if (a == "--direct-best-n-dp-fallback") cfg.direct_best_non_acgt_dp_fallback = true;
        else if (a == "--direct-best-resolution") {
            std::string value = next();
            if (value == "strict" || value == "strict-boundary") {
                cfg.direct_best_resolution = DirectBestResolutionMode::StrictBoundary;
            } else if (value == "overlap-span-dp" || value == "overlap") {
                cfg.direct_best_resolution = DirectBestResolutionMode::OverlapSpanDp;
            } else {
                std::cerr << "ERROR: --direct-best-resolution must be strict-boundary or overlap-span-dp\n";
                return 2;
            }
        }
        else if (a == "--overlap-dp-no-length-guard") {
            cfg.overlap_dp_length_guard = false;
        }
        else if (a == "--overlap-dp-bc1-min-obs-len") {
            cfg.overlap_dp_bc1_min_obs_len = std::stoi(next());
        }
        else if (a == "--overlap-dp-bc2-min-obs-len") {
            cfg.overlap_dp_bc2_min_obs_len = std::stoi(next());
        }
        else if (a == "--overlap-dp-bc1-max-len-delta") {
            cfg.overlap_dp_bc1_max_len_delta = std::stoi(next());
        }
        else if (a == "--overlap-dp-bc2-max-len-delta") {
            cfg.overlap_dp_bc2_max_len_delta = std::stoi(next());
        }
        else if (a == "--overlap-dp-anchor-gate-only") {
            cfg.overlap_dp_anchor_gate_only = true;
        }
        else if (a == "--overlap-dp-rank-anchor") {
            cfg.overlap_dp_anchor_gate_only = false;
        }
        else if (a == "--overlap-dp-target-scan") {
            cfg.overlap_dp_use_edit_hash = false;
        }
        else if (a == "--overlap-dp-edit-hash") {
            cfg.overlap_dp_use_edit_hash = true;
        }
        else if (a == "--overlap-dp-no-terminal-offset") {
            cfg.overlap_dp_offset_include_terminal_delta = false;
        }
        else if (a == "--overlap-dp-offset-includes-terminal-delta") {
            cfg.overlap_dp_offset_include_terminal_delta = true;
        }
        else if (a == "--score-mode") {
            std::string value = next();
            if (value == "legacy" || value == "edit") {
                cfg.score_mode = CandidateScoreMode::LegacyEdit;
            } else if (value == "affine") {
                cfg.score_mode = CandidateScoreMode::AffineGap;
            } else {
                std::cerr << "ERROR: --score-mode must be legacy or affine\n";
                return 2;
            }
        }
        else if (a == "--affine-sub-cost") cfg.affine.substitution_cost = std::stoi(next());
        else if (a == "--affine-gap-open") cfg.affine.gap_open_cost = std::stoi(next());
        else if (a == "--affine-gap-extend") cfg.affine.gap_extend_cost = std::stoi(next());
        else if (a == "--affine-n-cost") cfg.affine.n_cost = std::stoi(next());
        else if (a == "--emit-decode-annotations") cfg.emit_decode_annotations = true;
        else if (a == "--sr-oracle-tsv") cfg.sr_oracle_tsv = next();
        else if (a == "--candidate-audit-out") cfg.candidate_audit_out = next();
        else if (a == "--candidate-preserving-out") cfg.candidate_preserving_out = next();
        else if (a == "--latent-candidate-out") cfg.latent_candidate_out = next();
        else if (a == "--suppress-bc1-overlong-prefix-e1") {
            cfg.suppress_bc1_overlong_prefix_e1 = true;
        }
        else if (a == "--allow-bc2-short-terminal-e2") {
            cfg.suppress_bc2_short_terminal_e2 = false;
        }
        else if (a == "--audit-full-start-padding") cfg.audit_full_start_padding = std::stoi(next());
        else if (a == "--audit-bc2-start-padding") cfg.audit_bc2_start_padding = std::stoi(next());
        else if (a == "--audit-edit-max") cfg.audit_edit_max = std::stoi(next());
        else if (a == "--anchor-match-mode") mode = next();
        else if (a == "--anchor-signature-deletions") sig_del = std::stoi(next());
        else if (a == "--full-start-min") cfg.full_start_min = std::stoi(next());
        else if (a == "--full-start-max") cfg.full_start_max = std::stoi(next());
        else if (a == "--bc1-anchor-hamming") cfg.bc1_anchor_hamming = std::stoi(next());
        else if (a == "--bc2-anchor-hamming") cfg.bc2_anchor_hamming = std::stoi(next());
        else if (a == "--bc1-full-hamming") cfg.bc1_full_hamming = std::stoi(next());
        else if (a == "--bc2-full-hamming") cfg.bc2_full_hamming = std::stoi(next());
        else if (a == "--bc1-anchor-offset-window") cfg.bc1_anchor_offset_window = std::stoi(next());
        else if (a == "--bc2-offset-window") cfg.bc2_offset_window = std::stoi(next());
        else if (a == "--grid-rows") cfg.grid_rows = std::stoi(next());
        else if (a == "--grid-cols") cfg.grid_cols = std::stoi(next());
        else if (a == "--threads") threads = std::stoi(next());
        else if (a == "--non-acgt-reads-out") cfg.non_acgt_reads_out = next();
        else if (a == "--progress-log-every-batches") cfg.progress_log_every_batches = std::stol(next());
        else { std::cerr << "unknown arg: " << a << "\n"; return 2; }
    }
    if (r1_files.empty() || bc1_path.empty() || bc2_path.empty() ||
        (!suppress_read_output && out_path.empty()) ||
        (suppress_read_output && oligo_mutation_stats_out.empty())) {
        std::cerr << "usage: hd_r1_anchored_decode --r1-fastq F [--r1-fastq F...] --bc1-oligos F "
                     "--bc2-oligos F --out decode_reads.tsv [--cache-dir DIR] "
                     "[--oligo-mutation-stats-out TSV] [--suppress-read-output] "
                     "[--direct-tiered-decode] [--direct-tiered-fixed-decode] "
                     "[--direct-tiered-h2-decode] "
                     "[--tie-resolver frozen-support-consensus "
                     "--frozen-h0-prior-tsv TSV --frozen-oligo-stats-tsv TSV "
                     "--frozen-support-min-odds X] "
                     "[--direct-best-n-dp-fallback] [--non-acgt-reads-out TSV] "
                     "[--progress-log-every-batches N] "
                     "[--direct-best-resolution strict-boundary|overlap-span-dp] "
                     "[--overlap-dp-no-length-guard] "
                     "[--overlap-dp-bc1-min-obs-len N] "
                     "[--overlap-dp-bc2-min-obs-len N] "
                     "[--overlap-dp-bc1-max-len-delta N] "
                     "[--overlap-dp-bc2-max-len-delta N] "
                     "[--overlap-dp-anchor-gate-only] "
                     "[--overlap-dp-target-scan] [--overlap-dp-edit-hash] "
                     "[--overlap-dp-no-terminal-offset] "
                     "[--overlap-dp-offset-includes-terminal-delta] "
                     "[--score-mode legacy|affine --affine-sub-cost N "
                     "--affine-gap-open N --affine-gap-extend N --affine-n-cost N] "
                     "[--emit-decode-annotations] "
                     "[--suppress-bc1-overlong-prefix-e1] "
                     "[--allow-bc2-short-terminal-e2] "
                     "[--sr-oracle-tsv TSV --candidate-audit-out TSV] "
                     "[--candidate-preserving-out TSV] "
                     "[--latent-candidate-out TSV] [options]\n";
        return 2;
    }
    if (cfg.affine.substitution_cost < 0 || cfg.affine.gap_open_cost < 0 ||
        cfg.affine.gap_extend_cost < 0 || cfg.affine.n_cost < 0) {
        std::cerr << "ERROR: affine costs must be non-negative\n";
        return 2;
    }
    if (cfg.audit_full_start_padding < 0 || cfg.audit_bc2_start_padding < 0 ||
        cfg.audit_edit_max < 0) {
        std::cerr << "ERROR: audit padding/edit values must be non-negative\n";
        return 2;
    }
    if (!std::isfinite(cfg.frozen_support_min_odds) ||
        cfg.frozen_support_min_odds < 1.0) {
        std::cerr << "ERROR: --frozen-support-min-odds must be finite and >= 1\n";
        return 2;
    }
    if (cfg.frozen_support_consensus) {
        if (!cfg.direct_tiered_h2_decode) {
            std::cerr << "ERROR: frozen-support-consensus requires --direct-tiered-h2-decode\n";
            return 2;
        }
        if (frozen_h0_prior_path.empty() || frozen_oligo_stats_path.empty()) {
            std::cerr << "ERROR: frozen-support-consensus requires --frozen-h0-prior-tsv "
                         "and --frozen-oligo-stats-tsv\n";
            return 2;
        }
    } else if (!frozen_h0_prior_path.empty() || !frozen_oligo_stats_path.empty()) {
        std::cerr << "ERROR: frozen support TSVs require --tie-resolver "
                     "frozen-support-consensus\n";
        return 2;
    }
    if (cfg.overlap_dp_bc1_min_obs_len < 0 || cfg.overlap_dp_bc2_min_obs_len < 0 ||
        cfg.overlap_dp_bc1_max_len_delta < 0 || cfg.overlap_dp_bc2_max_len_delta < 0) {
        std::cerr << "ERROR: overlap DP length guard values must be non-negative\n";
        return 2;
    }
    if (cfg.sr_oracle_tsv.empty() != cfg.candidate_audit_out.empty()) {
        std::cerr << "ERROR: --sr-oracle-tsv and --candidate-audit-out must be supplied together\n";
        return 2;
    }
    if (cfg.direct_full_decode && !cfg.observed_cache_prebuild) {
        std::cerr << "ERROR: --direct-full-decode requires --observed-cache-prebuild\n";
        return 2;
    }
    const int standalone_direct_modes =
        static_cast<int>(cfg.direct_best_edit_decode) +
        static_cast<int>(cfg.direct_tiered_decode) +
        static_cast<int>(cfg.direct_tiered_fixed_decode) +
        static_cast<int>(cfg.direct_tiered_h2_decode);
    if (standalone_direct_modes > 1) {
        std::cerr << "ERROR: --direct-best-edit-decode, --direct-tiered-decode, and "
                     "--direct-tiered-fixed-decode/--direct-tiered-h2-decode are alternative "
                     "standalone decode paths; use only one\n";
        return 2;
    }
    if (standalone_direct_modes > 0 &&
        (cfg.observed_cache_experiment || cfg.observed_cache_prebuild || cfg.direct_full_decode)) {
        std::cerr << "ERROR: standalone direct decode paths must not be combined with "
                     "observed-cache or direct-full-decode modes\n";
        return 2;
    }
    if ((cfg.score_mode == CandidateScoreMode::AffineGap || cfg.emit_decode_annotations) &&
        cfg.observed_cache_prebuild) {
        std::cerr << "ERROR: affine scoring and decode annotations are implemented for the "
                     "standard anchored and direct-best decode paths, not observed-cache-prebuild\n";
        return 2;
    }
    if (cfg.score_mode == CandidateScoreMode::AffineGap &&
        (cfg.direct_tiered_decode || cfg.direct_tiered_fixed_decode ||
         cfg.direct_tiered_h2_decode)) {
        std::cerr << "ERROR: affine scoring is not implemented for direct-tiered decode paths\n";
        return 2;
    }
    if (!cfg.candidate_audit_out.empty() && cfg.observed_cache_prebuild) {
        std::cerr << "ERROR: candidate audit is implemented for the standard anchored and "
                     "direct-best decode paths, not observed-cache-prebuild\n";
        return 2;
    }
    if (!cfg.latent_candidate_out.empty() && !cfg.direct_tiered_h2_decode) {
        std::cerr << "ERROR: --latent-candidate-out requires --direct-tiered-h2-decode\n";
        return 2;
    }
    if (!ensure_directory(cache_dir)) return 1;

    load_lines({bc1_path}, cfg.bc1_oligos);
    load_lines({bc2_path}, cfg.bc2_oligos);
    cfg.bc1_len.resize(cfg.bc1_oligos.size());
    cfg.bc2_len.resize(cfg.bc2_oligos.size());
    for (size_t i = 0; i < cfg.bc1_oligos.size(); ++i) cfg.bc1_len[i] = static_cast<int>(cfg.bc1_oligos[i].size());
    for (size_t i = 0; i < cfg.bc2_oligos.size(); ++i) cfg.bc2_len[i] = static_cast<int>(cfg.bc2_oligos[i].size());

    if (cfg.frozen_support_consensus) {
        if (!load_frozen_support_counts(
                frozen_h0_prior_path, "exact_h0_read_count",
                cfg.bc1_oligos, cfg.bc2_oligos,
                cfg.bc1_frozen_h0, cfg.bc2_frozen_h0) ||
            !load_frozen_support_counts(
                frozen_oligo_stats_path, "exposure_h0_h1",
                cfg.bc1_oligos, cfg.bc2_oligos,
                cfg.bc1_frozen_exposure, cfg.bc2_frozen_exposure)) {
            return 2;
        }
    }

    // anchors: bc1 suffix(anchor_len), bc2 prefix(bc2_seg_len)
    cfg.bc1_anchors.resize(cfg.bc1_oligos.size());
    cfg.bc2_anchors.resize(cfg.bc2_oligos.size());
    for (size_t i = 0; i < cfg.bc1_oligos.size(); ++i)
        cfg.bc1_anchors[i] = cfg.bc1_oligos[i].substr(cfg.bc1_oligos[i].size() - cfg.anchor_len);
    for (size_t i = 0; i < cfg.bc2_oligos.size(); ++i)
        cfg.bc2_anchors[i] = cfg.bc2_oligos[i].substr(0, cfg.bc2_seg_len);

    std::unordered_map<std::string, OracleTarget> sr_oracle =
        load_sr_oracle_tsv(cfg.sr_oracle_tsv, cfg);

    std::vector<int> bc1_tiered_exact_lengths;
    std::vector<int> bc2_tiered_exact_lengths;
    std::vector<int> bc1_tiered_e1_query_lengths;
    std::vector<int> bc2_tiered_e1_query_lengths;
    std::vector<int> bc1_tiered_h2_query_lengths;
    std::vector<int> bc2_tiered_h2_query_lengths;
    int bc1_tiered_fixed_width = 0;
    int bc2_tiered_fixed_width = 0;

    if (!cfg.direct_best_edit_decode && !cfg.direct_tiered_decode &&
        !cfg.direct_tiered_fixed_decode && !cfg.direct_tiered_h2_decode) {
        cfg.bc1_lookup.build(cfg.bc1_anchors, cfg.bc1_anchor_hamming, mode, sig_del,
                             "bc1_anchor", cache_path_for(cache_dir, "bc1_anchor"));
        cfg.bc2_lookup.build(cfg.bc2_anchors, cfg.bc2_anchor_hamming, mode, sig_del,
                             "bc2_anchor", cache_path_for(cache_dir, "bc2_anchor"));
        cfg.bc1_full_lookup.build_hamming_universe(cfg.bc1_oligos, cfg.bc1_full_hamming,
                                                   PackedAnswerCache::ScoreMode::Edit, "bc1_full",
                                                   cache_path_for(cache_dir, "bc1_full"));
        cfg.bc2_full_lookup.build_hamming_universe(cfg.bc2_oligos, cfg.bc2_full_hamming,
                                                   PackedAnswerCache::ScoreMode::Edit, "bc2_full",
                                                   cache_path_for(cache_dir, "bc2_full"));
    }
    if (cfg.direct_best_edit_decode &&
        (cfg.direct_best_resolution == DirectBestResolutionMode::StrictBoundary ||
         cfg.overlap_dp_use_edit_hash)) {
        const std::vector<int> bc1_query_lengths(
            std::begin(kDirectBestBc1QueryLengths), std::end(kDirectBestBc1QueryLengths));
        const std::vector<int> bc2_query_lengths(
            std::begin(kDirectBestBc2QueryLengths), std::end(kDirectBestBc2QueryLengths));
        cfg.bc1_direct_best_edit_lookup.build_best_edit_universe_for_query_lengths(
            cfg.bc1_oligos, bc1_query_lengths, cfg.bc1_full_hamming,
            "bc1_direct_best_edit", TargetSliceMode::Full, false,
            cfg.direct_best_resolution == DirectBestResolutionMode::StrictBoundary);
        cfg.bc2_direct_best_edit_lookup.build_best_edit_universe_for_query_lengths(
            cfg.bc2_oligos, bc2_query_lengths, cfg.bc2_full_hamming,
            "bc2_direct_best_edit", TargetSliceMode::Full,
            cfg.direct_best_resolution == DirectBestResolutionMode::StrictBoundary,
            cfg.direct_best_resolution == DirectBestResolutionMode::StrictBoundary);
    }
    if (cfg.direct_tiered_decode) {
        bc1_tiered_exact_lengths = unique_sorted_lengths(cfg.bc1_len);
        bc2_tiered_exact_lengths = unique_sorted_lengths(cfg.bc2_len);
        bc1_tiered_e1_query_lengths = one_edit_query_lengths(bc1_tiered_exact_lengths);
        bc2_tiered_e1_query_lengths = one_edit_query_lengths(bc2_tiered_exact_lengths);
        cfg.bc1_tiered_exact_lookup.build_best_edit_universe_for_query_lengths(
            cfg.bc1_oligos, bc1_tiered_exact_lengths, 0,
            "bc1_tiered_exact", TargetSliceMode::Full, false, false);
        cfg.bc2_tiered_exact_lookup.build_best_edit_universe_for_query_lengths(
            cfg.bc2_oligos, bc2_tiered_exact_lengths, 0,
            "bc2_tiered_exact", TargetSliceMode::Full, false, false);
        cfg.bc1_tiered_e1_lookup.build_best_edit_universe_for_query_lengths(
            cfg.bc1_oligos, bc1_tiered_e1_query_lengths, 1,
            "bc1_tiered_e1", TargetSliceMode::Full, false, false);
        cfg.bc2_tiered_e1_lookup.build_best_edit_universe_for_query_lengths(
            cfg.bc2_oligos, bc2_tiered_e1_query_lengths, 1,
            "bc2_tiered_e1", TargetSliceMode::Full, false, false);
    }
    if (cfg.direct_tiered_fixed_decode) {
        bc1_tiered_exact_lengths = unique_sorted_lengths(cfg.bc1_len);
        bc2_tiered_exact_lengths = unique_sorted_lengths(cfg.bc2_len);
        if (bc1_tiered_exact_lengths.empty() || bc2_tiered_exact_lengths.empty()) {
            std::cerr << "ERROR: direct-tiered-fixed-decode requires non-empty BC oligo lists\n";
            return 2;
        }
        bc1_tiered_fixed_width = bc1_tiered_exact_lengths.back() + 1;
        bc2_tiered_fixed_width = bc2_tiered_exact_lengths.back() + 1;
        cfg.bc1_tiered_fixed_lookup.build_suffix_wildcard_edit_universe_for_fixed_width(
            cfg.bc1_oligos, bc1_tiered_fixed_width, 1, "bc1_tiered_fixed_e01");
        cfg.bc2_tiered_fixed_lookup.build_suffix_wildcard_edit_universe_for_fixed_width(
            cfg.bc2_oligos, bc2_tiered_fixed_width, 1, "bc2_tiered_fixed_e01");
    }
    if (cfg.direct_tiered_h2_decode) {
        bc1_tiered_exact_lengths = unique_sorted_lengths(cfg.bc1_len);
        bc2_tiered_exact_lengths = unique_sorted_lengths(cfg.bc2_len);
        bc1_tiered_h2_query_lengths = edit_query_lengths(bc1_tiered_exact_lengths, 2);
        bc2_tiered_h2_query_lengths = edit_query_lengths(bc2_tiered_exact_lengths, 2);
        cfg.bc1_tiered_h2_lookup.build_best_edit_universe_for_query_lengths(
            cfg.bc1_oligos, bc1_tiered_h2_query_lengths, 2,
            "bc1_tiered_h2", TargetSliceMode::Full, false, false);
        cfg.bc2_tiered_h2_lookup.build_best_edit_universe_for_query_lengths(
            cfg.bc2_oligos, bc2_tiered_h2_query_lengths, 2,
            "bc2_tiered_h2", TargetSliceMode::Full, false, false);
    }
    if (cfg.observed_cache_experiment) {
        if (cfg.observed_cache_use_e01) {
            cfg.bc1_anchor_e01.build_edit1_universe(cfg.bc1_anchors, "bc1_anchor");
            cfg.bc2_anchor_e01.build_edit1_universe(cfg.bc2_anchors, "bc2_anchor");
            cfg.bc1_full_e01.build_edit1_universe(cfg.bc1_oligos, "bc1_full");
            cfg.bc2_full_e01.build_edit1_universe(cfg.bc2_oligos, "bc2_full");
        }
        cfg.bc1_anchor_deletion2.build(cfg.bc1_anchors, cfg.bc1_anchor_hamming, "bc1_anchor");
        cfg.bc2_anchor_deletion2.build(cfg.bc2_anchors, cfg.bc2_anchor_hamming, "bc2_anchor");
        cfg.bc1_full_deletion2.build(cfg.bc1_oligos, 2, "bc1_full");
        cfg.bc2_full_deletion2.build(cfg.bc2_oligos, 2, "bc2_full");
        cfg.bc1_anchor_observed.init(cfg.bc1_anchors, cfg.bc1_anchor_e01, cfg.bc1_lookup.hamming_cache,
                                     cfg.bc1_anchor_deletion2, cfg.bc1_anchor_hamming,
                                     "bc1_anchor", cfg.observed_cache_use_e01);
        cfg.bc2_anchor_observed.init(cfg.bc2_anchors, cfg.bc2_anchor_e01, cfg.bc2_lookup.hamming_cache,
                                     cfg.bc2_anchor_deletion2, cfg.bc2_anchor_hamming,
                                     "bc2_anchor", cfg.observed_cache_use_e01);
        cfg.bc1_full_observed.init(cfg.bc1_oligos, cfg.bc1_full_e01, cfg.bc1_full_lookup,
                                   cfg.bc1_full_deletion2, cfg.bc1_full_hamming,
                                   "bc1_full", cfg.observed_cache_use_e01);
        cfg.bc2_full_observed.init(cfg.bc2_oligos, cfg.bc2_full_e01, cfg.bc2_full_lookup,
                                   cfg.bc2_full_deletion2, cfg.bc2_full_hamming,
                                   "bc2_full", cfg.observed_cache_use_e01);
    }
    std::cerr << "built decode indices (bc1=" << cfg.bc1_oligos.size()
              << " bc2=" << cfg.bc2_oligos.size() << ", mode="
              << (cfg.direct_tiered_h2_decode ? "direct-tiered-h2" :
                  (cfg.direct_tiered_fixed_decode ? "direct-tiered-fixed" :
                  (cfg.direct_tiered_decode ? "direct-tiered" :
                  (cfg.direct_best_edit_decode ? "direct-best-edit" : mode))
                  ))
              << ", score_mode=" << score_mode_name(cfg.score_mode)
              << ", direct_best_resolution="
              << direct_best_resolution_name(cfg.direct_best_resolution);
    if (cfg.direct_best_resolution == DirectBestResolutionMode::OverlapSpanDp) {
        std::cerr << ", overlap_dp_length_guard="
                  << (cfg.overlap_dp_length_guard ? "true" : "false")
                  << ", overlap_dp_bc1_min_obs_len="
                  << cfg.overlap_dp_bc1_min_obs_len
                  << ", overlap_dp_bc2_min_obs_len="
                  << cfg.overlap_dp_bc2_min_obs_len
                  << ", overlap_dp_bc1_max_len_delta="
                  << cfg.overlap_dp_bc1_max_len_delta
                  << ", overlap_dp_bc2_max_len_delta="
                  << cfg.overlap_dp_bc2_max_len_delta
                  << ", overlap_dp_anchor_gate_only="
                  << (cfg.overlap_dp_anchor_gate_only ? "true" : "false")
                  << ", overlap_dp_offset_include_terminal_delta="
                  << (cfg.overlap_dp_offset_include_terminal_delta ? "true" : "false")
                  << ", overlap_dp_use_edit_hash="
                  << (cfg.overlap_dp_use_edit_hash ? "true" : "false");
    }
    if (cfg.score_mode == CandidateScoreMode::AffineGap) {
        std::cerr << ", affine_sub=" << cfg.affine.substitution_cost
                  << ", affine_gap_open=" << cfg.affine.gap_open_cost
                  << ", affine_gap_extend=" << cfg.affine.gap_extend_cost
                  << ", affine_n=" << cfg.affine.n_cost;
    }
    if (cfg.frozen_support_consensus) {
        std::cerr << ", tie_resolver=frozen-support-consensus"
                  << ", frozen_support_min_odds=" << cfg.frozen_support_min_odds
                  << ", frozen_h0_prior=" << frozen_h0_prior_path
                  << ", frozen_oligo_stats=" << frozen_oligo_stats_path;
    }
    std::cerr << ")\n";
    if (cfg.direct_best_edit_decode && cfg.direct_best_non_acgt_dp_fallback) {
        std::cerr << "direct-best non-ACGT handling: null-skip counted; DP fallback disabled\n";
    }

    struct RunResult {
        long total = 0;
        long assigned = 0;
        long frozen_support_attempted = 0;
        long frozen_support_resolved = 0;
        long non_acgt_sequestered = 0;
        long non_acgt_dp_checked = 0;
        long non_acgt_null_reads = 0;
        long non_acgt_null_queries = 0;
        long prebuilt_fallbacks = 0;
        DpOnlyMemoStats dp_miss_stats;
        DecodeStats stats;
    };

    struct HalfOligoMutationCounters {
        std::vector<uint64_t> exposure_h0_h1;
        std::vector<uint64_t> h0_pair_reads;
        std::vector<uint64_t> h1_mutation_reads;
        std::vector<uint64_t> h1_substitution_reads;
        std::vector<uint64_t> h1_length_change_reads;

        void resize(size_t size) {
            exposure_h0_h1.resize(size, 0);
            h0_pair_reads.resize(size, 0);
            h1_mutation_reads.resize(size, 0);
            h1_substitution_reads.resize(size, 0);
            h1_length_change_reads.resize(size, 0);
        }

        void add(const HalfOligoMutationCounters& other) {
            for (size_t i = 0; i < exposure_h0_h1.size(); ++i) {
                exposure_h0_h1[i] += other.exposure_h0_h1[i];
                h0_pair_reads[i] += other.h0_pair_reads[i];
                h1_mutation_reads[i] += other.h1_mutation_reads[i];
                h1_substitution_reads[i] += other.h1_substitution_reads[i];
                h1_length_change_reads[i] += other.h1_length_change_reads[i];
            }
        }
    };

    struct OligoMutationCounters {
        HalfOligoMutationCounters bc1;
        HalfOligoMutationCounters bc2;
        uint64_t h0_reads = 0;
        uint64_t h1_reads = 0;
        uint64_t h1_substitution_reads = 0;
        uint64_t h1_length_change_reads = 0;

        void resize(size_t bc1_size, size_t bc2_size) {
            bc1.resize(bc1_size);
            bc2.resize(bc2_size);
        }

        void add(const OligoMutationCounters& other) {
            bc1.add(other.bc1);
            bc2.add(other.bc2);
            h0_reads += other.h0_reads;
            h1_reads += other.h1_reads;
            h1_substitution_reads += other.h1_substitution_reads;
            h1_length_change_reads += other.h1_length_change_reads;
        }
    };

    OligoMutationCounters oligo_mutation_stats;
    oligo_mutation_stats.resize(cfg.bc1_oligos.size(), cfg.bc2_oligos.size());

    auto run_decode_pass = [&](bool use_prebuilt, bool write_output, const char* label,
                               int pass_threads, long read_skip, long read_limit,
                               bool collect_oligo_mutation_stats) -> RunResult {
        auto pass_start = std::chrono::steady_clock::now();
        std::ofstream out;
        if (write_output) {
            out.open(out_path);
            out << "read_id\tuniverse_unique_assigned\trow2\tcol2\tfull_start";
            if (cfg.emit_decode_annotations) write_decode_annotation_header(out, cfg);
            out << "\n";
        }
        std::ofstream candidate_preserving_out;
        const bool write_candidate_preserving =
            write_output && !cfg.candidate_preserving_out.empty();
        if (write_candidate_preserving) {
            candidate_preserving_out.open(cfg.candidate_preserving_out);
            if (!candidate_preserving_out) {
                std::cerr << "ERROR: cannot open candidate-preserving output "
                          << cfg.candidate_preserving_out << "\n";
                std::exit(1);
            }
            candidate_preserving_out
                << "global_ordinal\tread_id\traw_umi\tcandidate_count\trow2\tcol2\tmin_tier\t"
                   "bc1_edit\tbc2_edit\tbc1_obs_len\tbc2_obs_len\tfull_start\t"
                   "log_sequence_likelihood\n";
        }
        std::ofstream candidate_audit_out;
        const bool write_candidate_audit =
            write_output && !cfg.candidate_audit_out.empty();
        if (write_candidate_audit) {
            candidate_audit_out.open(cfg.candidate_audit_out);
            if (!candidate_audit_out) {
                std::cerr << "ERROR: cannot open candidate audit output "
                          << cfg.candidate_audit_out << "\n";
                std::exit(1);
            }
            write_candidate_audit_header(candidate_audit_out);
            candidate_audit_out << "\n";
        }
        std::ofstream latent_candidate_out;
        const bool write_latent_candidates =
            write_output && !cfg.latent_candidate_out.empty();
        if (write_latent_candidates) {
            latent_candidate_out.open(cfg.latent_candidate_out);
            if (!latent_candidate_out) {
                std::cerr << "ERROR: cannot open latent-candidate output "
                          << cfg.latent_candidate_out << "\n";
                std::exit(1);
            }
            latent_candidate_out
                << "read_id\tbest_score0\tall_candidate_count\tall_candidate_traces\n";
        }
        std::ofstream non_acgt_out;
        const bool write_non_acgt =
            write_output && cfg.direct_best_edit_decode && !cfg.non_acgt_reads_out.empty();
        if (write_non_acgt) {
            non_acgt_out.open(cfg.non_acgt_reads_out);
            if (!non_acgt_out) {
                std::cerr << "ERROR: cannot open non-ACGT sidecar "
                          << cfg.non_acgt_reads_out << "\n";
                std::exit(1);
            }
            non_acgt_out << "read_id\tstatus\tseq\n";
        }

        const size_t BATCH = 1 << 21;  // ~2M reads/batch
        std::vector<Record> batch;
        batch.reserve(BATCH);
        std::atomic<long> total{0}, assigned{0};
        std::atomic<long> frozen_support_attempted{0}, frozen_support_resolved{0};
        std::atomic<long> non_acgt_sequestered{0};
        std::atomic<long> non_acgt_dp_checked{0};
        std::atomic<long> non_acgt_null_reads{0};
        std::atomic<long> non_acgt_null_queries{0};
        std::atomic<long> prebuilt_fallbacks{0};
        DecodeStats run_stats;
        DpOnlyMemoStats run_dp_miss_stats;
        long batch_index = 0;

        auto process_batch = [&](std::vector<Record>& recs) {
            const long batch_no = ++batch_index;
            auto batch_start = std::chrono::steady_clock::now();
            size_t n = recs.size();
            std::vector<std::string> lines(write_output ? n : 0);
            std::vector<std::string> candidate_audit_lines(write_candidate_audit ? n : 0);
            std::vector<std::string> candidate_preserving_lines(
                write_candidate_preserving ? n : 0);
            std::vector<std::string> latent_candidate_lines(
                write_latent_candidates ? n : 0);
            std::vector<char> non_acgt_flags(write_non_acgt ? n : 0, 0);
            int T = std::min<int>(pass_threads, static_cast<int>(n));
            if (T < 1) T = 1;
            std::vector<std::thread> pool;
            std::vector<DecodeStats> thread_stats(T);
            std::vector<DpOnlyMemoStats> thread_dp_miss_stats(T);
            std::vector<OligoMutationCounters> thread_oligo_mutation_stats(T);
            if (collect_oligo_mutation_stats) {
                for (auto& counters : thread_oligo_mutation_stats) {
                    counters.resize(cfg.bc1_oligos.size(), cfg.bc2_oligos.size());
                }
            }
            std::atomic<long> asg{0};
            std::atomic<long> batch_non_acgt_sequestered{0};
            std::atomic<long> batch_non_acgt_dp_checked{0};
            std::atomic<long> batch_non_acgt_null_reads{0};
            std::atomic<long> batch_non_acgt_null_queries{0};
            auto worker = [&](int tid) {
                std::vector<PackedHit> lk1, lk2, cache_hits, bc1_full_hits, bc2_full_hits;
                std::vector<FullCandidate> bc1_valid, bc2_valid;
                LookupScratch lookup_scratch;
                ObservedDecodeCacheLocals observed_locals;
                for (size_t i = tid; i < n; i += T) {
                    bool prebuilt_cache_miss = false;
                    Decoded d;
                    const OracleTarget* oracle = nullptr;
                    if (write_candidate_audit) {
                        auto oracle_it = sr_oracle.find(recs[i].read_id);
                        if (oracle_it != sr_oracle.end()) oracle = &oracle_it->second;
                    }
                    if (cfg.direct_tiered_h2_decode) {
                        d = decode_record_direct_tiered_h2(
                            recs[i].seq, cfg,
                            bc1_tiered_h2_query_lengths, bc2_tiered_h2_query_lengths,
                            oracle);
                    } else if (cfg.direct_tiered_fixed_decode) {
                        d = decode_record_direct_tiered_fixed(
                            recs[i].seq, cfg,
                            bc1_tiered_fixed_width, bc2_tiered_fixed_width,
                            oracle);
                    } else if (cfg.direct_tiered_decode) {
                        d = decode_record_direct_tiered(
                            recs[i].seq, cfg,
                            bc1_tiered_exact_lengths, bc2_tiered_exact_lengths,
                            bc1_tiered_e1_query_lengths, bc2_tiered_e1_query_lengths,
                            oracle);
                    } else if (cfg.direct_best_edit_decode) {
                        d = decode_record_direct_best_edit(recs[i].seq, cfg, observed_locals,
                                                           oracle);
                    } else if (use_prebuilt && cfg.direct_full_decode) {
                        d = decode_record_direct_prebuilt(recs[i].seq, cfg, observed_locals,
                                                          cfg.observed_cache_prebuild_fallback,
                                                          &prebuilt_cache_miss);
                    } else if (use_prebuilt) {
                        d = decode_record_prebuilt(recs[i].seq, cfg, bc1_valid, bc2_valid,
                                                   observed_locals,
                                                   cfg.observed_cache_prebuild_fallback,
                                                   &prebuilt_cache_miss);
                    } else {
                        d = decode_record(recs[i].seq, cfg, lk1, lk2,
                                          lookup_scratch, cache_hits,
                                          bc1_full_hits, bc2_full_hits,
                                          bc1_valid, bc2_valid, observed_locals,
                                          cfg.observed_cache_experiment,
                                          emit_stats ? &thread_stats[tid] : nullptr,
                                          oracle);
                    }
                    if (use_prebuilt && prebuilt_cache_miss) {
                        prebuilt_fallbacks.fetch_add(1, std::memory_order_relaxed);
                        d = decode_record(recs[i].seq, cfg, lk1, lk2,
                                          lookup_scratch, cache_hits,
                                          bc1_full_hits, bc2_full_hits,
                                          bc1_valid, bc2_valid, observed_locals,
                                          false,
                                          emit_stats ? &thread_stats[tid] : nullptr,
                                          oracle);
                    }
                    if (d.sequestered_non_acgt) {
                        non_acgt_sequestered.fetch_add(1, std::memory_order_relaxed);
                        batch_non_acgt_sequestered.fetch_add(1, std::memory_order_relaxed);
                        if (write_non_acgt) non_acgt_flags[i] = 'S';
                    }
                    if (d.non_acgt_dp_checked) {
                        non_acgt_dp_checked.fetch_add(1, std::memory_order_relaxed);
                        batch_non_acgt_dp_checked.fetch_add(1, std::memory_order_relaxed);
                        if (write_non_acgt) non_acgt_flags[i] = 'D';
                    }
                    if (d.non_acgt_null_queries > 0) {
                        non_acgt_null_reads.fetch_add(1, std::memory_order_relaxed);
                        non_acgt_null_queries.fetch_add(d.non_acgt_null_queries,
                                                        std::memory_order_relaxed);
                        batch_non_acgt_null_reads.fetch_add(1, std::memory_order_relaxed);
                        batch_non_acgt_null_queries.fetch_add(d.non_acgt_null_queries,
                                                              std::memory_order_relaxed);
                        if (write_non_acgt && !non_acgt_flags[i]) non_acgt_flags[i] = 'N';
                    }
                    if (d.frozen_support_attempted) {
                        frozen_support_attempted.fetch_add(1, std::memory_order_relaxed);
                    }
                    if (d.frozen_support_resolved) {
                        frozen_support_resolved.fetch_add(1, std::memory_order_relaxed);
                    }
                    if (write_output) {
                        std::string& line = lines[i];
                        line.reserve(recs[i].read_id.size() +
                                     (cfg.emit_decode_annotations ? 240 : 24));
                        line += recs[i].read_id;
                        if (d.assigned) {
                            line += "\tTrue\t";
                            line += std::to_string(d.row2);
                            line += '\t';
                            line += std::to_string(d.col2);
                            line += '\t';
                            line += std::to_string(d.full_start);
                        } else {
                            line += "\tFalse\t\t\t";
                        }
                        if (cfg.emit_decode_annotations) append_decode_annotations(line, d, cfg);
                    }
                    if (write_candidate_preserving) {
                        candidate_preserving_lines[i] =
                            make_candidate_preserving_lines(recs[i], d, cfg);
                    }
                    if (write_latent_candidates) {
                        latent_candidate_lines[i] =
                            make_latent_candidate_line(recs[i].read_id, d);
                    }
                    if (write_candidate_audit && oracle) {
                        const int current_edit_max =
                            std::max(cfg.bc1_full_hamming, cfg.bc2_full_hamming);
                        OracleTargetAudit current = audit_oracle_target_envelope(
                            recs[i].seq, cfg, *oracle, cfg.full_start_min,
                            cfg.full_start_max, 0, current_edit_max,
                            cfg.bc2_offset_window, true);
                        OracleTargetAudit current_no_anchor = audit_oracle_target_envelope(
                            recs[i].seq, cfg, *oracle, cfg.full_start_min,
                            cfg.full_start_max, 0, current_edit_max,
                            cfg.bc2_offset_window, false);
                        OracleTargetAudit broad = audit_oracle_target_envelope(
                            recs[i].seq, cfg, *oracle,
                            cfg.full_start_min - cfg.audit_full_start_padding,
                            cfg.full_start_max + cfg.audit_full_start_padding,
                            cfg.audit_bc2_start_padding, cfg.audit_edit_max,
                            cfg.audit_edit_max, true);
                        OracleTargetAudit broad_no_anchor = audit_oracle_target_envelope(
                            recs[i].seq, cfg, *oracle,
                            cfg.full_start_min - cfg.audit_full_start_padding,
                            cfg.full_start_max + cfg.audit_full_start_padding,
                            cfg.audit_bc2_start_padding, cfg.audit_edit_max,
                            cfg.audit_edit_max, false);
                        candidate_audit_lines[i] = make_candidate_audit_line(
                            recs[i].read_id, *oracle, d, current, current_no_anchor,
                            broad, broad_no_anchor, cfg);
                    }
                    if (d.assigned) asg.fetch_add(1, std::memory_order_relaxed);
                    if (collect_oligo_mutation_stats && d.assigned &&
                        !d.frozen_support_resolved &&
                        d.row2 >= 0 && d.row2 < static_cast<int>(cfg.bc2_oligos.size()) &&
                        d.col2 >= 0 && d.col2 < static_cast<int>(cfg.bc1_oligos.size())) {
                        const int total_edit = d.bc1_edit + d.bc2_edit;
                        if (d.bc1_edit >= 0 && d.bc2_edit >= 0 &&
                            (total_edit == 0 || total_edit == 1)) {
                            auto& counters = thread_oligo_mutation_stats[tid];
                            const size_t bc1_index = static_cast<size_t>(d.col2);
                            const size_t bc2_index = static_cast<size_t>(d.row2);
                            ++counters.bc1.exposure_h0_h1[bc1_index];
                            ++counters.bc2.exposure_h0_h1[bc2_index];
                            if (total_edit == 0) {
                                ++counters.h0_reads;
                                ++counters.bc1.h0_pair_reads[bc1_index];
                                ++counters.bc2.h0_pair_reads[bc2_index];
                            } else {
                                ++counters.h1_reads;
                                const bool mutation_is_bc1 = d.bc1_edit == 1;
                                auto& half = mutation_is_bc1 ? counters.bc1 : counters.bc2;
                                const size_t index = mutation_is_bc1 ? bc1_index : bc2_index;
                                const int observed_length = mutation_is_bc1 ? d.bc1_obs_len : d.bc2_obs_len;
                                const int target_length = mutation_is_bc1 ? cfg.bc1_len[bc1_index]
                                                                          : cfg.bc2_len[bc2_index];
                                ++half.h1_mutation_reads[index];
                                if (observed_length == target_length) {
                                    ++counters.h1_substitution_reads;
                                    ++half.h1_substitution_reads[index];
                                } else {
                                    ++counters.h1_length_change_reads;
                                    ++half.h1_length_change_reads[index];
                                }
                            }
                        }
                    }
                }
                if (cfg.observed_cache_experiment && !use_prebuilt) {
                    observed_locals.flush(cfg);
                }
                if (use_prebuilt) {
                    thread_dp_miss_stats[tid] = observed_locals.dp_miss_stats();
                }
            };
            for (int t = 0; t < T; ++t) pool.emplace_back(worker, t);
            for (auto& th : pool) th.join();
            if (emit_stats && !use_prebuilt) {
                for (auto& stats : thread_stats) run_stats.add(stats);
            }
            DpOnlyMemoStats batch_dp_miss_stats;
            if (use_prebuilt) {
                for (auto& stats : thread_dp_miss_stats) batch_dp_miss_stats.add(stats);
            }
            if (collect_oligo_mutation_stats) {
                for (const auto& counters : thread_oligo_mutation_stats) {
                    oligo_mutation_stats.add(counters);
                }
            }
            auto decode_done = std::chrono::steady_clock::now();
            if (write_output) {
                for (auto& l : lines) out << l << '\n';
            }
            if (write_candidate_audit) {
                for (auto& l : candidate_audit_lines) {
                    if (!l.empty()) candidate_audit_out << l << '\n';
                }
            }
            if (write_candidate_preserving) {
                for (auto& l : candidate_preserving_lines) {
                    if (!l.empty()) candidate_preserving_out << l;
                }
            }
            if (write_latent_candidates) {
                for (auto& l : latent_candidate_lines) latent_candidate_out << l << '\n';
            }
            if (write_non_acgt) {
                for (size_t i = 0; i < n; ++i) {
                    if (non_acgt_flags[i]) {
                        non_acgt_out << recs[i].read_id << '\t'
                                     << (non_acgt_flags[i] == 'D' ? "dp_checked" :
                                         (non_acgt_flags[i] == 'N' ? "null_skipped" : "sequestered"))
                                     << '\t' << recs[i].seq << '\n';
                    }
                }
            }
            auto write_done = std::chrono::steady_clock::now();
            const long batch_assigned = asg.load();
            const long total_after = total.fetch_add(static_cast<long>(n)) + static_cast<long>(n);
            const long assigned_after = assigned.fetch_add(batch_assigned) + batch_assigned;
            if (use_prebuilt) run_dp_miss_stats.add(batch_dp_miss_stats);
            if (cfg.progress_log_every_batches > 0 &&
                (batch_no % cfg.progress_log_every_batches) == 0) {
                const double decode_seconds =
                    std::chrono::duration<double>(decode_done - batch_start).count();
                const double write_seconds =
                    std::chrono::duration<double>(write_done - decode_done).count();
                const double elapsed_seconds =
                    std::chrono::duration<double>(write_done - pass_start).count();
                std::cerr << label
                          << " batch=" << batch_no
                          << " batch_reads=" << n
                          << " batch_unique_assigned=" << batch_assigned
                          << " batch_rate=" << (n ? double(batch_assigned) / double(n) : 0.0)
                          << " decode_seconds=" << decode_seconds
                          << " write_seconds=" << write_seconds
                          << " total_reads=" << total_after
                          << " total_unique_assigned=" << assigned_after
                          << " total_rate=" << (total_after ? double(assigned_after) / double(total_after) : 0.0)
                          << " elapsed_seconds=" << elapsed_seconds
                          << " reads_per_second=" << (elapsed_seconds ? double(total_after) / elapsed_seconds : 0.0)
                          << (batch_non_acgt_sequestered.load() > 0 ? " batch_non_acgt_sequestered=" : "")
                          << (batch_non_acgt_sequestered.load() > 0 ? std::to_string(batch_non_acgt_sequestered.load()) : "")
                          << (batch_non_acgt_dp_checked.load() > 0 ? " batch_non_acgt_dp_checked=" : "")
                          << (batch_non_acgt_dp_checked.load() > 0 ? std::to_string(batch_non_acgt_dp_checked.load()) : "")
                          << (batch_non_acgt_null_reads.load() > 0 ? " batch_non_acgt_null_reads=" : "")
                          << (batch_non_acgt_null_reads.load() > 0 ? std::to_string(batch_non_acgt_null_reads.load()) : "")
                          << (batch_non_acgt_null_queries.load() > 0 ? " batch_non_acgt_null_queries=" : "")
                          << (batch_non_acgt_null_queries.load() > 0 ? std::to_string(batch_non_acgt_null_queries.load()) : "")
                          << "\n";
            }
        };

        char buf[1 << 16];
        long records_seen = 0;
        long records_loaded = 0;
        bool reached_limit = false;
        for (auto& path : r1_files) {
            gzFile gz = gzopen(path.c_str(), "rb");
            if (!gz) {
                std::cerr << "ERROR: cannot open " << path << "\n";
                std::exit(1);
            }
            int line_no = 0;
            Record cur;
            while (gzgets(gz, buf, sizeof(buf))) {
                int m = line_no & 3;
                if (m == 0) cur.read_id = read_id_from_header(buf);
                else if (m == 1) {
                    cur.seq = buf;
                    while (!cur.seq.empty() &&
                           (cur.seq.back() == '\n' || cur.seq.back() == '\r')) {
                        cur.seq.pop_back();
                    }
                } else if (m == 3) {
                    cur.qual = buf;
                    while (!cur.qual.empty() &&
                           (cur.qual.back() == '\n' || cur.qual.back() == '\r')) {
                        cur.qual.pop_back();
                    }
                    if (records_seen < read_skip) {
                        ++records_seen;
                        ++line_no;
                        cur = Record();
                        continue;
                    }
                    if (read_limit > 0 && records_loaded >= read_limit) {
                        reached_limit = true;
                        break;
                    }
                    cur.global_ordinal = static_cast<uint64_t>(records_seen);
                    batch.push_back(std::move(cur));
                    cur = Record();
                    ++records_seen;
                    ++records_loaded;
                    if (batch.size() >= BATCH) {
                        process_batch(batch);
                        batch.clear();
                    }
                    if (read_limit > 0 && records_loaded >= read_limit) {
                        reached_limit = true;
                        break;
                    }
                }
                ++line_no;
            }
            gzclose(gz);
            if (reached_limit) break;
        }
        if (!batch.empty()) process_batch(batch);
        if (write_output) out.close();
        if (write_candidate_audit) candidate_audit_out.close();
        if (write_candidate_preserving) candidate_preserving_out.close();
        if (write_latent_candidates) latent_candidate_out.close();
        if (write_non_acgt) non_acgt_out.close();

        RunResult result;
        result.total = total.load();
        result.assigned = assigned.load();
        result.frozen_support_attempted = frozen_support_attempted.load();
        result.frozen_support_resolved = frozen_support_resolved.load();
        result.non_acgt_sequestered = non_acgt_sequestered.load();
        result.non_acgt_dp_checked = non_acgt_dp_checked.load();
        result.non_acgt_null_reads = non_acgt_null_reads.load();
        result.non_acgt_null_queries = non_acgt_null_queries.load();
        result.prebuilt_fallbacks = prebuilt_fallbacks.load();
        result.dp_miss_stats = run_dp_miss_stats;
        result.stats = run_stats;
        double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - pass_start).count();
        std::cerr << label << ": reads=" << result.total
                  << (read_skip > 0 ? " skipped_reads=" : "")
                  << (read_skip > 0 ? std::to_string(read_skip) : "")
                  << (read_limit > 0 ? " sample_limited=true" : "")
                  << " unique_assigned=" << result.assigned
                  << (cfg.frozen_support_consensus ? " frozen_support_attempted=" : "")
                  << (cfg.frozen_support_consensus
                          ? std::to_string(result.frozen_support_attempted) : "")
                  << (cfg.frozen_support_consensus ? " frozen_support_resolved=" : "")
                  << (cfg.frozen_support_consensus
                          ? std::to_string(result.frozen_support_resolved) : "")
                  << (result.non_acgt_sequestered > 0 ? " non_acgt_sequestered=" : "")
                  << (result.non_acgt_sequestered > 0 ? std::to_string(result.non_acgt_sequestered) : "")
                  << (result.non_acgt_dp_checked > 0 ? " non_acgt_dp_checked=" : "")
                  << (result.non_acgt_dp_checked > 0 ? std::to_string(result.non_acgt_dp_checked) : "")
                  << (result.non_acgt_null_reads > 0 ? " non_acgt_null_reads=" : "")
                  << (result.non_acgt_null_reads > 0 ? std::to_string(result.non_acgt_null_reads) : "")
                  << (result.non_acgt_null_queries > 0 ? " non_acgt_null_queries=" : "")
                  << (result.non_acgt_null_queries > 0 ? std::to_string(result.non_acgt_null_queries) : "")
                  << (result.prebuilt_fallbacks > 0 ? " prebuilt_fallback_reads=" : "")
                  << (result.prebuilt_fallbacks > 0 ? std::to_string(result.prebuilt_fallbacks) : "")
                  << (result.dp_miss_stats.lookup_calls > 0 ? " dp_miss_lookups=" : "")
                  << (result.dp_miss_stats.lookup_calls > 0 ? std::to_string(result.dp_miss_stats.lookup_calls) : "")
                  << (result.dp_miss_stats.lookup_calls > 0 ? " dp_miss_local_hits=" : "")
                  << (result.dp_miss_stats.lookup_calls > 0 ? std::to_string(result.dp_miss_stats.memo_hits) : "")
                  << (result.dp_miss_stats.lookup_calls > 0 ? " dp_miss_dp_calls=" : "")
                  << (result.dp_miss_stats.lookup_calls > 0 ? std::to_string(result.dp_miss_stats.dp_calls) : "")
                  << " rate=" << (result.total ? double(result.assigned) / result.total : 0.0)
                  << " seconds=" << elapsed
                  << "\n";
        return result;
    };

    RunResult result;
    if (cfg.observed_cache_prebuild) {
        RunResult build_result = run_decode_pass(false, false, "prebuild-pass", threads, 0,
                                                 cfg.observed_cache_prebuild_sample_reads, false);
        (void)build_result;
        cfg.bc1_anchor_observed.freeze_into(cfg.bc1_anchor_prebuilt);
        cfg.bc2_anchor_observed.freeze_into(cfg.bc2_anchor_prebuilt);
        cfg.bc1_full_observed.freeze_into(cfg.bc1_full_prebuilt);
        cfg.bc2_full_observed.freeze_into(cfg.bc2_full_prebuilt);
        result = run_decode_pass(true, !suppress_read_output, "decode-pass", threads,
                                 cfg.decode_skip_reads, cfg.decode_read_limit,
                                 !oligo_mutation_stats_out.empty());
    } else {
        result = run_decode_pass(false, !suppress_read_output, "decode", threads,
                                 cfg.decode_skip_reads, cfg.decode_read_limit,
                                 !oligo_mutation_stats_out.empty());
    }
    if (!oligo_mutation_stats_out.empty()) {
        std::ofstream mutation_out(oligo_mutation_stats_out);
        if (!mutation_out) {
            std::cerr << "ERROR: cannot open oligo mutation statistics output "
                      << oligo_mutation_stats_out << "\n";
            return 1;
        }
        mutation_out << "barcode_half\toligo_index\toligo_sequence\toligo_length"
                     << "\texposure_h0_h1\th0_pair_reads\th1_pair_reads"
                     << "\th1_mutation_reads\th1_substitution_reads"
                     << "\th1_length_change_reads\tempirical_h1_rate"
                     << "\tempirical_h1_substitution_rate\n";
        auto write_half = [&](const char* half_name,
                              const std::vector<std::string>& oligos,
                              const HalfOligoMutationCounters& counters) {
            for (size_t i = 0; i < oligos.size(); ++i) {
                const uint64_t exposure = counters.exposure_h0_h1[i];
                const uint64_t h0 = counters.h0_pair_reads[i];
                mutation_out << half_name << '\t' << i << '\t' << oligos[i] << '\t'
                             << oligos[i].size() << '\t' << exposure << '\t' << h0 << '\t'
                             << (exposure - h0) << '\t' << counters.h1_mutation_reads[i] << '\t'
                             << counters.h1_substitution_reads[i] << '\t'
                             << counters.h1_length_change_reads[i] << '\t';
                if (exposure > 0) {
                    mutation_out << (static_cast<double>(counters.h1_mutation_reads[i]) /
                                     static_cast<double>(exposure));
                }
                mutation_out << '\t';
                if (exposure > 0) {
                    mutation_out << (static_cast<double>(counters.h1_substitution_reads[i]) /
                                     static_cast<double>(exposure));
                }
                mutation_out << '\n';
            }
        };
        write_half("BC1", cfg.bc1_oligos, oligo_mutation_stats.bc1);
        write_half("BC2", cfg.bc2_oligos, oligo_mutation_stats.bc2);
        std::cerr << "oligo mutation stats: h0_reads=" << oligo_mutation_stats.h0_reads
                  << " h1_reads=" << oligo_mutation_stats.h1_reads
                  << " h1_substitution_reads=" << oligo_mutation_stats.h1_substitution_reads
                  << " h1_length_change_reads=" << oligo_mutation_stats.h1_length_change_reads
                  << " path=" << oligo_mutation_stats_out << "\n";
    }
    if (emit_stats) print_decode_stats(result.stats);
    if (cfg.observed_cache_experiment) {
        print_observed_cache_stats(cfg.bc1_anchor_observed);
        print_observed_cache_stats(cfg.bc2_anchor_observed);
        print_observed_cache_stats(cfg.bc1_full_observed);
        print_observed_cache_stats(cfg.bc2_full_observed);
    }
    if (cfg.observed_cache_prebuild) {
        print_prebuilt_observed_cache_stats(cfg.bc1_anchor_prebuilt);
        print_prebuilt_observed_cache_stats(cfg.bc2_anchor_prebuilt);
        print_prebuilt_observed_cache_stats(cfg.bc1_full_prebuilt);
        print_prebuilt_observed_cache_stats(cfg.bc2_full_prebuilt);
    }
    return 0;
}
#endif

namespace spatial_r1_decoder {

struct Decoder::Impl {
    InternalDecoderConfig config;
    std::vector<int> bc1QueryLengths;
    std::vector<int> bc2QueryLengths;

    explicit Impl(const spatial_r1_decoder::Config &input) {
        if (input.bc1OligosPath.empty() || input.bc2OligosPath.empty()) {
            throw std::invalid_argument("Visium HD R1 decoder requires both oligo files");
        }
        if (input.gridRows == 0 || input.gridColumns == 0) {
            throw std::invalid_argument("Visium HD R1 decoder grid dimensions must be positive");
        }
        if (input.fullStartMin < 0 || input.fullStartMax < input.fullStartMin) {
            throw std::invalid_argument("Visium HD R1 decoder start window is invalid");
        }

        load_lines({input.bc1OligosPath}, config.bc1_oligos);
        load_lines({input.bc2OligosPath}, config.bc2_oligos);
        if (config.bc1_oligos.size() != input.gridColumns ||
            config.bc2_oligos.size() != input.gridRows) {
            throw std::runtime_error("Visium HD oligo axes do not match the coordinate grid");
        }
        config.grid_rows = static_cast<int>(input.gridRows);
        config.grid_cols = static_cast<int>(input.gridColumns);
        config.full_start_min = input.fullStartMin;
        config.full_start_max = input.fullStartMax;
        config.direct_tiered_h2_decode = true;

        config.bc1_len.resize(config.bc1_oligos.size());
        config.bc2_len.resize(config.bc2_oligos.size());
        for (std::size_t i = 0; i < config.bc1_oligos.size(); ++i) {
            config.bc1_len[i] = static_cast<int>(config.bc1_oligos[i].size());
        }
        for (std::size_t i = 0; i < config.bc2_oligos.size(); ++i) {
            config.bc2_len[i] = static_cast<int>(config.bc2_oligos[i].size());
        }

        const std::vector<int> bc1Lengths = unique_sorted_lengths(config.bc1_len);
        const std::vector<int> bc2Lengths = unique_sorted_lengths(config.bc2_len);
        bc1QueryLengths = edit_query_lengths(bc1Lengths, 2);
        bc2QueryLengths = edit_query_lengths(bc2Lengths, 2);
        config.bc1_tiered_h2_lookup.build_best_edit_universe_for_query_lengths(
            config.bc1_oligos, bc1QueryLengths, 2, "bc1_tiered_h2",
            TargetSliceMode::Full, false, false);
        config.bc2_tiered_h2_lookup.build_best_edit_universe_for_query_lengths(
            config.bc2_oligos, bc2QueryLengths, 2, "bc2_tiered_h2",
            TargetSliceMode::Full, false, false);
    }
};

void ExactH0Counts::reset(std::size_t bc1Size, std::size_t bc2Size) {
    bc1.assign(bc1Size, 0);
    bc2.assign(bc2Size, 0);
    reads = 0;
}

void ExactH0Counts::add(const ExactH0Counts &other) {
    if (bc1.size() != other.bc1.size() || bc2.size() != other.bc2.size()) {
        throw std::invalid_argument("cannot merge incompatible Visium HD H0 axes");
    }
    for (std::size_t i = 0; i < bc1.size(); ++i) bc1[i] += other.bc1[i];
    for (std::size_t i = 0; i < bc2.size(); ++i) bc2[i] += other.bc2[i];
    reads += other.reads;
}

bool encodeRawUmi9(const char *sequence, std::size_t length,
                   std::uint32_t &packed) {
    if (sequence == nullptr || length < 9) return false;
    packed = 0;
    for (std::size_t i = 0; i < 9; ++i) {
        std::uint32_t value = 0;
        switch (static_cast<char>(std::toupper(static_cast<unsigned char>(sequence[i])))) {
            case 'A': value = 0; break;
            case 'C': value = 1; break;
            case 'G': value = 2; break;
            case 'T': value = 3; break;
            default: return false;
        }
        packed = (packed << 2) | value;
    }
    return true;
}

std::string decodeRawUmi9(std::uint32_t packed) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::string value(9, 'A');
    for (int i = 8; i >= 0; --i) {
        const std::uint32_t base = packed & 3u;
        value[static_cast<std::size_t>(i)] = bases[base];
        packed >>= 2;
    }
    return value;
}

std::uint32_t packAuditBits(std::uint8_t tier, std::uint8_t bc1Edit,
                            std::uint8_t bc2Edit, std::uint8_t bc1ObservedLength,
                            std::uint8_t bc2ObservedLength,
                            std::uint8_t fullStart) {
    if (tier > 7 || bc1Edit > 7 || bc2Edit > 7 || bc1ObservedLength > 63 ||
        bc2ObservedLength > 63 || fullStart > 63) {
        throw std::invalid_argument("Visium HD candidate audit field is out of range");
    }
    return static_cast<std::uint32_t>(tier)
        | (static_cast<std::uint32_t>(bc1Edit) << 3)
        | (static_cast<std::uint32_t>(bc2Edit) << 6)
        | (static_cast<std::uint32_t>(bc1ObservedLength) << 9)
        | (static_cast<std::uint32_t>(bc2ObservedLength) << 15)
        | (static_cast<std::uint32_t>(fullStart) << 21);
}

Decoder::Decoder(const spatial_r1_decoder::Config &config)
    : impl_(new Impl(config)) {}

Decoder::~Decoder() = default;

std::size_t Decoder::bc1Count() const { return impl_->config.bc1_oligos.size(); }
std::size_t Decoder::bc2Count() const { return impl_->config.bc2_oligos.size(); }
std::uint32_t Decoder::gridRows() const {
    return static_cast<std::uint32_t>(impl_->config.grid_rows);
}
std::uint32_t Decoder::gridColumns() const {
    return static_cast<std::uint32_t>(impl_->config.grid_cols);
}

bool Decoder::decode(const char *sequence, std::size_t sequenceLength,
                     const char *quality, std::size_t qualityLength,
                     Result &result, ExactH0Counts *h0, std::string &error) const {
    result = Result();
    error.clear();
    if (sequence == nullptr || quality == nullptr || sequenceLength != qualityLength) {
        error = "Visium HD R1 sequence/quality is missing or inconsistent";
        return false;
    }
    result.rawUmiValid = encodeRawUmi9(sequence, sequenceLength, result.rawUmi);
    if (sequenceLength >= 9) {
        for (std::size_t index = 0; index < 9; ++index) {
            if (static_cast<char>(std::toupper(
                    static_cast<unsigned char>(sequence[index]))) == 'N') {
                result.rawUmiHadN = true;
                break;
            }
        }
    }

    const std::string sequenceString(sequence, sequenceLength);
    const std::string qualityString(quality, qualityLength);
    const DecodeWindowClassification barcodeClassification =
        classify_direct_best_decode_window(sequenceString, impl_->config);
    result.barcodeHadN = barcodeClassification.n_count != 0;
    result.barcodeNCount =
        static_cast<std::uint8_t>(barcodeClassification.n_count);
    result.barcodeHadUnsupportedBase = barcodeClassification.unsupported;
    if (barcodeClassification.unsupported) return true;
    InternalDecoded decoded = decode_record_direct_tiered_h2(
        sequenceString, impl_->config, impl_->bc1QueryLengths,
        impl_->bc2QueryLengths, nullptr, &barcodeClassification);
    result.barcodeDpChecked = decoded.non_acgt_dp_checked;
    result.decoderAssigned = decoded.assigned;
    if (decoded.bc1_edit >= 0) result.bc1Edit = static_cast<std::uint8_t>(decoded.bc1_edit);
    if (decoded.bc2_edit >= 0) result.bc2Edit = static_cast<std::uint8_t>(decoded.bc2_edit);

    struct Choice {
        InternalCandidateTrace trace;
        double likelihood;
    };
    std::map<std::uint32_t, Choice> choices;
    for (const InternalCandidateTrace &trace : decoded.best_candidate_traces) {
        const int observedLength = trace.bc1_obs_len + trace.bc2_obs_len;
        if (trace.row2 < 0 || trace.col2 < 0 ||
            trace.row2 >= impl_->config.grid_rows ||
            trace.col2 >= impl_->config.grid_cols || trace.full_start < 0 ||
            observedLength <= 0 ||
            trace.full_start + observedLength > static_cast<int>(sequenceLength)) {
            error = "Visium HD decoder produced an invalid candidate span";
            return false;
        }
        const std::string observed = sequenceString.substr(
            static_cast<std::size_t>(trace.full_start),
            static_cast<std::size_t>(observedLength));
        const std::string observedQuality = qualityString.substr(
            static_cast<std::size_t>(trace.full_start),
            static_cast<std::size_t>(observedLength));
        const std::string candidate =
            impl_->config.bc1_oligos[static_cast<std::size_t>(trace.col2)] +
            impl_->config.bc2_oligos[static_cast<std::size_t>(trace.row2)];
        const double likelihood = raw_candidate_alignment_log_likelihood(
            observed, observedQuality, candidate);
        const std::uint32_t coordinate =
            static_cast<std::uint32_t>(trace.row2) *
                static_cast<std::uint32_t>(impl_->config.grid_cols) +
            static_cast<std::uint32_t>(trace.col2);
        const auto found = choices.find(coordinate);
        if (found == choices.end() || likelihood > found->second.likelihood ||
            (likelihood == found->second.likelihood &&
             std::tie(trace.full_start, trace.bc1_obs_len, trace.bc2_obs_len) <
                 std::tie(found->second.trace.full_start,
                          found->second.trace.bc1_obs_len,
                          found->second.trace.bc2_obs_len))) {
            choices[coordinate] = Choice{trace, likelihood};
        }
    }

    result.candidates.reserve(choices.size());
    for (const auto &entry : choices) {
        const InternalCandidateTrace &trace = entry.second.trace;
        Candidate candidate;
        candidate.coordinateIndex = entry.first;
        candidate.auditBits = packAuditBits(
            static_cast<std::uint8_t>(trace.bc1_edit + trace.bc2_edit),
            static_cast<std::uint8_t>(trace.bc1_edit),
            static_cast<std::uint8_t>(trace.bc2_edit),
            static_cast<std::uint8_t>(trace.bc1_obs_len),
            static_cast<std::uint8_t>(trace.bc2_obs_len),
            static_cast<std::uint8_t>(trace.full_start));
        candidate.logSequenceLikelihood = entry.second.likelihood;
        result.candidates.push_back(candidate);
    }

    if (h0 != nullptr && decoded.assigned && !decoded.frozen_support_resolved &&
        decoded.bc1_edit == 0 && decoded.bc2_edit == 0 &&
        decoded.col2 >= 0 && decoded.row2 >= 0 &&
        static_cast<std::size_t>(decoded.col2) < h0->bc1.size() &&
        static_cast<std::size_t>(decoded.row2) < h0->bc2.size()) {
        ++h0->bc1[static_cast<std::size_t>(decoded.col2)];
        ++h0->bc2[static_cast<std::size_t>(decoded.row2)];
        ++h0->reads;
    }
    return true;
}

} // namespace spatial_r1_decoder
