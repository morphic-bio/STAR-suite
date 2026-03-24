// Integration tests for applying clique corrections back to the cg_agg hash.
//
// Current intended behavior:
// - corrections are scoped to the final composite barcode surface
//   (CB + tagIdx + gene)
// - no cross-tag correction is allowed
// - within a group, every matching original key is updated when corrected

#include "UMICorrector.h"
#include "hash_shims_cpp_compat.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

static int g_pass = 0, g_fail = 0;

#define ASSERT_EQ(a, b, msg) do { \
    if ((a) != (b)) { \
        fprintf(stderr, "FAIL [%s:%d] %s: expected %llu, got %llu\n", \
                __FILE__, __LINE__, msg, (unsigned long long)(b), (unsigned long long)(a)); \
        g_fail++; return; \
    } \
} while(0)

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg); \
        g_fail++; return; \
    } \
} while(0)

static uint32_t packUMI(const char *seq) {
    uint32_t packed = 0;
    for (int i = 0; i < 12; ++i) {
        uint32_t base;
        switch (seq[i]) {
            case 'A': base = 0; break;
            case 'C': base = 1; break;
            case 'G': base = 2; break;
            case 'T': base = 3; break;
            default:  base = 0; break;
        }
        packed = (packed << 2) | base;
    }
    return packed;
}

// Reproduces the current logic from runCliqueCorrection and
// applyCliqueCorrectionsToHash, operating on a real khash.
static void test_fanout_multi_tag() {
    uint32_t cbIdx = 42;
    uint16_t geneIdx = 7;
    uint32_t winner_umi = packUMI("AAAAAAAAAAAA");
    uint32_t loser_umi  = packUMI("CAAAAAAAAAAA"); // Hamming-1 from winner

    // Build an inline hash with the loser UMI appearing under 3 different tags.
    // Only tag 1 has a local winner, so only tag 1 should be corrected.
    khash_t(cg_agg) *inlineHash = kh_init(cg_agg);
    int absent;

    // Winner: tag 1, count 100
    uint64_t key_winner_t1 = packCgAggKey(cbIdx, winner_umi, geneIdx, 1);
    khiter_t ki = kh_put(cg_agg, inlineHash, key_winner_t1, &absent);
    kh_val(inlineHash, ki) = 100;

    // Loser: tag 1, count 3
    uint64_t key_loser_t1 = packCgAggKey(cbIdx, loser_umi, geneIdx, 1);
    ki = kh_put(cg_agg, inlineHash, key_loser_t1, &absent);
    kh_val(inlineHash, ki) = 3;

    // Loser: tag 2, count 5
    uint64_t key_loser_t2 = packCgAggKey(cbIdx, loser_umi, geneIdx, 2);
    ki = kh_put(cg_agg, inlineHash, key_loser_t2, &absent);
    kh_val(inlineHash, ki) = 5;

    // Loser: tag 3, count 2
    uint64_t key_loser_t3 = packCgAggKey(cbIdx, loser_umi, geneIdx, 3);
    ki = kh_put(cg_agg, inlineHash, key_loser_t3, &absent);
    kh_val(inlineHash, ki) = 2;

    ASSERT_EQ(kh_size(inlineHash), 4u, "inline hash should have 4 entries");

    // --- Phase 1: Extract entries (same as runCliqueCorrection) ---
    struct FlatEntry {
        uint64_t groupKey;
        uint64_t cgAggKey;
        uint32_t umi24;
        uint32_t count;
    };
    std::vector<FlatEntry> entries;
    for (khiter_t it = kh_begin(inlineHash); it != kh_end(inlineHash); ++it) {
        if (!kh_exist(inlineHash, it)) continue;
        uint64_t key = kh_key(inlineHash, it);
        uint32_t count = kh_val(inlineHash, it);
        uint32_t cb, umi24;
        uint16_t gene;
        uint8_t tag;
        unpackCgAggKey(key, &cb, &umi24, &gene, &tag);
        uint64_t groupKey = (static_cast<uint64_t>(cb) << 24) |
                            (static_cast<uint64_t>(tag) << 16) |
                            static_cast<uint64_t>(gene);
        entries.push_back({groupKey, key, umi24, count});
    }
    std::sort(entries.begin(), entries.end(),
              [](const FlatEntry &a, const FlatEntry &b) { return a.groupKey < b.groupKey; });

    // --- Phase 2: Process groups (same as runCliqueCorrection) ---
    khash_t(cg_agg) *correctionHash = kh_init(cg_agg);
    UMIParams params(1, 2.0, 1000);
    uint32_t totalMerges = 0;

    size_t groupStart = 0;
    while (groupStart < entries.size()) {
        size_t groupEnd = groupStart + 1;
        while (groupEnd < entries.size() && entries[groupEnd].groupKey == entries[groupStart].groupKey) {
            ++groupEnd;
        }

        std::vector<UMICount> counts;
        std::vector<uint64_t> cgAggKeys;
        for (size_t i = groupStart; i < groupEnd; ++i) {
            counts.emplace_back(entries[i].umi24, entries[i].count);
            cgAggKeys.push_back(entries[i].cgAggKey);
        }

        UMICorrectionResult result = UMICorrector::correctClique(counts, params);
        totalMerges += result.merges;

        for (const auto &corr : result.urToUb) {
            if (corr.first == corr.second) continue;
            for (size_t i = 0; i < counts.size(); ++i) {
                if (counts[i].ur == corr.first) {
                    ki = kh_put(cg_agg, correctionHash, cgAggKeys[i], &absent);
                    kh_val(correctionHash, ki) = corr.second;
                }
            }
        }

        groupStart = groupEnd;
    }

    ASSERT_EQ(totalMerges, 1u, "should have 1 merge in total");
    ASSERT_EQ(kh_size(correctionHash), 1u,
              "correction hash must have 1 entry (only the same-tag loser)");

    khiter_t it;
    it = kh_get(cg_agg, correctionHash, key_loser_t1);
    ASSERT_TRUE(it != kh_end(correctionHash), "loser tag1 must be in correction hash");
    ASSERT_EQ(kh_val(correctionHash, it), winner_umi, "loser tag1 corrected to winner");

    it = kh_get(cg_agg, correctionHash, key_loser_t2);
    ASSERT_TRUE(it == kh_end(correctionHash), "loser tag2 must not be cross-tag corrected");

    it = kh_get(cg_agg, correctionHash, key_loser_t3);
    ASSERT_TRUE(it == kh_end(correctionHash), "loser tag3 must not be cross-tag corrected");

    // --- Phase 3: Apply corrections to inline hash (same as applyCliqueCorrectionsToHash) ---
    struct HashUpdate {
        uint64_t oldKey;
        uint64_t newKey;
        uint32_t count;
    };
    std::vector<HashUpdate> updates;

    for (khiter_t iter = kh_begin(inlineHash); iter != kh_end(inlineHash); ++iter) {
        if (!kh_exist(inlineHash, iter)) continue;
        uint64_t key = kh_key(inlineHash, iter);
        uint32_t count = kh_val(inlineHash, iter);
        khiter_t corrIt = kh_get(cg_agg, correctionHash, key);
        if (corrIt == kh_end(correctionHash)) continue;
        uint32_t correctedUmi = kh_val(correctionHash, corrIt);
        uint32_t cb;
        uint16_t gene;
        uint8_t tag;
        unpackCgAggKey(key, &cb, nullptr, &gene, &tag);
        uint64_t newKey = packCgAggKey(cb, correctedUmi, gene, tag);
        updates.push_back({key, newKey, count});
    }

    for (const auto &update : updates) {
        khiter_t old_iter = kh_get(cg_agg, inlineHash, update.oldKey);
        if (old_iter != kh_end(inlineHash))
            kh_del(cg_agg, inlineHash, old_iter);
        ki = kh_put(cg_agg, inlineHash, update.newKey, &absent);
        if (absent)
            kh_val(inlineHash, ki) = update.count;
        else
            kh_val(inlineHash, ki) += update.count;
    }

    // After correction: tag1 merges locally, other tags stay unchanged.
    // winner_umi tag1: original 100, plus loser tag1's 3 = 103
    uint64_t key_corrected_t1 = packCgAggKey(cbIdx, winner_umi, geneIdx, 1);
    it = kh_get(cg_agg, inlineHash, key_corrected_t1);
    ASSERT_TRUE(it != kh_end(inlineHash), "winner tag1 must exist after correction");
    ASSERT_EQ(kh_val(inlineHash, it), 103u, "winner tag1 count: 100 + 3 = 103");

    // loser tag1 should be gone after same-tag correction
    it = kh_get(cg_agg, inlineHash, key_loser_t1);
    ASSERT_TRUE(it == kh_end(inlineHash), "loser tag1 must be deleted");

    // Other tags remain unchanged: no cross-tag correction
    it = kh_get(cg_agg, inlineHash, key_loser_t2);
    ASSERT_TRUE(it != kh_end(inlineHash), "loser tag2 must remain");
    ASSERT_EQ(kh_val(inlineHash, it), 5u, "loser tag2 count remains unchanged");
    it = kh_get(cg_agg, inlineHash, key_loser_t3);
    ASSERT_TRUE(it != kh_end(inlineHash), "loser tag3 must remain");
    ASSERT_EQ(kh_val(inlineHash, it), 2u, "loser tag3 count remains unchanged");

    // Final entry count: winner_t1, loser_t2, loser_t3 = 3
    ASSERT_EQ(kh_size(inlineHash), 3u, "inline hash should have 3 entries after correction");

    kh_destroy(cg_agg, correctionHash);
    kh_destroy(cg_agg, inlineHash);

    g_pass++;
    printf("  PASS test_fanout_multi_tag\n");
}

// Regression test: grouping must remain tag-scoped. A winner on tag 1 must not
// pull in the same loser UMI from tag 2.
static void test_tag_scoped_grouping() {
    uint32_t cbIdx = 1;
    uint16_t geneIdx = 1;
    uint32_t winner_umi = packUMI("TTTTTTTTTTTT");
    uint32_t loser_umi  = packUMI("GTTTTTTTTTTT"); // Hamming-1

    struct FlatEntry {
        uint64_t groupKey;
        uint64_t cgAggKey;
        uint32_t umi24;
        uint32_t count;
    };

    std::vector<FlatEntry> entries = {
        {(static_cast<uint64_t>(cbIdx) << 24) | (static_cast<uint64_t>(1) << 16) | geneIdx,
         packCgAggKey(cbIdx, winner_umi, geneIdx, 1), winner_umi, 200},
        {(static_cast<uint64_t>(cbIdx) << 24) | (static_cast<uint64_t>(1) << 16) | geneIdx,
         packCgAggKey(cbIdx, loser_umi, geneIdx, 1), loser_umi, 4},
        {(static_cast<uint64_t>(cbIdx) << 24) | (static_cast<uint64_t>(2) << 16) | geneIdx,
         packCgAggKey(cbIdx, loser_umi, geneIdx, 2), loser_umi, 6},
    };
    std::sort(entries.begin(), entries.end(),
              [](const FlatEntry &a, const FlatEntry &b) { return a.groupKey < b.groupKey; });

    UMIParams params(1, 2.0, 1000);
    khash_t(cg_agg) *corrHash = kh_init(cg_agg);
    size_t groupStart = 0;
    uint32_t totalMerges = 0;
    while (groupStart < entries.size()) {
        size_t groupEnd = groupStart + 1;
        while (groupEnd < entries.size() && entries[groupEnd].groupKey == entries[groupStart].groupKey) {
            ++groupEnd;
        }

        std::vector<UMICount> counts;
        std::vector<uint64_t> cgAggKeys;
        for (size_t i = groupStart; i < groupEnd; ++i) {
            counts.emplace_back(entries[i].umi24, entries[i].count);
            cgAggKeys.push_back(entries[i].cgAggKey);
        }

        UMICorrectionResult result = UMICorrector::correctClique(counts, params);
        totalMerges += result.merges;
        for (const auto &corr : result.urToUb) {
            if (corr.first == corr.second) continue;
            for (size_t i = 0; i < counts.size(); ++i) {
                if (counts[i].ur == corr.first) {
                    int absent;
                    khiter_t ki = kh_put(cg_agg, corrHash, cgAggKeys[i], &absent);
                    kh_val(corrHash, ki) = corr.second;
                }
            }
        }

        groupStart = groupEnd;
    }

    ASSERT_EQ(totalMerges, 1u, "tag-scoped grouping: one local merge only");
    ASSERT_EQ(kh_size(corrHash), 1u,
              "tag-scoped grouping: only the tag1 loser should be corrected");

    kh_destroy(cg_agg, corrHash);
    g_pass++;
    printf("  PASS test_tag_scoped_grouping\n");
}

int main() {
    printf("=== Fan-out Integration Tests ===\n\n");

    test_fanout_multi_tag();
    test_tag_scoped_grouping();

    printf("\n=== Results: %d passed, %d failed ===\n", g_pass, g_fail);
    return g_fail > 0 ? 1 : 0;
}
