// Unit tests for UMICorrector::correctClique.
// Build: see Makefile in this directory.
// Run:   ./test_clique

#include "UMICorrector.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

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

// Pack a 12-base UMI string to 24-bit packed representation (2 bits/base).
// A=0, C=1, G=2, T=3
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

// Compute Hamming distance between two packed 12-base UMIs
static int hammingDist(uint32_t a, uint32_t b) {
    int dist = 0;
    for (int i = 0; i < 12; ++i) {
        if (((a >> (i*2)) & 3) != ((b >> (i*2)) & 3)) ++dist;
    }
    return dist;
}

// --- Test: singleton component (1 UMI, no neighbors) ---
static void test_singleton() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    std::vector<UMICount> counts = { {u1, 10} };
    UMIParams params(1, 2.0, 1000);

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.merges, 0u, "singleton: no merges");
    ASSERT_EQ(r.components, 1u, "singleton: 1 component");
    ASSERT_EQ(r.uniqueUmisInput, 1u, "singleton: 1 unique UMI input");
    ASSERT_EQ(r.uniqueUmisPostFilter, 1u, "singleton: 1 unique UMI post-filter");
    ASSERT_TRUE(r.urToUb.empty(), "singleton: no corrections");
    g_pass++;
    printf("  PASS test_singleton\n");
}

// --- Test: simple Hamming-1 merge (2 UMIs, 1 edit apart, clear winner) ---
static void test_hamming1_merge() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA"); // winner: high count
    uint32_t u2 = packUMI("CAAAAAAAAAAA"); // loser: 1 edit away, low count
    ASSERT_EQ(hammingDist(u1, u2), 1, "Hamming distance should be 1");

    std::vector<UMICount> counts = { {u1, 100}, {u2, 5} };
    UMIParams params(1, 2.0, 1000);

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.merges, 1u, "hamming1: 1 merge");
    ASSERT_EQ(r.components, 1u, "hamming1: 1 component");
    ASSERT_EQ(r.uniqueUmisInput, 2u, "hamming1: 2 unique UMIs input");
    ASSERT_TRUE(r.urToUb.count(u2) && r.urToUb.at(u2) == u1,
                "hamming1: u2 should map to u1");
    g_pass++;
    printf("  PASS test_hamming1_merge\n");
}

// --- Test: ratio threshold reject (2 UMIs close in count) ---
static void test_ratio_reject() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    uint32_t u2 = packUMI("CAAAAAAAAAAA");

    std::vector<UMICount> counts = { {u1, 10}, {u2, 8} };
    UMIParams params(1, 2.0, 1000); // ratio thresh = 2.0, but 10/8 = 1.25

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.merges, 0u, "ratio_reject: no merges (ratio below threshold)");
    ASSERT_EQ(r.components, 1u, "ratio_reject: 1 component");
    ASSERT_EQ(r.componentsBelowThreshold, 1u, "ratio_reject: 1 below-threshold component");
    ASSERT_TRUE(r.urToUb.empty(), "ratio_reject: no corrections");
    g_pass++;
    printf("  PASS test_ratio_reject\n");
}

// --- Test: ratio threshold accept (clear winner) ---
static void test_ratio_accept() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    uint32_t u2 = packUMI("CAAAAAAAAAAA");

    std::vector<UMICount> counts = { {u1, 100}, {u2, 5} };
    UMIParams params(1, 2.0, 1000); // 100/5 = 20.0, well above 2.0

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.merges, 1u, "ratio_accept: 1 merge");
    ASSERT_EQ(r.componentsBelowThreshold, 0u, "ratio_accept: 0 below-threshold");
    g_pass++;
    printf("  PASS test_ratio_accept\n");
}

// --- Test: min-count filter ---
static void test_min_count_filter() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    uint32_t u2 = packUMI("CAAAAAAAAAAA");

    std::vector<UMICount> counts = { {u1, 100}, {u2, 1} };
    UMIParams params(5, 2.0, 1000); // minCount=5, u2 has only 1 read

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.uniqueUmisInput, 2u, "min_count: 2 unique UMIs input");
    ASSERT_EQ(r.uniqueUmisPostFilter, 1u, "min_count: 1 unique UMI after filter");
    ASSERT_EQ(r.merges, 0u, "min_count: no merges (only 1 UMI passes filter)");
    ASSERT_EQ(r.components, 1u, "min_count: 1 component (the survivor)");
    g_pass++;
    printf("  PASS test_min_count_filter\n");
}

// --- Test: component cap behavior ---
// Oversized connected components should be rejected as a unit and should not be
// fragmented into multiple pseudo-components by the BFS traversal.
static void test_component_cap() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    uint32_t u2 = packUMI("CAAAAAAAAAAA");
    uint32_t u3 = packUMI("CCAAAAAAAAAA");
    ASSERT_EQ(hammingDist(u1, u2), 1, "u1-u2 Hamming distance");
    ASSERT_EQ(hammingDist(u2, u3), 1, "u2-u3 Hamming distance");

    std::vector<UMICount> counts = { {u1, 100}, {u2, 50}, {u3, 10} };
    UMIParams params(1, 2.0, 2); // maxComponentSize=2

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.components, 1u, "cap: one connected component");
    ASSERT_EQ(r.componentsCapped, 1u, "cap: oversized component rejected");
    ASSERT_EQ(r.merges, 0u, "cap: no merges for capped component");
    ASSERT_TRUE(r.urToUb.empty(), "cap: no corrections emitted for capped component");
    g_pass++;
    printf("  PASS test_component_cap\n");
}

// --- Test: two disconnected components ---
static void test_two_components() {
    // Component 1: u1, u2 (Hamming-1)
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    uint32_t u2 = packUMI("CAAAAAAAAAAA");
    // Component 2: u3, u4 (Hamming-1, but >1 edit from u1/u2)
    uint32_t u3 = packUMI("TTTTTTTTTTTT");
    uint32_t u4 = packUMI("GTTTTTTTTTTT");
    ASSERT_EQ(hammingDist(u1, u2), 1, "u1-u2 Hamming distance");
    ASSERT_EQ(hammingDist(u3, u4), 1, "u3-u4 Hamming distance");
    ASSERT_TRUE(hammingDist(u1, u3) > 1, "u1-u3 should be disconnected");

    std::vector<UMICount> counts = { {u1, 100}, {u2, 5}, {u3, 80}, {u4, 3} };
    UMIParams params(1, 2.0, 1000);

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.components, 2u, "two_comp: 2 components");
    ASSERT_EQ(r.merges, 2u, "two_comp: 2 merges (one per component)");
    ASSERT_TRUE(r.urToUb.count(u2) && r.urToUb.at(u2) == u1,
                "two_comp: u2 -> u1");
    ASSERT_TRUE(r.urToUb.count(u4) && r.urToUb.at(u4) == u3,
                "two_comp: u4 -> u3");
    g_pass++;
    printf("  PASS test_two_components\n");
}

// --- Test: tied top count rejects correction ---
static void test_tie_reject() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    uint32_t u2 = packUMI("CAAAAAAAAAAA");

    std::vector<UMICount> counts = { {u1, 50}, {u2, 50} };
    UMIParams params(1, 1.0, 1000); // even if ratioThresh would allow it, tie should reject

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.components, 1u, "tie: 1 component");
    ASSERT_EQ(r.merges, 0u, "tie: no merge on tied top count");
    ASSERT_TRUE(r.urToUb.empty(), "tie: no correction emitted");
    g_pass++;
    printf("  PASS test_tie_reject\n");
}

// --- Test: empty input ---
static void test_empty() {
    std::vector<UMICount> counts;
    UMIParams params(1, 2.0, 1000);

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.merges, 0u, "empty: no merges");
    ASSERT_EQ(r.components, 0u, "empty: no components");
    ASSERT_EQ(r.uniqueUmisInput, 0u, "empty: 0 UMIs");
    g_pass++;
    printf("  PASS test_empty\n");
}

// --- Test: duplicate UMI entries (same umi24, different tag variants) ---
static void test_duplicate_umi_aggregation() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");

    // Same UMI appearing twice (simulating different tagIdx entries)
    std::vector<UMICount> counts = { {u1, 10}, {u1, 20} };
    UMIParams params(1, 2.0, 1000);

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.uniqueUmisInput, 1u, "dup_agg: 1 unique UMI after aggregation");
    ASSERT_EQ(r.components, 1u, "dup_agg: 1 component");
    ASSERT_EQ(r.merges, 0u, "dup_agg: no merges (single UMI)");
    g_pass++;
    printf("  PASS test_duplicate_umi_aggregation\n");
}

// --- Test: Hamming-2 UMIs are NOT connected ---
static void test_hamming2_not_connected() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    uint32_t u2 = packUMI("CCAAAAAAAAAA"); // 2 edits from u1
    ASSERT_EQ(hammingDist(u1, u2), 2, "u1-u2 should be Hamming-2");

    std::vector<UMICount> counts = { {u1, 100}, {u2, 5} };
    UMIParams params(1, 2.0, 1000);

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.components, 2u, "hamming2: 2 separate components");
    ASSERT_EQ(r.merges, 0u, "hamming2: no merges");
    g_pass++;
    printf("  PASS test_hamming2_not_connected\n");
}

// --- Test: chain of 3 UMIs, all merge to winner ---
static void test_chain_merge() {
    uint32_t u1 = packUMI("AAAAAAAAAAAA");
    uint32_t u2 = packUMI("CAAAAAAAAAAA"); // 1 edit from u1
    uint32_t u3 = packUMI("CCAAAAAAAAAA"); // 1 edit from u2, 2 from u1

    std::vector<UMICount> counts = { {u1, 100}, {u2, 5}, {u3, 3} };
    UMIParams params(1, 2.0, 1000);

    UMICorrectionResult r = UMICorrector::correctClique(counts, params);

    ASSERT_EQ(r.components, 1u, "chain: 1 component (connected via u2)");
    ASSERT_EQ(r.merges, 2u, "chain: 2 merges");
    ASSERT_TRUE(r.urToUb.count(u2) && r.urToUb.at(u2) == u1, "chain: u2 -> u1");
    ASSERT_TRUE(r.urToUb.count(u3) && r.urToUb.at(u3) == u1, "chain: u3 -> u1");
    g_pass++;
    printf("  PASS test_chain_merge\n");
}

int main() {
    printf("=== Clique UMI Correction Unit Tests ===\n\n");

    test_singleton();
    test_hamming1_merge();
    test_ratio_reject();
    test_ratio_accept();
    test_min_count_filter();
    test_component_cap();
    test_two_components();
    test_tie_reject();
    test_empty();
    test_duplicate_umi_aggregation();
    test_hamming2_not_connected();
    test_chain_merge();

    printf("\n=== Results: %d passed, %d failed ===\n", g_pass, g_fail);
    return g_fail > 0 ? 1 : 0;
}
