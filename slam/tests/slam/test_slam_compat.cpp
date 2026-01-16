/**
 * Unit tests for SlamCompat helper functions.
 *
 * Tests the following helpers:
 *   - compatOverlapWeight: overlap-gene weighting
 *   - compatShouldCountPos: trim guard position filtering
 *   - compatIsIntronic: intronic classification
 *   - computeExonOverlap: exon overlap computation
 *
 * Compile:
 *   g++ -std=c++17 -I../../source -o test_slam_compat \
 *     test_slam_compat.cpp ../../source/SlamCompat.cpp
 *
 * Run:
 *   ./test_slam_compat
 */

#include "SlamCompat.h"
#include "Transcript.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <set>

static int g_passed = 0;
static int g_failed = 0;

static void check(bool ok, const std::string& label) {
    if (ok) {
        std::cout << "  PASS: " << label << "\n";
        ++g_passed;
    } else {
        std::cerr << "  FAIL: " << label << "\n";
        ++g_failed;
    }
}

static bool approxEqual(double a, double b, double tol = 1e-9) {
    return std::fabs(a - b) < tol;
}

// ============================================================================
// Test: compatOverlapWeight
// ============================================================================
void test_compatOverlapWeight() {
    std::cout << "\n=== Testing compatOverlapWeight ===\n";
    
    // Test with overlapWeight disabled
    {
        SlamCompatConfig cfg;
        cfg.overlapWeight = false;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        double w = compat.compatOverlapWeight(1.0, 3);
        check(approxEqual(w, 1.0), "disabled: weight unchanged");
    }
    
    // Test with overlapWeight enabled, single gene (no division)
    {
        SlamCompatConfig cfg;
        cfg.overlapWeight = true;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        double w = compat.compatOverlapWeight(1.0, 1);
        check(approxEqual(w, 1.0), "single gene: weight unchanged");
    }
    
    // Test with overlapWeight enabled, multiple genes
    {
        SlamCompatConfig cfg;
        cfg.overlapWeight = true;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        double w2 = compat.compatOverlapWeight(1.0, 2);
        check(approxEqual(w2, 0.5), "2 genes: weight = 0.5");
        
        double w5 = compat.compatOverlapWeight(1.0, 5);
        check(approxEqual(w5, 0.2), "5 genes: weight = 0.2");
        
        // Test with non-unit base weight
        double w3_base05 = compat.compatOverlapWeight(0.5, 3);
        check(approxEqual(w3_base05, 0.5 / 3.0), "3 genes, base=0.5: weight = 0.5/3");
    }
}

// ============================================================================
// Test: compatShouldCountPos
// ============================================================================
void test_compatShouldCountPos() {
    std::cout << "\n=== Testing compatShouldCountPos ===\n";
    
    // Test with no trimming
    {
        SlamCompatConfig cfg;
        cfg.trim5p = 0;
        cfg.trim3p = 0;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        check(compat.compatShouldCountPos(0, 100), "no trim: pos=0 allowed");
        check(compat.compatShouldCountPos(50, 100), "no trim: pos=50 allowed");
        check(compat.compatShouldCountPos(99, 100), "no trim: pos=99 allowed");
    }
    
    // Test with 5' trimming only
    {
        SlamCompatConfig cfg;
        cfg.trim5p = 5;
        cfg.trim3p = 0;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        check(!compat.compatShouldCountPos(0, 100), "trim5=5: pos=0 skipped");
        check(!compat.compatShouldCountPos(4, 100), "trim5=5: pos=4 skipped");
        check(compat.compatShouldCountPos(5, 100), "trim5=5: pos=5 allowed (boundary)");
        check(compat.compatShouldCountPos(50, 100), "trim5=5: pos=50 allowed");
        check(compat.compatShouldCountPos(99, 100), "trim5=5: pos=99 allowed");
    }
    
    // Test with 3' trimming only
    {
        SlamCompatConfig cfg;
        cfg.trim5p = 0;
        cfg.trim3p = 5;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        check(compat.compatShouldCountPos(0, 100), "trim3=5: pos=0 allowed");
        check(compat.compatShouldCountPos(50, 100), "trim3=5: pos=50 allowed");
        check(compat.compatShouldCountPos(94, 100), "trim3=5: pos=94 allowed (boundary)");
        check(!compat.compatShouldCountPos(95, 100), "trim3=5: pos=95 skipped");
        check(!compat.compatShouldCountPos(99, 100), "trim3=5: pos=99 skipped");
    }
    
    // Test with both 5' and 3' trimming
    {
        SlamCompatConfig cfg;
        cfg.trim5p = 10;
        cfg.trim3p = 10;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        check(!compat.compatShouldCountPos(0, 100), "trim5=10,trim3=10: pos=0 skipped");
        check(!compat.compatShouldCountPos(9, 100), "trim5=10,trim3=10: pos=9 skipped");
        check(compat.compatShouldCountPos(10, 100), "trim5=10,trim3=10: pos=10 allowed");
        check(compat.compatShouldCountPos(50, 100), "trim5=10,trim3=10: pos=50 allowed");
        check(compat.compatShouldCountPos(89, 100), "trim5=10,trim3=10: pos=89 allowed");
        check(!compat.compatShouldCountPos(90, 100), "trim5=10,trim3=10: pos=90 skipped");
        check(!compat.compatShouldCountPos(99, 100), "trim5=10,trim3=10: pos=99 skipped");
    }
    
    // Test underflow guard (trim >= mateLen)
    {
        SlamCompatConfig cfg;
        cfg.trim5p = 60;
        cfg.trim3p = 60;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        // Total trim exceeds mate length - all positions should be skipped
        check(!compat.compatShouldCountPos(0, 100), "trim>=len: pos=0 skipped");
        check(!compat.compatShouldCountPos(50, 100), "trim>=len: pos=50 skipped");
        check(!compat.compatShouldCountPos(99, 100), "trim>=len: pos=99 skipped");
    }
    
    // Test mate-2 style coordinates (after conversion from concatenated)
    {
        SlamCompatConfig cfg;
        cfg.trim5p = 3;
        cfg.trim3p = 3;
        std::vector<std::vector<uint32_t>> g2t;
        std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        // Simulate mate-2 with length 75
        uint32_t mate2Len = 75;
        check(!compat.compatShouldCountPos(0, mate2Len), "mate2: pos=0 skipped");
        check(!compat.compatShouldCountPos(2, mate2Len), "mate2: pos=2 skipped");
        check(compat.compatShouldCountPos(3, mate2Len), "mate2: pos=3 allowed");
        check(compat.compatShouldCountPos(71, mate2Len), "mate2: pos=71 allowed");
        check(!compat.compatShouldCountPos(72, mate2Len), "mate2: pos=72 skipped");
    }
}

// ============================================================================
// Test: compatIsIntronic
// ============================================================================
void test_compatIsIntronic() {
    std::cout << "\n=== Testing compatIsIntronic ===\n";
    
    // Build mock transcript data:
    // Gene 0 has 2 transcripts (tr0 and tr1), each with introns
    // Gene 1 has 1 transcript (tr2) with introns
    std::vector<std::vector<uint32_t>> geneToTranscripts = {
        {0, 1},  // Gene 0 -> transcripts 0, 1
        {2}      // Gene 1 -> transcript 2
    };
    
    // Intron intervals (genomic coordinates):
    // tr0: intron at [1050, 1099], [1200, 1299]
    // tr1: intron at [1040, 1089]
    // tr2: intron at [2000, 2100]
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> transcriptIntrons = {
        {{1050, 1099}, {1200, 1299}},  // tr0
        {{1040, 1089}},                  // tr1
        {{2000, 2100}}                   // tr2
    };
    
    // Test with intronic disabled
    {
        SlamCompatConfig cfg;
        cfg.intronic = false;
        auto g2t = geneToTranscripts;
        auto introns = transcriptIntrons;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 1060;
        aln.exons[0][EX_L] = 20;
        
        std::set<uint32_t> candidates = {0};
        std::set<uint32_t> outIntronic;
        bool result = compat.compatIsIntronic(aln, candidates, outIntronic);
        check(!result, "disabled: returns false");
        check(outIntronic.empty(), "disabled: no intronic genes");
    }
    
    // Test single-part gate (spliced read should fail)
    {
        SlamCompatConfig cfg;
        cfg.intronic = true;
        auto g2t = geneToTranscripts;
        auto introns = transcriptIntrons;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        Transcript aln;
        aln.nExons = 2;  // Spliced read
        aln.exons[0][EX_G] = 1000;
        aln.exons[0][EX_L] = 50;
        aln.exons[1][EX_G] = 1100;
        aln.exons[1][EX_L] = 50;
        
        std::set<uint32_t> candidates = {0};
        std::set<uint32_t> outIntronic;
        bool result = compat.compatIsIntronic(aln, candidates, outIntronic);
        check(!result, "spliced read: returns false");
        check(outIntronic.empty(), "spliced read: no intronic genes");
    }
    
    // Test alignment intersecting introns of >1 transcripts for gene 0
    {
        SlamCompatConfig cfg;
        cfg.intronic = true;
        auto g2t = geneToTranscripts;
        auto introns = transcriptIntrons;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        // Alignment at [1055, 1075] intersects:
        //   - tr0 intron [1050, 1099] ✓
        //   - tr1 intron [1040, 1089] ✓
        // Count = 2, so gene 0 should be intronic
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 1055;
        aln.exons[0][EX_L] = 21;  // [1055, 1075]
        
        std::set<uint32_t> candidates = {0};
        std::set<uint32_t> outIntronic;
        bool result = compat.compatIsIntronic(aln, candidates, outIntronic);
        check(result, "2 intron hits: returns true");
        check(outIntronic.count(0) == 1, "2 intron hits: gene 0 is intronic");
    }
    
    // Test alignment intersecting intron of only 1 transcript (should fail)
    {
        SlamCompatConfig cfg;
        cfg.intronic = true;
        auto g2t = geneToTranscripts;
        auto introns = transcriptIntrons;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        // Alignment at [1210, 1230] intersects only:
        //   - tr0 intron [1200, 1299] ✓
        //   - tr1 intron [1040, 1089] ✗ (no overlap)
        // Count = 1, not > 1, so gene 0 should NOT be intronic
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 1210;
        aln.exons[0][EX_L] = 21;  // [1210, 1230]
        
        std::set<uint32_t> candidates = {0};
        std::set<uint32_t> outIntronic;
        bool result = compat.compatIsIntronic(aln, candidates, outIntronic);
        check(!result, "1 intron hit: returns false");
        check(outIntronic.empty(), "1 intron hit: no intronic genes");
    }
    
    // Test alignment that doesn't intersect any intron
    {
        SlamCompatConfig cfg;
        cfg.intronic = true;
        auto g2t = geneToTranscripts;
        auto introns = transcriptIntrons;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        // Alignment at [500, 520] - no intron overlap
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 500;
        aln.exons[0][EX_L] = 21;
        
        std::set<uint32_t> candidates = {0, 1};
        std::set<uint32_t> outIntronic;
        bool result = compat.compatIsIntronic(aln, candidates, outIntronic);
        check(!result, "no intron overlap: returns false");
        check(outIntronic.empty(), "no intron overlap: no intronic genes");
    }
    
    // Test gene 1 with only 1 transcript (can never be intronic)
    {
        SlamCompatConfig cfg;
        cfg.intronic = true;
        auto g2t = geneToTranscripts;
        auto introns = transcriptIntrons;
        SlamCompat compat(cfg, std::move(g2t), std::move(introns));
        
        // Alignment at [2010, 2030] intersects tr2 intron [2000, 2100]
        // But gene 1 has only 1 transcript, so count = 1, not > 1
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 2010;
        aln.exons[0][EX_L] = 21;
        
        std::set<uint32_t> candidates = {1};
        std::set<uint32_t> outIntronic;
        bool result = compat.compatIsIntronic(aln, candidates, outIntronic);
        check(!result, "single-transcript gene: returns false");
        check(outIntronic.empty(), "single-transcript gene: no intronic genes");
    }
}

// ============================================================================
// Test: computeExonOverlap
// ============================================================================
void test_computeExonOverlap() {
    std::cout << "\n=== Testing computeExonOverlap ===\n";
    
    SlamCompatConfig cfg;
    std::vector<std::vector<uint32_t>> g2t;
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> introns;
    SlamCompat compat(cfg, std::move(g2t), std::move(introns));
    
    // Transcript exons (relative to trStart=1000):
    // Exon 0: [0, 49] -> genomic [1000, 1049]
    // Exon 1: [100, 149] -> genomic [1100, 1149]
    // Exon 2: [200, 299] -> genomic [1200, 1299]
    uint32_t trExSE[] = {0, 49, 100, 149, 200, 299};
    uint64_t trStart = 1000;
    uint16_t trExN = 3;
    
    // Test single-block alignment fully inside one exon
    {
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 1010;  // Start in exon 0
        aln.exons[0][EX_L] = 30;    // [1010, 1039]
        
        uint32_t overlap = compat.computeExonOverlap(aln, trStart, trExN, trExSE);
        check(overlap == 30, "fully inside exon: overlap = 30");
    }
    
    // Test single-block alignment spanning exon/intron boundary
    {
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 1040;  // Start in exon 0
        aln.exons[0][EX_L] = 30;    // [1040, 1069] - 10bp in exon, 20bp in intron
        
        uint32_t overlap = compat.computeExonOverlap(aln, trStart, trExN, trExSE);
        check(overlap == 10, "spans boundary: overlap = 10 (only exonic part)");
    }
    
    // Test single-block alignment fully in intron
    {
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 1060;  // In intron between exon 0 and 1
        aln.exons[0][EX_L] = 30;    // [1060, 1089]
        
        uint32_t overlap = compat.computeExonOverlap(aln, trStart, trExN, trExSE);
        check(overlap == 0, "fully in intron: overlap = 0");
    }
    
    // Test multi-block (spliced) alignment spanning two exons
    {
        Transcript aln;
        aln.nExons = 2;
        aln.exons[0][EX_G] = 1020;  // In exon 0: [1020, 1049] = 30bp
        aln.exons[0][EX_L] = 30;
        aln.exons[1][EX_G] = 1100;  // In exon 1: [1100, 1129] = 30bp
        aln.exons[1][EX_L] = 30;
        
        uint32_t overlap = compat.computeExonOverlap(aln, trStart, trExN, trExSE);
        check(overlap == 60, "spliced alignment: overlap = 60 (30 + 30)");
    }
    
    // Test alignment that partially overlaps multiple exons
    {
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 1130;  // Spans exon 1 end and into intron
        aln.exons[0][EX_L] = 100;   // [1130, 1229] - 20bp in exon 1, 51bp intron, 30bp in exon 2
        
        uint32_t overlap = compat.computeExonOverlap(aln, trStart, trExN, trExSE);
        check(overlap == 50, "partial multi-exon: overlap = 50 (20 + 30)");
    }
    
    // Test alignment outside all exons
    {
        Transcript aln;
        aln.nExons = 1;
        aln.exons[0][EX_G] = 2000;  // Far beyond transcript
        aln.exons[0][EX_L] = 50;
        
        uint32_t overlap = compat.computeExonOverlap(aln, trStart, trExN, trExSE);
        check(overlap == 0, "outside transcript: overlap = 0");
    }
}

// ============================================================================
// Main
// ============================================================================
int main() {
    std::cout << "SlamCompat Unit Tests\n";
    std::cout << "=====================\n";
    
    test_compatOverlapWeight();
    test_compatShouldCountPos();
    test_compatIsIntronic();
    test_computeExonOverlap();
    
    std::cout << "\n=====================\n";
    std::cout << "Results: " << g_passed << " passed, " << g_failed << " failed\n";
    
    if (g_failed == 0) {
        std::cout << "ALL TESTS PASSED\n";
        return 0;
    }
    return 1;
}
