/**
 * Unit tests for SNP mismatch threshold auto-estimation.
 *
 * Tests the knee/elbow detection algorithm used to automatically
 * determine the SNP mismatch fraction threshold.
 *
 * Compile:
 *   g++ -std=c++17 -I../../source -o test_snp_threshold \
 *     test_snp_threshold.cpp ../../source/SlamQuant.cpp ../../source/SlamSolver.cpp
 *
 * Run:
 *   ./test_snp_threshold
 *
 * Or use the helper script:
 *   tests/run_snp_threshold_test.sh
 */

#include "SlamQuant.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

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

static bool approxEqual(double a, double b, double tol = 0.05) {
    return std::fabs(a - b) < tol;
}

// ============================================================================
// Helper: simulate SNP observations with a given mismatch rate
// ============================================================================
static void simulatePositions(SlamQuant& sq, uint64_t startPos, uint64_t count,
                               uint32_t coverage, double mismatchRate,
                               std::mt19937& rng) {
    std::binomial_distribution<uint32_t> misDist(coverage, mismatchRate);
    for (uint64_t i = 0; i < count; ++i) {
        uint64_t pos = startPos + i;
        uint32_t mismatches = misDist(rng);
        // Record observations: coverage times, with mismatches times as mismatch
        for (uint32_t j = 0; j < coverage; ++j) {
            sq.recordSnpObservation(pos, j < mismatches);
        }
    }
}

// ============================================================================
// Test: Bimodal distribution (low background + high SNPs)
// Expected: knee around 0.2-0.3
// ============================================================================
void test_bimodal_distribution() {
    std::cout << "\n=== Testing bimodal distribution ===\n";
    
    // Create SlamQuant with auto-estimation (snpMismatchFrac <= 0)
    SlamQuant sq(10, true, -1.0);
    
    std::mt19937 rng(42);
    
    // Simulate ~8000 positions with low mismatch rate (~0.02)
    // These represent normal sequencing errors
    simulatePositions(sq, 0, 8000, 20, 0.02, rng);
    
    // Simulate ~2000 positions with high mismatch rate (~0.45)
    // These represent true SNPs
    simulatePositions(sq, 100000, 2000, 20, 0.45, rng);
    
    // Buffer some reads to trigger finalization
    for (uint32_t i = 0; i < 100; ++i) {
        std::vector<uint32_t> mismatchPos = {i};
        sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
    }
    
    SlamSnpBufferStats stats;
    sq.finalizeSnpMask(&stats);
    
    check(stats.mismatchFracMode == "auto", "bimodal: mode is auto");
    // With data at 0.02 and 0.45, knee detection finds threshold around 0.4-0.5
    // This correctly separates SNPs (0.45+) from noise (0.02)
    check(stats.mismatchFracUsed >= 0.10 && stats.mismatchFracUsed <= 0.60,
          "bimodal: threshold in valid range [0.10, 0.60]");
    check(stats.mismatchFracAuto > 0.0, "bimodal: auto value computed");
    check(stats.kneeBin > 0, "bimodal: knee bin detected");
    check(stats.eligibleSites >= 10000, "bimodal: eligible sites counted");
    
    std::cout << "    threshold=" << stats.mismatchFracUsed
              << " auto=" << stats.mismatchFracAuto
              << " knee_bin=" << stats.kneeBin
              << " eligible=" << stats.eligibleSites << "\n";
}

// ============================================================================
// Test: Insufficient data (fallback expected)
// ============================================================================
void test_insufficient_data() {
    std::cout << "\n=== Testing insufficient data ===\n";
    
    SlamQuant sq(10, true, -1.0);
    
    std::mt19937 rng(123);
    
    // Only 500 positions - below the 1000 threshold
    simulatePositions(sq, 0, 500, 15, 0.03, rng);
    
    // Buffer some reads
    for (uint32_t i = 0; i < 50; ++i) {
        std::vector<uint32_t> mismatchPos = {i};
        sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
    }
    
    SlamSnpBufferStats stats;
    sq.finalizeSnpMask(&stats);
    
    check(stats.mismatchFracMode == "auto_fallback", "insufficient: mode is auto_fallback");
    check(approxEqual(stats.mismatchFracUsed, 0.22, 0.01), "insufficient: fallback to 0.22");
    check(stats.mismatchFracAuto == 0.0, "insufficient: auto value is 0");
    check(stats.eligibleSites < 1000, "insufficient: eligible sites < 1000");
    
    std::cout << "    threshold=" << stats.mismatchFracUsed
              << " mode=" << stats.mismatchFracMode
              << " eligible=" << stats.eligibleSites << "\n";
}

// ============================================================================
// Test: Unimodal low distribution (weak knee -> fallback)
// ============================================================================
void test_unimodal_low() {
    std::cout << "\n=== Testing unimodal low distribution ===\n";
    
    SlamQuant sq(10, true, -1.0);
    
    std::mt19937 rng(456);
    
    // All positions have low mismatch rate - no clear knee expected
    simulatePositions(sq, 0, 5000, 20, 0.03, rng);
    
    // Buffer reads
    for (uint32_t i = 0; i < 50; ++i) {
        std::vector<uint32_t> mismatchPos = {i};
        sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
    }
    
    SlamSnpBufferStats stats;
    sq.finalizeSnpMask(&stats);
    
    // With unimodal low data, knee detection may succeed or fallback
    // Either way, the threshold should be reasonable
    check(stats.mismatchFracUsed >= 0.10 && stats.mismatchFracUsed <= 0.60,
          "unimodal_low: threshold in valid range");
    check(stats.mismatchFracMode == "auto" || stats.mismatchFracMode == "auto_fallback",
          "unimodal_low: mode is auto or auto_fallback");
    
    std::cout << "    threshold=" << stats.mismatchFracUsed
              << " mode=" << stats.mismatchFracMode
              << " eligible=" << stats.eligibleSites << "\n";
}

// ============================================================================
// Test: Unimodal high distribution (weak knee -> fallback)
// ============================================================================
void test_unimodal_high() {
    std::cout << "\n=== Testing unimodal high distribution ===\n";
    
    SlamQuant sq(10, true, -1.0);
    
    std::mt19937 rng(789);
    
    // All positions have high mismatch rate - extreme case
    simulatePositions(sq, 0, 3000, 20, 0.70, rng);
    
    // Buffer reads
    for (uint32_t i = 0; i < 50; ++i) {
        std::vector<uint32_t> mismatchPos = {i};
        sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
    }
    
    SlamSnpBufferStats stats;
    sq.finalizeSnpMask(&stats);
    
    // With all-high data, threshold should be clamped or fallback
    check(stats.mismatchFracUsed >= 0.10 && stats.mismatchFracUsed <= 0.60,
          "unimodal_high: threshold in valid range");
    
    std::cout << "    threshold=" << stats.mismatchFracUsed
              << " mode=" << stats.mismatchFracMode
              << " eligible=" << stats.eligibleSites << "\n";
}

// ============================================================================
// Test: Explicit threshold override
// ============================================================================
void test_explicit_threshold() {
    std::cout << "\n=== Testing explicit threshold ===\n";
    
    // Create with explicit threshold
    SlamQuant sq(10, true, 0.30);
    
    std::mt19937 rng(42);
    
    // Simulate positions (would normally trigger auto-estimation)
    simulatePositions(sq, 0, 5000, 20, 0.02, rng);
    simulatePositions(sq, 100000, 1000, 20, 0.50, rng);
    
    // Buffer reads
    for (uint32_t i = 0; i < 50; ++i) {
        std::vector<uint32_t> mismatchPos = {i};
        sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
    }
    
    SlamSnpBufferStats stats;
    sq.finalizeSnpMask(&stats);
    
    check(stats.mismatchFracMode == "explicit", "explicit: mode is explicit");
    check(approxEqual(stats.mismatchFracUsed, 0.30, 0.001), "explicit: threshold is 0.30");
    check(stats.mismatchFracAuto == 0.0, "explicit: auto value not computed");
    
    std::cout << "    threshold=" << stats.mismatchFracUsed
              << " mode=" << stats.mismatchFracMode << "\n";
}

// ============================================================================
// Test: Clamp behavior (knee below 0.10)
// ============================================================================
void test_clamp_low() {
    std::cout << "\n=== Testing clamp behavior (low) ===\n";
    
    SlamQuant sq(10, true, -1.0);
    
    std::mt19937 rng(111);
    
    // Create distribution with very low knee point
    // Many positions at very low rate, few at moderate rate
    simulatePositions(sq, 0, 9000, 25, 0.005, rng);  // Very low
    simulatePositions(sq, 100000, 1000, 25, 0.08, rng);  // Just above low
    
    // Buffer reads
    for (uint32_t i = 0; i < 50; ++i) {
        std::vector<uint32_t> mismatchPos = {i};
        sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
    }
    
    SlamSnpBufferStats stats;
    sq.finalizeSnpMask(&stats);
    
    // If knee is detected below 0.10, it should be clamped
    check(stats.mismatchFracUsed >= 0.10, "clamp_low: threshold >= 0.10");
    
    std::cout << "    threshold=" << stats.mismatchFracUsed
              << " mode=" << stats.mismatchFracMode
              << " knee_bin=" << stats.kneeBin << "\n";
}

// ============================================================================
// Test: Coverage filter (positions below kSnpMinCoverage ignored)
// ============================================================================
void test_coverage_filter() {
    std::cout << "\n=== Testing coverage filter ===\n";
    
    SlamQuant sq(10, true, -1.0);
    
    std::mt19937 rng(222);
    
    // Positions with coverage < 10 (should be ignored)
    for (uint64_t i = 0; i < 5000; ++i) {
        uint64_t pos = i;
        // Only 5 observations per position
        for (uint32_t j = 0; j < 5; ++j) {
            sq.recordSnpObservation(pos, j < 2);  // 40% mismatch
        }
    }
    
    // Positions with coverage >= 10 (should be counted)
    simulatePositions(sq, 100000, 2000, 15, 0.03, rng);
    
    // Buffer reads
    for (uint32_t i = 0; i < 50; ++i) {
        std::vector<uint32_t> mismatchPos = {i};
        sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
    }
    
    SlamSnpBufferStats stats;
    sq.finalizeSnpMask(&stats);
    
    // Only the 2000 high-coverage positions should count
    check(stats.eligibleSites >= 1500 && stats.eligibleSites <= 2500,
          "coverage_filter: only high-coverage positions counted");
    
    std::cout << "    eligible=" << stats.eligibleSites
              << " mode=" << stats.mismatchFracMode << "\n";
}

// ============================================================================
// Test: Determinism (same inputs yield same threshold)
// ============================================================================
void test_determinism() {
    std::cout << "\n=== Testing determinism ===\n";
    
    double threshold1, threshold2;
    
    // Run 1
    {
        SlamQuant sq(10, true, -1.0);
        std::mt19937 rng(42);
        simulatePositions(sq, 0, 5000, 20, 0.02, rng);
        simulatePositions(sq, 100000, 1000, 20, 0.45, rng);
        
        for (uint32_t i = 0; i < 50; ++i) {
            std::vector<uint32_t> mismatchPos = {i};
            sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
        }
        
        SlamSnpBufferStats stats;
        sq.finalizeSnpMask(&stats);
        threshold1 = stats.mismatchFracUsed;
    }
    
    // Run 2 (identical setup)
    {
        SlamQuant sq(10, true, -1.0);
        std::mt19937 rng(42);
        simulatePositions(sq, 0, 5000, 20, 0.02, rng);
        simulatePositions(sq, 100000, 1000, 20, 0.45, rng);
        
        for (uint32_t i = 0; i < 50; ++i) {
            std::vector<uint32_t> mismatchPos = {i};
            sq.bufferSnpRead(0, 40, mismatchPos, 1.0);
        }
        
        SlamSnpBufferStats stats;
        sq.finalizeSnpMask(&stats);
        threshold2 = stats.mismatchFracUsed;
    }
    
    check(threshold1 == threshold2, "determinism: same inputs yield same threshold");
    
    std::cout << "    run1=" << threshold1 << " run2=" << threshold2 << "\n";
}

// ============================================================================
// Test: SNP detection disabled (should not compute threshold)
// ============================================================================
void test_snp_detect_disabled() {
    std::cout << "\n=== Testing SNP detection disabled ===\n";
    
    // snpDetect = false
    SlamQuant sq(10, false, -1.0);
    
    // Try to record observations (should be ignored)
    for (uint64_t i = 0; i < 100; ++i) {
        sq.recordSnpObservation(i, true);
    }
    
    SlamSnpBufferStats stats;
    sq.finalizeSnpMask(&stats);
    
    // With SNP detection disabled, no threshold estimation
    check(stats.mismatchFracUsed == 0.0, "disabled: no threshold computed");
    check(stats.mismatchFracMode.empty(), "disabled: mode is empty");
    check(stats.eligibleSites == 0, "disabled: no eligible sites");
    
    std::cout << "    threshold=" << stats.mismatchFracUsed
              << " mode=" << stats.mismatchFracMode << "\n";
}

// ============================================================================
// Main
// ============================================================================
int main() {
    std::cout << "SNP Threshold Auto-Estimation Unit Tests\n";
    std::cout << "=========================================\n";
    
    test_bimodal_distribution();
    test_insufficient_data();
    test_unimodal_low();
    test_unimodal_high();
    test_explicit_threshold();
    test_clamp_low();
    test_coverage_filter();
    test_determinism();
    test_snp_detect_disabled();
    
    std::cout << "\n=========================================\n";
    std::cout << "Results: " << g_passed << " passed, " << g_failed << " failed\n";
    
    if (g_failed == 0) {
        std::cout << "ALL TESTS PASSED\n";
        return 0;
    }
    return 1;
}

