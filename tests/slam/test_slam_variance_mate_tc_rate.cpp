/**
 * Unit test: asymmetric per-mate T→C error window (computeGlobalTcErrorRatePerMate).
 *
 * Compile via tests/run_slam_unit_tests.sh
 */

#include "SlamVarianceAnalysis.h"

#include <cmath>
#include <iostream>

static int check(bool ok, const char* label) {
    if (!ok) {
        std::cerr << "FAIL: " << label << "\n";
        return 1;
    }
    std::cerr << "PASS: " << label << "\n";
    return 0;
}

int main() {
    int failed = 0;

    SlamVarianceAnalyzer an(1000, 1, 5, 3, 15);
    an.setSeparateMateHistograms(true);

    // One analyzed read; mate 0 and mate 1 each have one T at local position 3.
    an.recordRead();
    an.recordPositionMate(3, 0, 30, true, false);  // 1 T, 0 TC at pos 3 mate1
    an.recordPositionMate(3, 1, 30, true, true);   // 1 T, 1 TC at pos 3 mate2

    constexpr uint32_t mateLen = 10;
    int trim5p[2] = {4, 0};
    int trim3p[2] = {0, 0};

    auto tup = an.computeGlobalTcErrorRatePerMate(trim5p, trim3p, mateLen, mateLen);
    uint64_t t_total = std::get<0>(tup);
    uint64_t tc_total = std::get<1>(tup);
    double p_est = std::get<2>(tup);

    // Mate1 pos 3 excluded (3 < trim5p 4); mate2 pos 3 included → t=1, tc=1 → p=1
    failed += check(t_total == 1u && tc_total == 1u, "trim excludes mate1 only: t=1 tc=1");
    failed += check(std::abs(p_est - 1.0) < 1e-9, "p_est == 1 when only mate2 contributes");

    trim5p[0] = 0;
    trim5p[1] = 4;

    tup = an.computeGlobalTcErrorRatePerMate(trim5p, trim3p, mateLen, mateLen);
    t_total = std::get<0>(tup);
    tc_total = std::get<1>(tup);
    p_est = std::get<2>(tup);

    failed += check(t_total == 1u && tc_total == 0u, "trim excludes mate2 only: t=1 tc=0");
    failed += check(std::abs(p_est - 0.0) < 1e-9, "p_est == 0 when only mate1 contributes");

    if (failed) {
        std::cerr << failed << " CHECKS FAILED\n";
        return 1;
    }
    std::cerr << "ALL TESTS PASSED\n";
    return 0;
}
