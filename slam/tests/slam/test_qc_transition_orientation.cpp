/**
 * Unit test for transcript-oriented T->C aggregation in QC transition data.
 *
 * Compile:
 *   g++ -std=c++11 -O2 -I../../source \
 *     test_qc_transition_orientation.cpp ../../source/SlamQuant.cpp \
 *     -o test_qc_transition_orientation
 *
 * Run:
 *   ./test_qc_transition_orientation
 */

#include "SlamQuant.h"
#include <cmath>
#include <iostream>

static int check(bool ok, const std::string& label) {
    if (!ok) {
        std::cerr << "FAIL: " << label << "\n";
        return 1;
    }
    std::cerr << "PASS: " << label << "\n";
    return 0;
}

static bool approx(double a, double b, double tol = 1e-6) {
    return std::fabs(a - b) <= tol;
}

int main() {
    int failed = 0;
    SlamQuant slam(1, false);

    // ExonicSense position 0:
    // genomic T->C (plus strand) weight 10
    slam.addTransitionBase(SlamMismatchCategory::ExonicSense, 0, false, false, false, 3, 1, 10.0);
    // genomic A->G (minus strand transcript T->C) weight 5
    slam.addTransitionBase(SlamMismatchCategory::ExonicSense, 0, false, false, false, 0, 2, 5.0);
    // genomic T->A (control) weight 2
    slam.addTransitionBase(SlamMismatchCategory::ExonicSense, 0, false, false, false, 3, 0, 2.0);
    // genomic A->T (control on minus strand) weight 3
    slam.addTransitionBase(SlamMismatchCategory::ExonicSense, 0, false, false, false, 0, 3, 3.0);

    auto data = slam.getPositionTransitionData();
    failed += check(data.size() == 1, "transition data size == 1");
    auto it = data.find(0);
    failed += check(it != data.end(), "position 0 present");
    if (it != data.end()) {
        double cov_tc = std::get<0>(it->second);
        double mm_tc = std::get<1>(it->second);
        double cov_ta = std::get<2>(it->second);
        double mm_ta = std::get<3>(it->second);

        // Coverage is total transcript-T coverage (genomic T + genomic A).
        failed += check(approx(cov_tc, 20.0), "tc_cov == 20");
        failed += check(approx(cov_ta, 20.0), "ta_cov == 20");
        // Mismatch counts include T->C and A->G; T->A and A->T for control.
        failed += check(approx(mm_tc, 15.0), "tc_mm == 15");
        failed += check(approx(mm_ta, 5.0), "ta_mm == 5");
    }

    if (failed) {
        std::cerr << failed << " TEST(S) FAILED\n";
        return 1;
    }
    std::cerr << "ALL TESTS PASSED\n";
    return 0;
}
