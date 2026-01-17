/**
 * Basic unit checks for SlamSolver behavior.
 *
 * Compile (once SlamSolver exists):
 *   g++ -std=c++11 -I../source -o test_slam_solver \
 *     tests/slam/test_slam_solver.cpp ../source/SlamSolver.cpp
 *
 * Run:
 *   ./test_slam_solver
 *
 * Or use the helper script:
 *   tests/run_slam_solver_test.sh
 */

#include "SlamSolver.h"

#include <cmath>
#include <iostream>

static bool isFinite(double v) {
    return std::isfinite(v);
}

static int check(bool ok, const std::string& label) {
    if (!ok) {
        std::cerr << "FAIL: " << label << "\n";
        return 1;
    }
    return 0;
}

int main() {
    int failed = 0;

    SlamSolver solver(0.001, 0.05);

    // Low-conversion signal
    MismatchHistogram low;
    low[(40 << 8) | 0] = 500;
    low[(40 << 8) | 1] = 5;

    SlamResult low_res = solver.solve(low);
    failed += check(low_res.converged, "low: converged");
    failed += check(low_res.ntr >= 0.0 && low_res.ntr <= 1.0, "low: ntr range");
    failed += check(isFinite(low_res.ntr) && isFinite(low_res.sigma) && isFinite(low_res.log_likelihood),
                    "low: finite outputs");

    // Higher-conversion signal
    MismatchHistogram high;
    high[(40 << 8) | 0] = 100;
    high[(40 << 8) | 3] = 100;

    SlamResult high_res = solver.solve(high);
    failed += check(high_res.converged, "high: converged");
    failed += check(high_res.ntr >= 0.0 && high_res.ntr <= 1.0, "high: ntr range");
    failed += check(high_res.ntr >= low_res.ntr, "ordering: high >= low");

    if (failed == 0) {
        std::cout << "PASS: SlamSolver basic checks\n";
        return 0;
    }
    return 1;
}
