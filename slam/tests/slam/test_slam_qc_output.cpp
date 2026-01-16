/**
 * Smoke test for SLAM QC JSON output (comprehensive report).
 *
 * Compile:
 *   g++ -std=c++11 -O2 -I../../source \
 *     test_slam_qc_output.cpp ../../source/SlamQuant.cpp ../../source/SlamQcOutput.cpp \
 *     -o test_slam_qc_output
 *
 * Run:
 *   ./test_slam_qc_output
 */

#include "SlamQuant.h"
#include "SlamQcOutput.h"

#include <fstream>
#include <iostream>
#include <string>

static int check(bool ok, const std::string& label) {
    if (!ok) {
        std::cerr << "FAIL: " << label << "\n";
        return 1;
    }
    std::cerr << "PASS: " << label << "\n";
    return 0;
}

static bool fileContains(const std::string& path, const std::string& needle) {
    std::ifstream in(path.c_str());
    if (!in.good()) {
        return false;
    }
    std::string content((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
    return content.find(needle) != std::string::npos;
}

int main() {
    int failed = 0;
    const std::string outPath = "test_slam_qc_output.json";

    SlamQuant slam(1, false);
    slam.enableVarianceAnalysis(10, 1);
    slam.recordVarianceRead();
    slam.recordVariancePosition(0, 30, true, true);
    slam.recordVariancePosition(0, 30, true, false);

    // Add a couple of transitions so the JSON contains rates
    slam.addTransitionBase(SlamMismatchCategory::ExonicSense, 0, false, false, false, 3, 1, 1.0);
    slam.addTransitionBase(SlamMismatchCategory::ExonicSense, 0, false, false, false, 0, 2, 1.0);

    bool ok = writeSlamQcComprehensiveJson(slam, outPath, 2, 3, nullptr);
    failed += check(ok, "writeSlamQcComprehensiveJson returns true");

    failed += check(fileContains(outPath, "\"type\": \"comprehensive_qc\""),
                    "json contains type=comprehensive_qc");
    failed += check(fileContains(outPath, "\"position\": 1"),
                    "json contains position 1 (1-based)");
    failed += check(fileContains(outPath, "\"star_tc_rate\""),
                    "json contains star_tc_rate");
    failed += check(fileContains(outPath, "\"stddev_tc_rate\""),
                    "json contains stddev_tc_rate");

    if (failed) {
        std::cerr << failed << " TEST(S) FAILED\n";
        return 1;
    }
    std::cerr << "ALL TESTS PASSED\n";
    return 0;
}
