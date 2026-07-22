/**
 * Smoke test for SLAM QC JSON output (comprehensive report + legacy variance JSON).
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
#include "SlamVarianceAnalysis.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <unistd.h>

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
    char tmpdir[] = "/tmp/slam_qc_test_XXXXXX";
    if (mkdtemp(tmpdir) == nullptr) {
        std::cerr << "FAIL: mkdtemp\n";
        return 1;
    }
    const std::string outDir(tmpdir);
    const std::string compPath = outDir + "/comprehensive.json";
    const std::string curveOnlyPath = outDir + "/curve_only_comprehensive.json";
    const std::string perMatePath = outDir + "/per_mate.json";

    SlamQuant slam(1, false);
    slam.enableVarianceAnalysis(10, 1, 5, 3, 15, true);
    slam.recordVarianceRead();
    slam.recordVariancePosition(0, 30, true, true);
    slam.recordVariancePosition(0, 30, true, false);
    slam.recordVariancePosition(4, 0, 31, true, true);
    slam.recordVariancePosition(6, 1, 29, true, false);

    slam.addTransitionBase(SlamMismatchCategory::ExonicSense, 0, false, false, false, 3, 1, 1.0);
    slam.addTransitionBase(SlamMismatchCategory::ExonicSense, 0, false, false, false, 0, 2, 1.0);

    bool ok = writeSlamQcComprehensiveJson(slam, compPath, 2, 3, 7, 11, nullptr, nullptr, nullptr,
                                           10, 10);
    failed += check(ok, "writeSlamQcComprehensiveJson returns true");

    failed += check(fileContains(compPath, "\"type\": \"comprehensive_qc\""),
                    "comprehensive json contains type=comprehensive_qc");
    failed += check(fileContains(compPath, "\"position\": 1"),
                    "comprehensive json contains position 1 (1-based)");
    failed += check(fileContains(compPath, "\"star_tc_rate\""),
                    "comprehensive json contains star_tc_rate");
    failed += check(fileContains(compPath, "per_mate_separate"),
                    "comprehensive json reports variance_histogram_mode per_mate_separate");
    failed += check(fileContains(compPath, "\"trim5p_mate2\": 7"),
                    "comprehensive json emits trim5p_mate2");
    failed += check(fileContains(compPath, "\"trim3p_mate2\": 11"),
                    "comprehensive json emits trim3p_mate2");
    failed += check(fileContains(compPath, "\"mate_index\": 0"),
                    "comprehensive json contains mates[0]");
    failed += check(fileContains(compPath, "\"mate_index\": 1"),
                    "comprehensive json contains mates[1]");

    SlamQuant curveOnly(1, false);
    std::vector<double> curve0;
    curve0.push_back(0.10);
    curve0.push_back(std::numeric_limits<double>::quiet_NaN());
    curve0.push_back(0.30);
    std::vector<double> curve1;
    curve1.push_back(0.20);
    curve1.push_back(0.40);

    ok = writeSlamQcComprehensiveJson(curveOnly, curveOnlyPath, 1, 2, 3, 4, nullptr,
                                       &curve0, &curve1, 3, 2);
    failed += check(ok, "writeSlamQcComprehensiveJson curve-only returns true");
    failed += check(fileContains(curveOnlyPath, "per_mate_separate"),
                    "curve-only comprehensive json reports per-mate mode");
    failed += check(fileContains(curveOnlyPath, "\"stddev_tc_rate\": 0.100000"),
                    "curve-only comprehensive json emits mate1 curve samples");
    failed += check(fileContains(curveOnlyPath, "\"stddev_tc_rate\": 0.400000"),
                    "curve-only comprehensive json emits mate2 curve samples");
    failed += check(fileContains(curveOnlyPath, "\"position\": 5"),
                    "curve-only comprehensive json emits concatenated mate2 position");
    failed += check(fileContains(curveOnlyPath, "\"trim5p_mate2\": 3"),
                    "curve-only comprehensive json emits applied mate2 trim");

    SlamVarianceAnalyzer an;
    an.setSeparateMateHistograms(true);
    an.recordRead();
    an.recordPositionMate(5, 0, 30, true, false);
    an.recordPositionMate(12, 1, 28, true, true);

    ok = writeSlamQcJson(an, perMatePath, 0, "per_mate", 1, 2, 4, 8, 1u, nullptr, "", 0.0, 0.0, "");
    failed += check(ok, "writeSlamQcJson returns true");

    failed += check(fileContains(perMatePath, "per_mate_separate"),
                    "legacy qc json reports variance_histogram_mode");
    failed += check(fileContains(perMatePath, "\"trim5p_mate2\": 4"),
                    "legacy qc json emits trim5p_mate2");
    failed += check(fileContains(perMatePath, "\"trim3p_mate2\": 8"),
                    "legacy qc json emits trim3p_mate2");
    failed += check(fileContains(perMatePath, "\"mates\": ["),
                    "legacy qc json contains mates array");
    failed += check(fileContains(perMatePath, "\"mate_index\": 0"),
                    "legacy qc json contains mate_index 0");
    failed += check(fileContains(perMatePath, "\"mate_index\": 1"),
                    "legacy qc json contains mate_index 1");
    failed += check(fileContains(perMatePath, "\"position\": 5"),
                    "legacy qc json has mate-local position 5 (mate 1)");
    failed += check(fileContains(perMatePath, "\"position\": 12"),
                    "legacy qc json has mate-local position 12 (mate 2)");

    if (failed) {
        std::cerr << failed << " TEST(S) FAILED\n";
    } else {
        std::cerr << "ALL TESTS PASSED\n";
    }
    {
        const std::string cmd = "rm -rf " + outDir;
        if (std::system(cmd.c_str()) != 0) {
            std::cerr << "WARNING: could not remove temp dir " << outDir << "\n";
        }
    }
    return failed ? 1 : 0;
}
