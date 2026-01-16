#!/bin/bash
# Unified SLAM unit test runner
# Compiles and runs:
# - test_slam_snp_em
# - test_qc_transition_orientation
# - test_slam_qc_output

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TEST_DIR="$ROOT_DIR/tests/slam"
SOURCE_DIR="$ROOT_DIR/source"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++11 -O2}"
TMP_DIR="${TMP_DIR:-/tmp}"

echo "======================================================================"
echo "SLAM Unit Test Runner"
echo "======================================================================"
echo "Script dir: $SCRIPT_DIR"
echo "Root dir: $ROOT_DIR"
echo "Test dir: $TEST_DIR"
echo "Source dir: $SOURCE_DIR"
echo "Compiler: $CXX"
echo "Flags: $CXXFLAGS"
echo "======================================================================"
echo

# Track results
TESTS_PASSED=0
TESTS_FAILED=0
FAILED_TESTS=""

# Test 1: SLAM SNP EM
echo "[1/3] Building test_slam_snp_em..."
OUT_BIN_SNP_EM="$TMP_DIR/test_slam_snp_em"
if [[ ! -f "$TEST_DIR/test_slam_snp_em.cpp" ]]; then
    echo "FAIL: test_slam_snp_em.cpp not found at $TEST_DIR/test_slam_snp_em.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_slam_snp_em(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" -I"$SOURCE_DIR/libem" \
        "$TEST_DIR/test_slam_snp_em.cpp" "$SOURCE_DIR/libem/slam_snp_em.cpp" \
        -o "$OUT_BIN_SNP_EM" -lm 2>&1; then
        echo "Running test_slam_snp_em..."
        if "$OUT_BIN_SNP_EM" 2>&1; then
            echo "✓ test_slam_snp_em PASSED"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo "✗ test_slam_snp_em FAILED"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            FAILED_TESTS="$FAILED_TESTS test_slam_snp_em(run)"
        fi
    else
        echo "FAIL: Compilation failed for test_slam_snp_em"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        FAILED_TESTS="$FAILED_TESTS test_slam_snp_em(compile)"
    fi
fi
echo

# Test 2: QC Transition Orientation
echo "[2/3] Building test_qc_transition_orientation..."
OUT_BIN_QC_ORIENT="$TMP_DIR/test_qc_transition_orientation"
if [[ ! -f "$TEST_DIR/test_qc_transition_orientation.cpp" ]]; then
    echo "FAIL: test_qc_transition_orientation.cpp not found at $TEST_DIR/test_qc_transition_orientation.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_qc_transition_orientation(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" \
        "$TEST_DIR/test_qc_transition_orientation.cpp" "$SOURCE_DIR/SlamQuant.cpp" \
        "$SOURCE_DIR/SlamSolver.cpp" "$SOURCE_DIR/SlamVarianceAnalysis.cpp" "$SOURCE_DIR/SlamReadBuffer.cpp" "$SOURCE_DIR/SlamCompat.cpp" \
        -o "$OUT_BIN_QC_ORIENT" 2>&1; then
        echo "Running test_qc_transition_orientation..."
        if "$OUT_BIN_QC_ORIENT" 2>&1; then
            echo "✓ test_qc_transition_orientation PASSED"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo "✗ test_qc_transition_orientation FAILED"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            FAILED_TESTS="$FAILED_TESTS test_qc_transition_orientation(run)"
        fi
    else
        echo "FAIL: Compilation failed for test_qc_transition_orientation"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        FAILED_TESTS="$FAILED_TESTS test_qc_transition_orientation(compile)"
    fi
fi
echo

# Test 3: SLAM QC Output
echo "[3/3] Building test_slam_qc_output..."
OUT_BIN_QC_OUT="$TMP_DIR/test_slam_qc_output"
if [[ ! -f "$TEST_DIR/test_slam_qc_output.cpp" ]]; then
    echo "FAIL: test_slam_qc_output.cpp not found at $TEST_DIR/test_slam_qc_output.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_slam_qc_output(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" \
        "$TEST_DIR/test_slam_qc_output.cpp" "$SOURCE_DIR/SlamQuant.cpp" "$SOURCE_DIR/SlamQcOutput.cpp" \
        "$SOURCE_DIR/SlamSolver.cpp" "$SOURCE_DIR/SlamVarianceAnalysis.cpp" "$SOURCE_DIR/SlamReadBuffer.cpp" "$SOURCE_DIR/SlamCompat.cpp" \
        -o "$OUT_BIN_QC_OUT" 2>&1; then
        echo "Running test_slam_qc_output..."
        if "$OUT_BIN_QC_OUT" 2>&1; then
            echo "✓ test_slam_qc_output PASSED"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo "✗ test_slam_qc_output FAILED"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            FAILED_TESTS="$FAILED_TESTS test_slam_qc_output(run)"
        fi
    else
        echo "FAIL: Compilation failed for test_slam_qc_output"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        FAILED_TESTS="$FAILED_TESTS test_slam_qc_output(compile)"
    fi
fi
echo

# Summary
echo "======================================================================"
echo "Test Results Summary"
echo "======================================================================"
echo "Passed: $TESTS_PASSED"
echo "Failed: $TESTS_FAILED"
if [[ $TESTS_FAILED -gt 0 ]]; then
    echo "Failed tests:$FAILED_TESTS"
    exit 1
else
    echo "All tests passed!"
    exit 0
fi
