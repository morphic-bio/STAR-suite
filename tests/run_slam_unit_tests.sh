#!/bin/bash
# Unified SLAM unit test runner
# Compiles and runs:
# - test_slam_snp_em
# - test_slam_vb_overdisp
# - test_qc_transition_orientation
# - test_slam_variance_mate_tc_rate
# - test_slam_cb_output

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SLAM_DIR="$ROOT_DIR/slam"
CORE_DIR="$ROOT_DIR/core/legacy"
TEST_DIR="$ROOT_DIR/tests/slam"
SOURCE_DIR="$SLAM_DIR/source"
HTSLIB_DIR="$CORE_DIR/source/htslib"

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
echo "[1/6] Building test_slam_snp_em..."
OUT_BIN_SNP_EM="$TMP_DIR/test_slam_snp_em"
if [[ ! -f "$TEST_DIR/test_slam_snp_em.cpp" ]]; then
    echo "FAIL: test_slam_snp_em.cpp not found at $TEST_DIR/test_slam_snp_em.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_slam_snp_em(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" -I"$SOURCE_DIR/libem" -I"$CORE_DIR/source" \
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

# Test 2: SLAM VB Overdisp Solver
echo "[2/6] Building test_slam_vb_overdisp..."
OUT_BIN_VB="$TMP_DIR/test_slam_vb_overdisp"
if [[ ! -f "$TEST_DIR/test_slam_vb_overdisp.cpp" ]]; then
    echo "FAIL: test_slam_vb_overdisp.cpp not found at $TEST_DIR/test_slam_vb_overdisp.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_slam_vb_overdisp(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" -I"$SOURCE_DIR/libem" -I"$CORE_DIR/source" \
        "$TEST_DIR/test_slam_vb_overdisp.cpp" "$SOURCE_DIR/libem/slam_vb_overdisp.cpp" \
        -o "$OUT_BIN_VB" -lm 2>&1; then
        echo "Running test_slam_vb_overdisp..."
        if "$OUT_BIN_VB" 2>&1; then
            echo "✓ test_slam_vb_overdisp PASSED"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo "✗ test_slam_vb_overdisp FAILED"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            FAILED_TESTS="$FAILED_TESTS test_slam_vb_overdisp(run)"
        fi
    else
        echo "FAIL: Compilation failed for test_slam_vb_overdisp"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        FAILED_TESTS="$FAILED_TESTS test_slam_vb_overdisp(compile)"
    fi
fi
echo

# Test 3: QC Transition Orientation
echo "[3/6] Building test_qc_transition_orientation..."
OUT_BIN_QC_ORIENT="$TMP_DIR/test_qc_transition_orientation"
if [[ ! -f "$TEST_DIR/test_qc_transition_orientation.cpp" ]]; then
    echo "FAIL: test_qc_transition_orientation.cpp not found at $TEST_DIR/test_qc_transition_orientation.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_qc_transition_orientation(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" -I"$CORE_DIR/source" \
        "$TEST_DIR/test_qc_transition_orientation.cpp" "$SOURCE_DIR/SlamQuant.cpp" \
        "$SOURCE_DIR/SlamSolver.cpp" "$SOURCE_DIR/libem/slam_vb_overdisp.cpp" \
        "$SOURCE_DIR/SlamVarianceAnalysis.cpp" "$SOURCE_DIR/SlamReadBuffer.cpp" "$SOURCE_DIR/SlamCompat.cpp" \
        "$SOURCE_DIR/SlamDump.cpp" \
        -L"$HTSLIB_DIR" -lhts -lz -lssl -lcrypto -lpthread \
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

# Test 4: SLAM QC Output
echo "[4/6] Building test_slam_qc_output..."
OUT_BIN_QC_OUT="$TMP_DIR/test_slam_qc_output"
if [[ ! -f "$TEST_DIR/test_slam_qc_output.cpp" ]]; then
    echo "FAIL: test_slam_qc_output.cpp not found at $TEST_DIR/test_slam_qc_output.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_slam_qc_output(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" -I"$CORE_DIR/source" \
        "$TEST_DIR/test_slam_qc_output.cpp" "$SOURCE_DIR/SlamQuant.cpp" "$SOURCE_DIR/SlamQcOutput.cpp" \
        "$SOURCE_DIR/SlamSolver.cpp" "$SOURCE_DIR/libem/slam_vb_overdisp.cpp" \
        "$SOURCE_DIR/SlamVarianceAnalysis.cpp" "$SOURCE_DIR/SlamReadBuffer.cpp" "$SOURCE_DIR/SlamCompat.cpp" \
        "$SOURCE_DIR/SlamDump.cpp" \
        -L"$HTSLIB_DIR" -lhts -lz -lssl -lcrypto -lpthread \
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

# Test 5: per-mate T→C error window (asymmetric trims)
echo "[5/6] Building test_slam_variance_mate_tc_rate..."
OUT_BIN_MATE="$TMP_DIR/test_slam_variance_mate_tc_rate"
if [[ ! -f "$TEST_DIR/test_slam_variance_mate_tc_rate.cpp" ]]; then
    echo "FAIL: test_slam_variance_mate_tc_rate.cpp not found at $TEST_DIR/test_slam_variance_mate_tc_rate.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_slam_variance_mate_tc_rate(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" \
        "$TEST_DIR/test_slam_variance_mate_tc_rate.cpp" "$SOURCE_DIR/SlamVarianceAnalysis.cpp" \
        -o "$OUT_BIN_MATE" -lm 2>&1; then
        echo "Running test_slam_variance_mate_tc_rate..."
        if "$OUT_BIN_MATE" 2>&1; then
            echo "✓ test_slam_variance_mate_tc_rate PASSED"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo "✗ test_slam_variance_mate_tc_rate FAILED"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            FAILED_TESTS="$FAILED_TESTS test_slam_variance_mate_tc_rate(run)"
        fi
    else
        echo "FAIL: Compilation failed for test_slam_variance_mate_tc_rate"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        FAILED_TESTS="$FAILED_TESTS test_slam_variance_mate_tc_rate(compile)"
    fi
fi
echo

# Test 6: SLAM count-binomial output
echo "[6/6] Building test_slam_cb_output..."
OUT_BIN_CB="$TMP_DIR/test_slam_cb_output"
if [[ ! -f "$TEST_DIR/test_slam_cb_output.cpp" ]]; then
    echo "FAIL: test_slam_cb_output.cpp not found at $TEST_DIR/test_slam_cb_output.cpp"
    TESTS_FAILED=$((TESTS_FAILED + 1))
    FAILED_TESTS="$FAILED_TESTS test_slam_cb_output(missing)"
else
    if "$CXX" $CXXFLAGS -I"$SOURCE_DIR" -I"$SOURCE_DIR/libem" -I"$CORE_DIR/source" \
        "$TEST_DIR/test_slam_cb_output.cpp" "$SOURCE_DIR/SlamQuant.cpp" \
        "$SOURCE_DIR/SlamSolver.cpp" "$SOURCE_DIR/libem/slam_vb_overdisp.cpp" \
        "$SOURCE_DIR/SlamVarianceAnalysis.cpp" "$SOURCE_DIR/SlamReadBuffer.cpp" "$SOURCE_DIR/SlamCompat.cpp" \
        "$SOURCE_DIR/SlamDump.cpp" \
        -L"$HTSLIB_DIR" -lhts -lz -lssl -lcrypto -lpthread \
        -o "$OUT_BIN_CB" 2>&1; then
        echo "Running test_slam_cb_output..."
        if "$OUT_BIN_CB" 2>&1; then
            echo "✓ test_slam_cb_output PASSED"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo "✗ test_slam_cb_output FAILED"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            FAILED_TESTS="$FAILED_TESTS test_slam_cb_output(run)"
        fi
    else
        echo "FAIL: Compilation failed for test_slam_cb_output"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        FAILED_TESTS="$FAILED_TESTS test_slam_cb_output(compile)"
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
