#!/usr/bin/env bash
# Parameters Regression Test Runner
# Implements the test plan from plans/parms_tests_plan.md

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
STAR_SOURCE_DIR="${ROOT_DIR}/core/legacy/source"
RELEASE_BIN="${STAR_SOURCE_DIR}/STAR.release"
LSAN_SUPPRESSIONS="${ROOT_DIR}/tests/compat/lsan_suppressions.txt"
NORMALIZE_SCRIPT="${ROOT_DIR}/tests/compat/normalize_lsan_log.py"

# Source setup script to get environment variables
source "${SCRIPT_DIR}/setup_parms_tests.sh"

# Ensure STAR_BIN is set (setup script should export it, but set default if not)
export STAR_BIN="${STAR_BIN:-${STAR_SOURCE_DIR}/STAR}"

# Defaults
TEST_SET="standard"
BUILD_MODE="release"
KEEP_GOING=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --set)
            TEST_SET="$2"
            shift 2
            ;;
        --mode)
            BUILD_MODE="$2"
            shift 2
            ;;
        --keep-going)
            KEEP_GOING=true
            shift
            ;;
        -h|--help)
            cat << EOF
Usage: $0 [OPTIONS]

Run Parameters Regression Test Suite

OPTIONS:
    --set standard|extended    Test set to run (default: standard)
    --mode release|asan        Build mode (default: release)
    --keep-going               Continue on failure (default: stop on first failure)
    -h, --help                 Show this help message

NOTE: The setup script (setup_parms_tests.sh) is automatically sourced to set
      environment variables (STAR_BIN, REMOVE_Y_READS_BIN, TXIMPORT_BIN, etc.)
      that are exported to all test scripts.

EXAMPLES:
    $0                                    # Run standard set with release build
    $0 --set extended                    # Run extended set with release build
    $0 --mode asan                       # Run standard set with ASAN build
    $0 --set extended --mode asan        # Run extended set with ASAN build
    $0 --keep-going                      # Continue on failures
EOF
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            echo "Run '$0 --help' for usage information" >&2
            exit 1
            ;;
    esac
done

# Validate arguments
if [[ "$TEST_SET" != "standard" && "$TEST_SET" != "extended" ]]; then
    echo "ERROR: --set must be 'standard' or 'extended'" >&2
    exit 1
fi

if [[ "$BUILD_MODE" != "release" && "$BUILD_MODE" != "asan" ]]; then
    echo "ERROR: --mode must be 'release' or 'asan'" >&2
    exit 1
fi

# Create artifacts directory
TIMESTAMP=$(date +%Y%m%d)
ARTIFACTS_DIR="${ROOT_DIR}/plans/artifacts/parms_tests_${TIMESTAMP}"
mkdir -p "${ARTIFACTS_DIR}"
REPORT_FILE="${ARTIFACTS_DIR}/report.md"

# Initialize report
cat > "${REPORT_FILE}" << EOF
# Parameters Regression Test Report

Generated: $(date)

## Configuration
- Test Set: ${TEST_SET}
- Build Mode: ${BUILD_MODE}
- Keep Going: ${KEEP_GOING}
- STAR Binary: ${STAR_BIN}
- Artifacts Directory: ${ARTIFACTS_DIR}

## Test Results

EOF

# Track results
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0
SKIPPED_TESTS=0
FAILED_TEST_NAMES=()
SKIPPED_TEST_NAMES=()

# Function to log messages
log_info() {
    echo -e "${BLUE}INFO${NC}: $1" | tee -a "${REPORT_FILE}"
}

log_pass() {
    echo -e "${GREEN}PASS${NC}: $1" | tee -a "${REPORT_FILE}"
}

log_fail() {
    echo -e "${RED}FAIL${NC}: $1" | tee -a "${REPORT_FILE}"
}

log_warn() {
    echo -e "${YELLOW}WARN${NC}: $1" | tee -a "${REPORT_FILE}"
}

# Function to run a test script
run_test() {
    local test_script="$1"
    local test_name=$(basename "$test_script" .sh)
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    # Check if test script exists
    if [[ ! -f "$test_script" ]]; then
        log_warn "Test script not found: $test_script (skipping)"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        SKIPPED_TEST_NAMES+=("$test_name")
        echo "" >> "${REPORT_FILE}"
        echo "### $test_name" >> "${REPORT_FILE}"
        echo "- Status: SKIPPED (script not found)" >> "${REPORT_FILE}"
        echo "- Script: $test_script" >> "${REPORT_FILE}"
        echo "" >> "${REPORT_FILE}"
        return 0
    fi
    
    # Prefer executable scripts, but allow bash fallback.
    local runner=("$test_script")
    if [[ ! -x "$test_script" ]]; then
        log_warn "Test script not executable: $test_script (running via bash)"
        runner=(bash "$test_script")
    fi
    
    log_info "Running: $test_name"
    echo "" >> "${REPORT_FILE}"
    echo "### $test_name" >> "${REPORT_FILE}"
    echo "- Script: $test_script" >> "${REPORT_FILE}"
    
    # Create test-specific log directory
    local test_log_dir="${ARTIFACTS_DIR}/logs/${test_name}"
    mkdir -p "${test_log_dir}"
    local test_log="${test_log_dir}/test.log"
    local test_stderr="${test_log_dir}/stderr.log"
    
    # Run test with appropriate environment
    local exit_code=0
    local start_time=$(date +%s)
    
    if [[ "$BUILD_MODE" == "asan" ]]; then
        # Set ASAN/LSAN options
        export ASAN_OPTIONS="detect_leaks=1:abort_on_error=1:symbolize=1"
        export LSAN_OPTIONS="suppressions=${LSAN_SUPPRESSIONS}"
        
        # Run test and capture output
        if "${runner[@]}" > "${test_log}" 2> "${test_stderr}" || exit_code=$?; then
            : # Test passed
        fi
    else
        # Unset ASAN options for release mode to ensure clean environment
        unset ASAN_OPTIONS
        unset LSAN_OPTIONS
        
        # Run test and capture output (ensure clean environment)
        if env -u ASAN_OPTIONS -u LSAN_OPTIONS "${runner[@]}" > "${test_log}" 2> "${test_stderr}" || exit_code=$?; then
            : # Test passed
        fi
    fi
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    # Process ASAN/LSAN logs if in ASAN mode
    if [[ "$BUILD_MODE" == "asan" ]]; then
        local asan_log="${test_log_dir}/asan.log"
        local normalized_log="${test_log_dir}/normalized_leaks.txt"
        
        # Extract ASAN/LSAN output from both stdout and stderr
        {
            [[ -f "${test_log}" ]] && grep -E "(ERROR: LeakSanitizer|Direct leak|Indirect leak|SUMMARY: AddressSanitizer|AddressSanitizer:)" "${test_log}" || true
            [[ -f "${test_stderr}" ]] && grep -E "(ERROR: LeakSanitizer|Direct leak|Indirect leak|SUMMARY: AddressSanitizer|AddressSanitizer:)" "${test_stderr}" || true
        } > "${asan_log}" || true
        
        # Normalize leaks if normalize script exists and we have ASAN output
        if [[ -f "${NORMALIZE_SCRIPT}" && -s "${asan_log}" ]]; then
            python3 "${NORMALIZE_SCRIPT}" "${asan_log}" > "${normalized_log}" 2>/dev/null || true
            if [[ -s "${normalized_log}" ]]; then
                echo "- Normalized Leaks: ${normalized_log}" >> "${REPORT_FILE}"
            fi
        fi
        
        if [[ -s "${asan_log}" ]]; then
            echo "- ASAN Log: ${asan_log}" >> "${REPORT_FILE}"
        fi
    fi
    
    # Record result
    if [[ $exit_code -eq 0 ]]; then
        log_pass "$test_name (${duration}s)"
        PASSED_TESTS=$((PASSED_TESTS + 1))
        echo "- Status: PASS" >> "${REPORT_FILE}"
        echo "- Duration: ${duration}s" >> "${REPORT_FILE}"
        echo "- Log: ${test_log_dir}/" >> "${REPORT_FILE}"
    else
        log_fail "$test_name (exit code: $exit_code, ${duration}s)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        FAILED_TEST_NAMES+=("$test_name")
        echo "- Status: FAIL" >> "${REPORT_FILE}"
        echo "- Exit Code: $exit_code" >> "${REPORT_FILE}"
        echo "- Duration: ${duration}s" >> "${REPORT_FILE}"
        echo "- Log: ${test_log_dir}/" >> "${REPORT_FILE}"
        echo "- Error Output:" >> "${REPORT_FILE}"
        echo '```' >> "${REPORT_FILE}"
        tail -20 "${test_stderr}" >> "${REPORT_FILE}" || true
        echo '```' >> "${REPORT_FILE}"
        
        if [[ "$KEEP_GOING" == "false" ]]; then
            echo "" >> "${REPORT_FILE}"
            echo "## Summary" >> "${REPORT_FILE}"
            echo "- Total Tests: ${TOTAL_TESTS}" >> "${REPORT_FILE}"
            echo "- Passed: ${PASSED_TESTS}" >> "${REPORT_FILE}"
            echo "- Failed: ${FAILED_TESTS}" >> "${REPORT_FILE}"
            echo "- Skipped: ${SKIPPED_TESTS}" >> "${REPORT_FILE}"
            echo "" >> "${REPORT_FILE}"
            echo "**Stopped on first failure**" >> "${REPORT_FILE}"
            return 1
        fi
    fi
    
    echo "" >> "${REPORT_FILE}"
    return 0
}

# Build STAR binary
log_info "Building STAR (mode: ${BUILD_MODE})"
cd "${STAR_SOURCE_DIR}"

# Clean build to ensure no ASAN contamination in release mode
if [[ "$BUILD_MODE" == "release" ]]; then
    log_info "Cleaning previous build (ensuring no ASAN contamination)"
    make clean > /dev/null 2>&1 || true
fi

# Regenerate parametersDefault.xxd if needed (must be in source directory for correct variable names)
if [[ ! -f parametersDefault.xxd ]] || [[ parametersDefault -nt parametersDefault.xxd ]]; then
    log_info "Regenerating parametersDefault.xxd"
    # Must be in source directory and use relative filename to get correct variable names
    cd "${STAR_SOURCE_DIR}"
    xxd -i parametersDefault > parametersDefault.xxd
    cd - > /dev/null
fi

if [[ "$BUILD_MODE" == "asan" ]]; then
    log_info "Building with ASAN/LSAN enabled"
    make -j4 ASAN=1 STAR || {
        log_fail "STAR build failed (ASAN mode)"
        exit 1
    }
else
    log_info "Building release version"
    make -j4 STAR || {
        log_fail "STAR build failed (release mode)"
        exit 1
    }
    cp -f "${STAR_SOURCE_DIR}/STAR" "${RELEASE_BIN}"
    STAR_BIN="${RELEASE_BIN}"
fi

if [[ ! -f "${STAR_BIN}" ]]; then
    log_fail "STAR binary not found at ${STAR_BIN}"
    exit 1
fi

# Sanity check: ensure release build is not ASAN-instrumented
if [[ "$BUILD_MODE" == "release" ]]; then
    if ldd "${STAR_BIN}" 2>/dev/null | grep -qi "asan\\|lsan"; then
        log_fail "STAR binary links ASAN/LSAN but BUILD_MODE=release"
        exit 1
    fi
    if file "${STAR_BIN}" | grep -q "AddressSanitizer" || strings "${STAR_BIN}" | grep -q "__asan_" 2>/dev/null; then
        log_fail "STAR binary appears to be ASAN-instrumented but BUILD_MODE=release"
        exit 1
    fi
    export STAR_BIN
fi

log_pass "STAR build completed"

# Export all key variables so test scripts can access them
export STAR_BIN
export REMOVE_Y_READS_BIN
export TXIMPORT_BIN
export TXIMPORT_DIR
export SLAM_FIXTURE_ROOT
export FLEX_FIXTURE_ROOT
cd "${ROOT_DIR}"

# Export STAR_BIN for scripts that honor overrides.
export STAR_BIN

# Define test sets
declare -a STANDARD_TESTS=(
    # 1) Input parsing + guardrails
    "${SCRIPT_DIR}/test_readFilesIn_max_mates.sh"
    "${SCRIPT_DIR}/run_solo_smoke.sh"
    
    # 2) CR compatibility + CR multi
    "${SCRIPT_DIR}/run_cr_compat_integration_smoke.sh"
    "${SCRIPT_DIR}/cr_multi_smoke/run_gex_smoke.sh"
    "${SCRIPT_DIR}/cr_multi_smoke/run_cr_multi_smoke.sh"
    
    # 3) Flex pipelines
    "${SCRIPT_DIR}/flex_smoke/run_flex_smoke.sh"
    "${SCRIPT_DIR}/run_flex_inline_test.sh"
    "${SCRIPT_DIR}/run_flex_cbub_validation_test.sh"
    
    # 4) Y-chrom and FASTQ edge cases
    "${SCRIPT_DIR}/run_ychrom_fastq_se_test.sh"
    "${SCRIPT_DIR}/run_ychrom_bulk_pe_test.sh"
    "${SCRIPT_DIR}/run_ychrom_bam_split_test.sh"
    "${SCRIPT_DIR}/run_ychrom_regression_test.sh"
    "${SCRIPT_DIR}/run_remove_y_reads_test.sh"
    
    # 5) SLAM + transcriptVB + tximport
    "${SCRIPT_DIR}/run_slam_unit_tests.sh"
    "${SCRIPT_DIR}/run_slam_compat_test.sh"
    "${SCRIPT_DIR}/run_slam_fixture_parity.sh"
    "${SCRIPT_DIR}/transcriptvb/quick_test.sh"
    "${SCRIPT_DIR}/tximport/tximport_parity_test.sh"
)

if [[ "$BUILD_MODE" == "asan" ]]; then
    STANDARD_TESTS+=("${SCRIPT_DIR}/run_cr_compat_unit_tests.sh")
fi

declare -a EXTENDED_TESTS=(
    "${SCRIPT_DIR}/run_100K_regression_test.sh"
    "${SCRIPT_DIR}/run_downsample10M_alignment_and_gedi_goalset.sh"
    "${SCRIPT_DIR}/run_pe_preprod_validation.sh"
    "${SCRIPT_DIR}/run_flex_multisample_test.sh"
    "${SCRIPT_DIR}/run_flex_with_tags_test.sh"
    "${SCRIPT_DIR}/run_ychrom_flex_validation_test.sh"
    "${SCRIPT_DIR}/run_ychrom_spill_sorter.sh"
    "${SCRIPT_DIR}/run_slam_end_to_end.sh"
    "${SCRIPT_DIR}/run_slam_end_to_end_gedi_parity.sh"
    "${SCRIPT_DIR}/run_slam_end_to_end_vb_commonmask.sh"
    "${SCRIPT_DIR}/run_star_snp_mask_sweep.sh"
)

# Select test set
if [[ "$TEST_SET" == "standard" ]]; then
    TESTS=("${STANDARD_TESTS[@]}")
else
    TESTS=("${STANDARD_TESTS[@]}" "${EXTENDED_TESTS[@]}")
fi

log_info "Running ${TEST_SET} test set (${#TESTS[@]} tests)"
echo ""

# Run tests
for test_script in "${TESTS[@]}"; do
    if ! run_test "$test_script"; then
        # Test failed and KEEP_GOING is false
        break
    fi
done

# Final summary
echo "" >> "${REPORT_FILE}"
echo "## Summary" >> "${REPORT_FILE}"
echo "- Total Tests: ${TOTAL_TESTS}" >> "${REPORT_FILE}"
echo "- Passed: ${PASSED_TESTS}" >> "${REPORT_FILE}"
echo "- Failed: ${FAILED_TESTS}" >> "${REPORT_FILE}"
echo "- Skipped: ${SKIPPED_TESTS}" >> "${REPORT_FILE}"
echo "" >> "${REPORT_FILE}"

if [[ ${#FAILED_TEST_NAMES[@]} -gt 0 ]]; then
    echo "### Failed Tests" >> "${REPORT_FILE}"
    for test_name in "${FAILED_TEST_NAMES[@]}"; do
        echo "- ${test_name}" >> "${REPORT_FILE}"
    done
    echo "" >> "${REPORT_FILE}"
fi

if [[ ${#SKIPPED_TEST_NAMES[@]} -gt 0 ]]; then
    echo "### Skipped Tests" >> "${REPORT_FILE}"
    for test_name in "${SKIPPED_TEST_NAMES[@]}"; do
        echo "- ${test_name}" >> "${REPORT_FILE}"
    done
    echo "" >> "${REPORT_FILE}"
fi

# Print summary to console
echo ""
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Total Tests: ${TOTAL_TESTS}"
echo "Passed: ${PASSED_TESTS}"
echo "Failed: ${FAILED_TESTS}"
echo "Skipped: ${SKIPPED_TESTS}"
echo ""
echo "Report: ${REPORT_FILE}"
echo "Artifacts: ${ARTIFACTS_DIR}"
echo ""

# Exit with appropriate code
if [[ $FAILED_TESTS -gt 0 ]]; then
    exit 1
else
    exit 0
fi
