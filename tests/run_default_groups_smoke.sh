#!/usr/bin/env bash
# Smoke tests for --default-* parameter bundles
# This harness runs existing tests with default groups enabled via STAR_EXTRA_ARGS
#
# Usage: ./tests/run_default_groups_smoke.sh [--group GROUP] [--quick]
#        GROUP: bulk, scrna, crcompat, flex, slam, vbem, a375
#        --quick: run only the fastest tests for each group
#
# Set STAR_BIN to override the STAR binary path.
# Set SKIP_BUILD=1 to skip rebuilding STAR.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
STAR_BIN="${STAR_BIN:-$REPO_ROOT/core/legacy/source/STAR}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Track results
declare -A RESULTS
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0
SKIPPED_TESTS=0

# Parse arguments
RUN_GROUP=""
QUICK_MODE=0
while [[ $# -gt 0 ]]; do
    case $1 in
        --group)
            RUN_GROUP="$2"
            shift 2
            ;;
        --quick)
            QUICK_MODE=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

log_info() {
    echo -e "${GREEN}[INFO]${NC} $*"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $*"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $*"
}

# Check if a script exists and is executable
check_script() {
    local script="$1"
    if [[ -f "$script" && -x "$script" ]]; then
        return 0
    elif [[ -f "$script" ]]; then
        chmod +x "$script"
        return 0
    fi
    return 1
}

# Verify default bundle was applied by checking Log.out for the specific group header
# Stricter check: requires the exact "Applying --defaultX" header, not just any [default-group] line
verify_defaults_applied() {
    local group="$1"
    local log_file="$2"
    local start_time="$3"  # Unix timestamp when test started
    
    if [[ ! -f "$log_file" ]]; then
        return 1
    fi
    
    # Check if Log.out is stale (older than test start time)
    local log_mtime
    log_mtime=$(stat -c %Y "$log_file" 2>/dev/null || stat -f %m "$log_file" 2>/dev/null || echo "0")
    if [[ "$log_mtime" -lt "$start_time" ]]; then
        # Log.out is stale (predates this test run)
        return 2  # Return 2 to indicate stale log
    fi
    
    # Require the specific "Applying --defaultX" header for this group
    if grep -q "Applying --default${group}" "$log_file" 2>/dev/null; then
        return 0
    fi
    return 1
}

# Find Log.out in test output directory
find_log_out() {
    local test_log="$1"
    local log_out=""
    
    # Try to extract output prefix from test log
    local out_prefix
    out_prefix=$(grep -oP '(?<=outFileNamePrefix\s)[^\s]+' "$test_log" 2>/dev/null | head -1 || true)
    
    if [[ -n "$out_prefix" && -f "${out_prefix}Log.out" ]]; then
        log_out="${out_prefix}Log.out"
    fi
    
    # Common locations to check
    local common_dirs=(
        "$SCRIPT_DIR/solo_smoke/output"
        "$SCRIPT_DIR/flex_smoke/output/flex_on"
        "$SCRIPT_DIR/flex_smoke/output/flex_off"
    )
    
    if [[ -z "$log_out" ]]; then
        for dir in "${common_dirs[@]}"; do
            if [[ -f "$dir/Log.out" ]]; then
                log_out="$dir/Log.out"
                break
            fi
        done
    fi
    
    echo "$log_out"
}

# Run a test with a default group
run_test() {
    local group="$1"
    local test_script="$2"
    local description="$3"
    local skip_default_check="${4:-0}"  # Optional: skip default verification
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    if ! check_script "$test_script"; then
        log_warn "Skipping $description: script not found ($test_script)"
        RESULTS["$group:$description"]="SKIPPED (script not found)"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        return 0
    fi
    
    log_info "Running $description with --default${group^}..."
    
    # Export STAR_EXTRA_ARGS to inject the default group flag
    export STAR_EXTRA_ARGS="--default${group^} Yes"
    export STAR_BIN
    
    local start_time=$(date +%s)
    local test_passed=0
    if timeout 600 bash "$test_script" > /tmp/default_group_test_$$.log 2>&1; then
        test_passed=1
    fi
    
    local exit_code=$?
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    if [[ $test_passed -eq 1 ]]; then
        # Test script passed, now verify defaults were actually applied
        if [[ "$skip_default_check" -eq 0 ]]; then
            local log_out
            log_out=$(find_log_out /tmp/default_group_test_$$.log)
            
            if [[ -n "$log_out" ]]; then
                verify_defaults_applied "$group" "$log_out" "$start_time"
                local verify_result=$?
                
                if [[ $verify_result -eq 0 ]]; then
                    log_info "  PASSED ($duration s) - defaults verified in Log.out"
                    RESULTS["$group:$description"]="PASSED (${duration}s, defaults verified)"
                    PASSED_TESTS=$((PASSED_TESTS + 1))
                elif [[ $verify_result -eq 2 ]]; then
                    log_error "  FAILED - Log.out is stale (predates test run)"
                    RESULTS["$group:$description"]="FAILED (stale Log.out)"
                    FAILED_TESTS=$((FAILED_TESTS + 1))
                    echo "  Log.out location: $log_out"
                else
                    log_error "  FAILED - test passed but defaults NOT applied"
                    RESULTS["$group:$description"]="FAILED (defaults not applied)"
                    FAILED_TESTS=$((FAILED_TESTS + 1))
                    echo "  Log.out location: $log_out"
                    echo "  Expected to find: 'Applying --default${group^}' header"
                fi
            else
                # Can't find Log.out, assume pass but warn
                log_warn "  PASSED ($duration s) - could not verify defaults (Log.out not found)"
                RESULTS["$group:$description"]="PASSED (${duration}s, unverified)"
                PASSED_TESTS=$((PASSED_TESTS + 1))
            fi
        else
            log_info "  PASSED ($duration s)"
            RESULTS["$group:$description"]="PASSED (${duration}s)"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        fi
    else
        if [[ $exit_code -eq 124 ]]; then
            log_error "  TIMEOUT after 600s"
            RESULTS["$group:$description"]="TIMEOUT"
        else
            log_error "  FAILED (exit $exit_code, ${duration}s)"
            RESULTS["$group:$description"]="FAILED (exit $exit_code)"
        fi
        FAILED_TESTS=$((FAILED_TESTS + 1))
        # Show last 20 lines of output on failure
        echo "  Last 20 lines of output:"
        tail -20 /tmp/default_group_test_$$.log | sed 's/^/    /'
    fi
    
    unset STAR_EXTRA_ARGS
}

# Check fixtures for a test
check_fixtures() {
    local fixture_path="$1"
    if [[ -e "$fixture_path" ]]; then
        return 0
    fi
    return 1
}

# Build STAR if needed
build_star() {
    if [[ "${SKIP_BUILD:-0}" == "1" ]]; then
        log_info "Skipping STAR build (SKIP_BUILD=1)"
        return 0
    fi
    
    if [[ ! -f "$STAR_BIN" ]]; then
        log_info "Building STAR..."
        cd "$REPO_ROOT"
        make core -j$(nproc) || {
            log_error "Failed to build STAR"
            exit 1
        }
    fi
    
    log_info "Using STAR: $STAR_BIN"
}

# Test groups
test_bulk() {
    local group="Bulk"
    log_info "=== Testing --defaultBulk ==="
    
    # Primary: minimal fixture that relies entirely on default bundle
    if check_script "$SCRIPT_DIR/run_default_bundle_bulk_fixture.sh"; then
        # Minimal fixture has its own STAR_EXTRA_ARGS, so we run it directly
        log_info "Running minimal bulk fixture (relies on --defaultBulk)..."
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
        local start_time=$(date +%s)
        if timeout 300 bash "$SCRIPT_DIR/run_default_bundle_bulk_fixture.sh" > /tmp/default_group_test_$$.log 2>&1; then
            local end_time=$(date +%s)
            local duration=$((end_time - start_time))
            log_info "  PASSED ($duration s) - defaults verified by fixture"
            RESULTS["$group:minimal fixture"]="PASSED (${duration}s, verified)"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            local exit_code=$?
            local end_time=$(date +%s)
            local duration=$((end_time - start_time))
            log_error "  FAILED (exit $exit_code, ${duration}s)"
            RESULTS["$group:minimal fixture"]="FAILED (exit $exit_code)"
            FAILED_TESTS=$((FAILED_TESTS + 1))
            echo "  Last 20 lines of output:"
            tail -20 /tmp/default_group_test_$$.log | sed 's/^/    /'
        fi
    fi
    
    # Secondary: existing test script with injected defaults
    if [[ $QUICK_MODE -eq 0 ]]; then
        run_test "$group" "$SCRIPT_DIR/run_baseline_test.sh" "baseline bulk test"
        
        if check_fixtures "/storage/test_data/bulk"; then
            run_test "$group" "$SCRIPT_DIR/run_pe_preprod_validation.sh" "PE preprod validation"
        fi
    fi
}

test_scrna() {
    local group="CoreScrna"
    log_info "=== Testing --defaultCoreScrna ==="
    
    # Primary: minimal fixture that relies entirely on default bundle
    if check_script "$SCRIPT_DIR/run_default_bundle_scrna_fixture.sh"; then
        # Minimal fixture has its own STAR_EXTRA_ARGS, so we run it directly
        log_info "Running minimal scRNA fixture (relies on --defaultCoreScrna)..."
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
        local start_time=$(date +%s)
        if timeout 300 bash "$SCRIPT_DIR/run_default_bundle_scrna_fixture.sh" > /tmp/default_group_test_$$.log 2>&1; then
            local end_time=$(date +%s)
            local duration=$((end_time - start_time))
            log_info "  PASSED ($duration s) - defaults verified by fixture"
            RESULTS["$group:minimal fixture"]="PASSED (${duration}s, verified)"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            local exit_code=$?
            local end_time=$(date +%s)
            local duration=$((end_time - start_time))
            log_error "  FAILED (exit $exit_code, ${duration}s)"
            RESULTS["$group:minimal fixture"]="FAILED (exit $exit_code)"
            FAILED_TESTS=$((FAILED_TESTS + 1))
            echo "  Last 20 lines of output:"
            tail -20 /tmp/default_group_test_$$.log | sed 's/^/    /'
        fi
    fi
    
    # Secondary: existing test script with injected defaults
    if [[ $QUICK_MODE -eq 0 ]]; then
        run_test "$group" "$SCRIPT_DIR/run_solo_smoke.sh" "Solo scRNA smoke test"
    fi
}

test_crcompat() {
    local group="CrCompat"
    log_info "=== Testing --defaultCrCompat ==="
    
    # Check for A375 fixtures
    if check_fixtures "/storage/A375"; then
        run_test "$group" "$SCRIPT_DIR/test_cr_compat_crispr_calling.sh" "CR-compat CRISPR calling"
        
        if [[ $QUICK_MODE -eq 0 ]]; then
            run_test "$group" "$SCRIPT_DIR/run_cr_compat_integration_smoke.sh" "CR-compat integration"
        fi
    else
        log_warn "Skipping CR-compat tests: /storage/A375 fixtures not found"
        RESULTS["$group:CR-compat tests"]="SKIPPED (fixtures missing)"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
    fi
}

test_flex() {
    local group="Flex"
    log_info "=== Testing --defaultFlex ==="
    
    local flex_smoke="$SCRIPT_DIR/flex_smoke/run_flex_smoke.sh"
    if check_script "$flex_smoke"; then
        run_test "$group" "$flex_smoke" "Flex smoke test"
    else
        log_warn "Skipping Flex tests: $flex_smoke not found"
        RESULTS["$group:Flex smoke"]="SKIPPED (script not found)"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
    fi
}

test_slam() {
    local group="Slam"
    log_info "=== Testing --defaultSlam ==="
    
    # Find the smallest SLAM test
    local slam_tests=(
        "$SCRIPT_DIR/run_slam_quick_smoke.sh"
        "$SCRIPT_DIR/run_slam_smoke.sh"
        "$SCRIPT_DIR/slam/run_slam_test.sh"
    )
    
    local found_test=""
    for test in "${slam_tests[@]}"; do
        if check_script "$test"; then
            found_test="$test"
            break
        fi
    done
    
    if [[ -n "$found_test" ]]; then
        run_test "$group" "$found_test" "SLAM smoke test"
    else
        log_warn "Skipping SLAM tests: no slam test script found"
        RESULTS["$group:SLAM smoke"]="SKIPPED (script not found)"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
    fi
}

test_vbem() {
    local group="Vbem"
    log_info "=== Testing --defaultVbem ==="
    
    # Find a VBEM test
    local vbem_tests=(
        "$SCRIPT_DIR/transcriptvb/run_vbem_smoke.sh"
        "$REPO_ROOT/core/features/vbem/test/run_smoke.sh"
    )
    
    local found_test=""
    for test in "${vbem_tests[@]}"; do
        if check_script "$test"; then
            found_test="$test"
            break
        fi
    done
    
    if [[ -n "$found_test" ]]; then
        run_test "$group" "$found_test" "VBEM smoke test"
    else
        log_warn "Skipping VBEM tests: no vbem test script found"
        RESULTS["$group:VBEM smoke"]="SKIPPED (script not found)"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
    fi
}

test_a375parity() {
    local group="A375Parity"
    log_info "=== Testing --defaultA375Parity ==="
    
    if check_fixtures "/storage/A375"; then
        run_test "$group" "$SCRIPT_DIR/run_a375_gex_features_cr_parity.sh" "A375 GEX+features parity"
    else
        log_warn "Skipping A375 parity tests: /storage/A375 fixtures not found"
        RESULTS["$group:A375 parity"]="SKIPPED (fixtures missing)"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        TOTAL_TESTS=$((TOTAL_TESTS + 1))
    fi
}

test_cbub() {
    log_info "=== Testing CB/UB (uses parent pipeline defaults) ==="
    
    # CB/UB should work with the parent pipeline defaults
    # Test with --defaultCoreScrna since CB/UB is typically used with scRNA
    export STAR_EXTRA_ARGS="--defaultCoreScrna Yes"
    export STAR_BIN
    
    local group="CoreScrna"
    run_test "$group" "$SCRIPT_DIR/run_unsorted_cbub_smoke_test.sh" "CB/UB smoke test (with scRNA defaults)"
    
    unset STAR_EXTRA_ARGS
}

# Print summary
print_summary() {
    echo ""
    echo "=============================================="
    echo "Default Groups Smoke Test Summary"
    echo "=============================================="
    echo "Total: $TOTAL_TESTS | Passed: $PASSED_TESTS | Failed: $FAILED_TESTS | Skipped: $SKIPPED_TESTS"
    echo ""
    
    if [[ ${#RESULTS[@]} -gt 0 ]]; then
        echo "Detailed Results:"
        for key in "${!RESULTS[@]}"; do
            printf "  %-50s %s\n" "$key" "${RESULTS[$key]}"
        done
    fi
    
    echo ""
    if [[ $FAILED_TESTS -gt 0 ]]; then
        log_error "Some tests failed!"
        return 1
    elif [[ $PASSED_TESTS -eq 0 ]]; then
        log_warn "No tests passed (all skipped or none run)"
        return 0
    else
        log_info "All executed tests passed!"
        return 0
    fi
}

# Main
main() {
    log_info "Default Groups Smoke Test Harness"
    log_info "=================================="
    log_info "Repository: $REPO_ROOT"
    log_info "Quick mode: $([[ $QUICK_MODE -eq 1 ]] && echo 'yes' || echo 'no')"
    [[ -n "$RUN_GROUP" ]] && log_info "Running group: $RUN_GROUP"
    echo ""
    
    build_star
    
    if [[ -n "$RUN_GROUP" ]]; then
        # Run specific group
        case "$RUN_GROUP" in
            bulk)       test_bulk ;;
            scrna)      test_scrna ;;
            crcompat)   test_crcompat ;;
            flex)       test_flex ;;
            slam)       test_slam ;;
            vbem)       test_vbem ;;
            a375)       test_a375parity ;;
            cbub)       test_cbub ;;
            *)
                log_error "Unknown group: $RUN_GROUP"
                echo "Valid groups: bulk, scrna, crcompat, flex, slam, vbem, a375, cbub"
                exit 1
                ;;
        esac
    else
        # Run all groups
        test_bulk
        test_scrna
        test_crcompat
        test_flex
        test_slam
        test_vbem
        test_a375parity
        test_cbub
    fi
    
    print_summary
}

main "$@"
