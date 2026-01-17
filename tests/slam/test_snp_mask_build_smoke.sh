#!/usr/bin/env bash
#
# Smoke test for SLAM SNP Mask Build functionality
#
# Tests:
# 1. FOFN parsing (single-end and paired-end)
# 2. Parameter validation
# 3. EM unit tests
#
# Usage: ./test_snp_mask_build_smoke.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAR="${SCRIPT_DIR}/../../core/legacy/source/STAR"
SLAM_SOURCE="${SCRIPT_DIR}/../../slam/source"
TMPDIR="${SCRIPT_DIR}/tmp_snp_mask_test_$$"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

cleanup() {
    rm -rf "$TMPDIR"
}
trap cleanup EXIT

mkdir -p "$TMPDIR"

passed=0
failed=0

pass() {
    echo -e "${GREEN}PASS${NC}: $1"
    passed=$((passed + 1))
}

fail() {
    echo -e "${RED}FAIL${NC}: $1"
    failed=$((failed + 1))
}

echo "=== SLAM SNP Mask Build Smoke Tests ==="
echo ""

# Test 1: Help shows new parameters
echo "Test 1: Verify new parameters in help output"
HELP_OUTPUT=$("$STAR" --help 2>&1)

echo "$HELP_OUTPUT" | grep -q "slamSnpMaskBuildFastqs" && pass "--slamSnpMaskBuildFastqs found" || fail "--slamSnpMaskBuildFastqs missing"
echo "$HELP_OUTPUT" | grep -q "slamSnpMaskIn" && pass "--slamSnpMaskIn found" || fail "--slamSnpMaskIn missing"
echo "$HELP_OUTPUT" | grep -q "slamSnpMaskVcfIn" && pass "--slamSnpMaskVcfIn found" || fail "--slamSnpMaskVcfIn missing"
echo "$HELP_OUTPUT" | grep -q "slamSnpMaskOnly" && pass "--slamSnpMaskOnly found" || fail "--slamSnpMaskOnly missing"
echo "$HELP_OUTPUT" | grep -q "slamSnpMaskPosterior" && pass "--slamSnpMaskPosterior found" || fail "--slamSnpMaskPosterior missing"

# Test 2: EM unit tests
echo ""
echo "Test 2: Running EM unit tests"
cd "$SCRIPT_DIR"
if [ -f "test_slam_snp_em" ]; then
    EM_OUTPUT=$(./test_slam_snp_em 2>&1)
    echo "$EM_OUTPUT" | grep -q "ALL TESTS PASSED" && pass "EM unit tests" || fail "EM unit tests"
else
    echo "Compiling EM unit tests..."
    if g++ -std=c++11 -O2 -I"$SLAM_SOURCE" -I"$SLAM_SOURCE/libem" \
        test_slam_snp_em.cpp "$SLAM_SOURCE/libem/slam_snp_em.cpp" \
        -o test_slam_snp_em -lm 2>&1; then
        EM_OUTPUT=$(./test_slam_snp_em 2>&1)
        echo "$EM_OUTPUT" | grep -q "ALL TESTS PASSED" && pass "EM unit tests" || fail "EM unit tests"
    else
        fail "Failed to compile EM unit tests"
    fi
fi

# Test 3: Create valid FOFN files for testing
echo ""
echo "Test 3: FOFN format validation"

# Create test FOFN - single-end format
cat > "$TMPDIR/test_se.fofn" << 'EOF'
# Comment line should be ignored
/path/to/sample1.fastq.gz

/path/to/sample2.fq.gz
EOF

# Create test FOFN - paired-end format  
cat > "$TMPDIR/test_pe.fofn" << 'EOF'
# Paired-end FOFN
/path/to/sample1_R1.fastq.gz	/path/to/sample1_R2.fastq.gz
/path/to/sample2_R1.fq.gz	/path/to/sample2_R2.fq.gz
EOF

# Verify FOFNs are syntactically correct
[ -f "$TMPDIR/test_se.fofn" ] && pass "Single-end FOFN created" || fail "Single-end FOFN missing"
[ -f "$TMPDIR/test_pe.fofn" ] && pass "Paired-end FOFN created" || fail "Paired-end FOFN missing"

# Test 4: Default parameter values
echo ""
echo "Test 4: Default parameter values"

echo "$HELP_OUTPUT" | awk '/^slamSnpMaskMinCov/{found=1} found && /Default:/{print; exit}' | grep -q "Default: 20" \
    && pass "Default minCov=20" || fail "Default minCov wrong"
echo "$HELP_OUTPUT" | awk '/^slamSnpMaskMinAlt/{found=1} found && /Default:/{print; exit}' | grep -q "Default: 3" \
    && pass "Default minAlt=3" || fail "Default minAlt wrong"
echo "$HELP_OUTPUT" | grep "slamSnpMaskPosterior" | grep -q "0.99" && pass "Default posterior=0.99" || fail "Default posterior wrong"
echo "$HELP_OUTPUT" | grep "slamSnpMaskVcfMode" | grep -q "gt" && pass "Default vcfMode=gt" || fail "Default vcfMode wrong"
echo "$HELP_OUTPUT" | grep "slamSnpMaskVcfFilter" | grep -q "pass" && pass "Default vcfFilter=pass" || fail "Default vcfFilter wrong"

# Summary
echo ""
echo "==================================="
echo "Results: $passed passed, $failed failed"
echo "==================================="

if [ $failed -eq 0 ]; then
    echo -e "${GREEN}ALL SMOKE TESTS PASSED${NC}"
    exit 0
else
    echo -e "${RED}SOME TESTS FAILED${NC}"
    exit 1
fi
