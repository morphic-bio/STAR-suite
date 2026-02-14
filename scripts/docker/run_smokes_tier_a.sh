#!/bin/bash
# Run Tier A smoke tests in-container (self-contained; no /storage required)
# Usage: ./scripts/docker/run_smokes_tier_a.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BASE_IMAGE_TAG="${IMAGE_TAG:-biodepot/star-suite:latest}"
TEST_IMAGE_TAG="${TEST_IMAGE_TAG:-biodepot/star-suite:test-tier-a}"
MAKE_JOBS="${MAKE_JOBS:-16}"

cd "$REPO_ROOT"

# Ensure Tier A test image exists (built from suite-base).
if ! docker image inspect "$TEST_IMAGE_TAG" &>/dev/null; then
    echo "Tier A test image $TEST_IMAGE_TAG not found. Building it..."
    docker build \
        --target test-tier-a \
        -f docker/Dockerfile \
        -t "$TEST_IMAGE_TAG" \
        --build-arg MAKE_JOBS="$MAKE_JOBS" \
        .
fi

echo "Running Tier A smoke tests in container..."
echo "  Base image tag: $BASE_IMAGE_TAG"
echo "  Test image tag: $TEST_IMAGE_TAG"
echo ""

# STAR_BIN for tests that use source-relative paths; runtime has STAR in /usr/local/bin
export STAR_BIN="/usr/local/bin/STAR"

run_test() {
    local name="$1"
    local script="$2"
    echo "=== $name ==="
    if docker run --rm \
        -e STAR_BIN="$STAR_BIN" \
        -w /build \
        "$TEST_IMAGE_TAG" \
        bash -c "cd /build && $script"; then
        echo "PASS: $name"
        return 0
    else
        echo "FAIL: $name"
        return 1
    fi
}

failed=0

run_test "run_solo_smoke" "tests/run_solo_smoke.sh" || failed=1
run_test "test_snp_mask_build_smoke" "tests/slam/test_snp_mask_build_smoke.sh" || failed=1

echo ""
if [ $failed -eq 0 ]; then
    echo "All Tier A smoke tests passed."
else
    echo "Some Tier A smoke tests failed."
    exit 1
fi
