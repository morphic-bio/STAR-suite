#!/usr/bin/env bash
# Test: Phase 4 library-aware PF scheduler allocation logic.
# The harness mirrors the exact production code from PfMultiProcess.cpp,
# including the runtime clamp (max(1, ...)).
# Validates:
#   1. Single library gets full thread budget (backward compat)
#   2. Two equal libraries split budget proportionally
#   3. Two unequal libraries get proportional shares (larger gets more)
#   4. Three libraries, each gets at least 2 threads (min guarantee)
#   5. runThreadN=1: single library gets 1 thread (no clamp to 2)
#   6. runThreadN=2, two libraries: each gets 2 (min guarantee)
#   7. Strict conservation: sum == runThreadN when runThreadN >= numLibs*2
#   8. Dynamic permits disabled -> each gets full budget
#   9. Conservation for 10 threads, 4 equal libs
#  10. Conservation for 10 threads, 4 unequal libs
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0

pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: Phase 4 library-aware PF scheduler ==="

HARNESS="$WORK_DIR/test_scheduler"
cat > "$WORK_DIR/test_scheduler.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

struct LibrarySchedule {
    uint64_t estimatedWork = 0;
    int fileCount = 0;
    int threadBudget = 0;
};

// Exact mirror of PfMultiProcess.cpp scheduler + runtime clamp
static vector<int> computeSchedule(int runThreadN, const vector<uint64_t>& workEstimates,
                                   bool dynamicPermits) {
    size_t numLibs = workEstimates.size();
    vector<LibrarySchedule> schedules(numLibs);
    uint64_t totalWork = 0;
    for (size_t i = 0; i < numLibs; ++i) {
        uint64_t w = workEstimates[i];
        if (w == 0) w = 1;
        schedules[i].estimatedWork = w;
        totalWork += w;
    }

    if (numLibs > 1 && dynamicPermits) {
        const int minPerLib = 2;
        const int guaranteedTotal = static_cast<int>(numLibs) * minPerLib;

        if (runThreadN < guaranteedTotal) {
            for (size_t i = 0; i < numLibs; ++i) {
                schedules[i].threadBudget = minPerLib;
            }
        } else {
            const int surplus = runThreadN - guaranteedTotal;
            int floorSum = 0;
            vector<double> fractionalParts(numLibs);
            for (size_t i = 0; i < numLibs; ++i) {
                double exactExtra = static_cast<double>(schedules[i].estimatedWork)
                                  / static_cast<double>(std::max<uint64_t>(1, totalWork))
                                  * surplus;
                int floorExtra = static_cast<int>(exactExtra);
                schedules[i].threadBudget = minPerLib + floorExtra;
                fractionalParts[i] = exactExtra - floorExtra;
                floorSum += floorExtra;
            }
            int leftover = surplus - floorSum;
            while (leftover > 0) {
                size_t bestIdx = 0;
                for (size_t i = 1; i < numLibs; ++i) {
                    if (fractionalParts[i] > fractionalParts[bestIdx]) {
                        bestIdx = i;
                    }
                }
                schedules[bestIdx].threadBudget += 1;
                fractionalParts[bestIdx] = -1.0;
                --leftover;
            }
        }
    } else {
        for (size_t i = 0; i < numLibs; ++i) {
            schedules[i].threadBudget = runThreadN;
        }
    }

    // Apply runtime clamp: max(1, ...) — mirrors PfMultiProcess.cpp
    vector<int> result;
    for (size_t i = 0; i < numLibs; ++i) {
        result.push_back(std::max(1, schedules[i].threadBudget));
    }
    return result;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: test_scheduler <runThreadN> <dynamicPermits:0|1> <work1> [work2] ..." << endl;
        return 2;
    }
    int runThreadN = atoi(argv[1]);
    bool dynamicPermits = (atoi(argv[2]) != 0);
    vector<uint64_t> works;
    for (int i = 3; i < argc; ++i) {
        works.push_back(static_cast<uint64_t>(atoll(argv[i])));
    }
    vector<int> budgets = computeSchedule(runThreadN, works, dynamicPermits);
    for (size_t i = 0; i < budgets.size(); ++i) {
        if (i > 0) cout << ",";
        cout << budgets[i];
    }
    cout << endl;
    return 0;
}
HARNESS_EOF

g++ -std=c++11 -O2 "$WORK_DIR/test_scheduler.cpp" -o "$HARNESS" 2>&1
echo "  Harness built OK"

# Helper: parse comma-separated output, check sum and min
check_sum() {
    local out="$1" expected_sum="$2"
    IFS=',' read -ra B <<< "$out"
    local sum=0
    for b in "${B[@]}"; do sum=$((sum + b)); done
    [ "$sum" -eq "$expected_sum" ]
}

check_all_ge() {
    local out="$1" min_val="$2"
    IFS=',' read -ra B <<< "$out"
    for b in "${B[@]}"; do
        [ "$b" -ge "$min_val" ] || return 1
    done
    return 0
}

# --- Test 1: Single library gets full budget ---
echo ""
echo "--- Test 1: Single library gets full budget ---"
OUT=$("$HARNESS" 16 1 100000)
if [ "$OUT" = "16" ]; then
    pass "single library: budget=16"
else
    fail "single library: expected 16, got '$OUT'"
fi

# --- Test 2: Two equal libraries split proportionally ---
echo ""
echo "--- Test 2: Two equal libraries split proportionally ---"
OUT=$("$HARNESS" 16 1 100000 100000)
IFS=',' read -ra B <<< "$OUT"
if [ "${B[0]}" = "${B[1]}" ] && check_sum "$OUT" 16; then
    pass "two equal libs: ${B[0]}+${B[1]}=16"
else
    fail "two equal libs: got $OUT"
fi

# --- Test 3: Two unequal libraries, larger gets more ---
echo ""
echo "--- Test 3: Two unequal libraries, larger gets more ---"
OUT=$("$HARNESS" 16 1 900000 100000)
IFS=',' read -ra B <<< "$OUT"
if [ "${B[0]}" -gt "${B[1]}" ] && check_sum "$OUT" 16; then
    pass "unequal libs: ${B[0]}>${B[1]}, sum=16"
else
    fail "unequal libs: got $OUT"
fi

# --- Test 4: Three libraries, each >= 2 ---
echo ""
echo "--- Test 4: Three libraries, min 2 each ---"
OUT=$("$HARNESS" 8 1 100000 100000 100000)
if check_all_ge "$OUT" 2 && check_sum "$OUT" 8; then
    pass "three libs: all >=2, sum=8 ($OUT)"
else
    fail "three libs: got $OUT"
fi

# --- Test 5: runThreadN=1, single library: budget=1 (no clamp to 2) ---
echo ""
echo "--- Test 5: runThreadN=1, single library ---"
OUT=$("$HARNESS" 1 1 100000)
if [ "$OUT" = "1" ]; then
    pass "1-thread single lib: budget=1 (runtime clamp is max(1,...))"
else
    fail "1-thread single lib: expected 1, got '$OUT'"
fi

# --- Test 6: runThreadN=2, two libraries (min guarantee) ---
echo ""
echo "--- Test 6: runThreadN=2, two libraries ---"
OUT=$("$HARNESS" 2 1 100000 100000)
if check_all_ge "$OUT" 2; then
    pass "2-thread 2 libs: each >= 2 ($OUT)"
else
    fail "2-thread 2 libs: got $OUT"
fi

# --- Test 7: Strict conservation: 32 threads, 4 equal libs ---
echo ""
echo "--- Test 7: Strict conservation (32 threads, 4 equal libs) ---"
OUT=$("$HARNESS" 32 1 250000 250000 250000 250000)
if check_sum "$OUT" 32 && check_all_ge "$OUT" 2; then
    pass "32t/4libs: sum=32, all>=2 ($OUT)"
else
    fail "32t/4libs: got $OUT"
fi

# --- Test 8: Dynamic permits disabled -> each gets full budget ---
echo ""
echo "--- Test 8: Permits disabled ---"
OUT=$("$HARNESS" 16 0 100000 100000)
if [ "$OUT" = "16,16" ]; then
    pass "permits off: 16,16"
else
    fail "permits off: got $OUT"
fi

# --- Test 9: Conservation for 10 threads, 4 equal libs ---
echo ""
echo "--- Test 9: Conservation (10 threads, 4 equal libs) ---"
OUT=$("$HARNESS" 10 1 100000 100000 100000 100000)
if check_sum "$OUT" 10 && check_all_ge "$OUT" 2; then
    pass "10t/4equal: sum=10, all>=2 ($OUT)"
else
    fail "10t/4equal: got $OUT"
fi

# --- Test 10: Conservation for 10 threads, 4 unequal libs ---
echo ""
echo "--- Test 10: Conservation (10 threads, 4 unequal libs) ---"
OUT=$("$HARNESS" 10 1 700000 100000 100000 100000)
if check_sum "$OUT" 10 && check_all_ge "$OUT" 2; then
    pass "10t/4unequal: sum=10, all>=2 ($OUT)"
else
    fail "10t/4unequal: got $OUT"
fi

echo ""
echo "=========================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=========================================="

[ "$FAIL_COUNT" -eq 0 ]
