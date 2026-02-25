#!/usr/bin/env bash
# Test: pfMultiConfig star_chemistry column parsing and validation
#
# Validates that PfMultiConfig correctly parses the star_chemistry column,
# rejects invalid values, and that the parsed values propagate into
# LibraryEntry.starChemistry correctly.
#
# This is a compile-and-run test using a small C++ harness that links
# PfMultiConfig.o and exercises parseConfig().

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0

pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: star_chemistry column parsing ==="

# --- Build the test harness ---
HARNESS="$WORK_DIR/test_config_chem"
cat > "$WORK_DIR/test_config_chem.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
#include "PfMultiConfig.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: test_config_chem <config.csv> [expect_fail]" << endl;
        return 2;
    }
    bool expectFail = (argc >= 3 && string(argv[2]) == "expect_fail");
    try {
        PfMultiConfig::Config cfg = PfMultiConfig::parseConfig(argv[1]);
        if (expectFail) {
            cerr << "ERROR: expected parse failure but succeeded" << endl;
            return 1;
        }
        for (size_t i = 0; i < cfg.libraries.size(); ++i) {
            const auto& lib = cfg.libraries[i];
            cout << "LIB:" << i
                 << " fastqs=" << lib.fastqs
                 << " feature_types=" << lib.feature_types
                 << " star_chemistry=" << lib.starChemistry
                 << endl;
        }
        return 0;
    } catch (const std::exception& e) {
        if (expectFail) {
            cout << "EXPECTED_FAIL: " << e.what() << endl;
            return 0;
        }
        cerr << "UNEXPECTED ERROR: " << e.what() << endl;
        return 1;
    }
}
HARNESS_EOF

SOURCE_DIR="$REPO_ROOT/core/legacy/source"
g++ -std=c++11 -O2 -I"$SOURCE_DIR" \
    "$WORK_DIR/test_config_chem.cpp" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    -o "$HARNESS" 2>&1
if [ $? -ne 0 ]; then
    echo "FATAL: failed to build test harness"
    exit 1
fi
echo "  Harness built OK"

# --- Test 1: star_chemistry column with mixed values ---
echo ""
echo "--- Test 1: mixed star_chemistry values ---"
cat > "$WORK_DIR/config_mixed.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry
/path/to/mRNA,S1,Gene Expression,Gene Expression,TRU
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture,NXT
/path/to/lineage,S1,Custom,Custom,auto
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_mixed.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_chemistry=tru'; then
    pass "GEX row parsed as tru"
else
    fail "GEX row star_chemistry not tru"
fi

if echo "$OUTPUT" | grep -q 'LIB:1.*star_chemistry=nxt'; then
    pass "gRNA row parsed as nxt"
else
    fail "gRNA row star_chemistry not nxt"
fi

if echo "$OUTPUT" | grep -q 'LIB:2.*star_chemistry=auto'; then
    pass "lineage row parsed as auto"
else
    fail "lineage row star_chemistry not auto"
fi

# --- Test 2: no star_chemistry column (backward compatibility) ---
echo ""
echo "--- Test 2: no star_chemistry column ---"
cat > "$WORK_DIR/config_no_chem.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types
/path/to/mRNA,S1,Gene Expression,Gene Expression
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_no_chem.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_chemistry=$'; then
    pass "GEX row has empty star_chemistry (backward compat)"
else
    fail "GEX row should have empty star_chemistry"
fi

if echo "$OUTPUT" | grep -q 'LIB:1.*star_chemistry=$'; then
    pass "gRNA row has empty star_chemistry (backward compat)"
else
    fail "gRNA row should have empty star_chemistry"
fi

# --- Test 3: empty star_chemistry values ---
echo ""
echo "--- Test 3: empty star_chemistry values ---"
cat > "$WORK_DIR/config_empty.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry
/path/to/mRNA,S1,Gene Expression,Gene Expression,
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture,NXT
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_empty.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_chemistry=$'; then
    pass "empty column value → empty starChemistry"
else
    fail "empty column value should yield empty starChemistry"
fi

if echo "$OUTPUT" | grep -q 'LIB:1.*star_chemistry=nxt'; then
    pass "NXT column value → nxt"
else
    fail "NXT column value should yield nxt"
fi

# --- Test 4: invalid star_chemistry value ---
echo ""
echo "--- Test 4: invalid star_chemistry value ---"
cat > "$WORK_DIR/config_invalid.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry
/path/to/mRNA,S1,Gene Expression,Gene Expression,INVALID
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_invalid.csv" "expect_fail" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'EXPECTED_FAIL.*Invalid star_chemistry'; then
    pass "invalid value rejected with error"
else
    fail "invalid value should be rejected"
fi

# --- Test 5: case-insensitive column name ---
echo ""
echo "--- Test 5: case-insensitive column name (Star_Chemistry) ---"
cat > "$WORK_DIR/config_case.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,Star_Chemistry
/path/to/mRNA,S1,Gene Expression,Gene Expression,Tru
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_case.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_chemistry=tru'; then
    pass "case-insensitive column name and value"
else
    fail "should handle case-insensitive column name"
fi

# --- Test 6: alias 'starchemistry' (no underscore) ---
echo ""
echo "--- Test 6: alias starchemistry (no underscore) ---"
cat > "$WORK_DIR/config_alias.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,starchemistry
/path/to/mRNA,S1,Gene Expression,Gene Expression,nxt
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_alias.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_chemistry=nxt'; then
    pass "starchemistry alias works"
else
    fail "starchemistry alias should work"
fi

# --- Test 7: trailing field omitted (no trailing comma) ---
echo ""
echo "--- Test 7: trailing field omitted (no trailing comma) ---"
cat > "$WORK_DIR/config_trailing.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry
/path/to/mRNA,S1,Gene Expression,Gene Expression
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture,NXT
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_trailing.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*fastqs=/path/to/mRNA'; then
    pass "row with omitted trailing field is not dropped"
else
    fail "row with omitted trailing field should not be dropped"
fi

if echo "$OUTPUT" | grep -q 'LIB:0.*star_chemistry=$'; then
    pass "omitted trailing field → empty starChemistry"
else
    fail "omitted trailing field should yield empty starChemistry"
fi

if echo "$OUTPUT" | grep -q 'LIB:1.*star_chemistry=nxt'; then
    pass "second row with explicit value still parsed"
else
    fail "second row should still parse correctly"
fi

# --- Test 8: whitespace-padded star_chemistry value ---
echo ""
echo "--- Test 8: whitespace-padded star_chemistry value ---"
cat > "$WORK_DIR/config_ws.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry
/path/to/mRNA,S1,Gene Expression,Gene Expression, TRU 
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture,  nxt  
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_ws.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_chemistry=tru'; then
    pass "whitespace-padded TRU trimmed and accepted"
else
    fail "whitespace-padded TRU should be trimmed"
fi

if echo "$OUTPUT" | grep -q 'LIB:1.*star_chemistry=nxt'; then
    pass "whitespace-padded nxt trimmed and accepted"
else
    fail "whitespace-padded nxt should be trimmed"
fi

echo ""
echo "=========================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=========================================="

[ "$FAIL_COUNT" -eq 0 ]
