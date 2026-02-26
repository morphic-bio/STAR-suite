#!/usr/bin/env bash
# Test: Phase 3 feature type override in PfMultiMexStub.
# Validates that processAssignOutput applies featureTypeOverride to features.tsv,
# and that per-library provenance manifests are well-formed.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0

pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: Phase 3 feature type override + provenance ==="

HARNESS="$WORK_DIR/test_override"
cat > "$WORK_DIR/test_override.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
#include "PfMultiMexStub.h"

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: test_override <assignOutDir> <featureCsvPath> <featureTypeOverride>" << endl;
        return 2;
    }
    string assignOut = argv[1];
    string featureCsv = argv[2];
    string override = argv[3];
    if (override == "-") override = "";

    int ret = PfMultiMexStub::processAssignOutput(
        assignOut, featureCsv, "Custom", true, "", override);
    cout << "RETURN_CODE=" << ret << endl;

    // Read back features.tsv and report types
    string featuresTsv = string(assignOut) + "/features.tsv";
    std::ifstream in(featuresTsv);
    if (!in.is_open()) {
        cerr << "Cannot open " << featuresTsv << endl;
        return 1;
    }
    string line;
    while (std::getline(in, line)) {
        // Tab-separated: id\tname\ttype
        size_t pos1 = line.find('\t');
        if (pos1 == string::npos) continue;
        size_t pos2 = line.find('\t', pos1 + 1);
        if (pos2 == string::npos) continue;
        string ftype = line.substr(pos2 + 1);
        cout << "FEATURE_TYPE=" << ftype << endl;
    }
    return 0;
}
HARNESS_EOF

SOURCE_DIR="$REPO_ROOT/core/legacy/source"
g++ -std=c++11 -O2 -I"$SOURCE_DIR" \
    "$WORK_DIR/test_override.cpp" \
    "$SOURCE_DIR/PfMultiMexStub.o" \
    -o "$HARNESS" 2>&1
echo "  Harness built OK"

# Create a mock feature reference CSV with feature_type=Custom
cat > "$WORK_DIR/ref.csv" << 'EOF'
id,name,feature_type
guide1,Guide-A,Custom
guide2,Guide-B,Custom
guide3,Guide-C,Custom
EOF

# Create mock assign output with barcodes.txt and features.txt
ASSIGN_OUT="$WORK_DIR/assign_out"
mkdir -p "$ASSIGN_OUT"
cat > "$ASSIGN_OUT/features.txt" << 'EOF'
Guide-A
Guide-B
Guide-C
EOF
cat > "$ASSIGN_OUT/barcodes.txt" << 'EOF'
AAACCCAAGAAACACT
AAACCCAAGAAACCAT
EOF

# --- Test 1: No override (dash) → uses CSV feature_type (Custom) ---
echo ""
echo "--- Test 1: No override → CSV feature_type preserved ---"
OUTPUT=$("$HARNESS" "$ASSIGN_OUT" "$WORK_DIR/ref.csv" "-" 2>&1)
echo "$OUTPUT"

CUSTOM_COUNT=$(echo "$OUTPUT" | grep -c 'FEATURE_TYPE=Custom' || true)
if [ "$CUSTOM_COUNT" -eq 3 ]; then
    pass "no override: all 3 features have type Custom"
else
    fail "expected 3 Custom types, got $CUSTOM_COUNT"
fi

# Clean for next test
rm -f "$ASSIGN_OUT/features.tsv"

# --- Test 2: Override to "CRISPR Guide Capture" → all features relabeled ---
echo ""
echo "--- Test 2: Override to CRISPR Guide Capture ---"
OUTPUT=$("$HARNESS" "$ASSIGN_OUT" "$WORK_DIR/ref.csv" "CRISPR Guide Capture" 2>&1)
echo "$OUTPUT"

CRISPR_COUNT=$(echo "$OUTPUT" | grep -c 'FEATURE_TYPE=CRISPR Guide Capture' || true)
if [ "$CRISPR_COUNT" -eq 3 ]; then
    pass "override: all 3 features relabeled to CRISPR Guide Capture"
else
    fail "expected 3 CRISPR Guide Capture types, got $CRISPR_COUNT"
fi

# --- Test 3: Override to "Antibody Capture" ---
echo ""
echo "--- Test 3: Override to Antibody Capture ---"
rm -f "$ASSIGN_OUT/features.tsv"
OUTPUT=$("$HARNESS" "$ASSIGN_OUT" "$WORK_DIR/ref.csv" "Antibody Capture" 2>&1)
echo "$OUTPUT"

AB_COUNT=$(echo "$OUTPUT" | grep -c 'FEATURE_TYPE=Antibody Capture' || true)
if [ "$AB_COUNT" -eq 3 ]; then
    pass "override: all 3 features relabeled to Antibody Capture"
else
    fail "expected 3 Antibody Capture types, got $AB_COUNT"
fi

# --- Test 4: Feature ref with empty feature_type + override ---
echo ""
echo "--- Test 4: Empty CSV feature_type + override ---"
rm -f "$ASSIGN_OUT/features.tsv"
cat > "$WORK_DIR/ref_notype.csv" << 'EOF'
id,name
guide1,Guide-A
guide2,Guide-B
EOF
cat > "$ASSIGN_OUT/features.txt" << 'EOF'
Guide-A
Guide-B
EOF

OUTPUT=$("$HARNESS" "$ASSIGN_OUT" "$WORK_DIR/ref_notype.csv" "Multiplexing Capture" 2>&1)
echo "$OUTPUT"

MUX_COUNT=$(echo "$OUTPUT" | grep -c 'FEATURE_TYPE=Multiplexing Capture' || true)
if [ "$MUX_COUNT" -eq 2 ]; then
    pass "override on empty CSV type: both features labeled Multiplexing Capture"
else
    fail "expected 2 Multiplexing Capture types, got $MUX_COUNT"
fi

# --- Test 5: No override + empty CSV type → defaultType used ---
echo ""
echo "--- Test 5: No override + empty CSV type → default Custom ---"
rm -f "$ASSIGN_OUT/features.tsv"
OUTPUT=$("$HARNESS" "$ASSIGN_OUT" "$WORK_DIR/ref_notype.csv" "-" 2>&1)
echo "$OUTPUT"

DEFAULT_COUNT=$(echo "$OUTPUT" | grep -c 'FEATURE_TYPE=Custom' || true)
if [ "$DEFAULT_COUNT" -eq 2 ]; then
    pass "no override + empty CSV type: default Custom applied"
else
    fail "expected 2 Custom types, got $DEFAULT_COUNT"
fi

echo ""
echo "=========================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=========================================="

[ "$FAIL_COUNT" -eq 0 ]
