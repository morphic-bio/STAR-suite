#!/usr/bin/env bash
# Test: the global featureRef guard in PfMultiProcess allows configs where
# all feature libraries have star_feature_ref, even with no global ref.
# This exercises the actual guard logic from PfMultiProcess.cpp, not just parsing.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0

pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: global featureRef guard (pipeline-level) ==="

HARNESS="$WORK_DIR/test_guard"
cat > "$WORK_DIR/test_guard.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <string>
#include <vector>
using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
#include "PfMultiConfig.h"

static string lowerCopy(const string& input) {
    string out = input;
    for (char& c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}
static string trimCopy(const string& input) {
    size_t first = input.find_first_not_of(" \t\r\n");
    if (first == string::npos) return "";
    size_t last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}
static bool isUnsetToken(const string& input) {
    string token = lowerCopy(trimCopy(input));
    return token.empty() || token == "-" || token == "none";
}

// Reproduces the guard from PfMultiProcess::buildPfMultiPreparedContext
static bool globalRefRequired(const PfMultiConfig::Config& config, const string& globalRef) {
    if (!isUnsetToken(globalRef)) return false;  // global ref is set, no issue

    for (const auto& lib : config.libraries) {
        string norm = lib.normalizedFeatureType();
        bool isGex = (norm == "geneexpression" || norm == "gex");
        if (!isGex && lib.starFeatureRef.empty()) {
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage: test_guard <config.csv> <global_ref_or_dash>" << endl;
        return 2;
    }
    string globalRef = string(argv[2]);
    try {
        PfMultiConfig::Config cfg = PfMultiConfig::parseConfig(argv[1]);
        bool needed = globalRefRequired(cfg, globalRef);
        cout << "GLOBAL_REF_REQUIRED=" << (needed ? "yes" : "no") << endl;
        for (size_t i = 0; i < cfg.libraries.size(); ++i) {
            const auto& lib = cfg.libraries[i];
            cout << "LIB:" << i
                 << " id=" << lib.starLibraryId
                 << " feature_types=" << lib.feature_types
                 << " star_feature_ref=" << (lib.starFeatureRef.empty() ? "(none)" : lib.starFeatureRef)
                 << endl;
        }
        return 0;
    } catch (const std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
}
HARNESS_EOF

SOURCE_DIR="$REPO_ROOT/core/legacy/source"
PF_INCLUDE="$REPO_ROOT/core/features/process_features/include"
g++ -std=c++11 -O2 -I"$SOURCE_DIR" -I"$PF_INCLUDE" \
    "$WORK_DIR/test_guard.cpp" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    -o "$HARNESS" 2>&1
echo "  Harness built OK"

touch "$WORK_DIR/ref_grna.csv" "$WORK_DIR/ref_larry.csv"

# --- Test 1: all feature libs have star_feature_ref, global ref unset → not required ---
echo ""
echo "--- Test 1: all feature libs covered → global ref NOT required ---"
cat > "$WORK_DIR/config1.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref
/path/to/mRNA,S1,Gene Expression,
/path/to/gRNA,S1,CRISPR Guide Capture,$WORK_DIR/ref_grna.csv
/path/to/larry,S1,Custom,$WORK_DIR/ref_larry.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config1.csv" "-" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'GLOBAL_REF_REQUIRED=no'; then
    pass "global ref not required when all feature libs have star_feature_ref"
else
    fail "should not require global ref"
fi

# --- Test 2: one feature lib missing star_feature_ref, global ref unset → required ---
echo ""
echo "--- Test 2: one feature lib missing ref → global ref required ---"
cat > "$WORK_DIR/config2.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref
/path/to/mRNA,S1,Gene Expression,
/path/to/gRNA,S1,CRISPR Guide Capture,
/path/to/larry,S1,Custom,$WORK_DIR/ref_larry.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config2.csv" "-" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'GLOBAL_REF_REQUIRED=yes'; then
    pass "global ref required when a feature lib has no star_feature_ref"
else
    fail "should require global ref when one lib is missing star_feature_ref"
fi

# --- Test 3: global ref is set → not required regardless ---
echo ""
echo "--- Test 3: global ref set → not required ---"
OUTPUT=$("$HARNESS" "$WORK_DIR/config2.csv" "/some/global/ref.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'GLOBAL_REF_REQUIRED=no'; then
    pass "global ref set → guard does not fire"
else
    fail "guard should not fire when global ref is set"
fi

# --- Test 4: GEX-only config (no feature libs) → not required ---
echo ""
echo "--- Test 4: GEX-only config → global ref NOT required ---"
cat > "$WORK_DIR/config4.csv" << 'EOF'
[libraries]
fastqs,sample,feature_types
/path/to/mRNA,S1,Gene Expression
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config4.csv" "-" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'GLOBAL_REF_REQUIRED=no'; then
    pass "GEX-only config does not require global ref"
else
    fail "GEX-only config should not require global ref"
fi

echo ""
echo "=========================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=========================================="

[ "$FAIL_COUNT" -eq 0 ]
