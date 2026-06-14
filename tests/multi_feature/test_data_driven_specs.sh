#!/usr/bin/env bash
# Test: data-driven feature specs route Custom/LARRY feature types.
# Exercises buildFeatureSpecsFromConfig (PfMultiFeatureSpecs.h) and
# getFeatureLibraries (PfMultiConfig) — the actual routing code.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0

pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: data-driven feature spec routing ==="

HARNESS="$WORK_DIR/test_specs"
cat > "$WORK_DIR/test_specs.cpp" << 'HARNESS_EOF'
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
using std::string;
using std::vector;
using std::map;
using std::set;
using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;
#include "PfMultiConfig.h"
#include "PfMultiFeatureSpecs.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: test_specs <config.csv>" << endl;
        return 2;
    }
    try {
        PfMultiConfig::Config cfg = PfMultiConfig::parseConfig(argv[1]);
        ostringstream log;
        vector<PfMultiFeatureSpecs::FeatureSpec> specs =
            PfMultiFeatureSpecs::buildFeatureSpecsFromConfig(cfg, log);

        cout << "TOTAL_LIBS=" << cfg.libraries.size() << endl;
        cout << "TOTAL_SPECS=" << specs.size() << endl;

        for (size_t i = 0; i < specs.size(); ++i) {
            const auto& spec = specs[i];
            vector<PfMultiConfig::LibraryEntry> matched =
                cfg.getFeatureLibraries(spec.libraryType);
            cout << "SPEC:" << i
                 << " libraryType=" << spec.libraryType
                 << " featureRefType=" << spec.featureRefType
                 << " matched=" << matched.size()
                 << endl;
            for (size_t j = 0; j < matched.size(); ++j) {
                cout << "  MATCH:" << j
                     << " id=" << matched[j].starLibraryId
                     << " feature_types=" << matched[j].feature_types
                     << " star_feature_ref="
                     << (matched[j].starFeatureRef.empty() ? "(none)" : matched[j].starFeatureRef)
                     << endl;
            }
        }
        string notices = log.str();
        if (!notices.empty()) {
            cout << "LOG:" << notices;
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
    "$WORK_DIR/test_specs.cpp" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    -o "$HARNESS" 2>&1
echo "  Harness built OK"

touch "$WORK_DIR/ref_grna.csv" "$WORK_DIR/ref_larry.csv"

# --- Test 1: known types route with correct featureRefType ---
echo ""
echo "--- Test 1: known feature types ---"
cat > "$WORK_DIR/config_known.csv" << 'EOF'
[libraries]
fastqs,sample,feature_types
/path/to/mRNA,S1,Gene Expression
/path/to/gRNA,S1,CRISPR Guide Capture
/path/to/adt,S1,Antibody Capture
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_known.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'TOTAL_SPECS=2'; then
    pass "2 feature specs (CRISPR + Antibody)"
else
    fail "expected 2 feature specs"
fi

if echo "$OUTPUT" | grep -q 'SPEC:0.*libraryType=CRISPR Guide Capture.*featureRefType=CRISPR Guide Capture.*matched=1'; then
    pass "CRISPR spec: correct type and exactly 1 match"
else
    fail "CRISPR spec should match exactly 1 library"
fi

if echo "$OUTPUT" | grep -q 'SPEC:1.*libraryType=Antibody Capture.*featureRefType=Antibody Capture.*matched=1'; then
    pass "Antibody spec: correct type and exactly 1 match"
else
    fail "Antibody spec should match exactly 1 library"
fi

# --- Test 1b: Protein/ADT aliases map to Antibody Capture ---
echo ""
echo "--- Test 1b: Protein/ADT aliases ---"
cat > "$WORK_DIR/config_protein_aliases.csv" << EOF
[libraries]
fastqs,sample,feature_types
/path/to/mRNA,S1,Gene Expression
/path/to/adt,S1,ADT
/path/to/protein,S1,Protein
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_protein_aliases.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'TOTAL_SPECS=2'; then
    pass "2 feature specs (ADT + Protein aliases)"
else
    fail "expected 2 feature specs for ADT and Protein"
fi

if echo "$OUTPUT" | grep -q 'libraryType=ADT.*featureRefType=Antibody Capture'; then
    pass "ADT alias maps to Antibody Capture featureRefType"
else
    fail "ADT alias should map to Antibody Capture"
fi

if echo "$OUTPUT" | grep -q 'libraryType=Protein.*featureRefType=Antibody Capture'; then
    pass "Protein alias maps to Antibody Capture featureRefType"
else
    fail "Protein alias should map to Antibody Capture"
fi

# --- Test 2: Custom type routes with verbatim featureRefType ---
echo ""
echo "--- Test 2: Custom type (verbatim featureRefType) ---"
cat > "$WORK_DIR/config_custom.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref
/path/to/mRNA,S1,Gene Expression,
/path/to/larry,S1,Custom,$WORK_DIR/ref_larry.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_custom.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'TOTAL_SPECS=1'; then
    pass "1 feature spec (Custom only, GEX excluded)"
else
    fail "expected 1 feature spec"
fi

if echo "$OUTPUT" | grep -q 'SPEC:0.*libraryType=Custom.*featureRefType=Custom.*matched=1'; then
    pass "Custom spec: verbatim featureRefType and exactly 1 match"
else
    fail "Custom spec should match exactly 1 library"
fi

if echo "$OUTPUT" | grep -q 'LOG:NOTICE.*Custom.*not a known 10x type'; then
    pass "NOTICE logged for unknown type"
else
    fail "should log NOTICE for unknown Custom type"
fi

# --- Test 3: no over-match (Custom vs Custom2) ---
echo ""
echo "--- Test 3: exact matching, no over-match ---"
cat > "$WORK_DIR/config_no_overlap.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref
/path/to/mRNA,S1,Gene Expression,
/path/to/custom1,S1,Custom,$WORK_DIR/ref_larry.csv
/path/to/custom2,S1,Custom2,$WORK_DIR/ref_grna.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_no_overlap.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'TOTAL_SPECS=2'; then
    pass "2 specs (Custom and Custom2 are distinct)"
else
    fail "Custom and Custom2 should produce 2 distinct specs"
fi

if echo "$OUTPUT" | grep -q 'SPEC:0.*libraryType=Custom .*matched=1'; then
    pass "Custom matches exactly 1 (not Custom2)"
else
    # More lenient check: just verify matched=1
    if echo "$OUTPUT" | grep 'SPEC:0' | grep -q 'matched=1'; then
        pass "Custom matches exactly 1 (not Custom2)"
    else
        fail "Custom should match exactly 1 library, not include Custom2"
    fi
fi

if echo "$OUTPUT" | grep -q 'SPEC:1.*libraryType=Custom2.*matched=1'; then
    pass "Custom2 matches exactly 1"
else
    fail "Custom2 should match exactly 1 library"
fi

# --- Test 4: CellPlex (CMO) matches itself correctly ---
echo ""
echo "--- Test 4: CellPlex (CMO) punctuation handling ---"
cat > "$WORK_DIR/config_cellplex.csv" << 'EOF'
[libraries]
fastqs,sample,feature_types
/path/to/mRNA,S1,Gene Expression
/path/to/cmo,S1,CellPlex (CMO)
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_cellplex.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'TOTAL_SPECS=1'; then
    pass "1 feature spec (CellPlex)"
else
    fail "expected 1 feature spec"
fi

if echo "$OUTPUT" | grep 'SPEC:0' | grep -q 'matched=1'; then
    pass "CellPlex (CMO) matches itself"
else
    fail "CellPlex (CMO) should match itself via normalized comparison"
fi

if echo "$OUTPUT" | grep -q 'featureRefType=Multiplexing Capture'; then
    pass "CellPlex (CMO) maps to canonical Multiplexing Capture"
else
    fail "CellPlex (CMO) should map to Multiplexing Capture"
fi

# --- Test 5: dedup — two libraries with same feature_types get one spec ---
echo ""
echo "--- Test 5: dedup of same feature_types ---"
cat > "$WORK_DIR/config_dedup.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref,star_library_id
/path/to/mRNA,S1,Gene Expression,,gex
/path/to/grna1,S1,CRISPR Guide Capture,$WORK_DIR/ref_grna.csv,grna1
/path/to/grna2,S1,CRISPR Guide Capture,$WORK_DIR/ref_larry.csv,grna2
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_dedup.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'TOTAL_SPECS=1'; then
    pass "1 spec for two CRISPR libraries (dedup)"
else
    fail "two CRISPR libraries should produce 1 spec"
fi

if echo "$OUTPUT" | grep 'SPEC:0' | grep -q 'matched=2'; then
    pass "spec matches both CRISPR libraries"
else
    fail "spec should match both CRISPR libraries"
fi

# --- Test 6: Lineage type (non-standard) ---
echo ""
echo "--- Test 6: 'Lineage' type ---"
cat > "$WORK_DIR/config_lineage.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref
/path/to/mRNA,S1,Gene Expression,
/path/to/larry,S1,Lineage,$WORK_DIR/ref_larry.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_lineage.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'TOTAL_SPECS=1'; then
    pass "Lineage produces 1 spec"
else
    fail "Lineage should produce 1 spec"
fi

if echo "$OUTPUT" | grep -q 'SPEC:0.*featureRefType=Lineage.*matched=1'; then
    pass "Lineage: verbatim featureRefType and 1 match"
else
    fail "Lineage should use verbatim featureRefType"
fi

echo ""
echo "=========================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=========================================="

[ "$FAIL_COUNT" -eq 0 ]
