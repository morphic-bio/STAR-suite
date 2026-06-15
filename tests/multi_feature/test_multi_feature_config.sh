#!/usr/bin/env bash
# Test: pfMultiConfig star_feature_ref / star_library_id column parsing
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0

pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

echo "=== Test: multi-feature config columns ==="

# Build test harness (extended to output star_feature_ref and star_library_id)
HARNESS="$WORK_DIR/test_config_multi"
cat > "$WORK_DIR/test_config_multi.cpp" << 'HARNESS_EOF'
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
        cerr << "Usage: test_config_multi <config.csv> [expect_fail]" << endl;
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
                 << " star_feature_ref=" << lib.starFeatureRef
                 << " star_library_id=" << lib.starLibraryId
                 << " star_input_format=" << lib.starInputFormat
                 << " star_hash_demux=" << lib.starHashDemux
                 << " star_hash_feature_selector=" << lib.starHashFeatureSelector
                 << " star_hash_demux_method=" << lib.starHashDemuxMethod
                 << " star_hash_min_total=" << lib.starHashMinTotal
                 << " star_hash_min_top=" << lib.starHashMinTop
                 << " star_hash_min_ratio=" << lib.starHashMinRatio
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
PF_INCLUDE="$REPO_ROOT/core/features/process_features/include"
g++ -std=c++11 -O2 -I"$SOURCE_DIR" -I"$PF_INCLUDE" \
    "$WORK_DIR/test_config_multi.cpp" \
    "$SOURCE_DIR/PfMultiConfig.o" \
    -o "$HARNESS" 2>&1
if [ $? -ne 0 ]; then
    echo "FATAL: failed to build test harness"
    exit 1
fi
echo "  Harness built OK"

# Create dummy feature ref files for path validation
touch "$WORK_DIR/ref_grna.csv"
touch "$WORK_DIR/ref_larry.csv"

# --- Test 1: star_feature_ref per library ---
echo ""
echo "--- Test 1: star_feature_ref per library ---"
cat > "$WORK_DIR/config1.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref
/path/to/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,
/path/to/gRNA,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,$WORK_DIR/ref_grna.csv
/path/to/larry,DE_30KO,Custom,Custom,TRU,$WORK_DIR/ref_larry.csv

[feature]
ref,$WORK_DIR/ref_grna.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config1.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q "LIB:0.*star_feature_ref= "; then
    pass "GEX row has empty star_feature_ref"
else
    fail "GEX row should have empty star_feature_ref"
fi

if echo "$OUTPUT" | grep -q "LIB:1.*star_feature_ref=$WORK_DIR/ref_grna.csv"; then
    pass "gRNA row has per-library feature ref"
else
    fail "gRNA row should have per-library feature ref"
fi

if echo "$OUTPUT" | grep -q "LIB:2.*star_feature_ref=$WORK_DIR/ref_larry.csv"; then
    pass "LARRY row has per-library feature ref"
else
    fail "LARRY row should have per-library feature ref"
fi

# --- Test 2: star_library_id explicit ---
echo ""
echo "--- Test 2: explicit star_library_id ---"
cat > "$WORK_DIR/config2.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_library_id
/path/to/mRNA,S1,Gene Expression,Gene Expression,gex_main
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture,grna_pool
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config2.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_library_id=gex_main'; then
    pass "explicit star_library_id=gex_main"
else
    fail "should preserve explicit star_library_id"
fi

if echo "$OUTPUT" | grep -q 'LIB:1.*star_library_id=grna_pool'; then
    pass "explicit star_library_id=grna_pool"
else
    fail "should preserve explicit star_library_id for gRNA"
fi

# --- Test 3: star_library_id auto-generated when absent ---
echo ""
echo "--- Test 3: auto-generated star_library_id ---"
cat > "$WORK_DIR/config3.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types
/path/to/mRNA,S1,Gene Expression,Gene Expression
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config3.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_library_id=S1_Gene_Expression_0'; then
    pass "auto-generated ID for GEX: S1_Gene_Expression_0"
else
    fail "auto-generated ID should be S1_Gene_Expression_0"
fi

if echo "$OUTPUT" | grep -q 'LIB:1.*star_library_id=S1_CRISPR_Guide_Capture_1'; then
    pass "auto-generated ID for gRNA: S1_CRISPR_Guide_Capture_1"
else
    fail "auto-generated ID should be S1_CRISPR_Guide_Capture_1"
fi

# --- Test 4: duplicate star_library_id fails ---
echo ""
echo "--- Test 4: duplicate star_library_id rejected ---"
cat > "$WORK_DIR/config4.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,star_library_id
/path/to/mRNA,S1,Gene Expression,Gene Expression,same_id
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture,same_id
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config4.csv" "expect_fail" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'EXPECTED_FAIL.*Duplicate star_library_id'; then
    pass "duplicate star_library_id rejected"
else
    fail "duplicate star_library_id should be rejected"
fi

# --- Test 5: invalid star_feature_ref path fails ---
echo ""
echo "--- Test 5: nonexistent star_feature_ref rejected ---"
cat > "$WORK_DIR/config5.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types,star_feature_ref
/path/to/mRNA,S1,Gene Expression,Gene Expression,/nonexistent/path/ref.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config5.csv" "expect_fail" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'EXPECTED_FAIL.*star_feature_ref path does not exist'; then
    pass "nonexistent star_feature_ref rejected"
else
    fail "nonexistent star_feature_ref should be rejected"
fi

# --- Test 6: backward compat (no new columns) ---
echo ""
echo "--- Test 6: backward compatibility (no new columns) ---"
cat > "$WORK_DIR/config6.csv" << 'EOF'
[libraries]
fastqs,sample,library_type,feature_types
/path/to/mRNA,S1,Gene Expression,Gene Expression
/path/to/gRNA,S1,CRISPR Guide Capture,CRISPR Guide Capture
/path/to/larry,S1,Custom,Custom
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config6.csv" 2>&1)
echo "$OUTPUT"

NLIBS=$(echo "$OUTPUT" | grep -c '^LIB:')
if [ "$NLIBS" -eq 3 ]; then
    pass "all 3 libraries parsed without new columns"
else
    fail "expected 3 libraries, got $NLIBS"
fi

if echo "$OUTPUT" | grep -q 'LIB:0.*star_feature_ref= '; then
    pass "no star_feature_ref column → empty"
else
    fail "star_feature_ref should be empty when column absent"
fi

# --- Test 7: column aliases (starfeatureref, starlibraryid) ---
echo ""
echo "--- Test 7: column aliases ---"
cat > "$WORK_DIR/config7.csv" << EOF
[libraries]
fastqs,sample,feature_types,starfeatureref,starlibraryid
/path/to/mRNA,S1,Gene Expression,,my_gex
/path/to/gRNA,S1,CRISPR Guide Capture,$WORK_DIR/ref_grna.csv,my_grna
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config7.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q "LIB:1.*star_feature_ref=$WORK_DIR/ref_grna.csv"; then
    pass "starfeatureref alias works"
else
    fail "starfeatureref alias should work"
fi

if echo "$OUTPUT" | grep -q 'LIB:0.*star_library_id=my_gex'; then
    pass "starlibraryid alias works"
else
    fail "starlibraryid alias should work"
fi

# --- Test 8: trailing fields omitted with new columns ---
echo ""
echo "--- Test 8: trailing field padding with new columns ---"
cat > "$WORK_DIR/config8.csv" << 'EOF'
[libraries]
fastqs,sample,feature_types,star_feature_ref,star_library_id
/path/to/mRNA,S1,Gene Expression
/path/to/gRNA,S1,CRISPR Guide Capture
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config8.csv" 2>&1)
echo "$OUTPUT"

NLIBS=$(echo "$OUTPUT" | grep -c '^LIB:')
if [ "$NLIBS" -eq 2 ]; then
    pass "rows with omitted trailing fields not dropped ($NLIBS libs)"
else
    fail "expected 2 libraries with trailing padding, got $NLIBS"
fi

if echo "$OUTPUT" | grep -q 'LIB:0.*star_feature_ref=.*star_library_id='; then
    pass "padded trailing fields are empty"
else
    fail "padded trailing fields should be empty"
fi

# --- Test 9: all columns together ---
echo ""
echo "--- Test 9: full multi-feature config ---"
cat > "$WORK_DIR/config9.csv" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
/path/to/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,gex_de
/path/to/gRNA,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,$WORK_DIR/ref_grna.csv,grna_de
/path/to/larry,DE_30KO,Custom,Custom,TRU,$WORK_DIR/ref_larry.csv,larry_de

[feature]
ref,$WORK_DIR/ref_grna.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config9.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'LIB:0.*star_chemistry=tru.*star_feature_ref=.*star_library_id=gex_de'; then
    pass "GEX row: all fields correct"
else
    fail "GEX row should have all fields"
fi

if echo "$OUTPUT" | grep -q "LIB:1.*star_chemistry=nxt.*star_feature_ref=$WORK_DIR/ref_grna.csv.*star_library_id=grna_de"; then
    pass "gRNA row: all fields correct"
else
    fail "gRNA row should have all fields"
fi

if echo "$OUTPUT" | grep -q "LIB:2.*star_chemistry=tru.*star_feature_ref=$WORK_DIR/ref_larry.csv.*star_library_id=larry_de"; then
    pass "LARRY row: all fields correct"
else
    fail "LARRY row should have all fields"
fi

# --- Test 10: relative star_feature_ref resolved against config directory ---
echo ""
echo "--- Test 10: relative star_feature_ref resolution ---"
mkdir -p "$WORK_DIR/subdir"
touch "$WORK_DIR/subdir/my_ref.csv"
cat > "$WORK_DIR/subdir/config_rel.csv" << 'EOF'
[libraries]
fastqs,sample,feature_types,star_feature_ref
/path/to/mRNA,S1,Gene Expression,
/path/to/gRNA,S1,CRISPR Guide Capture,my_ref.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/subdir/config_rel.csv" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q "LIB:1.*star_feature_ref=$WORK_DIR/subdir/my_ref.csv"; then
    pass "relative path resolved against config directory"
else
    fail "relative star_feature_ref should resolve against config dir"
fi

# --- Test 11: all feature libs have star_feature_ref, no global ref needed ---
echo ""
echo "--- Test 11: no global ref when all libs have star_feature_ref ---"
cat > "$WORK_DIR/config_noref.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref
/path/to/mRNA,S1,Gene Expression,
/path/to/gRNA,S1,CRISPR Guide Capture,$WORK_DIR/ref_grna.csv
/path/to/larry,S1,Custom,$WORK_DIR/ref_larry.csv
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_noref.csv" 2>&1)
echo "$OUTPUT"

NLIBS=$(echo "$OUTPUT" | grep -c '^LIB:')
if [ "$NLIBS" -eq 3 ]; then
    pass "config parses OK without global [feature] ref"
else
    fail "config without global ref should parse when all feature libs have star_feature_ref"
fi

# --- Test 12: sanitized library ID collision rejected ---
echo ""
echo "--- Test 12: sanitized star_library_id collision rejected ---"
cat > "$WORK_DIR/config_sanitize_dup.csv" << 'EOF'
[libraries]
fastqs,sample,feature_types,star_library_id
/path/to/gRNA,S1,CRISPR Guide Capture,lib:a
/path/to/larry,S1,Custom,lib/a
EOF

OUTPUT=$("$HARNESS" "$WORK_DIR/config_sanitize_dup.csv" "expect_fail" 2>&1)
echo "$OUTPUT"

if echo "$OUTPUT" | grep -q 'EXPECTED_FAIL.*collides.*after path sanitization.*lib_a'; then
    pass "sanitized collision (lib:a vs lib/a → lib_a) rejected"
else
    fail "sanitized ID collision should be rejected"
fi

# --- Test 13: default star_input_format is fastq (empty) ---
echo ""
echo "--- Test 13: default star_input_format ---"
cat > "$WORK_DIR/config_fmt_default.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_feature_ref
/path/to/gRNA,S1,CRISPR Guide Capture,$WORK_DIR/ref_grna.csv
EOF
OUTPUT=$("$HARNESS" "$WORK_DIR/config_fmt_default.csv" 2>&1)
if echo "$OUTPUT" | grep -q 'star_input_format= star_hash'; then
    pass "default star_input_format is empty (fastq)"
else
    fail "default star_input_format should be empty"
fi

# --- Test 14: explicit fastq ---
echo ""
echo "--- Test 14: explicit star_input_format=fastq ---"
cat > "$WORK_DIR/config_fmt_fastq.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_input_format,star_feature_ref
/path/to/gRNA,S1,CRISPR Guide Capture,fastq,$WORK_DIR/ref_grna.csv
EOF
OUTPUT=$("$HARNESS" "$WORK_DIR/config_fmt_fastq.csv" 2>&1)
if echo "$OUTPUT" | grep -q 'star_input_format=fastq'; then
    pass "explicit star_input_format=fastq accepted"
else
    fail "explicit fastq should parse"
fi

# --- Test 15: table format ---
echo ""
echo "--- Test 15: star_input_format=table ---"
cat > "$WORK_DIR/hiv_counts.tsv" << 'EOF'
barcode	feature_id	count
GTATGTTCAGTAGCCT-1	HIV_DNA	5
GTATGTTCAGTAGCCT-1	HIV_RNA	2
EOF
cat > "$WORK_DIR/config_fmt_table.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_input_format,star_feature_ref,star_library_id
$WORK_DIR/hiv_counts.tsv,YW8,Custom,table,$WORK_DIR/ref_larry.csv,hiv_state_yw8
EOF
OUTPUT=$("$HARNESS" "$WORK_DIR/config_fmt_table.csv" 2>&1)
if echo "$OUTPUT" | grep -q 'star_input_format=table'; then
    pass "star_input_format=table accepted"
else
    fail "table format should parse"
fi

# --- Test 16: invalid star_input_format ---
echo ""
echo "--- Test 16: invalid star_input_format rejected ---"
cat > "$WORK_DIR/config_fmt_bad.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_input_format,star_feature_ref
/path/to/gRNA,S1,Custom,counts,$WORK_DIR/ref_larry.csv
EOF
OUTPUT=$("$HARNESS" "$WORK_DIR/config_fmt_bad.csv" "expect_fail" 2>&1)
if echo "$OUTPUT" | grep -q 'EXPECTED_FAIL.*Invalid star_input_format'; then
    pass "invalid star_input_format rejected"
else
    fail "invalid star_input_format should fail"
fi

# --- Test 17: table without star_feature_ref on non-GEX ---
echo ""
echo "--- Test 17: table non-GEX requires star_feature_ref ---"
cat > "$WORK_DIR/config_table_noref.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_input_format
$WORK_DIR/hiv_counts.tsv,YW8,Custom,table
EOF
OUTPUT=$("$HARNESS" "$WORK_DIR/config_table_noref.csv" "expect_fail" 2>&1)
if echo "$OUTPUT" | grep -q 'EXPECTED_FAIL.*star_feature_ref is required'; then
    pass "table non-GEX without star_feature_ref rejected"
else
    fail "table non-GEX should require star_feature_ref"
fi

# --- Test 18: star_hash_* columns ---
echo ""
echo "--- Test 18: star_hash_* demux columns ---"
cat > "$WORK_DIR/config_hash.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_library_id,star_hash_demux,star_hash_feature_selector,star_hash_demux_method,star_hash_min_total,star_hash_min_top,star_hash_min_ratio
/path/to/adt,S1,Protein,adt_s1,auto,id_prefix:hashtag,ratio,3,3,2.0

[feature]
ref,$WORK_DIR/ref_grna.csv
EOF
OUTPUT=$("$HARNESS" "$WORK_DIR/config_hash.csv" 2>&1)
if echo "$OUTPUT" | grep -q 'star_hash_demux=auto'; then
    pass "star_hash_demux parsed"
else
    fail "star_hash_demux should parse"
fi
if echo "$OUTPUT" | grep -q 'star_hash_feature_selector=id_prefix:hashtag'; then
    pass "star_hash_feature_selector parsed"
else
    fail "star_hash_feature_selector should parse"
fi
if echo "$OUTPUT" | grep -q 'star_hash_min_ratio=2'; then
    pass "star_hash_min_ratio parsed"
else
    fail "star_hash_min_ratio should parse"
fi

# --- Test 19: invalid star_hash_demux ---
echo ""
echo "--- Test 19: invalid star_hash_demux rejected ---"
cat > "$WORK_DIR/config_hash_bad.csv" << EOF
[libraries]
fastqs,sample,feature_types,star_library_id,star_hash_demux,star_feature_ref
/path/to/adt,S1,Protein,bad_hash,ye,$WORK_DIR/ref_grna.csv

[feature]
ref,$WORK_DIR/ref_grna.csv
EOF
OUTPUT=$("$HARNESS" "$WORK_DIR/config_hash_bad.csv" 2>&1 || true)
if echo "$OUTPUT" | grep -q 'Invalid star_hash_demux'; then
    pass "invalid star_hash_demux rejected"
else
    fail "invalid star_hash_demux should fail config validation"
fi

echo ""
echo "=========================================="
echo "Results: $PASS_COUNT passed, $FAIL_COUNT failed"
echo "=========================================="

[ "$FAIL_COUNT" -eq 0 ]
