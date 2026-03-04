#!/usr/bin/env bash
# End-to-end test for feature_barcodes --allow_union_whitelist
#
# Verifies that the preloaded-hash path in main.c performs union expansion
# when the flag is set, and that exact-only matching drops opposite-namespace
# barcodes when the flag is absent.
#
# Barcode namespace:
#   TRU: AAACCCAAGAAACCAT   pos7=A pos8=G
#   NXT: AAACCCATCAAACCAT   pos7=T pos8=C  (complement of A,G)
#
# The test creates a filtered-barcode file in NXT namespace, FASTQ reads
# with TRU-namespace barcodes, and a trivial feature list.  Without union
# expansion the TRU barcode will not match the NXT-only hash and the
# filtered/ output will be empty.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOOL="${SCRIPT_DIR}/../assignBarcodes"
TMPDIR_BASE=$(mktemp -d /tmp/test_union_wl_e2e.XXXXXX)
trap 'rm -rf "$TMPDIR_BASE"' EXIT

fail() { echo "FAIL: $1" >&2; exit 1; }

# ---- fixture constants ----
BC_TRU="AAACCCAAGAAACCAT"    # TRU form (pos7=A, pos8=G)
BC_NXT="AAACCCATCAAACCAT"    # NXT form (pos7=T, pos8=C)
UMI="AAAAAAAAAAAA"           # 12-base UMI
FEATURE_SEQ="ATCGATCGATCGATCG"
QUALS="IIIIIIIIIIIIIIIIIIIIIIIIIIII"

# ---- create fixtures ----
FIXTURE_DIR="$TMPDIR_BASE/fixtures"
mkdir -p "$FIXTURE_DIR"

echo "$BC_TRU" > "$FIXTURE_DIR/whitelist.txt"

cat > "$FIXTURE_DIR/features.csv" <<EOF
name,sequence,type
TestGuide,$FEATURE_SEQ,CRISPR Guide Capture
EOF

# Filtered barcodes: NXT form (opposite namespace from reads)
echo "$BC_NXT" > "$FIXTURE_DIR/filtered_barcodes.txt"

cat > "$FIXTURE_DIR/sample_R1_001.fastq" <<EOF
@read_0000
${BC_TRU}${UMI}
+
${QUALS}
@read_0001
${BC_TRU}${UMI}
+
${QUALS}
EOF

cat > "$FIXTURE_DIR/sample_R2_001.fastq" <<EOF
@read_0000
${FEATURE_SEQ}AAAAAAAAAAAA
+
${QUALS}
@read_0001
${FEATURE_SEQ}AAAAAAAAAAAA
+
${QUALS}
EOF

# ---- common args ----
COMMON_ARGS=(
    -w "$FIXTURE_DIR/whitelist.txt"
    -f "$FIXTURE_DIR/features.csv"
    --filtered_barcodes "$FIXTURE_DIR/filtered_barcodes.txt"
    --feature_constant_offset 0
    -b 16
    -u 12
    -t 1
    -S 1
)

# ====================================================================
# Test 1: WITHOUT --allow_union_whitelist
#   Filtered hash has NXT barcodes only (exact-only).
#   Reads carry TRU barcodes → no match → filtered/barcodes.txt empty.
# ====================================================================
echo "=== Test 1: without --allow_union_whitelist ==="

OUT1="$TMPDIR_BASE/out_no_union"
mkdir -p "$OUT1"

"$TOOL" "${COMMON_ARGS[@]}" -d "$OUT1" "$FIXTURE_DIR" >/dev/null 2>&1 || true

FILTERED_BC="$OUT1/fixtures/filtered/barcodes.txt"
if [ ! -f "$FILTERED_BC" ]; then
    fail "Test 1: filtered/barcodes.txt not produced"
fi
BC_COUNT=$(wc -l < "$FILTERED_BC" | tr -d ' ')
if [ "$BC_COUNT" -ne 0 ]; then
    fail "Test 1: expected 0 barcodes in filtered output, got $BC_COUNT"
fi
echo "  PASS: 0 barcodes in filtered output without --allow_union_whitelist"

# ====================================================================
# Test 2: WITH --allow_union_whitelist
#   Filtered hash is expanded to contain both NXT and TRU forms.
#   TRU barcodes from reads now match → filtered/barcodes.txt non-empty.
# ====================================================================
echo "=== Test 2: with --allow_union_whitelist ==="

OUT2="$TMPDIR_BASE/out_union"
mkdir -p "$OUT2"

"$TOOL" "${COMMON_ARGS[@]}" -d "$OUT2" --allow_union_whitelist "$FIXTURE_DIR" 2>"$OUT2/stderr.log"

FILTERED_BC2="$OUT2/fixtures/filtered/barcodes.txt"
if [ ! -f "$FILTERED_BC2" ]; then
    cat "$OUT2/stderr.log" >&2
    fail "Test 2: filtered/barcodes.txt not produced"
fi
BC_COUNT2=$(wc -l < "$FILTERED_BC2" | tr -d ' ')
if [ "$BC_COUNT2" -eq 0 ]; then
    cat "$OUT2/stderr.log" >&2
    fail "Test 2: expected >0 barcodes in filtered output, got 0"
fi
echo "  PASS: $BC_COUNT2 barcode(s) in filtered output with --allow_union_whitelist"

# ====================================================================
# Test 3: verify expansion message on stderr
# ====================================================================
echo "=== Test 3: expansion notice on stderr ==="
if ! grep -q "allow_union_whitelist active" "$OUT2/stderr.log"; then
    cat "$OUT2/stderr.log" >&2
    fail "Test 3: expected expansion notice on stderr"
fi
echo "  PASS: expansion notice present"

# ====================================================================
# Test 4: unfiltered output identical in both runs
#   Union flag should only affect the filtered path.
# ====================================================================
echo "=== Test 4: unfiltered output unchanged ==="
UNFILT1="$OUT1/fixtures/barcodes.txt"
UNFILT2="$OUT2/fixtures/barcodes.txt"
if [ ! -f "$UNFILT1" ] || [ ! -f "$UNFILT2" ]; then
    fail "Test 4: unfiltered barcodes.txt missing"
fi
if ! diff -q "$UNFILT1" "$UNFILT2" >/dev/null 2>&1; then
    fail "Test 4: unfiltered barcodes.txt differs between union and non-union runs"
fi
echo "  PASS: unfiltered barcodes.txt identical"

echo ""
echo "ALL TESTS PASSED"
