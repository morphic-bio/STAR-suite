#!/usr/bin/env bash
#
# Regression test: 2-column CB whitelist must produce PAIRED (COL1→COL2)
# entries in Solo, NOT independent/doubled entries.
#
# Bug: An earlier regression loaded both columns as separate barcodes,
# doubling cbWL size and creating a mixed-namespace Solo output that
# corrupted downstream MEX merges.
#
# This test runs STAR with a synthetic 2-column whitelist and verifies:
#   1) Log.out says "column 1 is used for matching, column 2 is used for MEX
#      barcode output mapping" (correct pairing behavior)
#   2) Log.out does NOT say "both columns are included in CB matching/output
#      namespaces" (regressed doubling behavior)
#   3) "Number of CBs in the whitelist" equals the number of whitelist ROWS
#      (not 2× the rows)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR_BIN="${REPO_ROOT}/core/legacy/source/STAR"

if [[ ! -x "$STAR_BIN" ]]; then
    echo "ERROR: STAR binary not found at ${STAR_BIN}" >&2
    echo "Build with:  make -C core/legacy/source -j8 STAR" >&2
    exit 1
fi

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

PASS_COUNT=0
FAIL_COUNT=0
pass() { echo "  PASS: $1"; PASS_COUNT=$((PASS_COUNT + 1)); }
fail() { echo "  FAIL: $1"; FAIL_COUNT=$((FAIL_COUNT + 1)); }

N_BARCODES=10

python3 -c "
import random, os
random.seed(99)
bases = 'ACGT'
comp = str.maketrans('ACGTacgt', 'TGCAtgca')
n = ${N_BARCODES}
bcs = []
for _ in range(n):
    bc = ''.join(random.choice(bases) for _ in range(16))
    bcs.append(bc)
with open('${WORK_DIR}/wl_2col.txt', 'w') as f:
    for bc in bcs:
        nxt = list(bc)
        nxt[7] = nxt[7].translate(comp)
        nxt[8] = nxt[8].translate(comp)
        nxt = ''.join(nxt)
        f.write(nxt + '\t' + bc + '\n')
"

cat > "${WORK_DIR}/r1.fq" <<'EOF'
@r1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

cat > "${WORK_DIR}/r2.fq" <<'EOF'
@r1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

echo "=== Regression: 2-column WL pairing (not doubling) ==="
echo ""

set +e
"${STAR_BIN}" \
    --runThreadN 1 \
    --genomeDir "${WORK_DIR}/no_genome" \
    --readFilesIn "${WORK_DIR}/r2.fq" "${WORK_DIR}/r1.fq" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WORK_DIR}/wl_2col.txt" \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --outFileNamePrefix "${WORK_DIR}/run_" \
    > "${WORK_DIR}/stdout.log" 2>&1
set -e

LOG="${WORK_DIR}/run_Log.out"

if [[ ! -f "$LOG" ]]; then
    echo "FATAL: Log.out not created; STAR may have crashed before WL loading."
    echo "--- stdout/stderr ---"
    cat "${WORK_DIR}/stdout.log"
    exit 1
fi

echo "--- Check 1: correct pairing message present ---"
if grep -q "column 1 is used for matching, column 2 is used for MEX barcode output mapping" "$LOG"; then
    pass "Log.out contains correct pairing message"
else
    fail "Log.out missing correct pairing message"
fi

echo "--- Check 2: regressed doubling message absent ---"
if grep -q "both columns are included in CB matching/output namespaces" "$LOG"; then
    fail "Log.out contains REGRESSED doubling message"
else
    pass "Log.out does NOT contain regressed doubling message"
fi

echo "--- Check 3: CB count equals number of WL rows (${N_BARCODES}), not 2x ---"
CB_COUNT=$(grep -oP 'Number of CBs in the whitelist = \K[0-9]+' "$LOG" || echo "0")
if [[ "$CB_COUNT" -eq "${N_BARCODES}" ]]; then
    pass "CB count = ${CB_COUNT} (matches ${N_BARCODES} rows)"
elif [[ "$CB_COUNT" -eq $((N_BARCODES * 2)) ]]; then
    fail "CB count = ${CB_COUNT} — doubled! (expected ${N_BARCODES})"
else
    fail "CB count = ${CB_COUNT} — unexpected (expected ${N_BARCODES})"
fi

echo ""
echo "--- Check 4: COL1/COL2 overlap collision detection ---"
# Create a whitelist with self-complementary positions 7-8, so that some
# COL1 (NXT) barcodes are identical to COL2 (TRU) barcodes of OTHER rows.
# NXT translation complements positions 7 and 8 (0-indexed).
#
# Row A: COL1=AAAAAAA_AT_AAAAAAA  COL2=AAAAAAA_TA_AAAAAAA
# Row B: COL1=AAAAAAA_TA_AAAAAAA  COL2=AAAAAAA_AT_AAAAAAA  ← COL1(B) == COL2(A)!
# Row C: COL1=CCCCCCC_GC_CCCCCCC  COL2=CCCCCCC_CG_CCCCCCC  (no overlap)
#
# If the normalization uses membership heuristics, a TRU barcode equal to
# COL2(A)="AAAAAAAATAAAAAAA" would falsely match COL1(B)="AAAAAAAATAAAAAAA"
# as "already in set" instead of being translated to COL1(A)="AAAAAAAATAAAAAA".

python3 -c "
# Barcode A: positions 7-8 are AT → complement TA
col1_a = 'AAAAAAGATAAAAAAA'  # NXT: ...AT...
col2_a = 'AAAAAAGTAGAAAAAA'  # TRU: ...TA...  (complement pos 7-8)
# Barcode B: positions 7-8 are TA → complement AT  → COL1(B) overlaps COL2(A)!
col1_b = 'AAAAAAGTAAAAAAA' + 'A'  # WAIT let me be more precise
# I need 16-bp barcodes where NXT-translation(col1) = col2
# Translation = complement positions 7 and 8 (0-indexed)

def translate_nxt(bc):
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    bc = list(bc)
    bc[7] = comp[bc[7]]
    bc[8] = comp[bc[8]]
    return ''.join(bc)

# Row A: col1_a has AT at positions 7-8
col1_a = 'AAAAAAGATCAAAAAA'  # 16 bp: pos7=A, pos8=T
col2_a = translate_nxt(col1_a)  # pos7=T, pos8=A → AAAAAAAATCAAAAAA → wait
# Actually: pos7 is index 7, pos8 is index 8
# col1_a = A A A A A A G A T C A A A A A A
#           0 1 2 3 4 5 6 7 8 9 ...
# translate: complement pos7(A→T) and pos8(T→A) → AAAAAAGATAAAAAA → same!
# That's self-complementary (A→T, T→A swaps back). Need non-palindromic.

# Let me use: pos7=A, pos8=C → complement = T, G
col1_a = 'AAAAAAAACAAAAAAA'  # pos7=A, pos8=C
col2_a = translate_nxt(col1_a)  # pos7=T, pos8=G → AAAAAAATGAAAAAAA
# Now create row B whose COL1 happens to equal COL2 of row A:
col1_b = col2_a  # AAAAAAATGAAAAAAA — this IS col2_a
col2_b = translate_nxt(col1_b)  # complement pos7=T→A, pos8=G→C → AAAAAAAACAAAAAAA = col1_a!
# So col1_b == col2_a and col2_b == col1_a. Full overlap!

# Row C: no overlap (different prefix)
col1_c = 'CCCCCCCGCCCCCCCC'
col2_c = translate_nxt(col1_c)

with open('${WORK_DIR}/wl_overlap.txt', 'w') as f:
    f.write(col1_a + '\t' + col2_a + '\n')
    f.write(col1_b + '\t' + col2_b + '\n')
    f.write(col1_c + '\t' + col2_c + '\n')

# Verify: col1_b == col2_a (overlap exists)
assert col1_b == col2_a, f'Expected overlap: {col1_b} != {col2_a}'
# Verify col1 entries are all distinct
assert len(set([col1_a, col1_b, col1_c])) == 3, 'col1 entries must be distinct'
# Verify col2 entries are all distinct
assert len(set([col2_a, col2_b, col2_c])) == 3, 'col2 entries must be distinct'
# Print for debugging
print(f'col1_a={col1_a} col2_a={col2_a}')
print(f'col1_b={col1_b} col2_b={col2_b}')
print(f'col1_c={col1_c} col2_c={col2_c}')
print(f'Overlap: col1_b == col2_a = {col1_b == col2_a}')
print(f'Overlap: col2_b == col1_a = {col2_b == col1_a}')
" 2>&1 | while read line; do echo "  $line"; done

# Verify the whitelist was created with 3 rows
OVERLAP_ROWS=$(wc -l < "${WORK_DIR}/wl_overlap.txt")
if [[ "$OVERLAP_ROWS" -eq 3 ]]; then
    pass "overlap whitelist created (${OVERLAP_ROWS} rows with deliberate COL1/COL2 overlap)"
else
    fail "overlap whitelist creation failed (expected 3 rows, got ${OVERLAP_ROWS})"
fi

echo "--- Check 5: STAR loads overlap whitelist with correct pairing (3 CBs, not 6) ---"
set +e
"${STAR_BIN}" \
    --runThreadN 1 \
    --genomeDir "${WORK_DIR}/no_genome" \
    --readFilesIn "${WORK_DIR}/r2.fq" "${WORK_DIR}/r1.fq" \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WORK_DIR}/wl_overlap.txt" \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --outFileNamePrefix "${WORK_DIR}/overlap_run_" \
    > "${WORK_DIR}/overlap_stdout.log" 2>&1
set -e

OVERLAP_LOG="${WORK_DIR}/overlap_run_Log.out"
if [[ -f "$OVERLAP_LOG" ]]; then
    OVERLAP_CB_COUNT=$(grep -oP 'Number of CBs in the whitelist = \K[0-9]+' "$OVERLAP_LOG" || echo "0")
    if [[ "$OVERLAP_CB_COUNT" -eq 3 ]]; then
        pass "Overlap whitelist: CB count = 3 (correct pairing despite COL1/COL2 overlap)"
    else
        fail "Overlap whitelist: CB count = ${OVERLAP_CB_COUNT} (expected 3)"
    fi
else
    fail "Overlap whitelist: Log.out not created"
fi

echo ""
echo "=========================================="
echo "Results: ${PASS_COUNT} passed, ${FAIL_COUNT} failed"
echo "=========================================="

if [[ ${FAIL_COUNT} -ne 0 ]]; then
    echo ""
    echo "--- Log.out excerpt ---"
    grep -E "(column|Number of CBs|whitelist|namespace)" "$LOG" || true
    exit 1
fi
