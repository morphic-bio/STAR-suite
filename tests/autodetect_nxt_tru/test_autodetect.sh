#!/usr/bin/env bash
# Test NXT/TRU auto-detection using standalone assignBarcodes with synthetic data.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FIXTURES="${SCRIPT_DIR}/fixtures"
ASSIGN="${SCRIPT_DIR}/../../core/features/process_features/assignBarcodes"

if [[ ! -x "$ASSIGN" ]]; then
    echo "ERROR: assignBarcodes not found at $ASSIGN"
    echo "Build with: make -C core/features/process_features tools"
    exit 1
fi

if [[ ! -d "$FIXTURES/reads_tru" ]]; then
    echo "Generating synthetic fixtures..."
    python3 "${SCRIPT_DIR}/generate_synthetic.py"
fi

PASS=0
FAIL=0

run_test() {
    local name="$1"
    local whitelist="$2"
    local reads_dir="$3"
    local expected_mode="$4"
    local outdir
    outdir=$(mktemp -d "/tmp/autodetect_test_${name}_XXXXXX")

    echo "--- Test: $name ---"
    echo "  whitelist: $(basename "$whitelist")"
    echo "  reads: $(basename "$reads_dir")"
    echo "  expected: $expected_mode"

    local log
    log=$("$ASSIGN" \
        -w "$whitelist" \
        -f "$FIXTURES/feature_ref.csv" \
        -d "$outdir" \
        --autodetect_chemistry 1 \
        --autodetect_chemistry_reads 500 \
        --autodetect_chemistry_min_hits 5 \
        "$reads_dir" 2>&1) || true

    local detected
    detected=$(echo "$log" | grep -o 'raw_frac=[0-9.]*' | head -1 || echo "")
    local notice
    notice=$(echo "$log" | grep 'NOTICE: chemistry auto-detect' || echo "")
    local warning
    warning=$(echo "$log" | grep 'WARNING: chemistry auto-detect' || echo "")

    echo "  log lines:"
    echo "$notice" | sed 's/^/    /'
    echo "$warning" | sed 's/^/    /'

    rm -rf "$outdir"

    # Require detection output lines to be present (not silent regression)
    if ! echo "$log" | grep -q 'NOTICE: chemistry auto-detect:'; then
        echo "  FAIL (no detection output found — detection may not have run)"
        FAIL=$((FAIL+1))
        return
    fi

    if echo "$log" | grep -q "result AMBIGUOUS" && [[ "$expected_mode" == "AMBIGUOUS" ]]; then
        echo "  PASS (AMBIGUOUS)"
        PASS=$((PASS+1))
        return
    fi

    if [[ -z "$detected" ]]; then
        echo "  FAIL (raw_frac not found in detection output)"
        FAIL=$((FAIL+1))
        return
    fi

    local raw_frac
    raw_frac=$(echo "$detected" | grep -oP '[0-9.]+' || echo "")
    if [[ -z "$raw_frac" ]]; then
        echo "  FAIL (could not parse raw_frac value)"
        FAIL=$((FAIL+1))
        return
    fi

    case "$expected_mode" in
        RAW_MATCH)
            if awk "BEGIN{exit !($raw_frac >= 0.8)}"; then
                echo "  PASS (raw_frac=$raw_frac >= 0.8)"
                PASS=$((PASS+1))
            else
                echo "  FAIL (raw_frac=$raw_frac, expected >= 0.8)"
                FAIL=$((FAIL+1))
            fi
            ;;
        TRANSLATED_MATCH)
            if awk "BEGIN{exit !($raw_frac <= 0.2)}"; then
                echo "  PASS (raw_frac=$raw_frac <= 0.2)"
                PASS=$((PASS+1))
            else
                echo "  FAIL (raw_frac=$raw_frac, expected <= 0.2)"
                FAIL=$((FAIL+1))
            fi
            ;;
        AMBIGUOUS)
            echo "  FAIL (expected AMBIGUOUS but got raw_frac=$raw_frac)"
            FAIL=$((FAIL+1))
            ;;
    esac
}

echo "=========================================="
echo "NXT/TRU Auto-Detection Tests"
echo "=========================================="
echo ""

# Case 1: TRU whitelist + TRU reads → RAW_MATCH
run_test "tru_wl_tru_reads" "$FIXTURES/whitelist_tru.txt" "$FIXTURES/reads_tru" "RAW_MATCH"
echo ""

# Case 2: TRU whitelist + NXT reads → TRANSLATED_MATCH
run_test "tru_wl_nxt_reads" "$FIXTURES/whitelist_tru.txt" "$FIXTURES/reads_nxt" "TRANSLATED_MATCH"
echo ""

# Case 3: NXT whitelist + NXT reads → RAW_MATCH
run_test "nxt_wl_nxt_reads" "$FIXTURES/whitelist_nxt.txt" "$FIXTURES/reads_nxt" "RAW_MATCH"
echo ""

# Case 4: NXT whitelist + TRU reads → TRANSLATED_MATCH
run_test "nxt_wl_tru_reads" "$FIXTURES/whitelist_nxt.txt" "$FIXTURES/reads_tru" "TRANSLATED_MATCH"
echo ""

# Case 5: TRU whitelist + mixed reads → AMBIGUOUS
run_test "tru_wl_mixed_reads" "$FIXTURES/whitelist_tru.txt" "$FIXTURES/reads_mixed" "AMBIGUOUS"
echo ""

# Case 6: TRU whitelist + tiny reads (drain finalization)
run_test "tru_wl_tiny_reads" "$FIXTURES/whitelist_tru.txt" "$FIXTURES/reads_tiny" "RAW_MATCH"
echo ""

# ==========================================
# Multi-library test 7: cross-chemistry detection isolation
# Same whitelist, different read chemistries. Validates that each
# pf_context detects independently. Library A (NXT reads, TRU whitelist)
# gets 0 barcode output because translate_NXT is not set — this is
# expected and explicitly checked.
# ==========================================
echo "--- Test: multi_cross_chemistry_detection ---"
echo "  Library A: CRISPR + NXT reads (TRU whitelist) → expect TRANSLATED_MATCH, 0 barcodes"
echo "  Library B: Lineage + TRU reads (TRU whitelist) → expect RAW_MATCH, >0 barcodes"

MULTI_OUTDIR_A=$(mktemp -d "/tmp/autodetect_multi_libA_XXXXXX")
MULTI_OUTDIR_B=$(mktemp -d "/tmp/autodetect_multi_libB_XXXXXX")

log_a=$("$ASSIGN" \
    -w "$FIXTURES/whitelist_tru.txt" \
    -f "$FIXTURES/feature_ref.csv" \
    -d "$MULTI_OUTDIR_A" \
    --autodetect_chemistry 1 \
    --autodetect_chemistry_reads 500 \
    --autodetect_chemistry_min_hits 5 \
    "$FIXTURES/reads_lib_a_nxt" 2>&1) || true

log_b=$("$ASSIGN" \
    -w "$FIXTURES/whitelist_tru.txt" \
    -f "$FIXTURES/feature_ref_lineage.csv" \
    -d "$MULTI_OUTDIR_B" \
    --autodetect_chemistry 1 \
    --autodetect_chemistry_reads 500 \
    --autodetect_chemistry_min_hits 5 \
    "$FIXTURES/reads_lib_b_tru" 2>&1) || true

frac_a=$(echo "$log_a" | grep -o 'raw_frac=[0-9.]*' | head -1 | grep -oP '[0-9.]+' || echo "")
frac_b=$(echo "$log_b" | grep -o 'raw_frac=[0-9.]*' | head -1 | grep -oP '[0-9.]+' || echo "")

echo "  Library A: $(echo "$log_a" | grep 'NOTICE: chemistry auto-detect:' | head -2 | sed 's/^/    /')"
echo "  Library B: $(echo "$log_b" | grep 'NOTICE: chemistry auto-detect:' | head -2 | sed 's/^/    /')"

multi7_pass=true

if [[ -n "$frac_a" ]] && awk "BEGIN{exit !($frac_a <= 0.2)}"; then
    echo "  Library A detect PASS (raw_frac=$frac_a → TRANSLATED_MATCH)"
else
    echo "  Library A detect FAIL (raw_frac=$frac_a, expected <= 0.2)"
    multi7_pass=false
fi

if [[ -n "$frac_b" ]] && awk "BEGIN{exit !($frac_b >= 0.8)}"; then
    echo "  Library B detect PASS (raw_frac=$frac_b → RAW_MATCH)"
else
    echo "  Library B detect FAIL (raw_frac=$frac_b, expected >= 0.8)"
    multi7_pass=false
fi

# Library A: NXT reads can't match TRU whitelist without translate_NXT → 0 output expected
barcodes_a=$(find "$MULTI_OUTDIR_A" -name 'barcodes.txt' -not -path '*/filtered/*' | head -1)
if [[ -n "$barcodes_a" && -f "$barcodes_a" ]]; then
    n_bc_a=$(wc -l < "$barcodes_a")
else
    n_bc_a=0
fi
if [[ "$n_bc_a" -eq 0 ]]; then
    echo "  Library A output PASS (0 barcodes — expected: NXT reads vs TRU whitelist, no translate_NXT)"
else
    echo "  Library A output FAIL ($n_bc_a barcodes — expected 0: NXT reads should not match TRU whitelist)"
    multi7_pass=false
fi

# Library B: TRU reads match TRU whitelist → must have non-zero output
barcodes_b=$(find "$MULTI_OUTDIR_B" -name 'barcodes.txt' -not -path '*/filtered/*' | head -1)
if [[ -n "$barcodes_b" && -f "$barcodes_b" ]]; then
    n_bc_b=$(wc -l < "$barcodes_b")
else
    n_bc_b=0
fi
if [[ "$n_bc_b" -gt 0 ]]; then
    echo "  Library B output PASS ($n_bc_b barcodes)"
else
    echo "  Library B output FAIL (0 barcodes — TRU reads should match TRU whitelist)"
    multi7_pass=false
fi

rm -rf "$MULTI_OUTDIR_A" "$MULTI_OUTDIR_B"

if $multi7_pass; then
    echo "  PASS"
    PASS=$((PASS+1))
else
    echo "  FAIL"
    FAIL=$((FAIL+1))
fi
echo ""

# ==========================================
# Multi-library test 8: matched-whitelist output correctness
# Each library uses a whitelist matching its read chemistry, so both
# produce non-zero output. Validates independent detection + output.
# ==========================================
echo "--- Test: multi_matched_whitelist_output ---"
echo "  Library A: CRISPR + NXT reads (NXT whitelist) → expect RAW_MATCH, >0 barcodes"
echo "  Library B: Lineage + TRU reads (TRU whitelist) → expect RAW_MATCH, >0 barcodes"

MATCHED_OUTDIR_A=$(mktemp -d "/tmp/autodetect_matched_libA_XXXXXX")
MATCHED_OUTDIR_B=$(mktemp -d "/tmp/autodetect_matched_libB_XXXXXX")

log_ma=$("$ASSIGN" \
    -w "$FIXTURES/whitelist_nxt.txt" \
    -f "$FIXTURES/feature_ref.csv" \
    -d "$MATCHED_OUTDIR_A" \
    --autodetect_chemistry 1 \
    --autodetect_chemistry_reads 500 \
    --autodetect_chemistry_min_hits 5 \
    "$FIXTURES/reads_lib_a_nxt" 2>&1) || true

log_mb=$("$ASSIGN" \
    -w "$FIXTURES/whitelist_tru.txt" \
    -f "$FIXTURES/feature_ref_lineage.csv" \
    -d "$MATCHED_OUTDIR_B" \
    --autodetect_chemistry 1 \
    --autodetect_chemistry_reads 500 \
    --autodetect_chemistry_min_hits 5 \
    "$FIXTURES/reads_lib_b_tru" 2>&1) || true

frac_ma=$(echo "$log_ma" | grep -o 'raw_frac=[0-9.]*' | head -1 | grep -oP '[0-9.]+' || echo "")
frac_mb=$(echo "$log_mb" | grep -o 'raw_frac=[0-9.]*' | head -1 | grep -oP '[0-9.]+' || echo "")

echo "  Library A: $(echo "$log_ma" | grep 'NOTICE: chemistry auto-detect:' | head -2 | sed 's/^/    /')"
echo "  Library B: $(echo "$log_mb" | grep 'NOTICE: chemistry auto-detect:' | head -2 | sed 's/^/    /')"

multi8_pass=true

# Both should be RAW_MATCH (reads match their own whitelist)
if [[ -n "$frac_ma" ]] && awk "BEGIN{exit !($frac_ma >= 0.8)}"; then
    echo "  Library A detect PASS (raw_frac=$frac_ma → RAW_MATCH)"
else
    echo "  Library A detect FAIL (raw_frac=$frac_ma, expected >= 0.8)"
    multi8_pass=false
fi

if [[ -n "$frac_mb" ]] && awk "BEGIN{exit !($frac_mb >= 0.8)}"; then
    echo "  Library B detect PASS (raw_frac=$frac_mb → RAW_MATCH)"
else
    echo "  Library B detect FAIL (raw_frac=$frac_mb, expected >= 0.8)"
    multi8_pass=false
fi

# Both must produce non-zero barcode output
barcodes_ma=$(find "$MATCHED_OUTDIR_A" -name 'barcodes.txt' -not -path '*/filtered/*' | head -1)
barcodes_mb=$(find "$MATCHED_OUTDIR_B" -name 'barcodes.txt' -not -path '*/filtered/*' | head -1)

n_bc_ma=0
n_bc_mb=0
if [[ -n "$barcodes_ma" && -f "$barcodes_ma" ]]; then n_bc_ma=$(wc -l < "$barcodes_ma"); fi
if [[ -n "$barcodes_mb" && -f "$barcodes_mb" ]]; then n_bc_mb=$(wc -l < "$barcodes_mb"); fi

if [[ "$n_bc_ma" -gt 0 ]]; then
    echo "  Library A output PASS ($n_bc_ma barcodes)"
else
    echo "  Library A output FAIL (0 barcodes — NXT reads should match NXT whitelist)"
    multi8_pass=false
fi

if [[ "$n_bc_mb" -gt 0 ]]; then
    echo "  Library B output PASS ($n_bc_mb barcodes)"
else
    echo "  Library B output FAIL (0 barcodes — TRU reads should match TRU whitelist)"
    multi8_pass=false
fi

# Output barcodes should be in different namespaces (NXT vs TRU)
if [[ "$n_bc_ma" -gt 0 && "$n_bc_mb" -gt 0 ]]; then
    shared=$(comm -12 <(sort "$barcodes_ma") <(sort "$barcodes_mb") | wc -l)
    if [[ "$shared" -eq 0 ]]; then
        echo "  Namespace cross-check PASS (0 shared — NXT and TRU outputs distinct)"
    else
        echo "  Namespace cross-check FAIL ($shared shared — NXT and TRU outputs should be fully distinct)"
        multi8_pass=false
    fi
fi

rm -rf "$MATCHED_OUTDIR_A" "$MATCHED_OUTDIR_B"

if $multi8_pass; then
    echo "  PASS"
    PASS=$((PASS+1))
else
    echo "  FAIL"
    FAIL=$((FAIL+1))
fi
echo ""

echo "=========================================="
echo "Results: $PASS passed, $FAIL failed"
echo "=========================================="

if [[ $FAIL -gt 0 ]]; then
    exit 1
fi
