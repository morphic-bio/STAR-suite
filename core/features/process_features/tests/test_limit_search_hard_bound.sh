#!/usr/bin/env bash
# Test: --limit_search is a hard bound; no out-of-window rescue.
# Also validates alias parity between --feature_limited_mode and
# --feature_limited_fallback (deprecated).
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSIGN_BIN="${ROOT_DIR}/assignBarcodes"

if [[ ! -x "${ASSIGN_BIN}" ]]; then
  echo "assignBarcodes not found or not executable: ${ASSIGN_BIN}" >&2
  echo "Build it with: make -C ${ROOT_DIR} tools" >&2
  exit 1
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

FASTQ_IN="${TMP_DIR}/fastqs_in_window"
FASTQ_OUT="${TMP_DIR}/fastqs_out_window"
mkdir -p "${FASTQ_IN}" "${FASTQ_OUT}"

WHITELIST="${TMP_DIR}/whitelist.txt"
FEATURE_REF="${TMP_DIR}/feature_ref.csv"

cat > "${WHITELIST}" <<'EOF'
AAAAAAAAAAAAAAAA
EOF

FEAT_SEQ="ACGTACGTACGTACGTACGT"  # 20 bp feature

cat > "${FEATURE_REF}" <<'EOF'
id,name,sequence,feature_type
feat1,Feat1,ACGTACGTACGTACGTACGT,CRISPR Guide Capture
EOF

# --- In-window read: feature starts at offset 0 ---
cat > "${FASTQ_IN}/sample_R1_001.fastq" <<'EOF'
@read_in
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
################################
EOF

cat > "${FASTQ_IN}/sample_R2_001.fastq" <<EOF
@read_in
${FEAT_SEQ}TTTTTTTTTTTTTTTTTTTT
+
########################################
EOF

# --- Out-of-window read: feature starts at offset 15, well outside window ---
cat > "${FASTQ_OUT}/sample_R1_001.fastq" <<'EOF'
@read_out
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
################################
EOF

cat > "${FASTQ_OUT}/sample_R2_001.fastq" <<EOF
@read_out
NNNNNNNNNNNNNNN${FEAT_SEQ}TTTTT
+
########################################
EOF

PASS_COUNT=0
FAIL_COUNT=0
pass() { echo "  PASS: $1"; (( PASS_COUNT++ )) || true; }
fail() { echo "  FAIL: $1" >&2; (( FAIL_COUNT++ )) || true; }

# Helper: extract total counts from matrix.mtx (skip comment lines and header)
matrix_total() {
  local mtx="$1"
  if [[ ! -f "${mtx}" ]]; then
    echo "0"
    return
  fi
  awk 'BEGIN {sum=0; hdr=0} /^%/ {next} {if(hdr==0){hdr=1;next} sum+=$3} END {printf "%.0f", sum}' "${mtx}"
}

# -----------------------------------------------------------------
# Test 1: in-window hit matches (limit_search=2, offset=0)
# -----------------------------------------------------------------
echo "=== Test 1: in-window hit matches ==="
for mode_name in in_window_full in_window_simple; do
  OUT_DIR="${TMP_DIR}/test1_${mode_name}"
  mkdir -p "${OUT_DIR}"
  "${ASSIGN_BIN}" \
    --whitelist "${WHITELIST}" \
    --featurelist "${FEATURE_REF}" \
    --directory "${OUT_DIR}" \
    --barcode_fastq_pattern "_R1_" \
    --forward_fastq_pattern "_R2_" \
    --barcode_length 16 \
    --umi_length 12 \
    --maxHammingDistance 1 \
    --feature_constant_offset 0 \
    --limit_search 2 \
    --feature_limited_mode "${mode_name}" \
    "${FASTQ_IN}" \
    > "${TMP_DIR}/test1_${mode_name}.log" 2>&1

  RESULT_DIR="${OUT_DIR}/$(basename "${FASTQ_IN}")"
  cnt=$(matrix_total "${RESULT_DIR}/matrix.mtx")
  if [[ "${cnt}" -ge 1 ]]; then
    pass "in-window hit matches (${mode_name}): count=${cnt}"
  else
    fail "in-window hit matches (${mode_name}): expected >=1, got ${cnt}"
  fi
done

# -----------------------------------------------------------------
# Test 2: out-of-window true hit remains unmatched
# -----------------------------------------------------------------
echo "=== Test 2: out-of-window hit unmatched ==="
for mode_name in in_window_full in_window_simple; do
  OUT_DIR="${TMP_DIR}/test2_${mode_name}"
  mkdir -p "${OUT_DIR}"
  "${ASSIGN_BIN}" \
    --whitelist "${WHITELIST}" \
    --featurelist "${FEATURE_REF}" \
    --directory "${OUT_DIR}" \
    --barcode_fastq_pattern "_R1_" \
    --forward_fastq_pattern "_R2_" \
    --barcode_length 16 \
    --umi_length 12 \
    --maxHammingDistance 1 \
    --feature_constant_offset 0 \
    --limit_search 2 \
    --feature_limited_mode "${mode_name}" \
    "${FASTQ_OUT}" \
    > "${TMP_DIR}/test2_${mode_name}.log" 2>&1

  RESULT_DIR="${OUT_DIR}/$(basename "${FASTQ_OUT}")"
  cnt=$(matrix_total "${RESULT_DIR}/matrix.mtx")
  if [[ "${cnt}" -eq 0 ]]; then
    pass "out-of-window hit unmatched (${mode_name}): count=${cnt}"
  else
    fail "out-of-window hit unmatched (${mode_name}): expected 0, got ${cnt}"
  fi
done

# -----------------------------------------------------------------
# Test 3: alias parity — --feature_limited_fallback produces same
#          output as --feature_limited_mode
# -----------------------------------------------------------------
echo "=== Test 3: alias parity ==="
for val_pair in "in_window_full:full" "in_window_simple:simple"; do
  IFS=: read -r new_val old_val <<< "${val_pair}"

  OUT_NEW="${TMP_DIR}/test3_new_${new_val}"
  OUT_OLD="${TMP_DIR}/test3_old_${old_val}"
  mkdir -p "${OUT_NEW}" "${OUT_OLD}"

  "${ASSIGN_BIN}" \
    --whitelist "${WHITELIST}" \
    --featurelist "${FEATURE_REF}" \
    --directory "${OUT_NEW}" \
    --barcode_fastq_pattern "_R1_" \
    --forward_fastq_pattern "_R2_" \
    --barcode_length 16 \
    --umi_length 12 \
    --maxHammingDistance 1 \
    --feature_constant_offset 0 \
    --limit_search 2 \
    --feature_limited_mode "${new_val}" \
    "${FASTQ_IN}" \
    > "${TMP_DIR}/test3_new_${new_val}.log" 2>&1

  "${ASSIGN_BIN}" \
    --whitelist "${WHITELIST}" \
    --featurelist "${FEATURE_REF}" \
    --directory "${OUT_OLD}" \
    --barcode_fastq_pattern "_R1_" \
    --forward_fastq_pattern "_R2_" \
    --barcode_length 16 \
    --umi_length 12 \
    --maxHammingDistance 1 \
    --feature_constant_offset 0 \
    --limit_search 2 \
    --feature_limited_fallback "${old_val}" \
    "${FASTQ_IN}" \
    > "${TMP_DIR}/test3_old_${old_val}.log" 2>&1

  NEW_RESULT="${OUT_NEW}/$(basename "${FASTQ_IN}")"
  OLD_RESULT="${OUT_OLD}/$(basename "${FASTQ_IN}")"

  cnt_new=$(matrix_total "${NEW_RESULT}/matrix.mtx")
  cnt_old=$(matrix_total "${OLD_RESULT}/matrix.mtx")

  if [[ "${cnt_new}" == "${cnt_old}" ]]; then
    pass "alias parity ${new_val}/${old_val}: both count=${cnt_new}"
  else
    fail "alias parity ${new_val}/${old_val}: new=${cnt_new} old=${cnt_old}"
  fi

  if diff -q "${NEW_RESULT}/matrix.mtx" "${OLD_RESULT}/matrix.mtx" > /dev/null 2>&1; then
    pass "alias parity ${new_val}/${old_val}: matrix.mtx identical"
  else
    fail "alias parity ${new_val}/${old_val}: matrix.mtx differs"
  fi
done

# -----------------------------------------------------------------
# Test 4: deprecation warning emitted for old flag
# -----------------------------------------------------------------
echo "=== Test 4: deprecation warning ==="
if grep -q "deprecated alias" "${TMP_DIR}/test3_old_full.log"; then
  pass "deprecation warning present for --feature_limited_fallback full"
else
  fail "deprecation warning missing for --feature_limited_fallback full"
fi
if grep -q "deprecated alias" "${TMP_DIR}/test3_old_simple.log"; then
  pass "deprecation warning present for --feature_limited_fallback simple"
else
  fail "deprecation warning missing for --feature_limited_fallback simple"
fi

# -----------------------------------------------------------------
# Summary
# -----------------------------------------------------------------
echo ""
echo "=== Summary: ${PASS_COUNT} passed, ${FAIL_COUNT} failed ==="
if [[ "${FAIL_COUNT}" -gt 0 ]]; then
  exit 1
fi
echo "PASS: limit_search hard-bound and alias parity verified."
