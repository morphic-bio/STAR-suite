#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSIGN_BIN="${ROOT_DIR}/assignBarcodes"

if [[ ! -x "${ASSIGN_BIN}" ]]; then
  echo "assignBarcodes not found or not executable: ${ASSIGN_BIN}" >&2
  echo "Build it with: make -C ${ROOT_DIR} assignBarcodes" >&2
  exit 1
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

FASTQ_DIR="${TMP_DIR}/fastqs"
OUT_DIR="${TMP_DIR}/out"
mkdir -p "${FASTQ_DIR}" "${OUT_DIR}"

WHITELIST="${TMP_DIR}/whitelist.txt"
FEATURE_REF="${TMP_DIR}/feature_ref.csv"

cat > "${WHITELIST}" <<'EOF'
ACGTACGTACGTACGT
EOF

cat > "${FEATURE_REF}" <<'EOF'
id,name,sequence,feature_type
feat1,Feat1,ACGTACGTACGTACGTACGT,CRISPR Guide Capture
EOF

cat > "${FASTQ_DIR}/sample_R1_001.fastq" <<'EOF'
@read1
ACGTACGTACGTACGTAAAAAAAAAAAA
+
############################
EOF

cat > "${FASTQ_DIR}/sample_R2_001.fastq" <<'EOF'
@read1
ACGTACGTACGTACGTACGA
+
####################
EOF

for search_threads in 1 4; do
  RUN_OUT="${OUT_DIR}/search_${search_threads}"
  mkdir -p "${RUN_OUT}"

  "${ASSIGN_BIN}" \
    --whitelist "${WHITELIST}" \
    --featurelist "${FEATURE_REF}" \
    --directory "${RUN_OUT}" \
    --barcode_fastq_pattern "_R1_" \
    --forward_fastq_pattern "_R2_" \
    --barcode_length 16 \
    --umi_length 12 \
    --maxHammingDistance 1 \
    --feature_n 1 \
    --barcode_n 1 \
    --stringency 1 \
    --min_counts 0 \
    --search_threads "${search_threads}" \
    "${FASTQ_DIR}" \
    > "${TMP_DIR}/assign_${search_threads}.log" 2>&1

  RESULT_DIR="${RUN_OUT}/$(basename "${FASTQ_DIR}")"
  MATRIX="${RESULT_DIR}/matrix.mtx"
  FEATURES="${RESULT_DIR}/feature_sequences.txt"

  if [[ ! -f "${MATRIX}" ]]; then
    echo "Missing matrix.mtx for search_threads=${search_threads}: ${MATRIX}" >&2
    exit 1
  fi
  if [[ ! -f "${FEATURES}" ]]; then
    echo "Missing feature_sequences.txt for search_threads=${search_threads}: ${FEATURES}" >&2
    exit 1
  fi

  total_counts="$(awk 'BEGIN {sum=0; header=0} /^%/ {next} {if (header==0) {header=1; next} sum+=$3} END {printf "%.0f", sum}' "${MATRIX}")"
  if [[ "${total_counts}" -lt 1 ]]; then
    echo "Expected at least 1 count in matrix.mtx for search_threads=${search_threads}, got ${total_counts}" >&2
    exit 1
  fi

  awk 'NR>1 && $3==1 {found=1} END {exit !found}' "${FEATURES}" || {
    echo "Expected a Hamming distance of 1 in feature_sequences.txt for search_threads=${search_threads}" >&2
    exit 1
  }
done

echo "PASS: maxHammingDistance is inclusive for feature matching."
