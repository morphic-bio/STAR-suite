#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${ROOT_DIR}/core/legacy/source/STAR"
SUPPRESSIONS="${ROOT_DIR}/tests/compat/lsan_suppressions.txt"

IDX_DIR="${ROOT_DIR}/tests/solo_smoke/ref/star_index"
R1="${ROOT_DIR}/tests/solo_smoke/fastq/R1.fastq"
R2="${ROOT_DIR}/tests/solo_smoke/fastq/R2.fastq"
WL="${ROOT_DIR}/tests/solo_smoke/whitelist.txt"

make -C "${ROOT_DIR}/core/legacy/source" -j ASAN=1 STAR

run_one() {
  local stub="${1}"
  local out_dir="/tmp/cr_compat_integration_stub_${stub}/"
  local log_file="${out_dir}/asan.log"
  rm -rf "${out_dir}"
  mkdir -p "${out_dir}"
  chmod 755 "${out_dir}" || true

  STAR_DISABLE_AMBIG_CB_RESOLVE="${stub}" \
  ASAN_OPTIONS="detect_leaks=1:exitcode=0" \
  LSAN_OPTIONS="detect_leaks=1:suppressions=${SUPPRESSIONS}" \
    "${STAR_BIN}" --runThreadN 2 \
      --genomeDir "${IDX_DIR}" \
      --readFilesIn "${R2}" "${R1}" \
      --outFileNamePrefix "${out_dir}" \
      --outSAMtype None \
      --soloType CB_UMI_Simple \
      --soloCBstart 1 --soloCBlen 16 \
      --soloUMIstart 17 --soloUMIlen 12 \
      --soloBarcodeReadLength 0 \
      --soloCBwhitelist "${WL}" \
      --soloFeatures Gene --soloCellFilter None \
    2>&1 | tee "${log_file}"
}

run_one 1
run_one 0

sig1="/tmp/cr_compat_integration_stub_1/lsan.sig"
sig0="/tmp/cr_compat_integration_stub_0/lsan.sig"
python3 "${ROOT_DIR}/tests/compat/normalize_lsan_log.py" "/tmp/cr_compat_integration_stub_1/asan.log" | sort > "${sig1}"
python3 "${ROOT_DIR}/tests/compat/normalize_lsan_log.py" "/tmp/cr_compat_integration_stub_0/asan.log" | sort > "${sig0}"

if ! diff -u "${sig1}" "${sig0}" >/tmp/cr_compat_integration_lsan.diff; then
  echo "FAIL: leak stack signatures differ between stub=1 and stub=0"
  echo "See: /tmp/cr_compat_integration_lsan.diff"
  exit 1
fi

echo "PASS: CR compat integration smoke (stub on/off)"
