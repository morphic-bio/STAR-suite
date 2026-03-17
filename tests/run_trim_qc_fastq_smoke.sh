#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${TRIM_QC_FASTQ_SMOKE_OUTDIR:-${REPO_ROOT}/tests/trim_qc_fastq_smoke_output}"
BIN="${TRIM_QC_FASTQ_BIN:-${REPO_ROOT}/core/legacy/source/trim_qc_fastq}"
INPUT_FASTQ="${REPO_ROOT}/tests/fixtures/trim_qc_fastq_tiny.fastq"
REPORT_PREFIX="${OUT_DIR}/tiny"

rm -rf "${OUT_DIR}"
mkdir -p "${OUT_DIR}"

if [[ ! -x "${BIN}" ]]; then
  make -C "${REPO_ROOT}/core/legacy/source" trim_qc_fastq
fi

"${BIN}" --input "${INPUT_FASTQ}" --report "${REPORT_PREFIX}" --stage fastq_smoke \
  > "${OUT_DIR}/summary.tsv"

python3 - "${REPORT_PREFIX}.trim_qc.json" <<'PY'
import json
import sys

path = sys.argv[1]
with open(path, "r", encoding="utf-8") as handle:
    data = json.load(handle)

assert data["stage"] == "fastq_smoke", data["stage"]
assert data["total_reads"] == 2, data["total_reads"]
assert data["read_mates"] == 1, data["read_mates"]

mate = data["mates"][0]
assert mate["reads"] == 2, mate["reads"]
assert mate["length_hist"][5] == 2, mate["length_hist"]
assert mate["gc_hist"][40] == 1, mate["gc_hist"][40]
assert mate["gc_hist"][60] == 1, mate["gc_hist"][60]

pos1 = mate["positions"][0]
assert pos1["base_counts"]["A"] == 1, pos1
assert pos1["base_counts"]["G"] == 1, pos1
assert abs(pos1["mean_qual"] - 21.0) < 1e-6, pos1["mean_qual"]
PY

[[ -s "${REPORT_PREFIX}.trim_qc.html" ]]

echo "PASS: trim_qc_fastq smoke"
