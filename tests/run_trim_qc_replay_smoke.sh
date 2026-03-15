#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
OUT_DIR="${TRIM_QC_REPLAY_SMOKE_OUTDIR:-${REPO_ROOT}/tests/trim_qc_replay_smoke_output}"
BIN="${TRIM_QC_REPLAY_BIN:-${REPO_ROOT}/core/legacy/source/trim_qc_replay}"
INPUT_SAM="${REPO_ROOT}/tests/fixtures/trim_qc_replay_tiny.sam"
REPORT_PREFIX="${OUT_DIR}/tiny"

mkdir -p "${OUT_DIR}"

if [[ ! -x "${BIN}" ]]; then
  make -C "${REPO_ROOT}/core/legacy/source" trim_qc_replay
fi

"${BIN}" \
  --input "${INPUT_SAM}" \
  --report "${REPORT_PREFIX}" \
  --stage smoke_replay > "${OUT_DIR}/summary.tsv"

python3 - "${REPORT_PREFIX}.trim_qc.json" <<'PY'
import json
import sys

path = sys.argv[1]
with open(path, "r", encoding="utf-8") as handle:
    data = json.load(handle)

assert data["stage"] == "smoke_replay", data["stage"]
assert data["read_mates"] == 2, data["read_mates"]
assert data["total_reads"] == 3, data["total_reads"]
assert data["quality_base"] == 33, data["quality_base"]
assert data["trim_stats"]["reads_processed"] == 0, data["trim_stats"]
assert data["mates"][0]["reads"] == 2, data["mates"][0]["reads"]
assert data["mates"][1]["reads"] == 1, data["mates"][1]["reads"]
assert data["mates"][0]["length_hist"][10] == 2, data["mates"][0]["length_hist"]
assert data["mates"][1]["length_hist"][10] == 1, data["mates"][1]["length_hist"]
assert data["mates"][0]["positions"][0]["base_counts"]["A"] == 1, data["mates"][0]["positions"][0]
assert data["mates"][0]["positions"][0]["base_counts"]["G"] == 1, data["mates"][0]["positions"][0]
assert abs(data["mates"][0]["positions"][0]["mean_qual"] - 41.0) < 1e-6, data["mates"][0]["positions"][0]
assert abs(data["mates"][1]["positions"][0]["mean_qual"] - 41.0) < 1e-6, data["mates"][1]["positions"][0]
PY

[[ -s "${REPORT_PREFIX}.trim_qc.html" ]]
[[ -s "${OUT_DIR}/summary.tsv" ]]

echo "PASS: trim_qc_replay smoke"
