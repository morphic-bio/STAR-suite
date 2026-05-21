#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
UNIT_BIN="${REPO_ROOT}/core/legacy/source/ocm_multi_unit_tests"
FIXTURE_ROOT="${REPO_ROOT}/tests/fixtures/ocm_multi_tiny"
RUN_DIR="${FIXTURE_ROOT}/run"
OUTS_DIR="${FIXTURE_ROOT}/outs"

if [[ ! -x "${UNIT_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/legacy/source" -j8 ocm-multi-unit-tests
fi

rm -rf "${OUTS_DIR}" "${RUN_DIR}/outs" "${FIXTURE_ROOT}/samples"
unset OCM_TEST_RUN_DIR OCM_TEST_CONFIG OCM_TEST_LOG
export OCM_TEST_FIXTURE_ROOT="${FIXTURE_ROOT}"
"${UNIT_BIN}" all

python3 - "${OUTS_DIR}" "${FIXTURE_ROOT}" <<'PY'
import gzip
import json
import sys
from pathlib import Path

outs = Path(sys.argv[1])
cells = json.loads((outs / "multi/multiplexing_analysis/cells_per_tag.json").read_text(encoding="utf-8"))
assert cells["OB1"] == ["AAACCCTGTAAGCGCG-1"]
assert cells["OB2"] == ["AAACCCGCAACTAGAC-1"]
assert set(cells["OB3"]) == set()
assert set(cells["OB4"]) == set()

filtered_csv = (outs / "per_sample_outs/GCM1-Day-4/count/sample_filtered_barcodes.csv").read_text(encoding="utf-8").strip().splitlines()
assert filtered_csv == ["GRCh38,AAACCCTGTAAGCGCG-1"]

raw_bc = outs / "multi/count/raw_feature_bc_matrix/barcodes.tsv.gz"
with gzip.open(raw_bc, "rt", encoding="utf-8") as handle:
    barcodes = [line.strip() for line in handle if line.strip()]
assert barcodes == [
    "AAAAAAAAAAAAAAAA-1",
    "AAACCAAAGCATTGAT-1",
    "AAACCATTCACCTGGG-1",
    "AAACCCGCAACTAGAC-1",
    "AAACCCTGTAAGCGCG-1",
]

raw_mtx = outs / "multi/count/raw_feature_bc_matrix/matrix.mtx.gz"
with gzip.open(raw_mtx, "rt", encoding="utf-8") as handle:
    raw_lines = [line.strip() for line in handle if line.strip() and not line.startswith("%")]
assert raw_lines[0] == "2 5 4", raw_lines[0]
assert raw_lines[1:] == ["1 5 1", "2 4 2", "1 3 3", "2 2 4"], raw_lines[1:]

union_dir = outs / "per_sample_outs/Union-Test/count/sample_raw_feature_bc_matrix/barcodes.tsv.gz"
with gzip.open(union_dir, "rt", encoding="utf-8") as handle:
    union_barcodes = [line.strip() for line in handle if line.strip()]
assert sorted(union_barcodes) == sorted(["AAACCCTGTAAGCGCG-1", "AAACCCGCAACTAGAC-1"])

fixture_root = Path(sys.argv[2])
sample_velo = fixture_root / "samples/GCM1-Day-4/run/outs/raw_velocyto_feature_bc_matrix"
for name in ("spliced.mtx.gz", "unspliced.mtx.gz", "ambiguous.mtx.gz"):
    assert (sample_velo / name).is_file(), f"missing downstream Velocyto layer {name}"

with gzip.open(sample_velo / "spliced.mtx.gz", "rt", encoding="utf-8") as handle:
    velo_lines = [line.strip() for line in handle if line.strip() and not line.startswith("%")]
_, velo_barcodes, velo_nnz = map(int, velo_lines[0].split())
_, _, velo_spliced = map(int, velo_lines[1].split())
assert velo_barcodes == 1 and velo_spliced == 5, (velo_barcodes, velo_spliced)

gex_sample = outs / "per_sample_outs/GCM1-Day-4/count/sample_raw_feature_bc_matrix"
with gzip.open(gex_sample / "matrix.mtx.gz", "rt", encoding="utf-8") as handle:
    gex_lines = [line.strip() for line in handle if line.strip() and not line.startswith("%")]
_, gex_barcodes, _ = map(int, gex_lines[0].split())
_, _, gex_count = map(int, gex_lines[1].split())
assert gex_barcodes == 1 and gex_count == 1, (gex_barcodes, gex_count)
assert velo_spliced == 5, velo_spliced
PY

echo "PASS: OCM tiny materializer structure checks"
