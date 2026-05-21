#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
UNIT_BIN="${REPO_ROOT}/core/legacy/source/ocm_multi_unit_tests"
FIXTURE_ROOT="${REPO_ROOT}/tests/fixtures/ocm_multi_velocyto_shuffle"
RUN_DIR="${FIXTURE_ROOT}/run"
OUTS_DIR="${FIXTURE_ROOT}/outs"

if [[ ! -x "${UNIT_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/legacy/source" -j8 ocm-multi-unit-tests
fi

rm -rf "${OUTS_DIR}" "${RUN_DIR}/outs" "${FIXTURE_ROOT}/samples"
unset OCM_TEST_RUN_DIR OCM_TEST_CONFIG OCM_TEST_LOG
export OCM_TEST_FIXTURE_ROOT="${FIXTURE_ROOT}"
"${UNIT_BIN}" materialize

python3 - "${FIXTURE_ROOT}" <<'PY'
import gzip
import sys
from pathlib import Path

fixture = Path(sys.argv[1])
sample_velo = fixture / "samples/GCM1-Day-4/run/outs/raw_velocyto_feature_bc_matrix"

with gzip.open(sample_velo / "barcodes.tsv.gz", "rt", encoding="utf-8") as handle:
    barcodes = [line.strip() for line in handle if line.strip()]
assert barcodes == ["AAACCCTGTAAGCGCG-1"], f"unexpected OB1 barcodes: {barcodes}"

# Velocyto raw barcodes are reversed vs GeneFull; OB1 is column 4 in pool Velocyto (count 10).
# Positional reuse of GeneFull column index 0 would incorrectly yield count 40 (OB4).
with gzip.open(sample_velo / "spliced.mtx.gz", "rt", encoding="utf-8") as handle:
    lines = [line.strip() for line in handle if line.strip() and not line.startswith("%")]
dims = lines[0].split()
assert len(dims) == 3
_, n_barcodes, nnz = map(int, dims)
assert n_barcodes == 1, n_barcodes
assert nnz == 1, nnz
gene_idx, cell_idx, count = map(int, lines[1].split())
assert gene_idx == 1 and cell_idx == 1 and count == 10, (gene_idx, cell_idx, count)

gex_dir = fixture / "outs/per_sample_outs/GCM1-Day-4/count/sample_raw_feature_bc_matrix"
with gzip.open(gex_dir / "matrix.mtx.gz", "rt", encoding="utf-8") as handle:
    gex_lines = [line.strip() for line in handle if line.strip() and not line.startswith("%")]
_, gex_barcodes, gex_nnz = map(int, gex_lines[0].split())
_, _, gex_count = map(int, gex_lines[1].split())
assert gex_barcodes == 1 and gex_count == 1, (gex_barcodes, gex_count)
PY

echo "PASS: OCM Velocyto barcode-order shuffle mapping"
