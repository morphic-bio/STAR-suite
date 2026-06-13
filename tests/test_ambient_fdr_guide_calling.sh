#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CALL_FEATURES="${ROOT_DIR}/core/features/process_features/call_features"
STAR_FEATURE_CALL="${ROOT_DIR}/core/legacy/source/star_feature_call"

if [[ ! -x "${CALL_FEATURES}" ]]; then
  make -C "${ROOT_DIR}/core/features/process_features" call_features
fi
if [[ ! -x "${STAR_FEATURE_CALL}" ]]; then
  make -C "${ROOT_DIR}/core/legacy/source" star_feature_call
fi

TMP_DIR="$(mktemp -d /tmp/star_ambient_fdr_test.XXXXXX)"
trap 'rm -rf "${TMP_DIR}"' EXIT

RAW_DIR="${TMP_DIR}/raw_feature_bc_matrix"
FILTERED_DIR="${TMP_DIR}/filtered_feature_bc_matrix"
OUT_DIR="${TMP_DIR}/crispr_analysis"
OUT_LOW="${TMP_DIR}/crispr_analysis_low"
OUT_GMM="${TMP_DIR}/crispr_analysis_gmm"
OUT_STAR="${TMP_DIR}/star_feature_call_out"
mkdir -p "${RAW_DIR}" "${FILTERED_DIR}"

python3 - "${RAW_DIR}" "${FILTERED_DIR}" <<'PY'
import sys
from pathlib import Path

raw = Path(sys.argv[1])
filtered = Path(sys.argv[2])

features = [
    ("geneX", "Gene X", "Gene Expression"),
    ("guideA", "guideA", "CRISPR Guide Capture"),
    ("guideB", "guideB", "CRISPR Guide Capture"),
    ("guideC", "guideC", "CRISPR Guide Capture"),
]
cells = [f"cell{i}-1" for i in range(1, 6)]
empties = [f"empty{i:02d}-1" for i in range(1, 21)]
raw_barcodes = cells + empties

def write_axes(directory: Path, barcodes):
    with (directory / "features.tsv").open("w") as fh:
        for row in features:
            fh.write("\t".join(row) + "\n")
    with (directory / "barcodes.tsv").open("w") as fh:
        for bc in barcodes:
            fh.write(bc + "\n")

def write_mtx(directory: Path, n_cols: int, entries):
    with (directory / "matrix.mtx").open("w") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write("% synthetic ambient-FDR guide caller test\n")
        fh.write(f"{len(features)} {n_cols} {len(entries)}\n")
        for row, col, value in entries:
            fh.write(f"{row} {col} {value}\n")

write_axes(raw, raw_barcodes)
write_axes(filtered, cells)

# Matrix Market rows are guides and columns are barcodes.
# Filtered cells: clear singlet, clear multiplet, ambient-like weak signal, no-call, low-count no-call.
filtered_entries = [
    (1, 1, 100),
    (1, 2, 100),
    (1, 3, 100),
    (1, 4, 100),
    (1, 5, 100),
    (2, 1, 25),
    (2, 2, 20),
    (3, 2, 19),
    (4, 3, 3),
    (2, 5, 1),
]
write_mtx(filtered, len(cells), filtered_entries)

raw_entries = list(filtered_entries)
# Empty droplets define ambient guide rates: guideA=20, guideB=20, guideC=60.
for i in range(10):
    col = len(cells) + i + 1
    raw_entries.extend([(1, col, 100), (2, col, 2), (3, col, 2), (4, col, 6)])
write_mtx(raw, len(raw_barcodes), raw_entries)
PY

"${CALL_FEATURES}" \
  --guide-caller ambient-fdr \
  --raw-mex-dir "${RAW_DIR}" \
  --guide-fdr 0.01 \
  --guide-fdr-min-umi 1 \
  "${FILTERED_DIR}" \
  "${OUT_DIR}"

"${CALL_FEATURES}" \
  --guide-caller ambient-fdr \
  --raw-mex-dir "${RAW_DIR}" \
  --guide-fdr 1e-20 \
  "${FILTERED_DIR}" \
  "${OUT_LOW}"

"${CALL_FEATURES}" \
  --guide-caller gmm \
  "${FILTERED_DIR}" \
  "${OUT_GMM}"

"${STAR_FEATURE_CALL}" \
  --call-only \
  --compat-perturb \
  --guide-caller both \
  --raw-mex-dir "${RAW_DIR}" \
  --mex-dir "${FILTERED_DIR}" \
  --output-dir "${OUT_STAR}"

python3 - "${OUT_DIR}" "${OUT_LOW}" "${OUT_GMM}" "${OUT_STAR}" <<'PY'
import csv
import json
import math
import sys
from pathlib import Path

out = Path(sys.argv[1])
out_low = Path(sys.argv[2])
out_gmm = Path(sys.argv[3])
out_star = Path(sys.argv[4])

rates = {}
with (out / "guide_ambient_rates.tsv").open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        rates[row["feature_id"]] = float(row["ambient_rate"])
assert set(rates) == {"guideA", "guideB", "guideC"}
assert abs(sum(rates.values()) - 1.0) < 1e-12
assert abs(rates["guideA"] - 0.2) < 1e-12
assert abs(rates["guideB"] - 0.2) < 1e-12
assert abs(rates["guideC"] - 0.6) < 1e-12

barcodes = (out / "guide_qvalues_barcodes.tsv").read_text().splitlines()
features = [line.split("\t")[0] for line in (out / "guide_qvalues_features.tsv").read_text().splitlines()]
assert barcodes == [f"cell{i}-1" for i in range(1, 6)]
assert not any(bc.startswith("empty") for bc in barcodes)
assert features == ["guideA", "guideB", "guideC"]

q_entries = {}
seen_dims = False
with (out / "guide_qvalues.mtx").open() as fh:
    for line in fh:
        if line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) == 3 and not seen_dims:
            n_rows, n_cols, nnz = map(int, parts)
            assert (n_rows, n_cols, nnz) == (5, 3, 5)
            seen_dims = True
            continue
        row, col, q = parts
        r = int(row)
        c = int(col)
        assert 1 <= r <= 5
        assert 1 <= c <= 3
        q_entries[(r - 1, c - 1)] = float(q)
assert len(q_entries) == 5
assert (3, 0) not in q_entries  # missing zero count is not materialized and implies q=1

observed = {
    (0, 0): 25,
    (1, 0): 20,
    (1, 1): 19,
    (2, 2): 3,
    (4, 0): 1,
}
depths = [25, 39, 3, 0, 1]
guide_rates = [0.2, 0.2, 0.6]

def poisson_sf(k, lam):
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    term = math.exp(-lam)
    cdf = term
    for i in range(1, k):
        term *= lam / i
        cdf += term
    return max(0.0, min(1.0, 1.0 - cdf))

pvals = []
for key, count in observed.items():
    cell, guide = key
    pvals.append((key, poisson_sf(count, guide_rates[guide] * depths[cell])))
pvals.sort(key=lambda item: item[1])
total_tests = 5 * 3
expected_q = {}
running = 1.0
for rank_index in range(len(pvals) - 1, -1, -1):
    key, p = pvals[rank_index]
    rank = rank_index + 1
    q = min(running, min(1.0, p * total_tests / rank))
    running = q
    expected_q[key] = q

for key, expected in expected_q.items():
    assert abs(q_entries[key] - expected) < 1e-9, (key, q_entries[key], expected)

calls = {}
with (out / "guide_fdr_calls_per_cell.csv").open() as fh:
    for row in csv.DictReader(fh):
        calls[row["cell_barcode"]] = row
assert calls["cell1-1"]["feature_call"] == "guideA"
assert calls["cell1-1"]["call_status"] == "singlet"
assert calls["cell1-1"]["min_called_umi"] == "25"
assert calls["cell1-1"]["max_called_umi"] == "25"
assert calls["cell2-1"]["feature_call"] == "guideA|guideB"
assert calls["cell2-1"]["call_status"] == "multiplet"
assert calls["cell2-1"]["min_called_umi"] == "19"
assert calls["cell2-1"]["max_called_umi"] == "20"
assert calls["cell3-1"]["feature_call"] == "None"
assert calls["cell4-1"]["feature_call"] == "None"
assert calls["cell5-1"]["feature_call"] == "None"

summary = json.loads((out / "guide_fdr_summary.json").read_text())
assert summary["status"] == "ok"
assert summary["total_tests"] == 15
assert summary["observed_filtered_entries"] == 5
assert summary["qvalue_entries"] == 5
assert summary["ambient_barcodes"] == 20
assert summary["cells_1_feature"] == 1
assert summary["cells_multi_feature"] == 1
assert summary["cells_no_call"] == 3

low_summary = json.loads((out_low / "guide_fdr_summary.json").read_text())
assert low_summary["assigned_cells"] < summary["assigned_cells"]

assert (out_gmm / "protospacer_calls_per_cell.csv").exists()
assert (out_gmm / "protospacer_calls_summary.csv").exists()
for gmm_file in [
    out_gmm / "protospacer_calls_per_cell.csv",
    out_gmm / "protospacer_calls_summary.csv",
    out_gmm / "protospacer_umi_thresholds.csv",
    out_gmm / "protospacer_umi_thresholds.json",
]:
    text = gmm_file.read_text()
    assert "geneX" not in text
    assert "Gene X" not in text
assert not (out_gmm / "ambient_fdr").exists()

star_crispr = out_star / "crispr_analysis"
assert (star_crispr / "protospacer_calls_per_cell.csv").exists()
assert "geneX" not in (star_crispr / "protospacer_calls_per_cell.csv").read_text()
assert (star_crispr / "ambient_fdr" / "guide_fdr_summary.json").exists()
assert (star_crispr / "ambient_fdr" / "guide_qvalues.mtx").exists()
PY

echo "ambient-FDR guide caller synthetic test passed"
