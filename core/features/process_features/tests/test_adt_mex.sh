#!/bin/bash
# test_adt_mex.sh - synthetic ADT/protein MEX output smoke test
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PF_DIR="${SCRIPT_DIR}/.."
ASSIGN="${PF_DIR}/assignBarcodes"
WORK_DIR="$(mktemp -d /tmp/adt_mex_test.XXXXXX)"
trap 'rm -rf "${WORK_DIR}"' EXIT

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'
pass() { echo -e "  ${GREEN}PASS${NC} $1"; }
fail() { echo -e "  ${RED}FAIL${NC} $1"; exit 1; }

INPUT="${WORK_DIR}/input"
OUT="${WORK_DIR}/output"
mkdir -p "${INPUT}"

cat > "${INPUT}/feature_ref.csv" <<'EOF'
id,name,sequence,feature_type,target_gene
CD29,CD29,ATCGATCGATCGATCG,Antibody Capture,ITGB1
CD46,CD46,TTAATTAATTAATTAA,Antibody Capture,CD46
EOF

cat > "${INPUT}/whitelist.txt" <<'EOF'
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCT
EOF

python3 - <<'PY' "${INPUT}"
import sys, os
input_dir = sys.argv[1]
features = {
    "CD29": "ATCGATCGATCGATCG",
    "CD46": "TTAATTAATTAATTAA",
}
barcodes = [
    "AAACCCAAGAAACCAT",
    "AAACCCAAGAAACCCA",
    "AAACCCAAGAAACCCT",
]
# Mix mirrors assignbarcodes_baseline: enough BC0 CD29 + BC1 CD46 reads for both barcodes.
reads = [
    (0, "CD29", "AAAAAAAAAAAA"),
    (0, "CD29", "AAAAAAAAAAAA"),  # duplicate UMI
    (0, "CD29", "AAAAAAAAAAAB"),
    (1, "CD46", "CCCCCCCCCCCC"),
    (1, "CD46", "CCCCCCCCCCCC"),  # duplicate UMI
    (1, "CD46", "CCCCCCCCCCD1"),
    (1, "CD46", "CCCCCCCCCCD2"),
    (2, "UNMATCHED", "DDDDDDDDDDDD"),
]
qual = lambda n: "I" * n
with open(os.path.join(input_dir, "sample_R1_001.fastq"), "w") as r1, \
     open(os.path.join(input_dir, "sample_R2_001.fastq"), "w") as r2:
    for i, (bc_idx, feat, umi) in enumerate(reads):
        bc = barcodes[bc_idx]
        seq = features.get(feat, "GGGGGGGGGGGGGGGG")
        r1.write(f"@read{i}\n{bc}{umi}\n+\n{qual(len(bc)+len(umi))}\n")
        r2.write(f"@read{i}\n{seq}{'A'*12}\n+\n{qual(len(seq)+12)}\n")
PY

if [ ! -x "${ASSIGN}" ]; then
    echo "Building assignBarcodes..."
    make -C "${PF_DIR}" assignBarcodes
fi

"${ASSIGN}" \
    -w "${INPUT}/whitelist.txt" \
    -f "${INPUT}/feature_ref.csv" \
    -d "${OUT}" \
    --output-mode adt_mex \
    --skip_empty_drops \
    --skip_qc_outputs \
    "${INPUT}" \
    -b 16 -u 12 > "${WORK_DIR}/run.log" 2>&1

SAMPLE_OUT="${OUT}/input"
for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz protein_quant_summary.json feature_reference.csv; do
    [ -f "${SAMPLE_OUT}/${f}" ] || fail "missing ${f}"
done
pass "MEX + provenance artifacts exist"

if ! zgrep -q "Antibody Capture" "${SAMPLE_OUT}/features.tsv.gz"; then
    fail "features.tsv missing Antibody Capture"
fi
pass "features.tsv has Antibody Capture"

python3 - <<'PY' "${SAMPLE_OUT}"
import gzip, sys
out = sys.argv[1]
with gzip.open(f"{out}/matrix.mtx.gz", "rt") as fh:
    lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("%")]
rows, cols, nnz = map(int, lines[0].split())
assert rows == 2, rows
assert cols == 2, cols
entries = {}
for ln in lines[1:]:
    r, c, v = ln.split()
    entries[(int(r), int(c))] = int(float(v))
expected = {(1, 1): 1, (2, 2): 1}
for key, val in expected.items():
    assert entries.get(key) == val, f"{key}: got {entries.get(key)} expected {val}"
print("matrix ok")
PY
pass "matrix dimensions and CD29/CD46 counts match expected"

if ! grep -q '"mode": "adt_mex"' "${SAMPLE_OUT}/protein_quant_summary.json"; then
    fail "protein_quant_summary.json missing mode=adt_mex"
fi
if ! grep -q 'feature_ref_path' "${SAMPLE_OUT}/protein_quant_summary.json"; then
    fail "protein_quant_summary.json missing feature_ref_path"
fi
if ! grep -q 'ITGB1' "${SAMPLE_OUT}/feature_reference.csv"; then
    fail "feature_reference.csv missing optional target_gene column"
fi
pass "provenance records feature ref and mode"

echo -e "${GREEN}ADT protein MEX test passed${NC}"
