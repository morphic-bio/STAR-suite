#!/bin/bash
# test_hash_demux_mex.sh - synthetic hash-only and mixed ADT+HTO MEX/demux smoke test
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PF_DIR="${SCRIPT_DIR}/.."
ASSIGN="${PF_DIR}/assignBarcodes"
WORK_DIR="$(mktemp -d /tmp/hash_demux_test.XXXXXX)"
trap 'rm -rf "${WORK_DIR}"' EXIT

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m'
pass() { echo -e "  ${GREEN}PASS${NC} $1"; }
fail() { echo -e "  ${RED}FAIL${NC} $1"; exit 1; }

INPUT="${WORK_DIR}/input"
mkdir -p "${INPUT}"

cat > "${INPUT}/whitelist.txt" <<'EOF'
AAACCCAAGAAACCAT
AAACCCAAGAAACCCA
AAACCCAAGAAACCCT
EOF

python3 - <<'PY' "${INPUT}"
import os, sys
input_dir = sys.argv[1]
barcodes = [
    "AAACCCAAGAAACCAT",
    "AAACCCAAGAAACCCA",
    "AAACCCAAGAAACCCT",
]
reads = [
    (0, "HTO1", "AAAAAAAAAAAA"),
    (0, "HTO1", "AAAAAAAAAAAC"),
    (1, "HTO2", "CCCCCCCCCCCC"),
    (1, "HTO2", "CCCCCCCCCCCA"),
    (2, "HTO1", "AAAAAAAAAAAG"),
    (2, "HTO2", "CCCCCCCCCCCT"),
]
features = {
    "HTO1": "ATCGATCGATCGATCG",
    "HTO2": "TTAATTAATTAATTAA",
}
qual = lambda n: "I" * n
with open(os.path.join(input_dir, "hash_R1_001.fastq"), "w") as r1, \
     open(os.path.join(input_dir, "hash_R2_001.fastq"), "w") as r2:
    for i, (bc_idx, feat, umi) in enumerate(reads):
        bc = barcodes[bc_idx]
        seq = features[feat]
        r1.write(f"@read{i}\n{bc}{umi}\n+\n{qual(len(bc)+len(umi))}\n")
        r2.write(f"@read{i}\n{seq}{'A'*12}\n+\n{qual(len(seq)+12)}\n")
PY

cat > "${INPUT}/hash_ref.csv" <<'EOF'
id,name,sequence,feature_type
HTO1,HTO1,ATCGATCGATCGATCG,Multiplexing Capture
HTO2,HTO2,TTAATTAATTAATTAA,Multiplexing Capture
EOF

if [ ! -x "${ASSIGN}" ]; then
    make -C "${PF_DIR}" assignBarcodes
fi

HASH_OUT="${WORK_DIR}/hash_only"
"${ASSIGN}" \
    -w "${INPUT}/whitelist.txt" \
    -f "${INPUT}/hash_ref.csv" \
    -d "${HASH_OUT}" \
    --output-mode adt_mex \
    --hash-demux yes \
    --hash-min-total 1 \
    --hash-min-top 1 \
    --hash-min-ratio 2.0 \
    --skip_empty_drops \
    --skip_qc_outputs \
    "${INPUT}" \
    -b 16 -u 12 > "${WORK_DIR}/hash_run.log" 2>&1

SAMPLE_OUT="${HASH_OUT}/input"
for f in hash/barcodes.tsv.gz hash/features.tsv.gz hash/matrix.mtx.gz \
         hash_demux_assignments.tsv singlet_barcodes.tsv doublet_barcodes.tsv \
         negative_barcodes.tsv hash_demux_summary.json hash_demux_command.txt; do
    [ -f "${SAMPLE_OUT}/${f}" ] || fail "missing ${f}"
done
pass "hash-only library emits hash MEX + demux contract"

if ! zgrep -q "Multiplexing Capture" "${SAMPLE_OUT}/hash/features.tsv.gz"; then
    fail "hash features.tsv missing Multiplexing Capture"
fi
pass "hash features.tsv uses Multiplexing Capture"

python3 - <<'PY' "${SAMPLE_OUT}/hash/matrix.mtx.gz"
import gzip, sys
with gzip.open(sys.argv[1], "rt") as fh:
    lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("%")]
rows, cols, nnz = map(int, lines[0].split())
entries = [ln for ln in lines[1:]]
assert rows == 2 and cols >= 2, (rows, cols)
assert nnz == len(entries) and nnz > 0, (nnz, len(entries))
print("hash matrix ok", nnz)
PY
pass "hash matrix.mtx.gz has coordinate entries"

python3 - <<'PY' "${SAMPLE_OUT}/hash_demux_summary.json"
import json, sys
summary = json.load(open(sys.argv[1], encoding="utf-8"))
assert summary["n_hash_features"] == 2, summary
total = summary["n_singlet"] + summary["n_doublet"] + summary["n_negative"]
assert total >= 2, summary
print("summary ok")
PY
pass "hash demux summary counts barcodes"

# Mixed ADT + HTO reference (hashtag rows first so assignment indexes them)
cat > "${INPUT}/mixed_ref.csv" <<'EOF'
id,name,sequence,feature_type
hashtag1,hashtag1,ATCGATCGATCGATCG,Antibody Capture
hashtag2,hashtag2,TTAATTAATTAATTAA,Antibody Capture
CD29,CD29,GGGGGGGGGGGGGGGG,Antibody Capture
EOF

MIXED_INPUT="${WORK_DIR}/mixed_input"
mkdir -p "${MIXED_INPUT}"
cp "${INPUT}/whitelist.txt" "${MIXED_INPUT}/"

python3 - <<'PY' "${MIXED_INPUT}"
import os, sys
input_dir = sys.argv[1]
barcodes = ["AAACCCAAGAAACCAT", "AAACCCAAGAAACCCA", "AAACCCAAGAAACCCT"]
reads = [
    (0, "hashtag1", "AAAAAAAAAAAA"),
    (0, "hashtag1", "AAAAAAAAAAAC"),
    (1, "hashtag2", "CCCCCCCCCCCC"),
    (1, "hashtag2", "CCCCCCCCCCCA"),
    (2, "CD29", "AAAAAAAAAAAA"),
    (2, "CD29", "AAAAAAAAAAAC"),
]
features = {
    "CD29": "GGGGGGGGGGGGGGGG",
    "hashtag1": "ATCGATCGATCGATCG",
    "hashtag2": "TTAATTAATTAATTAA",
}
qual = lambda n: "I" * n
with open(os.path.join(input_dir, "mixed_R1_001.fastq"), "w") as r1, \
     open(os.path.join(input_dir, "mixed_R2_001.fastq"), "w") as r2:
    for i, (bc_idx, feat, umi) in enumerate(reads):
        bc = barcodes[bc_idx]
        seq = features[feat]
        r1.write(f"@read{i}\n{bc}{umi}\n+\n{qual(len(bc)+len(umi))}\n")
        r2.write(f"@read{i}\n{seq}{'A'*12}\n+\n{qual(len(seq)+12)}\n")
PY

MIXED_OUT="${WORK_DIR}/mixed"
"${ASSIGN}" \
    -w "${INPUT}/whitelist.txt" \
    -f "${INPUT}/mixed_ref.csv" \
    -d "${MIXED_OUT}" \
    --output-mode adt_mex \
    --hash-demux yes \
    --hash-feature-selector id_prefix:hashtag \
    --hash-min-total 1 \
    --hash-min-top 1 \
    --hash-min-ratio 2.0 \
    --skip_empty_drops \
    --skip_qc_outputs \
    "${MIXED_INPUT}" \
    -b 16 -u 12 > "${WORK_DIR}/mixed_run.log" 2>&1

MIXED_SAMPLE="${MIXED_OUT}/mixed_input"
for d in hash protein; do
    for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
        [ -f "${MIXED_SAMPLE}/${d}/${f}" ] || fail "missing ${d}/${f}"
    done
done
pass "mixed ADT+HTO emits hash/ and protein/ MEX"

HASH_ROWS=$(zgrep -c . "${MIXED_SAMPLE}/hash/features.tsv.gz" || true)
PROT_ROWS=$(zgrep -c . "${MIXED_SAMPLE}/protein/features.tsv.gz" || true)
if [ "${HASH_ROWS}" -ne 2 ] || [ "${PROT_ROWS}" -ne 1 ]; then
    fail "mixed feature row split expected hash=2 protein=1, got hash=${HASH_ROWS} protein=${PROT_ROWS}"
fi
pass "mixed library splits 2 hash + 1 protein features"

python3 - <<'PY' "${MIXED_SAMPLE}/hash/matrix.mtx.gz" "${MIXED_SAMPLE}/protein/matrix.mtx.gz"
import gzip, sys

def check(path, min_rows, require_nnz):
    with gzip.open(path, "rt") as fh:
        lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("%")]
    rows, cols, nnz = map(int, lines[0].split())
    entries = [ln for ln in lines[1:]]
    assert rows >= min_rows and cols >= 1, (path, rows, cols)
    assert nnz == len(entries), (path, nnz, len(entries))
    if require_nnz:
        assert nnz > 0, (path, nnz)
    print(path, "ok", nnz)

check(sys.argv[1], 2, True)
check(sys.argv[2], 1, True)
PY
pass "mixed hash/protein matrix.mtx.gz headers match coordinate bodies"

# Explicit selector with zero matches must fail loudly.
BAD_OUT="${WORK_DIR}/bad_selector"
set +e
"${ASSIGN}" \
    -w "${INPUT}/whitelist.txt" \
    -f "${INPUT}/mixed_ref.csv" \
    -d "${BAD_OUT}" \
    --output-mode adt_mex \
    --hash-demux yes \
    --hash-feature-selector id_prefix:doesnotmatch \
    --skip_empty_drops \
    --skip_qc_outputs \
    "${MIXED_INPUT}" \
    -b 16 -u 12 > "${WORK_DIR}/bad_selector.log" 2>&1
BAD_RC=$?
set -e
if [ "${BAD_RC}" -eq 0 ]; then
    fail "bad hash selector should exit non-zero"
fi
if ! grep -q 'no hash rows matched' "${WORK_DIR}/bad_selector.log"; then
    fail "bad hash selector should log no hash rows matched"
fi
pass "explicit selector with zero hash rows fails loudly"

echo -e "${GREEN}Hash demux MEX test passed${NC}"
