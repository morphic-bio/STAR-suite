#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)"
HARNESS="${BGZF_INPUT_HARNESS_BIN:-${ROOT_DIR}/core/legacy/source/bgzf_input_harness}"
INPUT_DIR="${JAX_FASTQ_DIR:-/mnt/pikachu/JAX_sequences/JAX_scRNAseq01}"
THREADS="${THREADS:-32}"
LOG_DIR="${ROOT_DIR}/docs/benchmarks/bgzf_ingest_20260901"
SUMMARY="${LOG_DIR}/jax_bgzi_preindex.tsv"

[[ -x "${HARNESS}" ]] || {
    echo "ERROR: build ${HARNESS} with make -C core/legacy/source bgzf-input-harness" >&2
    exit 1
}

mapfile -t inputs < <(find "${INPUT_DIR}" -maxdepth 1 -type f \
    -name 'SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L*_R[12]_001.fastq.gz' | sort)
[[ "${#inputs[@]}" -eq 16 ]] || {
    echo "ERROR: expected 16 JAX SC2300771 FASTQs, found ${#inputs[@]}" >&2
    exit 1
}

printf 'file\trecord_count\tblock_count\tcache_status\tseconds\n' > "${SUMMARY}"
for input in "${inputs[@]}"; do
    start="$(date +%s)"
    result="$(${HARNESS} --mode index --input "${input}" --threads "${THREADS}" | \
        python3 -c 'import json,sys; value=json.load(sys.stdin); print("{}\t{}\t{}".format(value["record_count"], len(value["blocks"]), value["cache_status"]))')"
    end="$(date +%s)"
    printf '%s\t%s\t%s\n' "$(basename "${input}")" "${result}" "$((end - start))" \
        | tee -a "${SUMMARY}"
done

python3 - "${SUMMARY}" <<'PY'
import csv
import re
import sys

rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
by_lane = {}
for row in rows:
    match = re.search(r"_L(\d{3})_R([12])_", row["file"])
    if not match:
        raise SystemExit(f"cannot parse lane/mate from {row['file']}")
    by_lane.setdefault(match.group(1), {})[match.group(2)] = int(row["record_count"])
if len(by_lane) != 8:
    raise SystemExit(f"expected 8 lanes, found {len(by_lane)}")
for lane, mates in sorted(by_lane.items()):
    if mates.get("1") != mates.get("2"):
        raise SystemExit(f"lane {lane} mate record-count mismatch: {mates}")
print(f"PASS: {sum(mates['1'] for mates in by_lane.values())} paired records across 8 lanes")
PY
