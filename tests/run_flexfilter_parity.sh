#!/usr/bin/env bash
set -euo pipefail

# Flex EmptyDrops regression check against a known production MEX + summary.
# Builds a sample whitelist from flexfilter_debug.log if no whitelist is supplied.

FLEXFILTER_BIN=${FLEXFILTER_BIN:-/mnt/pikachu/STAR-suite/flex/tools/flexfilter/run_flexfilter_mex}
RAW_MEX=${RAW_MEX:-/mnt/pikachu/prod-12-28/SC2300771/Solo.out/Gene/raw}
GOLD_SUMMARY=${GOLD_SUMMARY:-/mnt/pikachu/prod-12-28/SC2300771/per_sample/flexfilter_summary.tsv}
DEBUG_LOG=${DEBUG_LOG:-/mnt/pikachu/prod-12-28/SC2300771/flex_debug/flexfilter_debug.log}
TOTAL_EXPECTED=${TOTAL_EXPECTED:-12000}
OUT_BASE=${OUT_BASE:-tests/flexfilter_parity_output_$(date +%Y%m%d_%H%M%S)}

if [[ ! -x "$FLEXFILTER_BIN" ]]; then
  echo "FAIL: FLEXFILTER_BIN not executable: $FLEXFILTER_BIN" >&2
  exit 1
fi
if [[ ! -d "$RAW_MEX" ]]; then
  echo "SKIP: RAW_MEX not found: $RAW_MEX" >&2
  exit 0
fi
if [[ ! -s "$GOLD_SUMMARY" ]]; then
  echo "SKIP: GOLD_SUMMARY not found: $GOLD_SUMMARY" >&2
  exit 0
fi

mkdir -p "$OUT_BASE"
OUT_BASE=$(realpath "$OUT_BASE")

SAMPLE_WHITELIST=${SAMPLE_WHITELIST:-"$OUT_BASE/sample_whitelist.tsv"}

if [[ ! -s "$SAMPLE_WHITELIST" ]]; then
  if [[ ! -s "$DEBUG_LOG" ]]; then
    echo "FAIL: No SAMPLE_WHITELIST provided and DEBUG_LOG missing: $DEBUG_LOG" >&2
    exit 1
  fi
  DEBUG_LOG="$DEBUG_LOG" python3 - <<'PY' >"$SAMPLE_WHITELIST"
import os
import re
from pathlib import Path

debug_log = Path(os.environ["DEBUG_LOG"])
text = debug_log.read_text()
pattern = re.compile(r"Processing tag: ([A-Z]+) \(sample: ([^)]+)\)")
seen = set()
for tag, sample in pattern.findall(text):
    if (sample, tag) in seen:
        continue
    seen.add((sample, tag))
    print(f"{sample}\t{tag}")
PY
fi

echo "[flexfilter_parity] FLEXFILTER_BIN=$FLEXFILTER_BIN"
echo "[flexfilter_parity] RAW_MEX=$RAW_MEX"
echo "[flexfilter_parity] GOLD_SUMMARY=$GOLD_SUMMARY"
echo "[flexfilter_parity] SAMPLE_WHITELIST=$SAMPLE_WHITELIST"
echo "[flexfilter_parity] TOTAL_EXPECTED=$TOTAL_EXPECTED"
echo "[flexfilter_parity] OUT_BASE=$OUT_BASE"

OUT_PREFIX="$OUT_BASE/output"
mkdir -p "$OUT_PREFIX"

"$FLEXFILTER_BIN" \
  --mex-dir "$RAW_MEX" \
  --output-prefix "$OUT_PREFIX" \
  --total-expected "$TOTAL_EXPECTED" \
  --sample-whitelist "$SAMPLE_WHITELIST" \
  >"$OUT_BASE/flexfilter.stdout.log" 2>"$OUT_BASE/flexfilter.stderr.log"

NEW_SUMMARY="$OUT_PREFIX/flexfilter_summary.tsv"
if [[ ! -s "$NEW_SUMMARY" ]]; then
  echo "FAIL: Missing summary: $NEW_SUMMARY" >&2
  exit 1
fi

GOLD_SUMMARY="$GOLD_SUMMARY" NEW_SUMMARY="$NEW_SUMMARY" OUT_BASE="$OUT_BASE" python3 - <<'PY'
import csv
import os
import pathlib
import sys

gold = pathlib.Path(os.environ["GOLD_SUMMARY"])
new = pathlib.Path(os.environ["NEW_SUMMARY"])
out_base = pathlib.Path(os.environ["OUT_BASE"])

def read_table(path):
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = {row["Sample"]: row for row in reader}
    return rows, reader.fieldnames

def normalize(row):
    row = dict(row)
    if "Simple_ED" not in row and "Simple" in row:
        row["Simple_ED"] = row["Simple"]
    if "ED_Pass" not in row and "Tail_Pass" in row:
        row["ED_Pass"] = row["Tail_Pass"]
    if "Occ_Rem" not in row and "Occ_Removed" in row:
        row["Occ_Rem"] = row["Occ_Removed"]
    return row

gold_rows, gold_cols = read_table(gold)
new_rows, new_cols = read_table(new)

gold_rows = {k: normalize(v) for k, v in gold_rows.items()}
new_rows = {k: normalize(v) for k, v in new_rows.items()}

missing = sorted(set(gold_rows) - set(new_rows))
extra = sorted(set(new_rows) - set(gold_rows))

diffs = []
check_cols = ["Expected","Retain","Simple_ED","Tail_Tested","ED_Pass","Occ_Rem","Final","Total_UMIs"]
for sample in sorted(gold_rows):
    if sample == "TOTAL":
        if sample not in new_rows:
            continue
    if sample not in new_rows:
        continue
    for col in check_cols:
        g = gold_rows[sample].get(col, "")
        n = new_rows[sample].get(col, "")
        if g != n:
            diffs.append((sample, col, g, n))

diff_path = out_base / "summary_diff.tsv"
with diff_path.open("w") as fh:
    fh.write("Sample\tColumn\tGold\tNew\n")
    for sample, col, g, n in diffs:
        fh.write(f"{sample}\t{col}\t{g}\t{n}\n")

meta_path = out_base / "summary_compare.txt"
meta_path.write_text(
    f"gold_samples={len(gold_rows)}\\n"
    f"new_samples={len(new_rows)}\\n"
    f"missing_samples={len(missing)}\\n"
    f"extra_samples={len(extra)}\\n"
    f"value_diffs={len(diffs)}\\n"
)

print(meta_path.read_text())
if missing:
    print("Missing samples:", ", ".join(missing))
if extra:
    print("Extra samples:", ", ".join(extra))
if diffs:
    print(f"Value diffs: {len(diffs)} (see {diff_path})")

if missing or extra or diffs:
    sys.exit(1)
PY

echo "[flexfilter_parity] PASS"
