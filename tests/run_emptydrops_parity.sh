#!/usr/bin/env bash
set -euo pipefail

# Compare legacy STAR EmptyDrops_CR vs libscrna backend on a reference raw MEX.
# Defaults target the A375 parity dataset; override with env vars if needed.

STAR_BIN=${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}
RAW_MEX=${RAW_MEX:-/storage/A375/star_gex_features_cr_parity_20260129_174602/Solo.out/Gene/raw}
OUT_BASE=${OUT_BASE:-tests/emptydrops_parity_output_$(date +%Y%m%d_%H%M%S)}
TOLERANCE=${TOLERANCE:-1}

if [[ ! -x "$STAR_BIN" ]]; then
  echo "FAIL: STAR_BIN not executable: $STAR_BIN" >&2
  exit 1
fi

if [[ ! -d "$RAW_MEX" ]]; then
  echo "SKIP: RAW_MEX not found: $RAW_MEX" >&2
  exit 0
fi

mkdir -p "$OUT_BASE/legacy" "$OUT_BASE/libscrna"
OUT_BASE=$(realpath "$OUT_BASE")

run_star() {
  local label=$1
  local legacy_flag=$2
  local outdir=$3
  local prefix="${outdir}/"

  mkdir -p "$outdir"
  "$STAR_BIN" --runMode soloCellFiltering "$RAW_MEX" "$prefix" \
    --soloCellFilter EmptyDrops_CR --soloEmptyDropsLegacy "$legacy_flag" \
    --outFileNamePrefix "$prefix" \
    >"$OUT_BASE/${label}.stdout.log" 2>"$OUT_BASE/${label}.stderr.log"

  if [[ ! -s "$outdir/barcodes.tsv" ]]; then
    echo "FAIL: missing or empty barcodes.tsv for $label at $outdir" >&2
    exit 1
  fi
}

echo "[emptydrops_parity] STAR_BIN=$STAR_BIN"
echo "[emptydrops_parity] RAW_MEX=$RAW_MEX"
echo "[emptydrops_parity] OUT_BASE=$OUT_BASE"
echo "[emptydrops_parity] TOLERANCE=$TOLERANCE"

run_star legacy yes "$OUT_BASE/legacy"
run_star libscrna no "$OUT_BASE/libscrna"

OUT_BASE="$OUT_BASE" python3 - <<'PY'
import os
import pathlib
import re
import sys

out_base = pathlib.Path(os.environ["OUT_BASE"])
legacy_dir = out_base / "legacy"
lib_dir = out_base / "libscrna"

legacy_barcodes = set((legacy_dir / "barcodes.tsv").read_text().splitlines())
lib_barcodes = set((lib_dir / "barcodes.tsv").read_text().splitlines())

inter = legacy_barcodes & lib_barcodes
union = legacy_barcodes | lib_barcodes
lib_only = sorted(lib_barcodes - legacy_barcodes)
legacy_only = sorted(legacy_barcodes - lib_barcodes)

def parse_simple(log_path):
    if not log_path.exists():
        return "NA"
    text = log_path.read_text()
    m = re.search(r"cellFiltering: simple: nUMImax=(\d+); nUMImin=(\d+); nCellsSimple=(\d+)", text)
    return ",".join(m.groups()) if m else "NA"

legacy_simple = parse_simple(legacy_dir / "Log.out")
lib_simple = parse_simple(lib_dir / "Log.out")

summary = out_base / "summary.tsv"
summary.write_text(
    "metric\tvalue\n"
    f"legacy_count\t{len(legacy_barcodes)}\n"
    f"libscrna_count\t{len(lib_barcodes)}\n"
    f"intersection\t{len(inter)}\n"
    f"union\t{len(union)}\n"
    f"legacy_only\t{len(legacy_only)}\n"
    f"libscrna_only\t{len(lib_only)}\n"
    f"legacy_simple(nUMImax,nUMImin,nCellsSimple)\t{legacy_simple}\n"
    f"libscrna_simple(nUMImax,nUMImin,nCellsSimple)\t{lib_simple}\n"
)

diff_path = out_base / "barcode_diff.tsv"
with diff_path.open("w") as fh:
    fh.write("source\tbarcode\n")
    for bc in legacy_only:
        fh.write(f"legacy_only\t{bc}\n")
    for bc in lib_only:
        fh.write(f"libscrna_only\t{bc}\n")

print(summary.read_text())
print(f"Wrote diff list: {diff_path}")
PY

# Enforce tolerance (lib-only + legacy-only)
diff_count=$(tail -n +2 "$OUT_BASE/barcode_diff.tsv" | wc -l | tr -d ' ')
if [[ "$diff_count" -gt "$TOLERANCE" ]]; then
  echo "FAIL: barcode diff count ($diff_count) exceeds tolerance ($TOLERANCE)" >&2
  exit 1
fi

echo "[emptydrops_parity] PASS: diff_count=$diff_count (<= $TOLERANCE)"
