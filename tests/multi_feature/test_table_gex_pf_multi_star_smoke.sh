#!/usr/bin/env bash
# Fast STAR-level E2E: GEX Solo + table-backed Custom pf-multi arm with production permits.
# Exercises PfMultiProcess assign, deferred filtering, provenance promotion, and MEX merge.
#
# Table barcodes are taken from a tiny GEX-only pilot (cached per fixture/read cap) so they
# overlap the EmptyDrops-called cells that feed outs/raw_feature_bc_matrix merge.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR_BIN="${STAR_BIN:-$REPO_ROOT/core/legacy/source/STAR}"

GENOME="${TABLE_GEX_GENOME:-/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star}"
FIXTURE="${TABLE_GEX_FIXTURE:-/tmp/msk_multi_fixture_ds}"
WL="${TABLE_GEX_WL:-/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt}"
THREADS="${TABLE_GEX_THREADS:-4}"
MAX_READS="${TABLE_GEX_MAX_READS:-5000}"
TABLE_BC_COUNT="${TABLE_GEX_TABLE_BC_COUNT:-5}"
OUT_DIR="${TABLE_GEX_OUT:-/tmp/table_gex_pf_multi_star_smoke_$(date +%Y%m%d_%H%M%S)}"

PASS=0
FAIL=0
pass() { echo "  PASS: $1"; PASS=$((PASS + 1)); }
fail() { echo "  FAIL: $1"; FAIL=$((FAIL + 1)); }

matrix_feature_sum() {
  local mex_dir="$1"
  shift
  python3 - "$mex_dir" "$@" <<'PY'
import gzip
import sys
from pathlib import Path

mex_dir = Path(sys.argv[1])
features = sys.argv[2:]

def open_maybe_gz(path: Path):
    if path.exists():
        return path.open("rt")
    gz = Path(str(path) + ".gz")
    if gz.exists():
        return gzip.open(gz, "rt")
    raise FileNotFoundError(path)

feat_path = mex_dir / "features.tsv"
bar_path = mex_dir / "barcodes.tsv"
mtx_path = mex_dir / "matrix.mtx"
for p in (feat_path, bar_path, mtx_path):
    if not p.exists() and not Path(str(p) + ".gz").exists():
        print("MISSING", p, file=sys.stderr)
        sys.exit(2)

with open_maybe_gz(feat_path) as f:
    feat_names = []
    for line in f:
        parts = line.rstrip("\n").split("\t")
        feat_names.append(parts[1] if len(parts) > 1 else parts[0])

wanted = {name: 0 for name in features}
with open_maybe_gz(mtx_path) as f:
    f.readline()
    f.readline()
    f.readline()
    for line in f:
        row_s, _col, val = line.split()
        idx = int(row_s) - 1
        if 0 <= idx < len(feat_names):
            name = feat_names[idx]
            if name in wanted:
                wanted[name] += int(val)

for name in features:
    print(f"{name}={wanted[name]}")
print(f"total={sum(wanted.values())}")
PY
}

echo "=== Table-backed GEX + Custom pf-multi STAR smoke ==="

if [[ ! -x "$STAR_BIN" ]]; then
  echo "ERROR: STAR binary not found at $STAR_BIN (run: make core)" >&2
  exit 1
fi

if [[ ! -d "$FIXTURE/mRNA" ]]; then
  echo "ERROR: fixture missing at $FIXTURE — run tests/multi_feature/create_fixture_downsampled.sh" >&2
  exit 1
fi
for req in "$GENOME/Genome" "$WL"; do
  [[ -r "$req" ]] || { echo "ERROR: missing $req" >&2; exit 1; }
done

GEX_R2=$(ls "$FIXTURE"/mRNA/*_R2_*.fastq.gz | head -n1)
GEX_R1=$(ls "$FIXTURE"/mRNA/*_R1_*.fastq.gz | head -n1)

mkdir -p "$OUT_DIR"
TABLE_DIR="$OUT_DIR/table_inputs"
mkdir -p "$TABLE_DIR"

FEATURE_REF="$TABLE_DIR/custom_feature_ref.csv"
TABLE_COUNTS="$TABLE_DIR/custom_counts.tsv"
BC_CACHE="$FIXTURE/.table_gex_filtered_barcodes_${MAX_READS}.txt"
PILOT_DIR="$OUT_DIR/gex_pilot"

cat > "$FEATURE_REF" << 'EOF'
id,name,sequence,feature_type
TABLE_A,TABLE_A,,Custom
TABLE_B,TABLE_B,,Custom
EOF

run_gex_star() {
  local prefix="$1"
  shift
  "$STAR_BIN" \
    --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$GENOME" \
    --readFilesIn "$GEX_R2" "$GEX_R1" \
    --readFilesCommand zcat \
    --readMapNumber "$MAX_READS" \
    --outFileNamePrefix "${prefix}/" \
    --outSAMtype None \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WL" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloFeatures GeneFull \
    --soloBarcodeReadLength 0 \
    --soloCellFilter EmptyDrops_CR \
    "$@"
}

if [[ ! -s "$BC_CACHE" ]]; then
  echo "  Pilot: GEX-only run to discover called-cell barcodes for table rows ..."
  rm -rf "$PILOT_DIR"
  mkdir -p "$PILOT_DIR"
  run_gex_star "$PILOT_DIR" \
    --defaultCrCompat yes \
    --crChemistry auto \
    2>&1 | tee "$PILOT_DIR/star_stdout.log"
  PILOT_BC="$PILOT_DIR/Solo.out/GeneFull/filtered/barcodes.tsv"
  if [[ ! -s "$PILOT_BC" ]]; then
    echo "ERROR: pilot missing filtered barcodes at $PILOT_BC" >&2
    exit 1
  fi
  head -n "$TABLE_BC_COUNT" "$PILOT_BC" > "$BC_CACHE"
  echo "  Cached ${TABLE_BC_COUNT} pilot barcodes -> $BC_CACHE"
fi

{
  echo -e "barcode\tfeature_id\tcount"
  while IFS= read -r bc; do
    [[ -n "$bc" ]] || continue
    echo -e "${bc}\tTABLE_A\t3"
    echo -e "${bc}\tTABLE_B\t1"
  done < "$BC_CACHE"
} > "$TABLE_COUNTS"

CONFIG="$OUT_DIR/pf_multi_config.csv"
cat > "$CONFIG" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_whitelist,star_feature_ref,star_library_id,star_input_format
$FIXTURE/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,$WL,,gex_de,fastq
$TABLE_COUNTS,DE_30KO,Custom,Custom,TRU,$WL,$FEATURE_REF,table_custom,table
EOF

echo "  STAR:    $STAR_BIN"
echo "  Genome:  $GENOME"
echo "  Fixture: $FIXTURE"
echo "  Output:  $OUT_DIR"
echo "  Table BC cache: $BC_CACHE"
echo ""

run_gex_star "$OUT_DIR" \
  --pfMultiConfig "$CONFIG" \
  --defaultCrCompat yes \
  --crChemistry auto \
  --dynamicThreadInterface 1 \
  --dynamicThreadTelemetry 1 \
  --crAssignConsumerThreads -1 \
  --crAssignSearchThreads 1 \
  2>&1 | tee "$OUT_DIR/star_stdout.log"

LOG="$OUT_DIR/Log.out"
PROV="$OUT_DIR/cr_assign/Custom/table_custom/pf_library_provenance.tsv"
API="$OUT_DIR/cr_assign/Custom/table_custom/table_feature_import.api_run.txt"
TABLE_MEX="$OUT_DIR/cr_assign/Custom/table_custom"
RAW_MEX="$OUT_DIR/outs/raw_feature_bc_matrix"
FILTERED_MEX="$OUT_DIR/outs/filtered_feature_bc_matrix"

grep -q 'table-backed feature library_id=table_custom' "$LOG" && pass "Log mentions table-backed library" || fail "missing table-backed log line"
[[ -f "$PROV" ]] && pass "table provenance exists" || fail "missing provenance"
grep -q '^table_rows_read' "$PROV" && pass "provenance has table_rows_read" || fail "missing table_rows_read in provenance"
[[ -f "$API" ]] && pass "table api_run telemetry exists" || fail "missing api_run"
FEAT_ACQ=$(grep -E '^dynamicPermitDelta\.feature\.acquires=' "$API" | cut -d= -f2 || echo 0)
TABLE_ACQ=$(grep -E '^table_feature_permit_acquires=' "$API" | cut -d= -f2 || echo 0)
if [[ "${TABLE_ACQ:-0}" =~ ^[0-9]+$ ]] && [[ "$TABLE_ACQ" -gt 0 ]]; then
  pass "table FEATURE permit acquires ($TABLE_ACQ)"
elif [[ "${FEAT_ACQ:-0}" =~ ^[0-9]+$ ]] && [[ "$FEAT_ACQ" -gt 0 ]]; then
  pass "telemetry FEATURE permit acquires ($FEAT_ACQ)"
else
  fail "no FEATURE permit acquires (table=$TABLE_ACQ telemetry=$FEAT_ACQ)"
fi

LIB_SUMS=$(matrix_feature_sum "$TABLE_MEX" TABLE_A TABLE_B) || LIB_SUMS=""
LIB_TOTAL=$(echo "$LIB_SUMS" | awk -F= '/^total=/{print $2}')
if [[ "${LIB_TOTAL:-0}" =~ ^[0-9]+$ ]] && [[ "$LIB_TOTAL" -gt 0 ]]; then
  pass "per-library table MEX has nonzero counts (total=$LIB_TOTAL)"
else
  fail "per-library table MEX total=$LIB_TOTAL"
fi

RAW_SUMS=$(matrix_feature_sum "$RAW_MEX" TABLE_A TABLE_B) || RAW_SUMS=""
RAW_TOTAL=$(echo "$RAW_SUMS" | awk -F= '/^total=/{print $2}')
if [[ "${RAW_TOTAL:-0}" =~ ^[0-9]+$ ]] && [[ "$RAW_TOTAL" -gt 0 ]]; then
  pass "merged raw_feature_bc_matrix has nonzero table counts (total=$RAW_TOTAL)"
else
  fail "merged raw_feature_bc_matrix table total=$RAW_TOTAL (mechanical merge hole)"
fi

FILT_SUMS=$(matrix_feature_sum "$FILTERED_MEX" TABLE_A TABLE_B) || FILT_SUMS=""
FILT_TOTAL=$(echo "$FILT_SUMS" | awk -F= '/^total=/{print $2}')
if [[ "${FILT_TOTAL:-0}" =~ ^[0-9]+$ ]] && [[ "$FILT_TOTAL" -gt 0 ]]; then
  pass "merged filtered_feature_bc_matrix has nonzero table counts (total=$FILT_TOTAL)"
else
  fail "merged filtered_feature_bc_matrix table total=$FILT_TOTAL"
fi

[[ -d "$FILTERED_MEX" ]] && pass "filtered_feature_bc_matrix exists" || fail "missing filtered MEX"
if [[ -f "$FILTERED_MEX/features.tsv.gz" ]]; then
  if zgrep -qF 'TABLE_A' "$FILTERED_MEX/features.tsv.gz" && zgrep -qF 'Gene Expression' "$FILTERED_MEX/features.tsv.gz"; then
    pass "filtered MEX merges GEX and table Custom features"
  else
    fail "filtered MEX missing expected feature types"
  fi
else
  fail "missing filtered features.tsv.gz"
fi

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ "$FAIL" -eq 0 ]] || exit 1
