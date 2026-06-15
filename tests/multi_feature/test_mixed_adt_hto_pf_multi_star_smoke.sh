#!/usr/bin/env bash
# STAR pf-multi E2E: GEX Solo + mixed ADT+HTO assignBarcodes hash demux merge.
# Asserts merged outs/raw_feature_bc_matrix and outs/filtered_feature_bc_matrix
# contain both Antibody Capture (protein/) and Multiplexing Capture (hash/) rows.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR_BIN="${STAR_BIN:-$REPO_ROOT/core/legacy/source/STAR}"

GENOME="${MIXED_ADT_HTO_GENOME:-/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star}"
FIXTURE="${MIXED_ADT_HTO_FIXTURE:-/tmp/msk_multi_fixture_ds}"
WL="${MIXED_ADT_HTO_WL:-/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt}"
THREADS="${MIXED_ADT_HTO_THREADS:-4}"
MAX_READS="${MIXED_ADT_HTO_MAX_READS:-5000}"
PILOT_BC_COUNT="${MIXED_ADT_HTO_PILOT_BC_COUNT:-3}"
OUT_DIR="${MIXED_ADT_HTO_OUT:-/tmp/mixed_adt_hto_pf_multi_star_smoke_$(date +%Y%m%d_%H%M%S)}"

PASS=0
FAIL=0
pass() { echo "  PASS: $1"; PASS=$((PASS + 1)); }
fail() { echo "  FAIL: $1"; FAIL=$((FAIL + 1)); }

mex_feature_type_counts() {
  local mex_dir="$1"
  python3 - "$mex_dir" <<'PY'
import gzip
import sys
from pathlib import Path

mex_dir = Path(sys.argv[1])

def open_maybe_gz(path: Path):
    if path.exists():
        return path.open("rt")
    gz = Path(str(path) + ".gz")
    if gz.exists():
        return gzip.open(gz, "rt")
    raise FileNotFoundError(path)

feat_path = mex_dir / "features.tsv"
mtx_path = mex_dir / "matrix.mtx"
for path in (feat_path, mtx_path):
    if not path.exists() and not Path(str(path) + ".gz").exists():
        print(f"MISSING {path}", file=sys.stderr)
        sys.exit(2)

row_types = []
row_names = []
with open_maybe_gz(feat_path) as handle:
    for line in handle:
        parts = line.rstrip("\n").split("\t")
        row_names.append(parts[1] if len(parts) > 1 else parts[0])
        row_types.append(parts[2] if len(parts) > 2 else "unknown")

type_counts = {}
name_counts = {name: 0 for name in row_names}
with open_maybe_gz(mtx_path) as handle:
    saw_shape = False
    for line in handle:
        if not line.strip() or line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) != 3:
            continue
        if not saw_shape:
            saw_shape = True
            continue
        row = int(parts[0]) - 1
        value = int(parts[2])
        if 0 <= row < len(row_types):
            ftype = row_types[row]
            type_counts[ftype] = type_counts.get(ftype, 0) + value
            name_counts[row_names[row]] += value

for ftype, count in sorted(type_counts.items()):
    print(f"type:{ftype}={count}")
for name, count in sorted(name_counts.items()):
    print(f"name:{name}={count}")
print(f"feature_rows={len(row_names)}")
PY
}

echo "=== Mixed ADT+HTO pf-multi STAR merge smoke ==="

if pgrep -x STAR >/dev/null 2>&1; then
  echo "ERROR: another STAR process is running; serialize benchmark runs" >&2
  exit 1
fi

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
ADT_DIR="$OUT_DIR/adt_fastqs"
PILOT_DIR="$OUT_DIR/gex_pilot"
BC_CACHE="$OUT_DIR/pilot_barcodes_${MAX_READS}.txt"
MIXED_REF="$OUT_DIR/mixed_adt_hto_ref.csv"
CONFIG="$OUT_DIR/pf_multi_config.csv"

cat > "$MIXED_REF" <<'EOF'
id,name,sequence,feature_type
hashtag1,hashtag1,ATCGATCGATCGATCG,Antibody Capture
hashtag2,hashtag2,TTAATTAATTAATTAA,Antibody Capture
CD29,CD29,GGGGGGGGGGGGGGGG,Antibody Capture
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
  echo "  Pilot: GEX-only run to discover filtered barcodes for feature FASTQs ..."
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
  head -n "$PILOT_BC_COUNT" "$PILOT_BC" > "$BC_CACHE"
  echo "  Cached ${PILOT_BC_COUNT} pilot barcodes -> $BC_CACHE"
fi

mkdir -p "$ADT_DIR"
python3 - "$BC_CACHE" "$ADT_DIR" <<'PY'
import os
import sys

cache_path, out_dir = sys.argv[1], sys.argv[2]
barcodes = []
with open(cache_path) as fh:
    for line in fh:
        bc = line.strip()
        if not bc:
            continue
        if "-" in bc:
            bc = bc.split("-", 1)[0]
        barcodes.append(bc)
if len(barcodes) < 3:
    raise SystemExit(f"need at least 3 pilot barcodes, got {len(barcodes)}")

features = {
    "hashtag1": "ATCGATCGATCGATCG",
    "hashtag2": "TTAATTAATTAATTAA",
    "CD29": "GGGGGGGGGGGGGGGG",
}
reads = [
    (0, "hashtag1", "AAAAAAAAAAAA"),
    (0, "hashtag1", "AAAAAAAAAAAC"),
    (1, "hashtag2", "CCCCCCCCCCCC"),
    (1, "hashtag2", "CCCCCCCCCCCA"),
    (2, "CD29", "AAAAAAAAAAAA"),
    (2, "CD29", "AAAAAAAAAAAC"),
]
qual = lambda n: "I" * n
r1_path = os.path.join(out_dir, "adt_R1_001.fastq")
r2_path = os.path.join(out_dir, "adt_R2_001.fastq")
with open(r1_path, "w") as r1, open(r2_path, "w") as r2:
    for i, (bc_idx, feat, umi) in enumerate(reads):
        bc = barcodes[bc_idx]
        seq = features[feat]
        r1.write(f"@read{i}\n{bc}{umi}\n+\n{qual(len(bc) + len(umi))}\n")
        r2.write(f"@read{i}\n{seq}{'A' * 12}\n+\n{qual(len(seq) + 12)}\n")
print("barcodes", barcodes[:3])
PY

cat > "$CONFIG" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_library_id,star_feature_ref,star_whitelist,star_hash_demux,star_hash_feature_selector,star_hash_min_total,star_hash_min_top,star_hash_min_ratio,star_input_format
$FIXTURE/mRNA,S1,Gene Expression,Gene Expression,TRU,gex_s1,,$WL,,,,,,fastq
$ADT_DIR,S1,Antibody Capture,Protein,TRU,adt_hto_s1,$MIXED_REF,$WL,auto,id_prefix:hashtag,1,1,2.0,fastq

[feature]
ref,$MIXED_REF
EOF

echo "  STAR:    $STAR_BIN"
echo "  Genome:  $GENOME"
echo "  Fixture: $FIXTURE"
echo "  Output:  $OUT_DIR"
echo ""

run_gex_star "$OUT_DIR" \
  --pfMultiConfig "$CONFIG" \
  --defaultCrCompat yes \
  --crChemistry auto \
  --dynamicThreadInterface 1 \
  --crAssignConsumerThreads -1 \
  --crAssignSearchThreads 1 \
  2>&1 | tee "$OUT_DIR/star_stdout.log"

LOG="$OUT_DIR/Log.out"
ASSIGN_ROOT="$OUT_DIR/cr_assign/Protein/adt_hto_s1"
RAW_MEX="$OUT_DIR/outs/raw_feature_bc_matrix"
FILTERED_MEX="$OUT_DIR/outs/filtered_feature_bc_matrix"

grep -q 'pf-multi merged feature breakdown' "$LOG" && pass "Log records pf-multi merge breakdown" || fail "missing merge breakdown in Log.out"

find_sample_dir() {
  local root="$1"
  if [[ -d "$root/hash" && -d "$root/protein" ]]; then
    echo "$root"
    return
  fi
  for sub in "$root"/*; do
    if [[ -d "$sub" ]]; then
      find_sample_dir "$sub"
      return
    fi
  done
}

ASSIGN_SAMPLE="$(find_sample_dir "$ASSIGN_ROOT" || true)"
if [[ -n "$ASSIGN_SAMPLE" ]]; then
  pass "assign output has hash/ and protein/ split MEX dirs"
else
  fail "missing hash/ and protein/ under $ASSIGN_ROOT"
fi

for mex_dir in "$RAW_MEX" "$FILTERED_MEX"; do
  label=$(basename "$mex_dir")
  if [[ ! -d "$mex_dir" ]]; then
    fail "missing $label"
    continue
  fi
  COUNTS=$(mex_feature_type_counts "$mex_dir") || COUNTS=""
  AB=$(echo "$COUNTS" | awk -F= '/^type:Antibody Capture=/{print $2}')
  MUX=$(echo "$COUNTS" | awk -F= '/^type:Multiplexing Capture=/{print $2}')
  H1=$(echo "$COUNTS" | awk -F= '/^name:hashtag1=/{print $2}')
  H2=$(echo "$COUNTS" | awk -F= '/^name:hashtag2=/{print $2}')
  CD29=$(echo "$COUNTS" | awk -F= '/^name:CD29=/{print $2}')
  if [[ "${AB:-0}" =~ ^[0-9]+$ ]] && [[ "$AB" -gt 0 ]]; then
    pass "$label has nonzero Antibody Capture counts ($AB)"
  else
    fail "$label Antibody Capture total=${AB:-missing}"
  fi
  if [[ "${MUX:-0}" =~ ^[0-9]+$ ]] && [[ "$MUX" -gt 0 ]]; then
    pass "$label has nonzero Multiplexing Capture counts ($MUX)"
  else
    fail "$label Multiplexing Capture total=${MUX:-missing}"
  fi
  if [[ "${H1:-0}" =~ ^[0-9]+$ ]] && [[ "$H1" -gt 0 ]] && [[ "${H2:-0}" =~ ^[0-9]+$ ]] && [[ "$H2" -gt 0 ]]; then
    pass "$label merges hash rows hashtag1/hashtag2"
  else
    fail "$label hash row totals hashtag1=${H1:-0} hashtag2=${H2:-0}"
  fi
  if [[ "${CD29:-0}" =~ ^[0-9]+$ ]] && [[ "$CD29" -gt 0 ]]; then
    pass "$label merges protein-only CD29 row ($CD29)"
  else
    fail "$label CD29 total=${CD29:-missing}"
  fi
  if zgrep -qF 'Antibody Capture' "$mex_dir/features.tsv.gz" \
      && zgrep -qF 'Multiplexing Capture' "$mex_dir/features.tsv.gz"; then
    pass "$label features.tsv.gz lists both capture types"
  else
    fail "$label features.tsv.gz missing Antibody or Multiplexing Capture"
  fi
done

echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ "$FAIL" -eq 0 ]] || exit 1
