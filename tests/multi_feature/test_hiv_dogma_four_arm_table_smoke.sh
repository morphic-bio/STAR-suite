#!/usr/bin/env bash
# HIV DOGMA four-arm downsample: GEX + ATAC + ADT FASTQ + HIV table-backed Custom in one STAR run.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR_BIN="${STAR_BIN:-$REPO_ROOT/core/legacy/source/STAR}"

ROOT="${HIV_DOGMA_ROOT:-/tmp/hiv_dogma_gse239916}"
FOUR_ARM_ROOT="${HIV_FOUR_ARM_ROOT:-$ROOT/star_four_arm_downsample}"
FASTQ_ROOT="$ROOT/star_trimodal_downsample/fastq"
GENOME="${HIV_DOGMA_GENOME:-/mnt/pikachu/autoindex_98_32/pe_index}"
GEX_WL="${HIV_DOGMA_GEX_WL:-/mnt/pikachu/GEX_whitelist/737K-arc-v1.txt}"
ATAC_WL="/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt"
ATAC2GEX="/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv"
CHROMAP_FASTA="/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
CHROMAP_IDX="/mnt/pikachu/catatac_gse288996/refs/GRCh38-arc.chromap.idx"
THREADS="${HIV_FOUR_ARM_THREADS:-16}"
MAX_READS="${HIV_FOUR_ARM_MAX_READS:-100000}"

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
wanted_names = sys.argv[2:]

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
for path in (feat_path, bar_path, mtx_path):
    if not path.exists() and not Path(str(path) + ".gz").exists():
        print(f"MISSING {path}", file=sys.stderr)
        sys.exit(2)

row_to_name = {}
with open_maybe_gz(feat_path) as handle:
    for idx, line in enumerate(handle):
        parts = line.rstrip("\n").split("\t")
        feature_id = parts[0] if parts else ""
        feature_name = parts[1] if len(parts) > 1 else feature_id
        for wanted in wanted_names:
            if feature_id == wanted or feature_name == wanted:
                row_to_name[idx] = wanted
                break

sums = {name: 0 for name in wanted_names}
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
        name = row_to_name.get(row)
        if name is not None:
            sums[name] += value

for name in wanted_names:
    print(f"{name}={sums[name]}")
print(f"total={sum(sums.values())}")
PY
}

echo "=== HIV DOGMA four-arm table-backed smoke ==="

if pgrep -x STAR >/dev/null 2>&1; then
  echo "ERROR: another STAR process is running; serialize benchmark runs" >&2
  exit 1
fi

bash "$SCRIPT_DIR/materialize_hiv_table_counts.sh" "$FOUR_ARM_ROOT" YW8

OUT_DIR="$FOUR_ARM_ROOT/star_run_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUT_DIR" "$FOUR_ARM_ROOT/chromap_tmp"
rm -rf "$FOUR_ARM_ROOT/star_tmp"

CONFIG="$FOUR_ARM_ROOT/pf_multi_config_four_arm.csv"
cat > "$CONFIG" << EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_library_id,star_feature_ref,star_whitelist,star_input_format
$ROOT/star_trimodal_downsample/gex,YW8,Gene Expression,Gene Expression,TRU,gex_yw8,,,fastq
$ROOT/star_trimodal_downsample/adt,YW8,Antibody Capture,Protein,TRU,adt_yw8,$ROOT/refs/YW8_viremic_protein_feature_ref.csv,$GEX_WL,fastq
$FOUR_ARM_ROOT/hiv_state_counts.tsv,YW8,Custom,Custom,TRU,hiv_state_yw8,$ROOT/refs/YW8_viremic_state_feature_ref.csv,$GEX_WL,table
EOF

for path in "$STAR_BIN" "$GENOME/Genome" "$FASTQ_ROOT/gex_R1.fastq.gz" "$FASTQ_ROOT/gex_R2.fastq.gz" \
  "$ROOT/star_trimodal_downsample/adt" \
  "$FASTQ_ROOT/atac_R1.fastq.gz" "$FASTQ_ROOT/atac_R3.fastq.gz" "$FASTQ_ROOT/atac_barcode.fastq.gz" \
  "$FOUR_ARM_ROOT/hiv_state_counts.tsv" "$ROOT/refs/YW8_viremic_state_feature_ref.csv" \
  "$ROOT/refs/YW8_viremic_protein_feature_ref.csv" "$GEX_WL" "$ATAC_WL" "$ATAC2GEX" "$CHROMAP_FASTA" "$CHROMAP_IDX"; do
  [[ -r "$path" ]] || { echo "ERROR: missing $path" >&2; exit 1; }
done

CMD=(
  "$STAR_BIN"
  --runThreadN "$THREADS"
  --genomeDir "$GENOME"
  --readFilesIn "$FASTQ_ROOT/gex_R2.fastq.gz" "$FASTQ_ROOT/gex_R1.fastq.gz"
  --readFilesCommand zcat
  --readMapNumber "$MAX_READS"
  --outFileNamePrefix "$OUT_DIR/"
  --outTmpDir "$FOUR_ARM_ROOT/star_tmp/"
  --outSAMtype None
  --clipAdapterType CellRanger4
  --clip3pPolyG yes
  --alignEndsType Local
  --chimSegmentMin 1000000
  --soloType CB_UMI_Simple
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
  --soloBarcodeReadLength 0
  --soloCBwhitelist "$GEX_WL"
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR
  --soloMultiMappers Unique
  --soloCellFilter EmptyDrops_CR
  --soloCbUbRequireTogether no
  --soloStrand Forward
  --soloFeatures GeneFull
  --soloCrGexFeature genefull
  --soloCrMultimapRescue yes
  --soloInlineHashMode no
  --pfMultiConfig "$CONFIG"
  --crWhitelist "$GEX_WL"
  --crChemistry TRU
  --crOutputChemistry TRU
  --crMinUmi 10
  --defaultCrCompat yes
  --dynamicThreadInterface 1
  --dynamicThreadTelemetry 1
  --crAssignConsumerThreads -1
  --crAssignSearchThreads 1
  --chromapAtacEnable 1
  --chromapAtacStartMode concurrent
  --chromapAtacReferenceFasta "$CHROMAP_FASTA"
  --chromapAtacIndex "$CHROMAP_IDX"
  --chromapAtacRead1 "$FASTQ_ROOT/atac_R1.fastq.gz"
  --chromapAtacRead2 "$FASTQ_ROOT/atac_R3.fastq.gz"
  --chromapAtacBarcode "$FASTQ_ROOT/atac_barcode.fastq.gz"
  --chromapAtacReadFormat "bc:8:23:-"
  --chromapAtacBarcodeWhitelist "$ATAC_WL"
  --chromapAtacBarcodeTranslate "$ATAC2GEX"
  --chromapAtacBarcodeTranslateFromFirst 1
  --chromapAtacOutputFormat BAM
  --chromapAtacOutputFragments "$OUT_DIR/atac_possorted.bam"
  --chromapAtacSecondaryFragments "$OUT_DIR/atac_fragments.bin"
  --chromapAtacSortBam 1
  --chromapAtacSummary "$OUT_DIR/chromap_summary.csv"
  --chromapAtacThreads "$THREADS"
  --chromapAtacLowMem 1
  --chromapAtacMacs3FragLowMem 1
  --chromapAtacTempDir "$FOUR_ARM_ROOT/chromap_tmp"
  --chromapAtacTn5ShiftMode classical
  --multiomeAtacPeakMexInline yes
  --multiomeAtacPeakMetricsTsv "$OUT_DIR/atac/atac_metrics.tsv"
  --multiomeAtacPeakMexOutDir "$OUT_DIR/atac/peak_mex"
  --multiomeAtacPeakNarrowPeak "$OUT_DIR/atac/atac_peaks.narrowPeak"
  --multiomeAtacPeakSummits "$OUT_DIR/atac/atac_summits.bed"
  --multiomeAtacPeakThreads "$THREADS"
)

printf 'Command:\n'; printf ' %q' "${CMD[@]}"; printf '\n\n'
"${CMD[@]}"

LOG="$OUT_DIR/Log.out"
PROV="$OUT_DIR/cr_assign/Custom/hiv_state_yw8/pf_library_provenance.tsv"
API="$OUT_DIR/cr_assign/Custom/hiv_state_yw8/table_feature_import.api_run.txt"
TABLE_MEX="$OUT_DIR/cr_assign/Custom/hiv_state_yw8"
RAW="$OUT_DIR/outs/raw_feature_bc_matrix"
FILTERED="$OUT_DIR/outs/filtered_feature_bc_matrix"

grep -q 'table-backed feature library_id=hiv_state_yw8' "$LOG" && pass "Log mentions table-backed HIV library" || fail "missing table-backed log"
[[ -f "$PROV" ]] && pass "HIV table provenance exists" || fail "missing provenance"
grep -q '^table_rows_read' "$PROV" && pass "provenance has table_rows_read" || fail "missing table_rows_read"
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

TABLE_SUMS=$(matrix_feature_sum "$TABLE_MEX" HIV_DNA HIV_RNA) || TABLE_SUMS=""
TABLE_TOTAL=$(echo "$TABLE_SUMS" | awk -F= '/^total=/{print $2}')
if [[ "${TABLE_TOTAL:-0}" =~ ^[0-9]+$ ]] && [[ "$TABLE_TOTAL" -gt 0 ]]; then
  pass "per-library HIV table MEX has nonzero counts (total=$TABLE_TOTAL)"
else
  fail "per-library HIV table MEX total=$TABLE_TOTAL"
fi

RAW_SUMS=$(matrix_feature_sum "$RAW" HIV_DNA HIV_RNA) || RAW_SUMS=""
RAW_TOTAL=$(echo "$RAW_SUMS" | awk -F= '/^total=/{print $2}')
if [[ "${RAW_TOTAL:-0}" =~ ^[0-9]+$ ]] && [[ "$RAW_TOTAL" -gt 0 ]]; then
  pass "merged raw_feature_bc_matrix has nonzero HIV counts (total=$RAW_TOTAL)"
else
  fail "merged raw_feature_bc_matrix HIV total=$RAW_TOTAL"
fi

FILT_SUMS=$(matrix_feature_sum "$FILTERED" HIV_DNA HIV_RNA) || FILT_SUMS=""
FILT_TOTAL=$(echo "$FILT_SUMS" | awk -F= '/^total=/{print $2}')
if [[ "${FILT_TOTAL:-0}" =~ ^[0-9]+$ ]] && [[ "$FILT_TOTAL" -gt 0 ]]; then
  pass "merged filtered_feature_bc_matrix has nonzero HIV counts (total=$FILT_TOTAL)"
else
  fail "merged filtered_feature_bc_matrix HIV total=$FILT_TOTAL"
fi

[[ -d "$FILTERED" ]] && pass "filtered_feature_bc_matrix exists" || fail "missing filtered MEX"
if [[ -f "$FILTERED/features.tsv.gz" ]]; then
  if zgrep -qF 'Custom' "$FILTERED/features.tsv.gz" && zgrep -qF 'Antibody Capture' "$FILTERED/features.tsv.gz"; then
    pass "filtered MEX has Custom and Antibody Capture features"
  else
    fail "filtered MEX missing expected feature types"
  fi
else
  fail "filtered features.tsv.gz missing"
fi

NBC=$(zgrep -c . "$FILTERED/barcodes.tsv.gz" 2>/dev/null || echo 0)
NFEAT=$(zgrep -c . "$FILTERED/features.tsv.gz" 2>/dev/null || echo 0)
echo ""
echo "Output: $OUT_DIR"
echo "Matrix: ${NFEAT} features x ${NBC} barcodes (filtered)"
echo "Provenance: $PROV"
echo "Permit telemetry: $API (feature.acquires=$FEAT_ACQ)"
echo ""
echo "Results: $PASS passed, $FAIL failed"
[[ "$FAIL" -eq 0 ]]
