#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
STAR="$REPO_ROOT/core/legacy/source/STAR"
FIXTURE="${1:-/tmp/msk_multi_fixture_ds_autodetect_smoke}"
OUT_DIR="${2:-$REPO_ROOT/tests/mixed_chemistry_autodetect_smoke_output_$(date +%Y%m%d_%H%M%S)}"
THREADS="${MIXED_CHEM_AUTODETECT_THREADS:-4}"
READS="${MIXED_CHEM_AUTODETECT_READS:-10000}"
LARRY_FEATURES="${MIXED_CHEM_AUTODETECT_LARRY_FEATURES:-500}"
GENOME="${MIXED_CHEM_GENOME:-/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-crstar/star}"
WHITELIST="${MIXED_CHEM_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt}"
FIXTURE_MAKER="$REPO_ROOT/tests/multi_feature/create_fixture_downsampled.sh"
VALIDATOR="$REPO_ROOT/tests/multi_feature/validate_mixed_filtered_merge.py"

if [[ ! -x "$STAR" ]]; then
  echo "FATAL: STAR binary not found at $STAR" >&2
  exit 1
fi
if [[ ! -x "$FIXTURE_MAKER" ]]; then
  echo "FATAL: fixture maker missing at $FIXTURE_MAKER" >&2
  exit 1
fi
if [[ ! -x "$VALIDATOR" ]]; then
  echo "FATAL: validator missing at $VALIDATOR" >&2
  exit 1
fi
if [[ ! -d "$GENOME" ]]; then
  echo "FATAL: genome index not found at $GENOME" >&2
  exit 1
fi
if [[ ! -f "$WHITELIST" ]]; then
  echo "FATAL: whitelist not found at $WHITELIST" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

if [[ ! -d "$FIXTURE/mRNA" || ! -d "$FIXTURE/PolyIII" || ! -d "$FIXTURE/LARRY" ]]; then
  bash "$FIXTURE_MAKER" "$FIXTURE" "$READS" "$LARRY_FEATURES"
fi

GRNA_REF="$FIXTURE/ref_feature_geneBC.csv"
LARRY_REF="$FIXTURE/ref_feature_larryBC.csv"
GEX_R1=$(ls "$FIXTURE"/mRNA/*_R1_*.fastq.gz | head -n1)
GEX_R2=$(ls "$FIXTURE"/mRNA/*_R2_*.fastq.gz | head -n1)
LOG="$OUT_DIR/Log.out"
STDOUT_LOG="$OUT_DIR/star_stdout.log"
SUMMARY="$OUT_DIR/SUMMARY.txt"
REPORT="$OUT_DIR/filtered_barcode_report.tsv"
CONFIG="$OUT_DIR/multi_config.csv"

cat > "$CONFIG" <<EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_feature_ref,star_library_id
$FIXTURE/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,gex_de
$FIXTURE/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,auto,$GRNA_REF,grna_de
$FIXTURE/LARRY,DE_30KO,Custom,Custom,TRU,$LARRY_REF,larry_de
EOF

"$STAR" \
  --runMode alignReads \
  --runThreadN "$THREADS" \
  --genomeDir "$GENOME" \
  --readFilesIn "$GEX_R2" "$GEX_R1" \
  --readFilesCommand zcat \
  --pfMultiConfig "$CONFIG" \
  --defaultCrCompat yes \
  --crChemistry auto \
  --outFileNamePrefix "$OUT_DIR/" \
  --outSAMtype BAM Unsorted \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$WHITELIST" \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloFeatures Gene GeneFull \
  --soloBarcodeReadLength 0 \
  --soloCellFilter EmptyDrops_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloMultiMappers Rescue \
  --soloCbUbRequireTogether no \
  --soloCrGexFeature GeneFull \
  --soloCrMultimapRescue yes \
  > "$STDOUT_LOG" 2>&1

pass=0
fail=0
pass_msg() { echo "PASS\t$1" | tee -a "$SUMMARY"; pass=$((pass+1)); }
fail_msg() { echo "FAIL\t$1" | tee -a "$SUMMARY"; fail=$((fail+1)); }
check_file() {
  local label="$1"
  local path="$2"
  if [[ -s "$path" ]]; then pass_msg "$label"; else fail_msg "$label"; fi
}
check_grep() {
  local label="$1"
  local pattern="$2"
  local path="$3"
  if grep -q "$pattern" "$path"; then pass_msg "$label"; else fail_msg "$label"; fi
}

: > "$SUMMARY"
echo "OUT_DIR\t$OUT_DIR" >> "$SUMMARY"
echo "FIXTURE\t$FIXTURE" >> "$SUMMARY"
echo "THREADS\t$THREADS" >> "$SUMMARY"

SOLO_FILTERED="$OUT_DIR/Solo.out/GeneFull/filtered/barcodes.tsv"
MERGED_FILTERED="$OUT_DIR/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
GRNA_FILTERED=$(find "$OUT_DIR/cr_assign/CRISPR_Guide_Capture/grna_de" -path '*/filtered/barcodes.tsv' | head -n1)
LARRY_FILTERED=$(find "$OUT_DIR/cr_assign/Custom/larry_de" -path '*/filtered/barcodes.tsv' | head -n1)
AUTODETECT_PROBE=$(find "$OUT_DIR" -name '.autodetect_probe' -print -quit)
RAW_FRAC=$(grep -o 'raw_frac=[0-9.]*' "$STDOUT_LOG" | head -n1 | cut -d= -f2)

check_file "Solo GeneFull filtered barcodes exist" "$SOLO_FILTERED"
check_file "Merged filtered_feature_bc_matrix barcodes exist" "$MERGED_FILTERED"
if [[ -n "$GRNA_FILTERED" ]]; then pass_msg "grna filtered barcodes exist"; else fail_msg "grna filtered barcodes exist"; fi
if [[ -n "$LARRY_FILTERED" ]]; then pass_msg "larry filtered barcodes exist"; else fail_msg "larry filtered barcodes exist"; fi
check_grep "Log records grna auto-detect eligibility" 'library grna_de .*resolved=auto.*auto-detect eligible' "$LOG"
check_grep "Log records grna TRANSLATED_MATCH composition" 'NOTICE: auto-detect TRANSLATED_MATCH for CRISPR Guide Capture → effectiveReadNamespace=NXT, assignmentWhitelistNamespace=TRU' "$LOG"
check_grep "Log records larry explicit TRU skip" 'NOTICE: Custom star_chemistry=tru → effectiveChem=TRU (auto-detect skipped)' "$LOG"
check_grep "Log summary records TRANSLATED_MATCH for grna" $'grna_de\tCRISPR Guide Capture\tDE_30KO\tNXT\tTRANSLATED_MATCH\tOK' "$LOG"
check_grep "STDOUT records chemistry auto-detect raw_frac" 'NOTICE: chemistry auto-detect: raw_frac=' "$STDOUT_LOG"
if [[ -n "$RAW_FRAC" ]] && awk "BEGIN{exit !($RAW_FRAC <= 0.2)}"; then
  pass_msg "Auto-detect raw_frac indicates TRANSLATED_MATCH"
else
  fail_msg "Auto-detect raw_frac indicates TRANSLATED_MATCH"
fi
if ! grep -q 'WARNING: chemistry auto-detect: ambiguous' "$STDOUT_LOG"; then
  pass_msg "Auto-detect did not warn ambiguous"
else
  fail_msg "Auto-detect did not warn ambiguous"
fi
if [[ -z "$AUTODETECT_PROBE" ]]; then
  pass_msg "Auto-detect probe directory cleaned"
else
  fail_msg "Auto-detect probe directory cleaned"
fi

if [[ -s "$SOLO_FILTERED" && -s "$MERGED_FILTERED" && -n "$GRNA_FILTERED" && -n "$LARRY_FILTERED" ]]; then
  if python3 "$VALIDATOR" \
      --solo-barcodes "$SOLO_FILTERED" \
      --merged-barcodes "$MERGED_FILTERED" \
      --subset-barcodes "$GRNA_FILTERED::NXT" \
      --subset-barcodes "$LARRY_FILTERED::TRU" \
      --log "$LOG" \
      --report "$REPORT"; then
    pass_msg "Filtered barcode invariants hold"
  else
    fail_msg "Filtered barcode invariants hold"
  fi
else
  fail_msg "Filtered barcode invariants hold"
fi

echo "RESULTS\t${pass}_passed\t${fail}_failed" >> "$SUMMARY"
cat "$SUMMARY"

[[ "$fail" -eq 0 ]]
