#!/usr/bin/env bash
set -euo pipefail



PROBE_OFFSET=${PROBE_OFFSET:-68}
OUTDIR_BASE=${OUTDIR:-/storage/JAX_demux}
INPUT_DIR=${INPUT_DIR:-/mnt/pikachu/7-7/JAX_scRNAseq01}
PROBE_BARCODES=${PROBE_BARCODES:-/mnt/pikachu/process_features/tables/probe-barcodes-fixed-rna-profiling-rna.txt}
SAMPLE_MAP=${SAMPLE_MAP:-/mnt/pikachu/process_features/tables/probe-barcode-to-sample-mapping.txt}
MAX_RECORDS=${MAX_RECORDS:-0}

HASH_OUTDIR="${OUTDIR_BASE}_hash"
DIRECT_OUTDIR="${OUTDIR_BASE}_direct"

echo "Using probe_offset=$PROBE_OFFSET"
echo "Using test dir=$INPUT_DIR"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "ERROR: Test input directory not found: $INPUT_DIR" >&2
  exit 1
fi

if [[ ! -f "$PROBE_BARCODES" ]]; then
  echo "ERROR: Probe barcode table not found: $PROBE_BARCODES" >&2
  exit 1
fi
if [[ ! -f "$SAMPLE_MAP" ]]; then
  echo "ERROR: Sample mapping table not found: $SAMPLE_MAP" >&2
  exit 1
fi

# Build
make -k | cat

# Helper: cat .gz
catz() {
  if command -v zcat >/dev/null 2>&1; then zcat "$@"; else gunzip -c "$@"; fi
}

count_reads_in_dir() {
  local dir="$1"; shift || true
  local pattern="$1"; shift || true
  local total=0
  # shellcheck disable=SC2044
  for f in $(find "$dir" -type f -name "$pattern" | sort); do
    local n
    n=$(catz "$f" | wc -l)
    # divide by 4 for FASTQ
    total=$(( total + n/4 ))
  done
  echo "$total"
}

echo "Counting input records..."
input_total=$(count_reads_in_dir "$INPUT_DIR" "*_R1_*.fastq.gz")
echo "Input total records (sum over R1): $input_total"

echo "Cleaning output directories..."
rm -rf "$HASH_OUTDIR" "$DIRECT_OUTDIR"
mkdir -p "$HASH_OUTDIR" "$DIRECT_OUTDIR"

run_demux() {
  local outdir="$1"; shift
  local direct_flag="$1"; shift # "" or "--direct_search"
  echo "Running demux_fastq to $outdir ${direct_flag:+($direct_flag)}..."
  set -x
  start_ns=$(date +%s%N)
  ./demux_fastq \
    --threads 2 \
    --probe_barcodes "$PROBE_BARCODES" \
    --sample_map "$SAMPLE_MAP" \
    --outdir "$outdir" \
    --probe_read R2 \
    --probe_offset "$PROBE_OFFSET" \
    ${MAX_RECORDS:+--max_records "$MAX_RECORDS"} \
    "$INPUT_DIR" \
    ${direct_flag:+"$direct_flag"}
  end_ns=$(date +%s%N)
  set +x
  elapsed_ms=$(( (end_ns - start_ns) / 1000000 ))
  echo "Elapsed: ${elapsed_ms} ms for outdir=$outdir ${direct_flag:+($direct_flag)}"
}


run_demux "$DIRECT_OUTDIR" "--direct_search"

echo "Counting output records..."
hash_total=$(count_reads_in_dir "$HASH_OUTDIR" "*_R1_*.fastq.gz")
direct_total=$(count_reads_in_dir "$DIRECT_OUTDIR" "*_R1_*.fastq.gz")
echo "Output total (hash):   $hash_total"
echo "Output total (direct): $direct_total"

if [[ "$input_total" -eq "$hash_total" ]]; then
  echo "OK: Input and output totals match for hash run."
else
  echo "WARNING: Input and output totals differ for hash run." >&2
fi

if [[ "$input_total" -eq "$direct_total" ]]; then
  echo "OK: Input and output totals match for direct run."
else
  echo "WARNING: Input and output totals differ for direct run." >&2
fi

echo "Spot-checking per-file-set consistency..."
# Verify that for each input set, R1/R2/R3 counts (if present) are consistent
check_set_consistency() {
  local prefix="$1"
  local r1="$2" r2="$3" r3="$4"
  local c1 c2 c3
  c1=$(catz "$r1" | wc -l); c1=$((c1/4))
  if [[ -f "$r2" ]]; then c2=$(catz "$r2" | wc -l); c2=$((c2/4)); else c2=$c1; fi
  if [[ -f "$r3" ]]; then c3=$(catz "$r3" | wc -l); c3=$((c3/4)); else c3=$c1; fi
  if [[ "$c1" -ne "$c2" || "$c1" -ne "$c3" ]]; then
    echo "Set mismatch for $prefix: R1=$c1 R2=$c2 R3=$c3" >&2
  fi
}

# shellcheck disable=SC2044
for r1 in $(find "$INPUT_DIR" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" | sort); do
  base="${r1##*/}"
  pref="${base%%_R1_*}"
  r2="$INPUT_DIR/${base/_R1_/_R2_}"
  r3="$INPUT_DIR/${base/_R1_/_R3_}"
  if [[ -f "$r2" || -f "$r3" ]]; then
    check_set_consistency "$pref" "$r1" "$r2" "$r3"
  fi
done

echo "Done."


