#!/usr/bin/env bash
set -euo pipefail

# Simple test runner for demux_fastq against test_files/downsampled
# Defaults: probe_offset=68, probe_read=R2

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
cd "$repo_root"

PROBE_OFFSET=${PROBE_OFFSET:-68}
OUTDIR_BASE=${OUTDIR:-demux_test_out}
TEST_DIR=${TEST_DIR:-/mnt/pikachu/process_features/test_files/downsampled}
PROBE_BARCODES=${PROBE_BARCODES:-/mnt/pikachu/process_features/tables/probe-barcodes-fixed-rna-profiling-rna.txt}
SAMPLE_MAP=${SAMPLE_MAP:-/mnt/pikachu/process_features/tables/probe-barcode-to-sample-mapping.txt}
MAX_RECORDS=${MAX_RECORDS:-0}

HASH_OUTDIR="${OUTDIR_BASE}_hash"
DIRECT_OUTDIR="${OUTDIR_BASE}_direct"

echo "Using probe_offset=$PROBE_OFFSET"
echo "Using test dir=$TEST_DIR"

if [[ ! -d "$TEST_DIR" ]]; then
  echo "ERROR: Test input directory not found: $TEST_DIR" >&2
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
input_total=$(count_reads_in_dir "$TEST_DIR" "*_R1_*.fastq.gz")
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
    --threads 4 \
    --probe_barcodes "$PROBE_BARCODES" \
    --sample_map "$SAMPLE_MAP" \
    --outdir "$outdir" \
    --probe_read R2 \
    --probe_offset "$PROBE_OFFSET" \
    ${MAX_RECORDS:+--max_records "$MAX_RECORDS"} \
    "$TEST_DIR" \
    ${direct_flag:+"$direct_flag"}
  end_ns=$(date +%s%N)
  set +x
  elapsed_ms=$(( (end_ns - start_ns) / 1000000 ))
  echo "Elapsed: ${elapsed_ms} ms for outdir=$outdir ${direct_flag:+($direct_flag)}"
}

# Run without direct_search (hash path)
run_demux "$HASH_OUTDIR" ""

# Run with direct_search enabled
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
for r1 in $(find "$TEST_DIR" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" | sort); do
  base="${r1##*/}"
  pref="${base%%_R1_*}"
  r2="$TEST_DIR/${base/_R1_/_R2_}"
  r3="$TEST_DIR/${base/_R1_/_R3_}"
  if [[ -f "$r2" || -f "$r3" ]]; then
    check_set_consistency "$pref" "$r1" "$r2" "$r3"
  fi
done

# Test demux_bam with CBUB file
echo "Testing demux_bam with CBUB file..."
TAGS_DIR="/mnt/pikachu/process_features/test_files/tags"
if [[ -f "$TAGS_DIR/Aligned.out.bam" && -f "$TAGS_DIR/Aligned.out.cb_ub.bin" && -f "$TAGS_DIR/Aligned.out.tags.bam" ]]; then
    echo "Found CBUB test files, running comparison test..."
    
    # Create temporary directories
    LEGACY_OUT="tmp/legacy"
    CBUB_OUT="tmp/cbub"
    rm -rf tmp
    mkdir -p "$LEGACY_OUT" "$CBUB_OUT"
    
    # Test with legacy mode (using tagged BAM)
    echo "Running legacy mode with Aligned.out.tags.bam..."
    ./demux_bam --bam "$TAGS_DIR/Aligned.out.tags.bam" --outdir "$LEGACY_OUT" --sample_probes "$PROBE_BARCODES" || {
        echo "Legacy mode failed, skipping CBUB test"
    }
    
    # Test with CBUB mode
    echo "Running CBUB mode with Aligned.out.bam + Aligned.out.cb_ub.bin..."
    ./demux_bam --bam "$TAGS_DIR/Aligned.out.bam" --outdir "$CBUB_OUT" --sample_probes "$PROBE_BARCODES" --CBUB_file "$TAGS_DIR/Aligned.out.cb_ub.bin" --whitelist "/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt" || {
        echo "CBUB mode failed"
        exit 1
    }
    
    # Compare outputs
    echo "Comparing outputs..."
    files_to_compare=("matrix.mtx" "barcodes.tsv" "features.tsv" "stats.txt")
    all_match=true
    
    for file in "${files_to_compare[@]}"; do
        if [[ -f "$LEGACY_OUT/$file" && -f "$CBUB_OUT/$file" ]]; then
            if cmp -s "$LEGACY_OUT/$file" "$CBUB_OUT/$file"; then
                echo "✓ $file matches"
            else
                echo "✗ $file differs"
                all_match=false
            fi
        else
            echo "? $file missing in one or both outputs"
            all_match=false
        fi
    done
    
    if $all_match; then
        echo "✓ All output files match between legacy and CBUB modes"
    else
        echo "✗ Some output files differ between modes"
        exit 1
    fi
    
    # Cleanup
    rm -rf tmp
    echo "CBUB test completed successfully"
else
    echo "CBUB test files not found, skipping CBUB test"
fi

echo "Done."


