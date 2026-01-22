#!/usr/bin/env bash
# downsample_fastq.sh
# Usage: ./downsample_fastq.sh <nReads> [dir]
# Truncates each FASTQ(.gz) in the directory to the first nReads (1 read = 4 lines)
# Outputs go to ./downsampled with the same basenames and extensions.

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <nReads> [dir]" >&2
  exit 1
fi

nReads="$1"
[[ "$nReads" =~ ^[0-9]+$ ]] || { echo "nReads must be a non-negative integer"; exit 1; }
dir="${2:-.}"

outdir="$dir/downsampled"
mkdir -p "$outdir"

# Number of lines to keep (4 lines per read)
lines=$(( nReads * 4 ))

# Find FASTQ files (plain or gzipped) in the top level of the directory
while IFS= read -r -d '' f; do
  base="$(basename "$f")"
  dest="$outdir/$base"

  # Decide how to read and write based on extension
  if [[ "$f" == *.gz ]]; then
    # Input is gzipped; keep output gzipped
    if (( lines > 0 )); then
      # Avoid set -o pipefail exiting on gzip SIGPIPE by masking gzip status
      ( gzip -dc -- "$f" || true ) | head -n "$lines" | gzip -c > "$dest"
    else
      # nReads = 0 -> empty gz file
      : | gzip -c > "$dest"
    fi
  else
    # Plain FASTQ; keep output plain
    if (( lines > 0 )); then
      head -n "$lines" -- "$f" > "$dest"
    else
      : > "$dest"
    fi
  fi

  echo "Wrote: $dest"
done < <(find "$dir" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) -print0)

echo "Done. Outputs in: $outdir"
