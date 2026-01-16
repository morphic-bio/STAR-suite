#!/bin/bash

# filter_probes.sh
# Remove DEPRECATED probes from probe set CSV and generate manifest
# Usage: ./filter_probes.sh --probe-set <input.csv> --output-dir <dir> [--log <file>] [--keep-temp]

set -euo pipefail

PROBE_SET=""
OUTPUT_DIR=""
LOG_FILE=""
KEEP_TEMP=false
VERBOSE=true

usage() {
    cat << EOF
Usage: $0 --probe-set <input.csv> --output-dir <dir> [OPTIONS]

Remove DEPRECATED probes from probe set CSV and generate manifest.

Required arguments:
  --probe-set PATH    Input probe set CSV file
  --output-dir PATH   Output directory for filtered results

Optional arguments:
  --log PATH          Log file path (default: stderr)
  --keep-temp         Keep temporary files for debugging
  --quiet             Suppress informational messages
  -h, --help          Show this help message

Output files:
  <output-dir>/filtered_probe_set.csv                Filtered probe CSV
  <output-dir>/metadata/filtered_probe_manifest.json Provenance manifest
EOF
    exit "${1:-0}"
}

log() {
    if [ "$VERBOSE" = true ]; then
        if [ -n "$LOG_FILE" ]; then
            echo "[filter_probes] $*" | tee -a "$LOG_FILE" >&2
        else
            echo "[filter_probes] $*" >&2
        fi
    fi
}

error_exit() {
    echo "[filter_probes] ERROR: $*" >&2
    exit 1
}

compute_sha256() {
    local file="$1"
    if command -v sha256sum &> /dev/null; then
        sha256sum "$file" | awk '{print $1}'
    elif command -v shasum &> /dev/null; then
        shasum -a 256 "$file" | awk '{print $1}'
    else
        error_exit "Neither sha256sum nor shasum found. Cannot compute checksums."
    fi
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --probe-set) PROBE_SET="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --log) LOG_FILE="$2"; shift 2 ;;
        --keep-temp) KEEP_TEMP=true; shift ;;
        --quiet) VERBOSE=false; shift ;;
        -h|--help) usage 0 ;;
        *) echo "Error: Unknown argument '$1'" >&2; usage 1 ;;
    esac
done

if [ -z "$PROBE_SET" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments" >&2
    usage 1
fi

[ -f "$PROBE_SET" ] || error_exit "Input probe set file not found: $PROBE_SET"
mkdir -p "$OUTPUT_DIR/metadata"

PROBE_SET_ABS=$(readlink -f "$PROBE_SET")

log "Starting probe filtering"
log "Input: $PROBE_SET_ABS"
log "Output directory: $OUTPUT_DIR"

FILTERED_CSV="${OUTPUT_DIR}/filtered_probe_set.csv"
MANIFEST_JSON="${OUTPUT_DIR}/metadata/filtered_probe_manifest.json"
TEMP_DATA="${OUTPUT_DIR}/.temp_filtered_data.$$"

set +o pipefail
HEADER=$(grep -v "^#" "$PROBE_SET_ABS" | head -n 1)
set -o pipefail

echo "$HEADER" | grep -q "gene_id" || error_exit "Required column 'gene_id' not found"
echo "$HEADER" | grep -q "probe_id" || error_exit "Required column 'probe_id' not found"

PROBE_ID_COL=$(echo "$HEADER" | tr ',' '\n' | grep -n "probe_id" | cut -d: -f1)
[ -n "$PROBE_ID_COL" ] || error_exit "Could not determine probe_id column index"

log "probe_id is column $PROBE_ID_COL"

log "Extracting data rows..."
grep -v "^#" "$PROBE_SET_ABS" > "$TEMP_DATA.all_rows"
tail -n +2 "$TEMP_DATA.all_rows" > "$TEMP_DATA"
rm "$TEMP_DATA.all_rows"
[ -s "$TEMP_DATA" ] || error_exit "Failed to extract data rows"

TOTAL_PROBES=$(wc -l < "$TEMP_DATA")
log "Total probes in input: $TOTAL_PROBES"

log "Writing header..."
grep "^#" "$PROBE_SET_ABS" > "$FILTERED_CSV"
echo "$HEADER" >> "$FILTERED_CSV"

log "Filtering DEPRECATED probes..."
awk -F',' -v col="$PROBE_ID_COL" '
{
    probe_id = $col
    if (toupper(probe_id) !~ /DEPRECATED/) {
        print $0
    }
}' "$TEMP_DATA" > "$TEMP_DATA.filtered"

cat "$TEMP_DATA.filtered" >> "$FILTERED_CSV"

FILTERED_PROBES=$(wc -l < "$TEMP_DATA.filtered")
REMOVED_PROBES=$((TOTAL_PROBES - FILTERED_PROBES))

log "Filtered probes retained: $FILTERED_PROBES"
log "Probes removed (DEPRECATED): $REMOVED_PROBES"

[ "$FILTERED_PROBES" -gt 0 ] || error_exit "Filtered output is empty"

INPUT_SHA256=$(compute_sha256 "$PROBE_SET_ABS")
OUTPUT_SHA256=$(compute_sha256 "$FILTERED_CSV")

log "Input SHA256: $INPUT_SHA256"
log "Output SHA256: $OUTPUT_SHA256"

cat > "$MANIFEST_JSON" << EOF
{
  "manifest_version": "1.0",
  "generated_at": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "filter_script": "$0",
  "source_files": {
    "probe_set": {
      "path": "$PROBE_SET_ABS",
      "sha256": "$INPUT_SHA256"
    }
  },
  "filtering": {
    "method": "exclude_deprecated",
    "pattern": "DEPRECATED (case-insensitive)",
    "total_probes_input": $TOTAL_PROBES,
    "total_probes_output": $FILTERED_PROBES,
    "probes_removed": $REMOVED_PROBES,
    "removal_percentage": $(awk "BEGIN {printf \"%.2f\", ($REMOVED_PROBES/$TOTAL_PROBES)*100}")
  },
  "output_files": {
    "filtered_probe_set": {
      "path": "$(readlink -f "$FILTERED_CSV")",
      "sha256": "$OUTPUT_SHA256"
    }
  }
}
EOF

log "Manifest written to: $MANIFEST_JSON"

if [ "$KEEP_TEMP" = false ]; then
    rm -f "$TEMP_DATA" "$TEMP_DATA.filtered"
fi

cat << EOF

============================================
Probe Filtering Summary
============================================
Input file:      $PROBE_SET_ABS
Output file:     $FILTERED_CSV
Manifest:        $MANIFEST_JSON

Total probes:    $TOTAL_PROBES
Retained:        $FILTERED_PROBES
Removed:         $REMOVED_PROBES ($(awk "BEGIN {printf \"%.2f\", ($REMOVED_PROBES/$TOTAL_PROBES)*100}")%)

Filter criterion: Exclude probes with "DEPRECATED" in probe_id
============================================

EOF

exit 0
