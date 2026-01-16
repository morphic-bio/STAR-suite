#!/bin/bash

# filter_probes_to_gtf.sh
# Filter 50bp gene probes to match a target GTF and build hybrid reference.
#
# Filtering rules:
#   1. 50bp A/C/G/T only (fail on invalid)
#   2. Skip DEPRECATED probes (case-insensitive)
#   3. Keep only probes whose gene_id exists in the target GTF
#   4. Stable sort by gene_id then probe_id
#
# Outputs:
#   - filtered_probe_set.csv        Filtered probes (sorted)
#   - probes_only.fa                Probe-only FASTA
#   - probes_only.gtf               Probe-only GTF (no header)
#   - genome.filtered.fa            Hybrid FASTA (base + probes)
#   - genes.filtered.gtf            Hybrid GTF (base + probes)
#   - probe_genes_exons.bed         BED file for probes
#   - probe_list.txt                Unique gene IDs
#   - metadata/reference_manifest.json
#
# Usage: ./filter_probes_to_gtf.sh --probe-set <csv> --gtf <gtf[.gz]> --base-fasta <fa> --output-dir <dir>

set -euo pipefail

PROBE_SET=""
GTF_PATH=""
BASE_FASTA=""
OUTPUT_DIR=""
ENFORCE_LENGTH=50
VERBOSE=true
FORCE=false

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Filter 50bp gene probes to match a target GTF and build hybrid reference.

Required arguments:
  --probe-set PATH     Input 50bp gene probe CSV file
  --gtf PATH           Target GTF file (.gtf or .gtf.gz)
  --base-fasta PATH    Base genome FASTA file
  --output-dir PATH    Output directory

Optional arguments:
  --enforce-length N   Expected probe length (default: 50)
  --force              Force regeneration even if outputs exist
  --quiet              Suppress informational messages
  -h, --help           Show this help message

Output files:
  <output-dir>/filtered_probe_set.csv
  <output-dir>/probes_only.fa
  <output-dir>/probes_only.gtf
  <output-dir>/genome.filtered.fa      (hybrid = base + probes)
  <output-dir>/genes.filtered.gtf      (hybrid = base + probes)
  <output-dir>/probe_genes_exons.bed
  <output-dir>/probe_list.txt
  <output-dir>/metadata/reference_manifest.json
EOF
    exit "${1:-0}"
}

log() {
    if [ "$VERBOSE" = true ]; then
        echo "[filter_probes_to_gtf] $*" >&2
    fi
}

error_exit() {
    echo "[filter_probes_to_gtf] ERROR: $*" >&2
    exit 1
}

compute_sha256() {
    local file="$1"
    if command -v sha256sum &> /dev/null; then
        sha256sum "$file" | awk '{print $1}'
    elif command -v shasum &> /dev/null; then
        shasum -a 256 "$file" | awk '{print $1}'
    else
        echo "N/A"
    fi
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --probe-set) PROBE_SET="$2"; shift 2 ;;
        --gtf) GTF_PATH="$2"; shift 2 ;;
        --base-fasta) BASE_FASTA="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --enforce-length) ENFORCE_LENGTH="$2"; shift 2 ;;
        --force) FORCE=true; shift ;;
        --quiet) VERBOSE=false; shift ;;
        -h|--help) usage 0 ;;
        *) echo "Error: Unknown argument '$1'" >&2; usage 1 ;;
    esac
done

# Validate required arguments
if [ -z "$PROBE_SET" ] || [ -z "$GTF_PATH" ] || [ -z "$BASE_FASTA" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments" >&2
    usage 1
fi

[ -f "$PROBE_SET" ] || error_exit "Probe set file not found: $PROBE_SET"
[ -f "$GTF_PATH" ] || error_exit "GTF file not found: $GTF_PATH"
[ -f "$BASE_FASTA" ] || error_exit "Base FASTA file not found: $BASE_FASTA"

mkdir -p "$OUTPUT_DIR/metadata"

PROBE_SET_ABS=$(readlink -f "$PROBE_SET")
GTF_PATH_ABS=$(readlink -f "$GTF_PATH")
BASE_FASTA_ABS=$(readlink -f "$BASE_FASTA")

log "========================================"
log "Filter Probes to GTF + Build Hybrid Ref"
log "========================================"
log "Probe set: $PROBE_SET_ABS"
log "GTF: $GTF_PATH_ABS"
log "Base FASTA: $BASE_FASTA_ABS"
log "Output directory: $OUTPUT_DIR"
log "Enforce probe length: $ENFORCE_LENGTH bp"
log ""

# Output files
FILTERED_CSV="$OUTPUT_DIR/filtered_probe_set.csv"
PROBES_ONLY_FA="$OUTPUT_DIR/probes_only.fa"
PROBES_ONLY_GTF="$OUTPUT_DIR/probes_only.gtf"
HYBRID_FA="$OUTPUT_DIR/genome.filtered.fa"
HYBRID_GTF="$OUTPUT_DIR/genes.filtered.gtf"
PROBE_BED="$OUTPUT_DIR/probe_genes_exons.bed"
PROBE_LIST="$OUTPUT_DIR/probe_list.txt"
MANIFEST="$OUTPUT_DIR/metadata/reference_manifest.json"

# Temp files
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

GTF_GENE_IDS="$TEMP_DIR/gtf_gene_ids.txt"
TEMP_PROBES="$TEMP_DIR/temp_probes.txt"
BASE_CONTIGS="$TEMP_DIR/base_contigs.txt"

# ============================================================================
# Step 1: Extract gene_ids from GTF
# ============================================================================
log "Step 1/7: Extracting gene_ids from GTF..."

# Handle .gz files
if [[ "$GTF_PATH_ABS" == *.gz ]]; then
    GTF_CAT="zcat"
else
    GTF_CAT="cat"
fi

# Extract unique gene_ids from GTF (attribute format: gene_id "ENSG...")
$GTF_CAT "$GTF_PATH_ABS" | grep -v "^#" | \
    grep -oP 'gene_id\s+"[^"]+"' | \
    sed 's/gene_id\s*"\([^"]*\)"/\1/' | \
    sort -u > "$GTF_GENE_IDS"

GTF_GENE_COUNT=$(wc -l < "$GTF_GENE_IDS")
log "  Found $GTF_GENE_COUNT unique gene_ids in GTF"

# ============================================================================
# Step 2: Extract contig names from base FASTA (for collision check)
# ============================================================================
log "Step 2/7: Extracting contig names from base FASTA..."

grep "^>" "$BASE_FASTA_ABS" | sed 's/^>//' | cut -d' ' -f1 > "$BASE_CONTIGS"
BASE_CONTIG_COUNT=$(wc -l < "$BASE_CONTIGS")
log "  Found $BASE_CONTIG_COUNT contigs in base FASTA"

# ============================================================================
# Step 3: Filter probes with combined criteria
# ============================================================================
log "Step 3/7: Filtering probes..."

# Read probe CSV header (avoid pipefail)
set +o pipefail
HEADER=$(grep -v "^#" "$PROBE_SET_ABS" | head -n 1)
set -o pipefail

if [ -z "$HEADER" ]; then
    error_exit "Could not read header from probe CSV"
fi

if ! echo "$HEADER" | grep -q "gene_id"; then
    error_exit "Required column 'gene_id' not found"
fi
if ! echo "$HEADER" | grep -q "probe_seq"; then
    error_exit "Required column 'probe_seq' not found"
fi
if ! echo "$HEADER" | grep -q "probe_id"; then
    error_exit "Required column 'probe_id' not found"
fi

# Determine column indices (1-based) - avoid pipefail
set +o pipefail
GENE_ID_COL=$(echo "$HEADER" | tr ',' '\n' | grep -n "^gene_id$" | cut -d: -f1)
PROBE_SEQ_COL=$(echo "$HEADER" | tr ',' '\n' | grep -n "^probe_seq$" | cut -d: -f1)
PROBE_ID_COL=$(echo "$HEADER" | tr ',' '\n' | grep -n "^probe_id$" | cut -d: -f1)
set -o pipefail

if [ -z "$GENE_ID_COL" ] || [ -z "$PROBE_SEQ_COL" ] || [ -z "$PROBE_ID_COL" ]; then
    error_exit "Could not determine column indices"
fi

log "  Column indices: gene_id=$GENE_ID_COL, probe_seq=$PROBE_SEQ_COL, probe_id=$PROBE_ID_COL"

# Count input probes
TOTAL_INPUT=$(grep -v "^#" "$PROBE_SET_ABS" | tail -n +2 | wc -l)
log "  Total input probes: $TOTAL_INPUT"

# Filter probes with awk
# Outputs: gene_id,probe_seq,probe_id,... (same format as input)
# Also tracks rejection reasons

# First, skip header and extract data rows (avoid pipefail issues)
set +o pipefail
grep -v "^#" "$PROBE_SET_ABS" | tail -n +2 > "$TEMP_DIR/probe_data.csv" || true
set -o pipefail

if [ ! -s "$TEMP_DIR/probe_data.csv" ]; then
    error_exit "No data rows found in probe CSV"
fi

awk -F',' -v gene_col="$GENE_ID_COL" -v seq_col="$PROBE_SEQ_COL" -v id_col="$PROBE_ID_COL" \
    -v enforce_len="$ENFORCE_LENGTH" -v gtf_file="$GTF_GENE_IDS" \
    -v stats_file="$TEMP_DIR/filter_stats.txt" -v error_file="$TEMP_DIR/filter_errors.txt" '
BEGIN {
    # Load GTF gene_ids into set
    while ((getline line < gtf_file) > 0) {
        gtf_genes[line] = 1
    }
    close(gtf_file)
    deprecated = 0
    no_match = 0
    invalid_seq = 0
    valid = 0
    had_error = 0
}
{
    gene_id = $gene_col
    probe_seq = $seq_col
    probe_id = $id_col
    
    # Check DEPRECATED (case-insensitive)
    if (toupper(probe_id) ~ /DEPRECATED/) {
        deprecated++
        next
    }
    
    # Check sequence length
    if (length(probe_seq) != enforce_len) {
        print "ERROR: line " (NR+1) ": probe_seq is " length(probe_seq) "bp, expected " enforce_len "bp" > error_file
        invalid_seq++
        had_error = 1
        exit 1
    }
    
    # Check sequence characters (A/C/G/T only)
    if (probe_seq !~ /^[ACGTacgt]+$/) {
        print "ERROR: line " (NR+1) ": probe_seq contains invalid characters (only A/C/G/T allowed)" > error_file
        invalid_seq++
        had_error = 1
        exit 1
    }
    
    # Check gene_id in GTF
    if (!(gene_id in gtf_genes)) {
        no_match++
        next
    }
    
    # Passed all filters - output the line
    valid++
    print $0
}
END {
    print deprecated ":" no_match ":" invalid_seq ":" valid > stats_file
}
' "$TEMP_DIR/probe_data.csv" > "$TEMP_PROBES"

# Check for errors in filtering
if [ -f "$TEMP_DIR/filter_errors.txt" ] && [ -s "$TEMP_DIR/filter_errors.txt" ]; then
    cat "$TEMP_DIR/filter_errors.txt" >&2
    error_exit "Probe validation failed"
fi

# Parse stats
if [ ! -f "$TEMP_DIR/filter_stats.txt" ]; then
    error_exit "Filter stats file not created"
fi

STATS_LINE=$(cat "$TEMP_DIR/filter_stats.txt")
DROPPED_DEPRECATED=$(echo "$STATS_LINE" | cut -d: -f1)
DROPPED_NO_MATCH=$(echo "$STATS_LINE" | cut -d: -f2)
DROPPED_INVALID=$(echo "$STATS_LINE" | cut -d: -f3)
TOTAL_OUTPUT=$(echo "$STATS_LINE" | cut -d: -f4)

log "  Dropped (DEPRECATED): $DROPPED_DEPRECATED"
log "  Dropped (no GTF match): $DROPPED_NO_MATCH"
log "  Dropped (invalid seq): $DROPPED_INVALID"
log "  Retained: $TOTAL_OUTPUT"

if [ "$TOTAL_OUTPUT" -eq 0 ]; then
    error_exit "0 probes matched GTF genes; check gene_id format"
fi

# ============================================================================
# Step 4: Stable sort by gene_id then probe_id
# ============================================================================
log "Step 4/7: Sorting probes (gene_id, probe_id)..."

# Sort with stable sort (-s)
sort -t',' -k"$GENE_ID_COL","$GENE_ID_COL" -k"$PROBE_ID_COL","$PROBE_ID_COL" -s "$TEMP_PROBES" > "$TEMP_DIR/sorted_probes.txt"

# Write filtered CSV with header
grep "^#" "$PROBE_SET_ABS" > "$FILTERED_CSV" || true
echo "$HEADER" >> "$FILTERED_CSV"
cat "$TEMP_DIR/sorted_probes.txt" >> "$FILTERED_CSV"

log "  Wrote: $FILTERED_CSV"

# ============================================================================
# Step 5: Generate probe-only FASTA and GTF + check collisions
# ============================================================================
log "Step 5/7: Generating probe-only FASTA and GTF..."

> "$PROBES_ONLY_FA"
> "$PROBES_ONLY_GTF"

# Track synthetic contigs for manifest
SYNTHETIC_CONTIGS_JSON="$TEMP_DIR/synthetic_contigs.json"
echo "[" > "$SYNTHETIC_CONTIGS_JSON"
FIRST_CONTIG=true

while IFS=',' read -r line; do
    # Parse CSV fields
    gene_id=$(echo "$line" | cut -d',' -f"$GENE_ID_COL")
    probe_seq=$(echo "$line" | cut -d',' -f"$PROBE_SEQ_COL")
    probe_id=$(echo "$line" | cut -d',' -f"$PROBE_ID_COL")
    region=$(echo "$line" | cut -d',' -f5)  # region is column 5
    
    # Extract gene_name from probe_id (format: ENSG|GENE|hash)
    gene_name=$(echo "$probe_id" | cut -d'|' -f2)
    seq_length=${#probe_seq}
    
    # Check contig collision with base FASTA
    if grep -Fxq "$probe_id" "$BASE_CONTIGS"; then
        error_exit "Synthetic contig '$probe_id' collides with existing contig in base FASTA"
    fi
    
    # Write FASTA entry
    echo ">$probe_id" >> "$PROBES_ONLY_FA"
    echo "$probe_seq" >> "$PROBES_ONLY_FA"
    
    # Write GTF entry (single exon spanning entire sequence)
    printf "%s\tFLEX\texon\t1\t%d\t.\t+\t.\tgene_id \"%s\"; gene_name \"%s\"; transcript_id \"%s\"; probe_id \"%s\"; region \"%s\";\n" \
        "$probe_id" "$seq_length" "$gene_name" "$gene_name" "$probe_id" "$probe_id" "$region" >> "$PROBES_ONLY_GTF"
    
    # Add to synthetic contigs JSON
    if [ "$FIRST_CONTIG" = true ]; then
        FIRST_CONTIG=false
    else
        echo "," >> "$SYNTHETIC_CONTIGS_JSON"
    fi
    printf '  {"name": "%s", "length": %d}' "$probe_id" "$seq_length" >> "$SYNTHETIC_CONTIGS_JSON"
    
done < "$TEMP_DIR/sorted_probes.txt"

echo "" >> "$SYNTHETIC_CONTIGS_JSON"
echo "]" >> "$SYNTHETIC_CONTIGS_JSON"

log "  Wrote: $PROBES_ONLY_FA"
log "  Wrote: $PROBES_ONLY_GTF"

# ============================================================================
# Step 6: Generate hybrid FASTA and GTF
# ============================================================================
log "Step 6/7: Creating hybrid reference (base + probes)..."

# Copy base and append probes
cp "$BASE_FASTA_ABS" "$HYBRID_FA"
cat "$PROBES_ONLY_FA" >> "$HYBRID_FA"
log "  Wrote: $HYBRID_FA"

# Handle .gz GTF for base
if [[ "$GTF_PATH_ABS" == *.gz ]]; then
    zcat "$GTF_PATH_ABS" > "$HYBRID_GTF"
else
    cp "$GTF_PATH_ABS" "$HYBRID_GTF"
fi
cat "$PROBES_ONLY_GTF" >> "$HYBRID_GTF"
log "  Wrote: $HYBRID_GTF"

# ============================================================================
# Step 6b: Generate probe BED and probe list
# ============================================================================
log "  Generating probe BED and gene list..."

# Generate BED file (probe coordinates on synthetic contigs)
> "$PROBE_BED"
while IFS=',' read -r line; do
    probe_id=$(echo "$line" | cut -d',' -f"$PROBE_ID_COL")
    probe_seq=$(echo "$line" | cut -d',' -f"$PROBE_SEQ_COL")
    seq_length=${#probe_seq}
    gene_name=$(echo "$probe_id" | cut -d'|' -f2)
    
    # BED format: chrom start end name score strand
    printf "%s\t0\t%d\t%s\t0\t+\n" "$probe_id" "$seq_length" "$gene_name" >> "$PROBE_BED"
done < "$TEMP_DIR/sorted_probes.txt"
log "  Wrote: $PROBE_BED"

# Generate probe list (unique gene IDs, sorted)
cut -d',' -f"$GENE_ID_COL" "$TEMP_DIR/sorted_probes.txt" | sort -u > "$PROBE_LIST"
UNIQUE_GENES=$(wc -l < "$PROBE_LIST")
log "  Wrote: $PROBE_LIST ($UNIQUE_GENES unique genes)"

# ============================================================================
# Step 7: Generate manifest
# ============================================================================
log "Step 7/7: Generating manifest..."

PROBE_CSV_SHA=$(compute_sha256 "$PROBE_SET_ABS")
GTF_SHA=$(compute_sha256 "$GTF_PATH_ABS")
BASE_FASTA_SHA=$(compute_sha256 "$BASE_FASTA_ABS")
FILTERED_CSV_SHA=$(compute_sha256 "$FILTERED_CSV")
PROBES_FA_SHA=$(compute_sha256 "$PROBES_ONLY_FA")
PROBES_GTF_SHA=$(compute_sha256 "$PROBES_ONLY_GTF")
HYBRID_FA_SHA=$(compute_sha256 "$HYBRID_FA")
HYBRID_GTF_SHA=$(compute_sha256 "$HYBRID_GTF")
PROBE_BED_SHA=$(compute_sha256 "$PROBE_BED")
PROBE_LIST_SHA=$(compute_sha256 "$PROBE_LIST")

cat > "$MANIFEST" << EOF
{
  "manifest_version": "1.0",
  "generated_at": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "generator": "$0",
  "inputs": {
    "probe_csv": {
      "path": "$PROBE_SET_ABS",
      "sha256": "$PROBE_CSV_SHA"
    },
    "gtf": {
      "path": "$GTF_PATH_ABS",
      "sha256": "$GTF_SHA"
    },
    "base_fasta": {
      "path": "$BASE_FASTA_ABS",
      "sha256": "$BASE_FASTA_SHA"
    }
  },
  "filtering": {
    "total_input_probes": $TOTAL_INPUT,
    "dropped_deprecated": $DROPPED_DEPRECATED,
    "dropped_no_gtf_match": $DROPPED_NO_MATCH,
    "dropped_invalid_seq": $DROPPED_INVALID,
    "total_output_probes": $TOTAL_OUTPUT,
    "unique_genes": $UNIQUE_GENES,
    "enforce_probe_length": $ENFORCE_LENGTH,
    "ordering": "stable_sort_gene_id_probe_id"
  },
  "synthetic_contigs": $(cat "$SYNTHETIC_CONTIGS_JSON"),
  "outputs": {
    "filtered_csv": {
      "path": "$(readlink -f "$FILTERED_CSV")",
      "sha256": "$FILTERED_CSV_SHA"
    },
    "probes_fasta": {
      "path": "$(readlink -f "$PROBES_ONLY_FA")",
      "sha256": "$PROBES_FA_SHA"
    },
    "probes_gtf": {
      "path": "$(readlink -f "$PROBES_ONLY_GTF")",
      "sha256": "$PROBES_GTF_SHA"
    },
    "hybrid_fasta": {
      "path": "$(readlink -f "$HYBRID_FA")",
      "sha256": "$HYBRID_FA_SHA"
    },
    "hybrid_gtf": {
      "path": "$(readlink -f "$HYBRID_GTF")",
      "sha256": "$HYBRID_GTF_SHA"
    },
    "probe_bed": {
      "path": "$(readlink -f "$PROBE_BED")",
      "sha256": "$PROBE_BED_SHA"
    },
    "probe_list": {
      "path": "$(readlink -f "$PROBE_LIST")",
      "sha256": "$PROBE_LIST_SHA"
    }
  }
}
EOF

log "  Wrote: $MANIFEST"

# ============================================================================
# Summary
# ============================================================================
log ""
log "========================================"
log "Summary"
log "========================================"
log "Input probes:     $TOTAL_INPUT"
log "Output probes:    $TOTAL_OUTPUT"
log "Unique genes:     $UNIQUE_GENES"
log "Dropped:"
log "  - DEPRECATED:   $DROPPED_DEPRECATED"
log "  - No GTF match: $DROPPED_NO_MATCH"
log "  - Invalid seq:  $DROPPED_INVALID"
log ""
log "Outputs:"
log "  $FILTERED_CSV"
log "  $PROBES_ONLY_FA"
log "  $PROBES_ONLY_GTF"
log "  $HYBRID_FA (use for genomeGenerate)"
log "  $HYBRID_GTF (use for genomeGenerate)"
log "  $PROBE_BED"
log "  $PROBE_LIST"
log "  $MANIFEST"
log "========================================"

exit 0

