#!/bin/bash

# build_filtered_reference.sh
# Complete workflow to build filtered reference from probe set (removes DEPRECATED probes)
# Generates: filtered probe CSV, FASTA, GTF, probe BED, probe list, and manifests
# Usage: ./build_filtered_reference.sh --probe-set <csv> --base-fasta <fa> --base-gtf <gtf> --work-dir <dir> [OPTIONS]
#
# NOTE: This script produces a HYBRID reference by concatenating probe sequences to the
# base genome. The output files (genome.filtered.fa, genes.filtered.gtf) contain:
#   1. The original base genome chromosomes/annotations
#   2. Probe sequences appended as additional "chromosomes"
# You must use these hybrid files (NOT the probes_only.fa/gtf alone) when building
# the STAR index for proper alignment to both genomic regions and probe targets.

set -euo pipefail

# Default values
PROBE_SET=""
BASE_FASTA=""
BASE_GTF=""
WORK_DIR="${FLEX_WORKDIR:-.}"
SKIP_FILTER=false
REUSE_FILTERED=""
FORCE=false
VERBOSE=true

# Derived paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FILTER_SCRIPT="$SCRIPT_DIR/filter_probes.sh"
MAKE_FLEX_SCRIPT="$SCRIPT_DIR/make_flex_reference.sh"
MK_PROBE_BED_SCRIPT="$SCRIPT_DIR/mk_probe_bed.py"

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Build complete filtered reference (removes DEPRECATED probes) including:
  - Filtered probe CSV
  - Hybrid FASTA (base genome + filtered probes)
  - Hybrid GTF (base genome + filtered probes)
  - Probe BED file and gene list
  - Provenance manifests

Required arguments:
  --probe-set PATH      Input probe set CSV file
  --base-fasta PATH     Base genome FASTA (will be concatenated with probe sequences)
  --base-gtf PATH       Base genome GTF (will be concatenated with probe GTF)
  --work-dir PATH       Working directory for outputs (default: current dir or \$FLEX_WORKDIR)

Optional arguments:
  --skip-filter         Skip filtering step (use existing probe set as-is)
  --reuse-filtered DIR  Reuse existing filtered artifacts from DIR
  --force               Force regeneration even if outputs exist
  --quiet               Suppress informational messages
  -h, --help            Show this help message
EOF
    exit "${1:-0}"
}

log() {
    if [ "$VERBOSE" = true ]; then
        echo "[build_filtered_reference] $*" >&2
    fi
}

error_exit() {
    echo "[build_filtered_reference] ERROR: $*" >&2
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --probe-set)
            PROBE_SET="$2"; shift 2 ;;
        --base-fasta)
            BASE_FASTA="$2"; shift 2 ;;
        --base-gtf)
            BASE_GTF="$2"; shift 2 ;;
        --work-dir)
            WORK_DIR="$2"; shift 2 ;;
        --skip-filter)
            SKIP_FILTER=true; shift ;;
        --reuse-filtered)
            REUSE_FILTERED="$2"; shift 2 ;;
        --force)
            FORCE=true; shift ;;
        --quiet)
            VERBOSE=false; shift ;;
        -h|--help)
            usage 0 ;;
        *)
            echo "Error: Unknown argument '$1'" >&2
            usage 1 ;;
    esac
done

# Validate required arguments
if [ -z "$PROBE_SET" ] || [ -z "$BASE_FASTA" ] || [ -z "$BASE_GTF" ]; then
    echo "Error: Missing required arguments (need --probe-set, --base-fasta, --base-gtf)" >&2
    usage 1
fi

[ -f "$PROBE_SET" ] || error_exit "Probe set file not found: $PROBE_SET"
[ -f "$BASE_FASTA" ] || error_exit "Base FASTA file not found: $BASE_FASTA"
[ -f "$BASE_GTF" ] || error_exit "Base GTF file not found: $BASE_GTF"
[ -x "$FILTER_SCRIPT" ] || error_exit "Filter script not found or not executable: $FILTER_SCRIPT"
[ -x "$MAKE_FLEX_SCRIPT" ] || error_exit "make_flex_reference.sh not found or not executable: $MAKE_FLEX_SCRIPT"
[ -f "$MK_PROBE_BED_SCRIPT" ] || error_exit "mk_probe_bed.py not found: $MK_PROBE_BED_SCRIPT"

OUTPUT_DIR="$WORK_DIR/filtered_reference"
METADATA_DIR="$OUTPUT_DIR/metadata"
mkdir -p "$OUTPUT_DIR" "$METADATA_DIR"

log "========================================"
log "Building Filtered Hybrid Reference"
log "========================================"
log "Probe set: $PROBE_SET"
log "Base FASTA: $BASE_FASTA"
log "Base GTF: $BASE_GTF"
log "Output directory: $OUTPUT_DIR"
log ""

PROBE_SET_ABS=$(readlink -f "$PROBE_SET")
BASE_FASTA_ABS=$(readlink -f "$BASE_FASTA")
BASE_GTF_ABS=$(readlink -f "$BASE_GTF")

# Stage 1: Filter probes (or reuse existing)
FILTERED_PROBE_CSV="$OUTPUT_DIR/filtered_probe_set.csv"

if [ -n "$REUSE_FILTERED" ] && [ -f "$REUSE_FILTERED/filtered_probe_set.csv" ]; then
    log "Reusing filtered probe set from: $REUSE_FILTERED"
    cp "$REUSE_FILTERED/filtered_probe_set.csv" "$FILTERED_PROBE_CSV"
    if [ -f "$REUSE_FILTERED/metadata/filtered_probe_manifest.json" ]; then
        cp "$REUSE_FILTERED/metadata/filtered_probe_manifest.json" "$METADATA_DIR/"
    fi
elif [ "$SKIP_FILTER" = true ]; then
    log "Skipping filter, using probe set as-is"
    cp "$PROBE_SET_ABS" "$FILTERED_PROBE_CSV"
elif [ "$FORCE" = true ] || [ ! -f "$FILTERED_PROBE_CSV" ]; then
    log "Step 1/5: Filtering DEPRECATED probes..."
    "$FILTER_SCRIPT" \
        --probe-set "$PROBE_SET_ABS" \
        --output-dir "$OUTPUT_DIR" \
        $([ "$VERBOSE" = false ] && echo "--quiet")
    log "Filtering complete"
else
    log "Using existing filtered probe set: $FILTERED_PROBE_CSV"
fi

FILTERED_PROBE_COUNT=$(grep -v "^#" "$FILTERED_PROBE_CSV" | tail -n +2 | wc -l)
log "Filtered probe count: $FILTERED_PROBE_COUNT"
log ""

# Stage 2: Generate probe FASTA and GTF from filtered probes
PROBE_ONLY_FASTA="$OUTPUT_DIR/probes_only.fa"
PROBE_ONLY_GTF="$OUTPUT_DIR/probes_only.gtf"
FILTERED_FASTA="$OUTPUT_DIR/genome.filtered.fa"
FILTERED_GTF="$OUTPUT_DIR/genes.filtered.gtf"

if [ "$FORCE" = true ] || [ ! -f "$PROBE_ONLY_FASTA" ] || [ ! -f "$PROBE_ONLY_GTF" ]; then
    log "Step 2a/5: Generating FASTA and GTF from filtered probes..."
    "$MAKE_FLEX_SCRIPT" \
        "$FILTERED_PROBE_CSV" \
        "$OUTPUT_DIR/probes_only" \
        --no-gtf-header
    log "Generated probe-only files (no GTF header for appending)"
else
    log "Using existing probe-only FASTA and GTF"
fi

NEED_HYBRID_REBUILD=false
if [ "$FORCE" = true ] || [ ! -f "$FILTERED_FASTA" ] || [ ! -f "$FILTERED_GTF" ]; then
    NEED_HYBRID_REBUILD=true
elif [ "$PROBE_ONLY_FASTA" -nt "$FILTERED_FASTA" ] || [ "$PROBE_ONLY_GTF" -nt "$FILTERED_GTF" ]; then
    log "Probe files are newer than hybrid files - rebuilding hybrid reference"
    NEED_HYBRID_REBUILD=true
fi

if [ "$NEED_HYBRID_REBUILD" = true ]; then
    log "Step 2b/5: Creating hybrid reference (copy base + append probes)..."
    log "  Copying base genome FASTA (~GB scale)..."
    cp "$BASE_FASTA_ABS" "$FILTERED_FASTA"
    log "  Appending filtered probe sequences..."
    cat "$PROBE_ONLY_FASTA" >> "$FILTERED_FASTA"
    log "Created hybrid FASTA: $FILTERED_FASTA"

    log "  Copying base genome GTF..."
    cp "$BASE_GTF_ABS" "$FILTERED_GTF"
    log "  Appending filtered probe annotations..."
    cat "$PROBE_ONLY_GTF" >> "$FILTERED_GTF"
    log "Created hybrid GTF: $FILTERED_GTF"
else
    log "Using existing hybrid FASTA and GTF"
fi

FASTA_COUNT=$(grep -c "^>" "$FILTERED_FASTA" || echo "0")
GTF_COUNT=$(grep -c "^[^#]" "$FILTERED_GTF" || echo "0")
log "FASTA entries: $FASTA_COUNT"
log "GTF entries: $GTF_COUNT"
log ""

# Stage 3: Generate probe BED file
PROBE_BED="$OUTPUT_DIR/probe_genes_exons.bed"

if [ "$FORCE" = true ] || [ ! -f "$PROBE_BED" ]; then
    log "Step 3/5: Generating probe BED file..."
    python3 "$MK_PROBE_BED_SCRIPT" \
        --probeset "$FILTERED_PROBE_CSV" \
        --gtf "$BASE_GTF_ABS" \
        --out "$PROBE_BED"
    log "Generated: $PROBE_BED"
else
    log "Using existing probe BED: $PROBE_BED"
fi

BED_COUNT=$(wc -l < "$PROBE_BED" || echo "0")
log "BED entries: $BED_COUNT"
log ""

# Stage 4: Generate probe list (gene IDs)
PROBE_LIST="$OUTPUT_DIR/probe_list.txt"

if [ "$FORCE" = true ] || [ ! -f "$PROBE_LIST" ]; then
    log "Step 4/5: Generating probe gene list..."
    python3 - <<'PYEOF' "$FILTERED_PROBE_CSV" "$PROBE_LIST"
import sys
import csv

probeset_file = sys.argv[1]
output_file = sys.argv[2]

gene_ids = set()
with open(probeset_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('#') or not line:
            continue
        if line.startswith('gene_id,'):
            continue
        parts = line.split(',')
        if len(parts) > 0:
            gene_id = parts[0].strip()
            if gene_id:
                gene_ids.add(gene_id)

with open(output_file, 'w') as f:
    for gene_id in sorted(gene_ids):
        f.write(f"{gene_id}\n")

print(f"Generated probe list with {len(gene_ids)} unique genes", file=sys.stderr)
PYEOF
    log "Generated: $PROBE_LIST"
else
    log "Using existing probe list: $PROBE_LIST"
fi

GENE_COUNT=$(wc -l < "$PROBE_LIST" || echo "0")
log "Unique genes: $GENE_COUNT"
log ""

# Stage 5: Generate reference manifest
REFERENCE_MANIFEST="$METADATA_DIR/reference_manifest.json"

log "Step 5/5: Generating reference manifest..."

compute_sha256() {
    if command -v sha256sum &> /dev/null; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum &> /dev/null; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        echo "N/A"
    fi
}

FASTA_SHA256=$(compute_sha256 "$FILTERED_FASTA")
GTF_SHA256=$(compute_sha256 "$FILTERED_GTF")
BED_SHA256=$(compute_sha256 "$PROBE_BED")
LIST_SHA256=$(compute_sha256 "$PROBE_LIST")

cat > "$REFERENCE_MANIFEST" << EOF
{
  "manifest_version": "1.0",
  "generated_at": "$(date -u +"%Y-%m-%dT%H:%M:%SZ")",
  "build_script": "$0",
  "command_line": "$0 $@",
  "reference_type": "hybrid_genome_plus_probes",
  "source_files": {
    "probe_set": {
      "path": "$PROBE_SET_ABS",
      "filtered": $([ "$SKIP_FILTER" = true ] && echo "false" || echo "true")
    },
    "base_fasta": {
      "path": "$BASE_FASTA_ABS",
      "sha256": "$(compute_sha256 "$BASE_FASTA_ABS")"
    },
    "base_gtf": {
      "path": "$BASE_GTF_ABS",
      "sha256": "$(compute_sha256 "$BASE_GTF_ABS")"
    }
  },
  "output_files": {
    "filtered_probe_set": {
      "path": "$(readlink -f "$FILTERED_PROBE_CSV")",
      "sha256": "$(compute_sha256 "$FILTERED_PROBE_CSV")",
      "probe_count": $FILTERED_PROBE_COUNT
    },
    "fasta": {
      "path": "$(readlink -f "$FILTERED_FASTA")",
      "sha256": "$FASTA_SHA256",
      "entry_count": $FASTA_COUNT
    },
    "gtf": {
      "path": "$(readlink -f "$FILTERED_GTF")",
      "sha256": "$GTF_SHA256",
      "entry_count": $GTF_COUNT
    },
    "probe_bed": {
      "path": "$(readlink -f "$PROBE_BED")",
      "sha256": "$BED_SHA256",
      "entry_count": $BED_COUNT
    },
    "probe_list": {
      "path": "$(readlink -f "$PROBE_LIST")",
      "sha256": "$LIST_SHA256",
      "gene_count": $GENE_COUNT
    }
  },
  "statistics": {
    "total_probes": $FILTERED_PROBE_COUNT,
    "unique_genes": $GENE_COUNT,
    "fasta_entries": $FASTA_COUNT,
    "gtf_entries": $GTF_COUNT,
    "bed_entries": $BED_COUNT
  }
}
EOF

log "Generated: $REFERENCE_MANIFEST"
log ""

cat << EOF

============================================
Filtered Reference Build Complete
============================================
Output directory: $OUTPUT_DIR

Files generated:
  - filtered_probe_set.csv       ($FILTERED_PROBE_COUNT probes)
  - genome.filtered.fa           ($FASTA_COUNT entries)
  - genes.filtered.gtf           ($GTF_COUNT entries)
  - probe_genes_exons.bed        ($BED_COUNT regions)
  - probe_list.txt               ($GENE_COUNT genes)

Manifests:
  - metadata/filtered_probe_manifest.json
  - metadata/reference_manifest.json

Statistics:
  - Total probes:   $FILTERED_PROBE_COUNT
  - Unique genes:   $GENE_COUNT
  - BED regions:    $BED_COUNT

Next steps:
  1. Review manifests in $METADATA_DIR
  2. Run STAR index generation (see scripts/make_filtered_star_index.sh)
  3. Integrate into pipeline workflow
============================================

EOF

exit 0
