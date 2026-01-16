#!/bin/bash

# make_filtered_star_index.sh
# Generate STAR index from filtered reference (FASTA + GTF with DEPRECATED probes removed)
# Usage: ./make_filtered_star_index.sh --filtered-reference <dir> --output-dir <dir> [OPTIONS]

set -euo pipefail

FILTERED_REF_DIR=""
OUTPUT_DIR=""
NUM_THREADS="${STAR_THREADS:-24}"
SA_INDEX_N_BASES="${STAR_SA_INDEX:-14}"
CHR_BIN_N_BITS="${STAR_CHR_BIN:-18}"
STAR_BIN="${STAR_BIN:-STAR}"
SKIP_IF_EXISTS=false
GENOME_LIMIT="${STAR_GENOME_LIMIT:-67108864000}"

usage() {
    cat << EOF
Usage: $0 --filtered-reference <dir> --output-dir <dir> [OPTIONS]

Generate STAR genome index from filtered reference files (DEPRECATED probes removed).

Required arguments:
  --filtered-reference PATH   Directory containing filtered_reference/ artifacts
                              (output from build_filtered_reference.sh)
  --output-dir PATH           Output directory for STAR index

Optional arguments:
  --threads N                 Number of threads (default: 24 or \$STAR_THREADS)
  --sa-index-n-bases N        genomeSAindexNbases parameter (default: 14 or \$STAR_SA_INDEX)
  --chr-bin-n-bits N          genomeChrBinNbits parameter (default: 18 or \$STAR_CHR_BIN)
  --star-bin PATH             Path to STAR binary (default: STAR or \$STAR_BIN)
  --genome-limit N            limitGenomeGenerateRAM (default: 67108864000 or \$STAR_GENOME_LIMIT)
  --skip-if-exists            Skip index generation if output directory exists
  -h, --help                  Show this help message
EOF
    exit "${1:-0}"
}

log() { echo "[make_filtered_star_index] $*" >&2; }
error_exit() { echo "[make_filtered_star_index] ERROR: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --filtered-reference) FILTERED_REF_DIR="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) NUM_THREADS="$2"; shift 2 ;;
        --sa-index-n-bases) SA_INDEX_N_BASES="$2"; shift 2 ;;
        --chr-bin-n-bits) CHR_BIN_N_BITS="$2"; shift 2 ;;
        --star-bin) STAR_BIN="$2"; shift 2 ;;
        --genome-limit) GENOME_LIMIT="$2"; shift 2 ;;
        --skip-if-exists) SKIP_IF_EXISTS=true; shift ;;
        -h|--help) usage 0 ;;
        *) echo "Error: Unknown argument '$1'" >&2; usage 1 ;;
    esac
done

[ -n "$FILTERED_REF_DIR" ] || usage 1
[ -n "$OUTPUT_DIR" ] || usage 1
[ -d "$FILTERED_REF_DIR" ] || error_exit "Filtered reference directory not found: $FILTERED_REF_DIR"

FASTA_FILE="$FILTERED_REF_DIR/genome.filtered.fa"
GTF_FILE="$FILTERED_REF_DIR/genes.filtered.gtf"
MANIFEST_FILE="$FILTERED_REF_DIR/metadata/reference_manifest.json"

[ -f "$FASTA_FILE" ] || error_exit "FASTA file not found: $FASTA_FILE"
[ -f "$GTF_FILE" ] || error_exit "GTF file not found: $GTF_FILE"

if ! command -v "$STAR_BIN" >/dev/null 2>&1; then
    error_exit "STAR binary not found: $STAR_BIN"
fi

STAR_VERSION=$("$STAR_BIN" --version 2>/dev/null || echo "unknown")

log "========================================"
log "STAR Index Generation (Filtered Reference)"
log "========================================"
log "Filtered reference: $FILTERED_REF_DIR"
log "Output directory: $OUTPUT_DIR"
log "FASTA: $FASTA_FILE"
log "GTF: $GTF_FILE"
log ""
log "STAR binary: $STAR_BIN"
log "STAR version: $STAR_VERSION"
log "Threads: $NUM_THREADS"
log "SA index N bases: $SA_INDEX_N_BASES"
log "Chr bin N bits: $CHR_BIN_N_BITS"
log "Genome limit: $GENOME_LIMIT"
log ""

if [ -d "$OUTPUT_DIR" ] && [ "$SKIP_IF_EXISTS" = true ]; then
    if [ -f "$OUTPUT_DIR/SA" ] && [ -f "$OUTPUT_DIR/Genome" ]; then
        log "Output directory exists and contains index files, skipping generation"
        log "Output: $OUTPUT_DIR"
        exit 0
    fi
fi

FASTA_ABS=$(readlink -f "$FASTA_FILE")
GTF_ABS=$(readlink -f "$GTF_FILE")
OUTPUT_ABS=$(readlink -f "$OUTPUT_DIR" 2>/dev/null || echo "$OUTPUT_DIR")

mkdir -p "$OUTPUT_DIR"

"$STAR_BIN" --runMode genomeGenerate \
    --genomeDir "$OUTPUT_ABS" \
    --runThreadN "$NUM_THREADS" \
    --genomeFastaFiles "$FASTA_ABS" \
    --sjdbGTFfile "$GTF_ABS" \
    --genomeSAindexNbases "$SA_INDEX_N_BASES" \
    --genomeChrBinNbits "$CHR_BIN_N_BITS" \
    --limitGenomeGenerateRAM "$GENOME_LIMIT"

log "STAR index generation complete: $OUTPUT_ABS"
log "Copied manifest: ${MANIFEST_FILE:-none}"

if [ -f "$MANIFEST_FILE" ]; then
    cp "$MANIFEST_FILE" "$OUTPUT_DIR/manifest_copy.json"
fi
