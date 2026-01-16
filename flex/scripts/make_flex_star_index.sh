#!/bin/bash

# make_flex_star_index.sh
# Build a STAR index that includes flex probe pseudo-chromosomes.
# Steps:
#   1) Generate probe-only FASTA/GTF from the 10x probe CSV
#   2) Concatenate base genome FASTA/GTF with probe FASTA/GTF
#   3) Run STAR --runMode genomeGenerate on the combined reference
#
# Usage:
#   ./make_flex_star_index.sh \
#     --probe-csv probes.csv \
#     --base-fasta genome.fa \
#     --base-gtf genes.gtf \
#     --out-dir ./flex_ref \
#     [--threads 24] [--sa-index-n-bases 14] [--chr-bin-n-bits 18] [--star-bin STAR] [--genome-limit 67108864000]
#
# Outputs (under --out-dir):
#   probes_only.fa / probes_only.gtf
#   genome.flex.fa / genome.flex.gtf   (base + probes)
#   star_index/                         (STAR genomeGenerate output)

set -euo pipefail

PROBE_CSV=""
BASE_FASTA=""
BASE_GTF=""
OUT_DIR=""
THREADS="${STAR_THREADS:-24}"
SA_INDEX_N_BASES="${STAR_SA_INDEX:-14}"
CHR_BIN_N_BITS="${STAR_CHR_BIN:-18}"
STAR_BIN="${STAR_BIN:-STAR}"
GENOME_LIMIT="${STAR_GENOME_LIMIT:-67108864000}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAKE_FLEX_REF="${SCRIPT_DIR}/make_flex_reference.sh"

usage() {
    cat << EOF
Usage: $0 --probe-csv <file> --base-fasta <file> --base-gtf <file> --out-dir <dir> [OPTIONS]

Required:
  --probe-csv FILE     10x probe CSV (contains probe sequences and metadata)
  --base-fasta FILE    Base genome FASTA
  --base-gtf FILE      Base genome GTF
  --out-dir DIR        Output directory for combined reference + STAR index

Options:
  --threads N          STAR threads (default: 24 or \$STAR_THREADS)
  --sa-index-n-bases N genomeSAindexNbases (default: 14 or \$STAR_SA_INDEX)
  --chr-bin-n-bits N   genomeChrBinNbits (default: 18 or \$STAR_CHR_BIN)
  --star-bin PATH      STAR binary to use (default: STAR or \$STAR_BIN)
  --genome-limit N     limitGenomeGenerateRAM (default: 67108864000 or \$STAR_GENOME_LIMIT)
  -h, --help           Show this message
EOF
    exit "${1:-0}"
}

log() {
    echo "[make_flex_star_index] $*" >&2
}

error_exit() {
    echo "[make_flex_star_index] ERROR: $*" >&2
    exit 1
}

# Parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --probe-csv) PROBE_CSV="$2"; shift 2 ;;
        --base-fasta) BASE_FASTA="$2"; shift 2 ;;
        --base-gtf) BASE_GTF="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --sa-index-n-bases) SA_INDEX_N_BASES="$2"; shift 2 ;;
        --chr-bin-n-bits) CHR_BIN_N_BITS="$2"; shift 2 ;;
        --star-bin) STAR_BIN="$2"; shift 2 ;;
        --genome-limit) GENOME_LIMIT="$2"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

# Validate
if [ -z "$PROBE_CSV" ] || [ -z "$BASE_FASTA" ] || [ -z "$BASE_GTF" ] || [ -z "$OUT_DIR" ]; then
    usage 1
fi

[ -f "$PROBE_CSV" ] || error_exit "Probe CSV not found: $PROBE_CSV"
[ -f "$BASE_FASTA" ] || error_exit "Base FASTA not found: $BASE_FASTA"
[ -f "$BASE_GTF" ] || error_exit "Base GTF not found: $BASE_GTF"
[ -x "$MAKE_FLEX_REF" ] || error_exit "make_flex_reference.sh missing or not executable: $MAKE_FLEX_REF"

if ! command -v "$STAR_BIN" >/dev/null 2>&1; then
    error_exit "STAR binary not found: $STAR_BIN"
fi

mkdir -p "$OUT_DIR"

PROBE_PREFIX="$OUT_DIR/probes_only"
COMBINED_PREFIX="$OUT_DIR/genome.flex"
STAR_INDEX_DIR="$OUT_DIR/star_index"

log "Generating probe pseudo-chromosomes..."
"$MAKE_FLEX_REF" "$PROBE_CSV" "$PROBE_PREFIX"

log "Combining base genome with probe pseudo-chromosomes..."
cat "$BASE_FASTA" "$PROBE_PREFIX.fa" > "${COMBINED_PREFIX}.fa"

# Preserve any headers in base GTF; probe GTF already has header unless suppressed
cat "$BASE_GTF" "$PROBE_PREFIX.gtf" > "${COMBINED_PREFIX}.gtf"

log "Running STAR genomeGenerate..."
mkdir -p "$STAR_INDEX_DIR"
"$STAR_BIN" --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --runThreadN "$THREADS" \
    --genomeFastaFiles "${COMBINED_PREFIX}.fa" \
    --sjdbGTFfile "${COMBINED_PREFIX}.gtf" \
    --genomeSAindexNbases "$SA_INDEX_N_BASES" \
    --genomeChrBinNbits "$CHR_BIN_N_BITS" \
    --limitGenomeGenerateRAM "$GENOME_LIMIT"

log "Done. Outputs:"
log "  Probe FASTA/GTF : ${PROBE_PREFIX}.fa / ${PROBE_PREFIX}.gtf"
log "  Combined FASTA  : ${COMBINED_PREFIX}.fa"
log "  Combined GTF    : ${COMBINED_PREFIX}.gtf"
log "  STAR index      : ${STAR_INDEX_DIR}"
