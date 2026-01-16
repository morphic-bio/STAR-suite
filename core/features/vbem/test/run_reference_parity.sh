#!/bin/bash
#
# run_reference_parity.sh
# Stepwise parity harness for CellRanger-style formatting and optional Flex probe artifacts.
#
# Example:
#   ./test/run_reference_parity.sh \
#     --input-fasta /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
#     --input-gtf /path/to/gencode.v32.primary_assembly.annotation.gtf \
#     --reference-dir /mnt/pikachu/patrick_bwb_mnt/processing/genome
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DEFAULT_INPUT_FASTA="/mnt/pikachu/patrick_bwb_mnt/processing/morphic-jax/JAX_RNAseq_ExtraEmbryonic/downloads/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DEFAULT_INPUT_GTF="/mnt/pikachu/patrick_bwb_mnt/processing/morphic-jax/JAX_RNAseq_ExtraEmbryonic/downloads/gencode.v32.primary_assembly.annotation.gtf"
DEFAULT_REFERENCE_DIR="/mnt/pikachu/patrick_bwb_mnt/processing/genome"

INPUT_FASTA=""
INPUT_GTF=""
REFERENCE_DIR=""
PROBE_CSV=""
FLEX_REFERENCE_DIR=""
WORK_DIR="$REPO_ROOT/test/tmp_reference_parity"
FORCE=false
DEEP_DIFF=false
NO_BUILD=false

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Stepwise parity harness for CellRanger-style formatting (and optional Flex probe artifacts).

Required inputs (defaults provided if available on this machine):
  --input-fasta PATH        Raw genome FASTA (pre-format)
  --input-gtf PATH          Raw GTF (pre-format, with versioned gene_id)
  --reference-dir PATH      Reference directory containing genome.fa and genome.gtf

Optional:
  --probe-csv PATH          Flex probe CSV to generate probe artifacts
  --flex-reference-dir PATH Reference flex_probe_artifacts or genomeDir to compare against
  --work-dir PATH           Output working directory (default: $WORK_DIR)
  --force                   Re-run formatting even if outputs exist
  --deep-diff               Extract contig/gene_id lists for deeper mismatch hints (slow)
  --no-build                Skip building helper binaries if missing
  -h, --help                Show this help
EOF
    exit "${1:-0}"
}

error_exit() {
    echo "[reference_parity] ERROR: $*" >&2
    exit 1
}

log() {
    echo "[reference_parity] $*" >&2
}

hash_file() {
    local f="$1"
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$f" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$f" | awk '{print $1}'
    else
        echo "N/A"
    fi
}

compare_files() {
    local label="$1"
    local ref="$2"
    local out="$3"
    local kind="$4"

    if cmp -s "$ref" "$out"; then
        log "✓ ${label} matches"
        return 0
    fi

    log "✗ ${label} differs"
    log "  ref: $ref"
    log "  out: $out"
    log "  ref_sha256: $(hash_file "$ref")"
    log "  out_sha256: $(hash_file "$out")"
    log "  ref_size: $(stat -c '%s' "$ref" 2>/dev/null || wc -c < "$ref") bytes"
    log "  out_size: $(stat -c '%s' "$out" 2>/dev/null || wc -c < "$out") bytes"

    if [ "$DEEP_DIFF" = true ]; then
        if [ "$kind" = "fasta" ]; then
            deep_diff_fasta "$ref" "$out"
        elif [ "$kind" = "gtf" ]; then
            deep_diff_gtf "$ref" "$out"
        fi
    fi
    return 1
}

deep_diff_fasta() {
    local ref="$1"
    local out="$2"
    local tmp_dir="$WORK_DIR/deep_diff_fasta"
    mkdir -p "$tmp_dir"
    local ref_contigs="$tmp_dir/ref_contigs.txt"
    local out_contigs="$tmp_dir/out_contigs.txt"

    log "  Deep diff (FASTA): extracting contig lists..."
    grep "^>" "$ref" | cut -d' ' -f1 | sed 's/^>//' | sort -u > "$ref_contigs"
    grep "^>" "$out" | cut -d' ' -f1 | sed 's/^>//' | sort -u > "$out_contigs"
    log "  ref_contigs: $(wc -l < "$ref_contigs")"
    log "  out_contigs: $(wc -l < "$out_contigs")"

    if ! comm -3 "$ref_contigs" "$out_contigs" | head -n 20; then
        true
    fi
}

deep_diff_gtf() {
    local ref="$1"
    local out="$2"
    local tmp_dir="$WORK_DIR/deep_diff_gtf"
    mkdir -p "$tmp_dir"
    local ref_genes="$tmp_dir/ref_genes.txt"
    local out_genes="$tmp_dir/out_genes.txt"

    log "  Deep diff (GTF): extracting gene_id lists (slow on large files)..."
    awk '
        $0 !~ /^#/ {
            if (match($0, /gene_id "[^"]+"/)) {
                v = substr($0, RSTART + 9, RLENGTH - 10);
                print v;
            }
        }
    ' "$ref" | sort -u > "$ref_genes"

    awk '
        $0 !~ /^#/ {
            if (match($0, /gene_id "[^"]+"/)) {
                v = substr($0, RSTART + 9, RLENGTH - 10);
                print v;
            }
        }
    ' "$out" | sort -u > "$out_genes"

    log "  ref_gene_ids: $(wc -l < "$ref_genes")"
    log "  out_gene_ids: $(wc -l < "$out_genes")"
    log "  sample missing in output (first 20):"
    comm -23 "$ref_genes" "$out_genes" | head -n 20 || true
    log "  sample extra in output (first 20):"
    comm -13 "$ref_genes" "$out_genes" | head -n 20 || true
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --input-fasta) INPUT_FASTA="$2"; shift 2 ;;
        --input-gtf) INPUT_GTF="$2"; shift 2 ;;
        --reference-dir) REFERENCE_DIR="$2"; shift 2 ;;
        --probe-csv) PROBE_CSV="$2"; shift 2 ;;
        --flex-reference-dir) FLEX_REFERENCE_DIR="$2"; shift 2 ;;
        --work-dir) WORK_DIR="$2"; shift 2 ;;
        --force) FORCE=true; shift ;;
        --deep-diff) DEEP_DIFF=true; shift ;;
        --no-build) NO_BUILD=true; shift ;;
        -h|--help) usage 0 ;;
        *) echo "Unknown argument: $1" >&2; usage 1 ;;
    esac
done

if [ -z "$INPUT_FASTA" ]; then
    INPUT_FASTA="$DEFAULT_INPUT_FASTA"
fi
if [ -z "$INPUT_GTF" ]; then
    INPUT_GTF="$DEFAULT_INPUT_GTF"
fi
if [ -z "$REFERENCE_DIR" ]; then
    REFERENCE_DIR="$DEFAULT_REFERENCE_DIR"
fi

[ -f "$INPUT_FASTA" ] || error_exit "Input FASTA not found: $INPUT_FASTA"
[ -f "$INPUT_GTF" ] || error_exit "Input GTF not found: $INPUT_GTF"
[ -d "$REFERENCE_DIR" ] || error_exit "Reference dir not found: $REFERENCE_DIR"

REF_FASTA="$REFERENCE_DIR/genome.fa"
REF_GTF="$REFERENCE_DIR/genome.gtf"
[ -f "$REF_FASTA" ] || error_exit "Reference FASTA not found: $REF_FASTA"
[ -f "$REF_GTF" ] || error_exit "Reference GTF not found: $REF_GTF"

CELLRANGER_TOOL="$REPO_ROOT/test/test_cellranger_format"
FLEX_TOOL="$REPO_ROOT/test/test_flex_probe_index"

if [ "$NO_BUILD" = false ]; then
    if [ ! -x "$CELLRANGER_TOOL" ]; then
        log "Building test_cellranger_format..."
        (cd "$REPO_ROOT/source" && make test_cellranger_format)
    fi
    if [ -n "$PROBE_CSV" ] && [ ! -x "$FLEX_TOOL" ]; then
        log "Building test_flex_probe_index..."
        (cd "$REPO_ROOT/source" && make test_flex_probe_index)
    fi
fi

mkdir -p "$WORK_DIR/cellranger_ref"
OUT_FASTA="$WORK_DIR/cellranger_ref/genome.fa"
OUT_GTF="$WORK_DIR/cellranger_ref/genes.gtf"

if [ "$FORCE" = false ] && [ -s "$OUT_GTF" ]; then
    ref_size="$(stat -c '%s' "$REF_GTF" 2>/dev/null || wc -c < "$REF_GTF")"
    out_size="$(stat -c '%s' "$OUT_GTF" 2>/dev/null || wc -c < "$OUT_GTF")"
    if [ "$out_size" -lt "$ref_size" ]; then
        log "WARNING: Existing formatted GTF is smaller than reference (possible partial run)."
        log "         Consider re-running with --force."
    fi
fi

if [ "$FORCE" = true ] || [ ! -s "$OUT_FASTA" ] || [ ! -s "$OUT_GTF" ]; then
    log "Running CellRangerFormatter..."
    "$CELLRANGER_TOOL" "$INPUT_FASTA" "$INPUT_GTF" "$OUT_FASTA" "$OUT_GTF"
else
    log "Skipping CellRangerFormatter (outputs exist). Use --force to re-run."
fi

log "Comparing CellRanger-formatted outputs to reference..."
compare_files "FASTA" "$REF_FASTA" "$OUT_FASTA" "fasta" || true
compare_files "GTF" "$REF_GTF" "$OUT_GTF" "gtf" || true

if [ -n "$PROBE_CSV" ]; then
    [ -f "$PROBE_CSV" ] || error_exit "Probe CSV not found: $PROBE_CSV"
    FLEX_OUT_DIR="$WORK_DIR/flex_probe_artifacts"
    mkdir -p "$FLEX_OUT_DIR"

    log "Running FlexProbeIndex..."
    "$FLEX_TOOL" "$PROBE_CSV" "$OUT_GTF" "$OUT_FASTA" "$FLEX_OUT_DIR"

    if [ -n "$FLEX_REFERENCE_DIR" ]; then
        if [ -d "$FLEX_REFERENCE_DIR/flex_probe_artifacts" ]; then
            FLEX_REFERENCE_DIR="$FLEX_REFERENCE_DIR/flex_probe_artifacts"
        fi
        [ -d "$FLEX_REFERENCE_DIR" ] || error_exit "Flex reference dir not found: $FLEX_REFERENCE_DIR"

        log "Comparing Flex probe artifacts..."
        compare_files "filtered_probe_set.csv" "$FLEX_REFERENCE_DIR/filtered_probe_set.csv" "$FLEX_OUT_DIR/filtered_probe_set.csv" "text" || true
        compare_files "probes_only.fa" "$FLEX_REFERENCE_DIR/probes_only.fa" "$FLEX_OUT_DIR/probes_only.fa" "fasta" || true
        compare_files "probes_only.gtf" "$FLEX_REFERENCE_DIR/probes_only.gtf" "$FLEX_OUT_DIR/probes_only.gtf" "gtf" || true
        compare_files "genome.filtered.fa" "$FLEX_REFERENCE_DIR/genome.filtered.fa" "$FLEX_OUT_DIR/genome.filtered.fa" "fasta" || true
        compare_files "genes.filtered.gtf" "$FLEX_REFERENCE_DIR/genes.filtered.gtf" "$FLEX_OUT_DIR/genes.filtered.gtf" "gtf" || true
        compare_files "probe_list.txt" "$FLEX_REFERENCE_DIR/probe_list.txt" "$FLEX_OUT_DIR/probe_list.txt" "text" || true
    fi
fi

log "Done."
