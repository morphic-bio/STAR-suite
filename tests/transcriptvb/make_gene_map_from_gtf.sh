#!/bin/bash
# make_gene_map_from_gtf.sh
# Generate a Salmon gene map (transcript_id<TAB>gene_id) from a GTF.
#
# Usage:
#   ./make_gene_map_from_gtf.sh --gtf genes.gtf --out gene_map.tsv [--feature transcript]

set -euo pipefail

GTF=""
OUT=""
FEATURE="transcript"

usage() {
    cat <<'EOF'
Usage: make_gene_map_from_gtf.sh --gtf FILE --out FILE [--feature NAME]

Required:
  --gtf FILE         Input GTF (optionally .gz)
  --out FILE         Output gene map (tab-delimited)

Optional:
  --feature NAME     Feature type to parse (default: transcript)
  -h, --help         Show this help message
EOF
    exit "${1:-0}"
}

log() { echo "[make_gene_map] $*" >&2; }
error_exit() { echo "[make_gene_map] ERROR: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --gtf) GTF="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        --feature) FEATURE="$2"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) echo "Error: Unknown argument '$1'" >&2; usage 1 ;;
    esac
done

[ -n "$GTF" ] || usage 1
[ -n "$OUT" ] || usage 1
[ -f "$GTF" ] || error_exit "GTF not found: $GTF"

read_cmd=("cat")
if [[ "$GTF" == *.gz ]]; then
    if command -v zcat >/dev/null 2>&1; then
        read_cmd=("zcat")
    elif command -v gzip >/dev/null 2>&1; then
        read_cmd=("gzip" "-dc")
    else
        error_exit "Need zcat or gzip to read: $GTF"
    fi
fi

mkdir -p "$(dirname "$OUT")"
tmp_out="${OUT}.tmp.$$"

log "Input GTF: $GTF"
log "Feature type: $FEATURE"
log "Output: $OUT"

if ! "${read_cmd[@]}" "$GTF" | awk -F'\t' -v feature="$FEATURE" '
    BEGIN { count=0 }
    $0 ~ /^#/ { next }
    $3 != feature { next }
    {
        if (match($9, /transcript_id "([^"]+)"/, t) && match($9, /gene_id "([^"]+)"/, g)) {
            tx=t[1];
            gene=g[1];
            if (!(tx in seen)) {
                print tx "\t" gene;
                seen[tx]=1;
                count++;
            }
        }
    }
    END { if (count==0) { exit 2 } }
' > "$tmp_out"; then
    rm -f "$tmp_out"
    error_exit "Failed to generate gene map (no '${FEATURE}' entries with transcript_id/gene_id?)"
fi

if [ ! -s "$tmp_out" ]; then
    rm -f "$tmp_out"
    error_exit "Gene map is empty: $OUT"
fi

mv "$tmp_out" "$OUT"
log "Gene map written: $OUT"
log "Transcripts mapped: $(wc -l < "$OUT")"
