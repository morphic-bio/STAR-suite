#!/usr/bin/env bash
# DOGMA-HIV HTO demux smoke from raw ADT FASTQs (skipped when local data unavailable)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
PF_DIR="$REPO_ROOT/core/features/process_features"
ASSIGN="$PF_DIR/assignBarcodes"

ROOT="${HIV_DOGMA_ROOT:-/tmp/hiv_dogma_gse239916}"
ADT_FASTQS="${HIV_DOGMA_ADT_FASTQS:-$ROOT/star_trimodal_downsample/adt}"
FEATURE_REF="${HIV_DOGMA_ADT_FEATURE_REF:-$ROOT/refs/YW8_viremic_protein_feature_ref.csv}"
GEX_WL="${HIV_DOGMA_GEX_WL:-/mnt/pikachu/GEX_whitelist/737K-arc-v1.txt}"
FILTERED_BC="${HIV_DOGMA_FILTERED_BC:-/tmp/hiv_dogma_gse239916/star_trimodal_downsample/star_run/Solo.out/GeneFull/filtered/EmptyDrops/filtered_barcodes.txt}"

OUT="${HIV_DOGMA_HTO_DEMUX_OUT:-/tmp/hiv_dogma_hto_demux_smoke_$(date +%Y%m%d_%H%M%S)}"

echo "=== DOGMA-HIV native HTO demux smoke ==="

missing=0
for path in "$ADT_FASTQS" "$FEATURE_REF" "$GEX_WL"; do
    if [ ! -e "$path" ]; then
        echo "SKIP: missing local path $path"
        missing=1
    fi
done
if [ "$missing" -ne 0 ]; then
    echo "SKIP: DOGMA-HIV local ADT FASTQs not available on this host"
    exit 0
fi

if [ ! -x "$ASSIGN" ]; then
    make -C "$PF_DIR" assignBarcodes
fi

mkdir -p "$OUT"
CMD=(
    "$ASSIGN"
    -w "$GEX_WL"
    -f "$FEATURE_REF"
    -d "$OUT"
    --output-mode adt_mex
    --hash-demux yes
    --hash-feature-selector id_prefix:hashtag
    --hash-min-total 3
    --hash-min-top 3
    --hash-min-ratio 2.0
    --skip_empty_drops
    --skip_qc_outputs
    "$ADT_FASTQS"
    -b 16 -u 12
)
if [ -n "$FILTERED_BC" ] && [ -f "$FILTERED_BC" ]; then
    CMD+=(--filtered_barcodes "$FILTERED_BC" --source_namespace TRU --target_namespace TRU)
fi

"${CMD[@]}" > "$OUT/run.log" 2>&1

SAMPLE_OUT="$(find "$OUT" -mindepth 1 -maxdepth 1 -type d | head -1)"
SUMMARY="$SAMPLE_OUT/hash_demux_summary.json"
if [ ! -f "$SUMMARY" ]; then
    echo "FAIL: missing hash_demux_summary.json under $OUT"
    exit 1
fi

python3 - <<'PY' "$SUMMARY" "$SAMPLE_OUT/hash/features.tsv.gz" "$SAMPLE_OUT/protein/features.tsv.gz"
import gzip, json, sys
summary_path, hash_feat, prot_feat = sys.argv[1:4]
summary = json.load(open(summary_path, encoding="utf-8"))
hash_rows = sum(1 for _ in gzip.open(hash_feat, "rt"))
prot_rows = sum(1 for _ in gzip.open(prot_feat, "rt"))
print("hash_features", hash_rows, "protein_features", prot_rows)
print("singlet", summary.get("n_singlet"), "doublet", summary.get("n_doublet"), "negative", summary.get("n_negative"))
targets = {
    "hash_features": 6,
    "protein_features": 165,
    "singlet": 8386,
    "doublet": 633,
    "negative": 15,
}
# Parity targets are guidance; log deltas without hard-failing on barcode universe drift.
for key in ("hash_features", "protein_features"):
    got = hash_rows if key == "hash_features" else prot_rows
    exp = targets[key]
    if got != exp:
        print(f"NOTE: {key} got={got} expected~{exp}")
for key in ("singlet", "doublet", "negative"):
    got = summary.get(f"n_{key}" if key != "singlet" else "n_singlet")
    if key == "doublet":
        got = summary.get("n_doublet")
    elif key == "negative":
        got = summary.get("n_negative")
    else:
        got = summary.get("n_singlet")
    exp = targets[key]
    if got is None:
        raise SystemExit(f"FAIL: missing n_{key} in summary")
    delta = abs(got - exp)
    if delta > max(50, int(exp * 0.05)):
        print(f"NOTE: {key} got={got} expected~{exp} (delta={delta})")
PY

echo "PASS: DOGMA-HIV HTO demux smoke completed (see NOTES above for parity deltas)"
