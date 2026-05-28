#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PF_BIN="${CBQ_PF_ADAPTER_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/cbq_pf_adapter_harness}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_pf_adapter_smoke}"
FIXTURE_DIR="$ROOT_DIR/core/features/process_features/tests/fixtures/assignbarcodes_baseline/input"

skip() {
    echo "SKIP: $*"
    exit 0
}

if [[ -n "${BQTOOLS:-}" ]]; then
    if [[ -x "$BQTOOLS" ]]; then
        BQTOOLS_BIN="$BQTOOLS"
    elif command -v "$BQTOOLS" >/dev/null 2>&1; then
        BQTOOLS_BIN="$(command -v "$BQTOOLS")"
    else
        skip "BQTOOLS is set but not executable or resolvable: $BQTOOLS"
    fi
elif [[ -x /tmp/star_suite_bqtools/bin/bqtools ]]; then
    BQTOOLS_BIN="/tmp/star_suite_bqtools/bin/bqtools"
else
    command -v bqtools >/dev/null 2>&1 || skip "bqtools not found in PATH; set BQTOOLS=/path/to/bqtools"
    BQTOOLS_BIN="$(command -v bqtools)"
fi

if [[ ! -x "$PF_BIN" ]]; then
    make -C "$ROOT_DIR/core/legacy/source" cbq-pf-adapter-harness
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{inputs,fastq,cbq}

FEATURES="$FIXTURE_DIR/features.csv"
WHITELIST="$FIXTURE_DIR/whitelist.txt"
R1_SRC="$OUT_ROOT/inputs/sample_R1_001.fastq"
R2_SRC="$OUT_ROOT/inputs/sample_R2_001.fastq"
R1_GZ="$OUT_ROOT/inputs/sample_R1_001.fastq.gz"
R2_GZ="$OUT_ROOT/inputs/sample_R2_001.fastq.gz"
CBQ="$OUT_ROOT/inputs/sample_pair.cbq"

python3 - "$R1_SRC" "$R2_SRC" <<'PY'
import sys
from pathlib import Path

r1_path = Path(sys.argv[1])
r2_path = Path(sys.argv[2])
barcodes = [
    "AAACCCAAGAAACCAT",
    "AAACCCAAGAAACCCA",
    "AAACCCAAGAAACCCT",
]
features = {
    "FeatureA": "ATCGATCGATCGATCG",
    "FeatureB": "GCTAGCTAGCTAGCTA",
    "FeatureC": "TTAATTAATTAATTAA",
}
reads = [
    (barcodes[0], "AAAAAAAAAAAA", "FeatureA"),
    (barcodes[0], "AAAAAAAAAAAA", "FeatureA"),
    (barcodes[0], "AAAAAAAACCCC", "FeatureA"),
    (barcodes[0], "CCCCCCCCCCCC", "FeatureB"),
    (barcodes[1], "GGGGGGGGGGGG", "FeatureC"),
    (barcodes[1], "TTTTTTTTTTTT", "FeatureC"),
    (barcodes[2], "ACGTACGTACGT", "FeatureB"),
]
with r1_path.open("wt", encoding="ascii") as r1, r2_path.open("wt", encoding="ascii") as r2:
    for idx, (barcode, umi, feature_name) in enumerate(reads):
        name = f"read_{idx:04d}"
        r1_seq = barcode + umi
        r2_seq = features[feature_name] + "AAAAAAAAAAAA"
        r1.write(f"@{name}\n{r1_seq}\n+\n{'I' * len(r1_seq)}\n")
        r2.write(f"@{name}\n{r2_seq}\n+\n{'I' * len(r2_seq)}\n")
PY

gzip -c "$R1_SRC" > "$R1_GZ"
gzip -c "$R2_SRC" > "$R2_GZ"

encode_pair_cbq() {
    local r1="$1"
    local r2="$2"
    local out="$3"

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" --mode cbq -o "$out" -T 2 \
        > "$OUT_ROOT/encode_pair.stdout" 2> "$OUT_ROOT/encode_pair.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" --mode cbq -T 2 \
        >> "$OUT_ROOT/encode_pair.stdout" 2>> "$OUT_ROOT/encode_pair.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" -T 2 \
        >> "$OUT_ROOT/encode_pair.stdout" 2>> "$OUT_ROOT/encode_pair.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" \
        >> "$OUT_ROOT/encode_pair.stdout" 2>> "$OUT_ROOT/encode_pair.stderr"
    [[ -s "$out" ]]
}

if ! encode_pair_cbq "$R1_SRC" "$R2_SRC" "$CBQ"; then
    echo "ERROR: failed to encode paired FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

"$PF_BIN" \
    --mode fastq \
    --barcodeFastq "$R1_GZ" \
    --featureFastq "$R2_GZ" \
    --whitelist "$WHITELIST" \
    --featureRef "$FEATURES" \
    --outputDir "$OUT_ROOT/fastq" \
    --sampleName sample \
    --barcodeLength 16 \
    --umiLength 12 \
    > "$OUT_ROOT/fastq/stdout.tsv" \
    2> "$OUT_ROOT/fastq/stderr.log"

"$PF_BIN" \
    --mode cbq \
    --readFilesIn "$CBQ" \
    --whitelist "$WHITELIST" \
    --featureRef "$FEATURES" \
    --outputDir "$OUT_ROOT/cbq" \
    --sampleName sample \
    --barcodeLength 16 \
    --umiLength 12 \
    > "$OUT_ROOT/cbq/stdout.tsv" \
    2> "$OUT_ROOT/cbq/stderr.log"

compare_file() {
    local rel="$1"
    cmp -s "$OUT_ROOT/fastq/sample/$rel" "$OUT_ROOT/cbq/sample/$rel" || {
        echo "ERROR: process_features CBQ adapter parity failed for $rel" >&2
        diff -u "$OUT_ROOT/fastq/sample/$rel" "$OUT_ROOT/cbq/sample/$rel" >&2 || true
        exit 1
    }
}

compare_file barcodes.txt
compare_file features.txt
compare_file feature_sequences.txt
compare_file matrix.mtx
compare_file feature_per_cell.csv
compare_file deduped_counts_histograms.txt

echo "PASS: process_features CBQ adapter smoke completed at $OUT_ROOT"
