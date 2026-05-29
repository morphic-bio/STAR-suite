#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTX_BIN="${FASTX_INPUT_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/fastx_input_harness}"
CBQ_BIN="${CBQ_READER_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/cbq_reader_harness}"
ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN:-$ROOT_DIR/core/legacy/source/cbq_ordered_encoder}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_ordered_encoder_smoke}"

if [[ ! -x "$FASTX_BIN" || ! -x "$CBQ_BIN" || ! -x "$ENCODER_BIN" ]]; then
    make -C "$ROOT_DIR/core/legacy/source" fastx-input-harness cbq-reader-harness cbq-ordered-encoder
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{inputs,dumps,logs}

cat > "$OUT_ROOT/inputs/lane1_R1.fastq" <<'FASTQ'
@readA/1 1:N:0:ACGT
ACGTACGTNN
+
IIIIIIIIII
@readB/1 1:Y:0:ACGT
TTGGAACCAA
+
HHHHHHHHHH
@readC/1 1:N:0:ACGT
NNNNACGTACGTACGTACGTACGTACGTNN
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@readD/1 1:N:0:ACGT
ACACACACACACACAC
+
FFFFFFFFFFFFFFFF
FASTQ

cat > "$OUT_ROOT/inputs/lane1_R2.fastq" <<'FASTQ'
@readA/2 2:N:0:ACGT
TGCATGCANN
+
IIIIIIIIII
@readB/2 2:Y:0:ACGT
AACCTTGGTT
+
HHHHHHHHHH
@readC/2 2:N:0:ACGT
TGCATGCATGCATGCATGCATGCATGCANNNN
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@readD/2 2:N:0:ACGT
TGTGTGTGTGTGTGTG
+
FFFFFFFFFFFFFFFF
FASTQ

cat > "$OUT_ROOT/inputs/single.fastq" <<'FASTQ'
@singleA 1:N:0:TTTT
ACGTNNACGT
+
IIIIIIIIII
@singleB 1:Y:0:TTTT
TTTTCCCCAAAA
+
HHHHHHHHHHHH
FASTQ

gzip -c "$OUT_ROOT/inputs/lane1_R1.fastq" > "$OUT_ROOT/inputs/lane1_R1.fastq.gz"
gzip -c "$OUT_ROOT/inputs/lane1_R2.fastq" > "$OUT_ROOT/inputs/lane1_R2.fastq.gz"

PAIR_CBQ="$OUT_ROOT/inputs/lane1_pair.cbq"
PAIR_GZ_CBQ="$OUT_ROOT/inputs/lane1_pair_gz.cbq"
SINGLE_CBQ="$OUT_ROOT/inputs/single.cbq"

"$ENCODER_BIN" \
    --readFilesIn "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" \
    --outFile "$PAIR_CBQ" \
    --blockSize 90 \
    > "$OUT_ROOT/logs/encode_pair.stdout" 2> "$OUT_ROOT/logs/encode_pair.stderr"

"$FASTX_BIN" --readNameSeparator / \
    --readFilesIn "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" \
    > "$OUT_ROOT/dumps/source_fastq_pair.tsv"
"$CBQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --readFilesIn "$PAIR_CBQ" \
    > "$OUT_ROOT/dumps/cbq_ordered_pair.tsv"
cmp -s "$OUT_ROOT/dumps/source_fastq_pair.tsv" "$OUT_ROOT/dumps/cbq_ordered_pair.tsv"

"$ENCODER_BIN" \
    --readFilesIn "$OUT_ROOT/inputs/lane1_R1.fastq.gz" "$OUT_ROOT/inputs/lane1_R2.fastq.gz" \
    --outFile "$PAIR_GZ_CBQ" \
    --blockSize 90 \
    > "$OUT_ROOT/logs/encode_pair_gz.stdout" 2> "$OUT_ROOT/logs/encode_pair_gz.stderr"

"$FASTX_BIN" --readNameSeparator / \
    --readFilesIn "$OUT_ROOT/inputs/lane1_R1.fastq.gz" "$OUT_ROOT/inputs/lane1_R2.fastq.gz" \
    > "$OUT_ROOT/dumps/source_fastq_pair_gz.tsv"
"$CBQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --readFilesIn "$PAIR_GZ_CBQ" \
    > "$OUT_ROOT/dumps/cbq_ordered_pair_gz.tsv"
cmp -s "$OUT_ROOT/dumps/source_fastq_pair_gz.tsv" "$OUT_ROOT/dumps/cbq_ordered_pair_gz.tsv"

"$ENCODER_BIN" \
    --readFilesIn "$OUT_ROOT/inputs/single.fastq" \
    --outFile "$SINGLE_CBQ" \
    --blockSize 64 \
    > "$OUT_ROOT/logs/encode_single.stdout" 2> "$OUT_ROOT/logs/encode_single.stderr"

"$FASTX_BIN" --readNameSeparator / \
    --readFilesIn "$OUT_ROOT/inputs/single.fastq" \
    > "$OUT_ROOT/dumps/source_fastq_single.tsv"
"$CBQ_BIN" --readNameSeparator / \
    --mateCount 1 \
    --readFilesIn "$SINGLE_CBQ" \
    > "$OUT_ROOT/dumps/cbq_ordered_single.tsv"
cmp -s "$OUT_ROOT/dumps/source_fastq_single.tsv" "$OUT_ROOT/dumps/cbq_ordered_single.tsv"

if [[ -n "${BQTOOLS:-}" ]]; then
    if [[ -x "$BQTOOLS" ]]; then
        BQTOOLS_BIN="$BQTOOLS"
    elif command -v "$BQTOOLS" >/dev/null 2>&1; then
        BQTOOLS_BIN="$(command -v "$BQTOOLS")"
    else
        BQTOOLS_BIN=""
    fi
elif command -v bqtools >/dev/null 2>&1; then
    BQTOOLS_BIN="$(command -v bqtools)"
elif [[ -x /tmp/star_suite_bqtools/bin/bqtools ]]; then
    BQTOOLS_BIN="/tmp/star_suite_bqtools/bin/bqtools"
else
    BQTOOLS_BIN=""
fi

if [[ -n "$BQTOOLS_BIN" ]]; then
    "$BQTOOLS_BIN" info "$PAIR_CBQ" > "$OUT_ROOT/dumps/bqtools_pair_info.txt"
    mkdir -p "$OUT_ROOT/bqtools_decode"
    "$BQTOOLS_BIN" decode "$PAIR_CBQ" \
        --prefix "$OUT_ROOT/bqtools_decode/pair" \
        -f q \
        -T 1 \
        > "$OUT_ROOT/logs/bqtools_decode.stdout" 2> "$OUT_ROOT/logs/bqtools_decode.stderr"
    cmp -s "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/bqtools_decode/pair_R1.fq"
    cmp -s "$OUT_ROOT/inputs/lane1_R2.fastq" "$OUT_ROOT/bqtools_decode/pair_R2.fq"
fi

echo "PASS: ordered CBQ encoder smoke output at $OUT_ROOT"
