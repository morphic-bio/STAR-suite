#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTX_BIN="${FASTX_INPUT_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/fastx_input_harness}"
BINSEQ_BIN="${BINSEQ_PROBE_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/binseq_probe_harness}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_binseq_probe_smoke}"

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
else
    command -v bqtools >/dev/null 2>&1 || skip "bqtools not found in PATH; set BQTOOLS=/path/to/bqtools"
    BQTOOLS_BIN="$(command -v bqtools)"
fi

if [[ ! -x "$FASTX_BIN" || ! -x "$BINSEQ_BIN" ]]; then
    make -C "$ROOT_DIR/core/legacy/source" fastx-input-harness binseq-probe-harness
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{inputs,dumps,decoded,probe}

cat > "$OUT_ROOT/inputs/lane1_R1.fastq" <<'FASTQ'
@readA/1 1:N:0:ACGT
ACGTACGTNN
+
IIIIIIIIII
@readB/1 1:Y:0:ACGT
TTGGAACCAA
+
HHHHHHHHHH
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
FASTQ

CBQ="$OUT_ROOT/inputs/lane1_pair.cbq"

encode_pair_cbq() {
    local r1="$1"
    local r2="$2"
    local out="$3"

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" --mode cbq -o "$out" -T 2 \
        > "$OUT_ROOT/encode.stdout" 2> "$OUT_ROOT/encode.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" --mode cbq -T 2 \
        >> "$OUT_ROOT/encode.stdout" 2>> "$OUT_ROOT/encode.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" -T 2 \
        >> "$OUT_ROOT/encode.stdout" 2>> "$OUT_ROOT/encode.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" \
        >> "$OUT_ROOT/encode.stdout" 2>> "$OUT_ROOT/encode.stderr"
    [[ -s "$out" ]]
}

if ! encode_pair_cbq "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" "$CBQ"; then
    echo "ERROR: failed to encode paired FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode.stderr" >&2
    exit 1
fi

"$BQTOOLS_BIN" info "$CBQ" > "$OUT_ROOT/info.txt" 2> "$OUT_ROOT/info.stderr"

"$FASTX_BIN" --readNameSeparator / \
    --readFilesIn "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" \
    > "$OUT_ROOT/dumps/source_fastq.tsv"

"$BINSEQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --bqtools "$BQTOOLS_BIN" \
    --workDir "$OUT_ROOT/probe/direct" \
    --readFilesIn "$CBQ" \
    > "$OUT_ROOT/dumps/binseq_probe.tsv"

diff -u "$OUT_ROOT/dumps/source_fastq.tsv" "$OUT_ROOT/dumps/binseq_probe.tsv"

printf '%s\t-\tID:lane1\n' "$CBQ" > "$OUT_ROOT/manifest.tsv"
"$BINSEQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --bqtools "$BQTOOLS_BIN" \
    --workDir "$OUT_ROOT/probe/manifest" \
    --readFilesManifest "$OUT_ROOT/manifest.tsv" \
    > "$OUT_ROOT/dumps/binseq_probe_manifest.tsv"

diff -u "$OUT_ROOT/dumps/source_fastq.tsv" "$OUT_ROOT/dumps/binseq_probe_manifest.tsv"

echo "PASS: BINSEQ probe smoke completed at $OUT_ROOT"
