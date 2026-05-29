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

cp "$ROOT_DIR/tests/fixtures/trim_qc_fastq_tiny.fastq" "$OUT_ROOT/inputs/single.fastq"

CBQ="$OUT_ROOT/inputs/lane1_pair.cbq"
CBQ_UNCOMPRESSED="$OUT_ROOT/inputs/lane1_pair_uncompressed.cbq"
SE_CBQ="$OUT_ROOT/inputs/single.cbq"
SE_CBQ_UNCOMPRESSED="$OUT_ROOT/inputs/single_uncompressed.cbq"

encode_pair_cbq() {
    local r1="$1"
    local r2="$2"
    local out="$3"
    local level="${4:-}"
    local level_args=()
    [[ -n "$level" ]] && level_args=(-l "$level")

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" --mode cbq "${level_args[@]}" -o "$out" -T 2 \
        > "$OUT_ROOT/encode_pair.stdout" 2> "$OUT_ROOT/encode_pair.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" --mode cbq "${level_args[@]}" -T 2 \
        >> "$OUT_ROOT/encode_pair.stdout" 2>> "$OUT_ROOT/encode_pair.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" "${level_args[@]}" -T 2 \
        >> "$OUT_ROOT/encode_pair.stdout" 2>> "$OUT_ROOT/encode_pair.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    "$BQTOOLS_BIN" encode "$r1" "$r2" -o "$out" "${level_args[@]}" \
        >> "$OUT_ROOT/encode_pair.stdout" 2>> "$OUT_ROOT/encode_pair.stderr"
    [[ -s "$out" ]]
}

encode_single_cbq() {
    local r1="$1"
    local out="$2"
    local level="${3:-}"
    local level_args=()
    [[ -n "$level" ]] && level_args=(-l "$level")

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" --mode cbq "${level_args[@]}" -o "$out" -T 2 \
        > "$OUT_ROOT/encode_single.stdout" 2> "$OUT_ROOT/encode_single.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" -o "$out" --mode cbq "${level_args[@]}" -T 2 \
        >> "$OUT_ROOT/encode_single.stdout" 2>> "$OUT_ROOT/encode_single.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" -o "$out" "${level_args[@]}" -T 2 \
        >> "$OUT_ROOT/encode_single.stdout" 2>> "$OUT_ROOT/encode_single.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    "$BQTOOLS_BIN" encode "$r1" -o "$out" "${level_args[@]}" \
        >> "$OUT_ROOT/encode_single.stdout" 2>> "$OUT_ROOT/encode_single.stderr"
    [[ -s "$out" ]]
}

compare_contract_tsv_unordered() {
    python3 - "$1" "$2" <<'PY'
from collections import Counter
import sys

expected_path, observed_path = sys.argv[1:3]

def load(path):
    rows = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 8:
                raise SystemExit(f"{path}:{line_number}: expected 8 TSV fields, observed {len(fields)}")
            # Field 1 is read_ordinal. It is emitted-stream order, not a
            # source-order guarantee, so parity ignores it.
            rows.append(tuple([fields[0]] + fields[2:]))
    return Counter(rows)

expected = load(expected_path)
observed = load(observed_path)
if expected != observed:
    expected_only = list((expected - observed).elements())[:5]
    observed_only = list((observed - expected).elements())[:5]
    raise SystemExit(
        "order-independent input contract parity failed\n"
        f"expected_only={expected_only}\n"
        f"observed_only={observed_only}"
    )
PY
}

if ! encode_pair_cbq "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" "$CBQ"; then
    echo "ERROR: failed to encode paired FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

if ! encode_pair_cbq "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" "$CBQ_UNCOMPRESSED" 0; then
    echo "ERROR: failed to encode paired FASTQ to uncompressed CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

if ! encode_single_cbq "$OUT_ROOT/inputs/single.fastq" "$SE_CBQ"; then
    echo "ERROR: failed to encode single-end FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_single.stderr" >&2
    exit 1
fi

if ! encode_single_cbq "$OUT_ROOT/inputs/single.fastq" "$SE_CBQ_UNCOMPRESSED" 0; then
    echo "ERROR: failed to encode single-end FASTQ to uncompressed CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_single.stderr" >&2
    exit 1
fi

"$BQTOOLS_BIN" info "$CBQ" > "$OUT_ROOT/pair_info.txt" 2> "$OUT_ROOT/pair_info.stderr"
"$BQTOOLS_BIN" info "$CBQ_UNCOMPRESSED" > "$OUT_ROOT/pair_uncompressed_info.txt" 2> "$OUT_ROOT/pair_uncompressed_info.stderr"
"$BQTOOLS_BIN" info "$SE_CBQ" > "$OUT_ROOT/single_info.txt" 2> "$OUT_ROOT/single_info.stderr"
"$BQTOOLS_BIN" info "$SE_CBQ_UNCOMPRESSED" > "$OUT_ROOT/single_uncompressed_info.txt" 2> "$OUT_ROOT/single_uncompressed_info.stderr"

"$FASTX_BIN" --readNameSeparator / \
    --readFilesIn "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" \
    > "$OUT_ROOT/dumps/source_fastq_pair.tsv"

"$BINSEQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --bqtools "$BQTOOLS_BIN" \
    --workDir "$OUT_ROOT/probe/direct" \
    --readFilesIn "$CBQ" \
    > "$OUT_ROOT/dumps/binseq_probe_pair.tsv"

compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_pair.tsv" "$OUT_ROOT/dumps/binseq_probe_pair.tsv"

"$BINSEQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --bqtools "$BQTOOLS_BIN" \
    --workDir "$OUT_ROOT/probe/direct_uncompressed" \
    --readFilesIn "$CBQ_UNCOMPRESSED" \
    > "$OUT_ROOT/dumps/binseq_probe_pair_uncompressed.tsv"

compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_pair.tsv" "$OUT_ROOT/dumps/binseq_probe_pair_uncompressed.tsv"

printf '%s\t-\tID:lane1\n' "$CBQ" > "$OUT_ROOT/manifest.tsv"
"$BINSEQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --bqtools "$BQTOOLS_BIN" \
    --workDir "$OUT_ROOT/probe/manifest" \
    --readFilesManifest "$OUT_ROOT/manifest.tsv" \
    > "$OUT_ROOT/dumps/binseq_probe_manifest.tsv"

compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_pair.tsv" "$OUT_ROOT/dumps/binseq_probe_manifest.tsv"

"$FASTX_BIN" --readNameSeparator / \
    --readFilesIn "$OUT_ROOT/inputs/single.fastq" \
    > "$OUT_ROOT/dumps/source_fastq_single.tsv"

"$BINSEQ_BIN" --readNameSeparator / \
    --mateCount 1 \
    --bqtools "$BQTOOLS_BIN" \
    --workDir "$OUT_ROOT/probe/single" \
    --readFilesIn "$SE_CBQ" \
    > "$OUT_ROOT/dumps/binseq_probe_single.tsv"

compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_single.tsv" "$OUT_ROOT/dumps/binseq_probe_single.tsv"

"$BINSEQ_BIN" --readNameSeparator / \
    --mateCount 1 \
    --bqtools "$BQTOOLS_BIN" \
    --workDir "$OUT_ROOT/probe/single_uncompressed" \
    --readFilesIn "$SE_CBQ_UNCOMPRESSED" \
    > "$OUT_ROOT/dumps/binseq_probe_single_uncompressed.tsv"

compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_single.tsv" "$OUT_ROOT/dumps/binseq_probe_single_uncompressed.tsv"

echo "PASS: BINSEQ probe smoke completed at $OUT_ROOT"
