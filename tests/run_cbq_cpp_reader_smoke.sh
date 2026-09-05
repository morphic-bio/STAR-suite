#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTX_BIN="${FASTX_INPUT_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/fastx_input_harness}"
CBQ_BIN="${CBQ_READER_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/cbq_reader_harness}"
CBQ_STAR_BIN="${CBQ_STAR_ADAPTER_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/cbq_star_adapter_harness}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_cpp_reader_smoke}"

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

if [[ ! -x "$FASTX_BIN" || ! -x "$CBQ_BIN" || ! -x "$CBQ_STAR_BIN" ]]; then
    make -C "$ROOT_DIR/core/legacy/source" fastx-input-harness cbq-reader-harness cbq-star-adapter-harness
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{inputs,dumps}

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
FASTQ

cp "$ROOT_DIR/tests/fixtures/trim_qc_fastq_tiny.fastq" "$OUT_ROOT/inputs/single.fastq"

CBQ="$OUT_ROOT/inputs/lane1_pair.cbq"
CBQ_LEVEL0="$OUT_ROOT/inputs/lane1_pair_level0.cbq"
SE_CBQ="$OUT_ROOT/inputs/single.cbq"
SE_CBQ_LEVEL0="$OUT_ROOT/inputs/single_level0.cbq"

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

if ! encode_pair_cbq "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" "$CBQ_LEVEL0" 0; then
    echo "ERROR: failed to encode paired FASTQ to level-0 CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

if ! encode_single_cbq "$OUT_ROOT/inputs/single.fastq" "$SE_CBQ"; then
    echo "ERROR: failed to encode single-end FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_single.stderr" >&2
    exit 1
fi

if ! encode_single_cbq "$OUT_ROOT/inputs/single.fastq" "$SE_CBQ_LEVEL0" 0; then
    echo "ERROR: failed to encode single-end FASTQ to level-0 CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_single.stderr" >&2
    exit 1
fi

"$FASTX_BIN" --readNameSeparator / \
    --readFilesIn "$OUT_ROOT/inputs/lane1_R1.fastq" "$OUT_ROOT/inputs/lane1_R2.fastq" \
    > "$OUT_ROOT/dumps/source_fastq_pair.tsv"

"$CBQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --verify-packed-windows \
    --readFilesIn "$CBQ" \
    > "$OUT_ROOT/dumps/cbq_reader_pair.tsv"
compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_pair.tsv" "$OUT_ROOT/dumps/cbq_reader_pair.tsv"

"$CBQ_STAR_BIN" --readNameSeparator / \
    --mateCount 2 \
    --readFilesIn "$CBQ" \
    --mode direct \
    > "$OUT_ROOT/dumps/cbq_star_pair.direct.bin"
"$CBQ_STAR_BIN" --readNameSeparator / \
    --mateCount 2 \
    --readFilesIn "$CBQ" \
    --mode reference \
    > "$OUT_ROOT/dumps/cbq_star_pair.reference.bin"
cmp -s "$OUT_ROOT/dumps/cbq_star_pair.direct.bin" "$OUT_ROOT/dumps/cbq_star_pair.reference.bin"

"$CBQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --verify-packed-windows \
    --readFilesIn "$CBQ_LEVEL0" \
    > "$OUT_ROOT/dumps/cbq_reader_pair_level0.tsv"
compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_pair.tsv" "$OUT_ROOT/dumps/cbq_reader_pair_level0.tsv"

"$CBQ_STAR_BIN" --readNameSeparator / \
    --mateCount 2 \
    --readFilesIn "$CBQ_LEVEL0" \
    --mode direct \
    > "$OUT_ROOT/dumps/cbq_star_pair_level0.direct.bin"
"$CBQ_STAR_BIN" --readNameSeparator / \
    --mateCount 2 \
    --readFilesIn "$CBQ_LEVEL0" \
    --mode reference \
    > "$OUT_ROOT/dumps/cbq_star_pair_level0.reference.bin"
cmp -s "$OUT_ROOT/dumps/cbq_star_pair_level0.direct.bin" "$OUT_ROOT/dumps/cbq_star_pair_level0.reference.bin"

printf '%s\t-\tID:lane1\n' "$CBQ" > "$OUT_ROOT/manifest.tsv"
"$CBQ_BIN" --readNameSeparator / \
    --mateCount 2 \
    --readFilesManifest "$OUT_ROOT/manifest.tsv" \
    > "$OUT_ROOT/dumps/cbq_reader_manifest.tsv"
compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_pair.tsv" "$OUT_ROOT/dumps/cbq_reader_manifest.tsv"

"$FASTX_BIN" --readNameSeparator / \
    --readFilesIn "$OUT_ROOT/inputs/single.fastq" \
    > "$OUT_ROOT/dumps/source_fastq_single.tsv"

"$CBQ_BIN" --readNameSeparator / \
    --mateCount 1 \
    --verify-packed-windows \
    --readFilesIn "$SE_CBQ" \
    > "$OUT_ROOT/dumps/cbq_reader_single.tsv"
compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_single.tsv" "$OUT_ROOT/dumps/cbq_reader_single.tsv"

"$CBQ_STAR_BIN" --readNameSeparator / \
    --mateCount 1 \
    --readFilesIn "$SE_CBQ" \
    --mode direct \
    > "$OUT_ROOT/dumps/cbq_star_single.direct.bin"
"$CBQ_STAR_BIN" --readNameSeparator / \
    --mateCount 1 \
    --readFilesIn "$SE_CBQ" \
    --mode reference \
    > "$OUT_ROOT/dumps/cbq_star_single.reference.bin"
cmp -s "$OUT_ROOT/dumps/cbq_star_single.direct.bin" "$OUT_ROOT/dumps/cbq_star_single.reference.bin"

"$CBQ_BIN" --readNameSeparator / \
    --mateCount 1 \
    --verify-packed-windows \
    --readFilesIn "$SE_CBQ_LEVEL0" \
    > "$OUT_ROOT/dumps/cbq_reader_single_level0.tsv"
compare_contract_tsv_unordered "$OUT_ROOT/dumps/source_fastq_single.tsv" "$OUT_ROOT/dumps/cbq_reader_single_level0.tsv"

"$CBQ_STAR_BIN" --readNameSeparator / \
    --mateCount 1 \
    --readFilesIn "$SE_CBQ_LEVEL0" \
    --mode direct \
    > "$OUT_ROOT/dumps/cbq_star_single_level0.direct.bin"
"$CBQ_STAR_BIN" --readNameSeparator / \
    --mateCount 1 \
    --readFilesIn "$SE_CBQ_LEVEL0" \
    --mode reference \
    > "$OUT_ROOT/dumps/cbq_star_single_level0.reference.bin"
cmp -s "$OUT_ROOT/dumps/cbq_star_single_level0.direct.bin" "$OUT_ROOT/dumps/cbq_star_single_level0.reference.bin"

echo "PASS: native CBQ C++ reader smoke completed at $OUT_ROOT"
