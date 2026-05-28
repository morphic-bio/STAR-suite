#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ADAPTER_BIN="${CBQ_CHROMAP_ADAPTER_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/cbq_chromap_adapter_harness}"
CHROMAP_BIN="${CHROMAP_BIN:-/mnt/pikachu/Chromap-suite/chromap}"
RUN_CHROMAP_MAPPING_SMOKE="${RUN_CHROMAP_MAPPING_SMOKE:-auto}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_chromap_adapter_smoke}"

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

if [[ ! -x "$ADAPTER_BIN" ]]; then
    make -C "$ROOT_DIR/core/legacy/source" cbq-chromap-adapter-harness
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{inputs,adapted,chromap}

python3 - "$OUT_ROOT" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(74219)
alphabet = "ACGT"
genome = "".join(rng.choice(alphabet) for _ in range(3000))

def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

barcodes = [
    "ACGTACGTACGTACGT",
    "TGCATGCATGCATGCA",
    "GGGGAAAACCCCTTTT",
]
fragments = [
    ("read001", 200, 320, barcodes[0], "I"),
    ("read002", 600, 735, barcodes[1], "H"),
    ("read003", 1200, 1360, barcodes[2], "J"),
]

(out / "inputs" / "genome.fa").write_text(">chrSynthetic\n" + genome + "\n", encoding="ascii")
(out / "inputs" / "whitelist.txt").write_text("\n".join(barcodes) + "\n", encoding="ascii")

with (out / "inputs" / "reads_R1.fastq").open("wt", encoding="ascii") as r1, \
     (out / "inputs" / "reads_R2.fastq").open("wt", encoding="ascii") as r2, \
     (out / "inputs" / "barcodes.fastq").open("wt", encoding="ascii") as bc:
    for name, start, end, barcode, qchar in fragments:
        seq1 = genome[start:start + 60]
        seq2 = rc(genome[end - 60:end])
        r1.write(f"@{name}/1 1:N:0:ACGT\n{seq1}\n+\n{qchar * len(seq1)}\n")
        r2.write(f"@{name}/2 2:N:0:ACGT\n{seq2}\n+\n{qchar * len(seq2)}\n")
        bc.write(f"@{name}/3 3:N:0:ACGT\n{barcode}\n+\n{qchar * len(barcode)}\n")
PY

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

encode_single_cbq() {
    local r1="$1"
    local out="$2"

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" --mode cbq -o "$out" -T 2 \
        > "$OUT_ROOT/encode_single.stdout" 2> "$OUT_ROOT/encode_single.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" -o "$out" --mode cbq -T 2 \
        >> "$OUT_ROOT/encode_single.stdout" 2>> "$OUT_ROOT/encode_single.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    if "$BQTOOLS_BIN" encode "$r1" -o "$out" -T 2 \
        >> "$OUT_ROOT/encode_single.stdout" 2>> "$OUT_ROOT/encode_single.stderr"; then
        [[ -s "$out" ]] && return 0
    fi

    rm -f "$out"
    "$BQTOOLS_BIN" encode "$r1" -o "$out" \
        >> "$OUT_ROOT/encode_single.stdout" 2>> "$OUT_ROOT/encode_single.stderr"
    [[ -s "$out" ]]
}

READ_PAIR_CBQ="$OUT_ROOT/inputs/reads_pair.cbq"
BARCODE_CBQ="$OUT_ROOT/inputs/barcodes.cbq"

if ! encode_pair_cbq "$OUT_ROOT/inputs/reads_R1.fastq" "$OUT_ROOT/inputs/reads_R2.fastq" "$READ_PAIR_CBQ"; then
    echo "ERROR: failed to encode read pair FASTQs to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

if ! encode_single_cbq "$OUT_ROOT/inputs/barcodes.fastq" "$BARCODE_CBQ"; then
    echo "ERROR: failed to encode barcode FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_single.stderr" >&2
    exit 1
fi

"$ADAPTER_BIN" \
    --readPairCbq "$READ_PAIR_CBQ" \
    --barcodeCbq "$BARCODE_CBQ" \
    --outputDir "$OUT_ROOT/adapted" \
    --sampleName smoke \
    > "$OUT_ROOT/adapted/stdout.tsv" \
    2> "$OUT_ROOT/adapted/stderr.log"

compare_fastq_payload() {
    python3 - "$1" "$2" <<'PY'
import sys

source_path, observed_path = sys.argv[1:3]

def payload(path):
    rows = []
    with open(path, "rt", encoding="ascii") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip("\n")
            plus = handle.readline()
            qual = handle.readline().rstrip("\n")
            if not qual or not header.startswith("@") or not plus.startswith("+"):
                raise SystemExit(f"invalid FASTQ record in {path}")
            rows.append((seq, qual))
    return rows

source = payload(source_path)
observed = payload(observed_path)
if source != observed:
    raise SystemExit(f"FASTQ payload mismatch: {source_path} vs {observed_path}")
PY
}

compare_fastq_payload "$OUT_ROOT/inputs/reads_R1.fastq" "$OUT_ROOT/adapted/smoke.chromap.R1.fastq"
compare_fastq_payload "$OUT_ROOT/inputs/reads_R2.fastq" "$OUT_ROOT/adapted/smoke.chromap.R2.fastq"
compare_fastq_payload "$OUT_ROOT/inputs/barcodes.fastq" "$OUT_ROOT/adapted/smoke.chromap.barcode.fastq"

if [[ "$RUN_CHROMAP_MAPPING_SMOKE" == "1" || ( "$RUN_CHROMAP_MAPPING_SMOKE" == "auto" && -x "$CHROMAP_BIN" ) ]]; then
    "$CHROMAP_BIN" \
        -i \
        -r "$OUT_ROOT/inputs/genome.fa" \
        -o "$OUT_ROOT/chromap/genome.index" \
        -k 11 \
        -w 5 \
        > "$OUT_ROOT/chromap/index.stdout" \
        2> "$OUT_ROOT/chromap/index.stderr"

    "$CHROMAP_BIN" \
        --preset atac \
        -r "$OUT_ROOT/inputs/genome.fa" \
        -x "$OUT_ROOT/chromap/genome.index" \
        -1 "$OUT_ROOT/adapted/smoke.chromap.R1.fastq" \
        -2 "$OUT_ROOT/adapted/smoke.chromap.R2.fastq" \
        -b "$OUT_ROOT/adapted/smoke.chromap.barcode.fastq" \
        --barcode-whitelist "$OUT_ROOT/inputs/whitelist.txt" \
        --read-format bc:0:-1 \
        -o "$OUT_ROOT/chromap/fragments.bed" \
        --BED \
        -t 1 \
        --min-read-length 20 \
        -q 0 \
        > "$OUT_ROOT/chromap/map.stdout" \
        2> "$OUT_ROOT/chromap/map.stderr"

    [[ -s "$OUT_ROOT/chromap/fragments.bed" ]] || {
        echo "ERROR: Chromap mapping smoke produced no fragments" >&2
        exit 1
    }
fi

echo "PASS: Chromap CBQ adapter smoke completed at $OUT_ROOT"
