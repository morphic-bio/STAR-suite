#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${STAR_BIN:-$ROOT_DIR/core/legacy/source/STAR}"
CBQ_BIN="${CBQ_READER_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/cbq_reader_harness}"
ENCODER_BIN="${CBQ_ORDERED_ENCODER_BIN:-$ROOT_DIR/core/legacy/source/cbq_ordered_encoder}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_ynoy_smoke}"

if [[ ! -x "$STAR_BIN" || ! -x "$CBQ_BIN" || ! -x "$ENCODER_BIN" ]]; then
    make -C "$ROOT_DIR/core/legacy/source" -j"${THREADS:-4}" STAR cbq-reader-harness cbq-ordered-encoder
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{inputs,genome,dumps,logs}

python3 - "$OUT_ROOT" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(76029)
alphabet = "ACGT"
chr1 = "".join(rng.choice(alphabet) for _ in range(5000))
chry = "".join(rng.choice(alphabet) for _ in range(5000))

def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

(out / "inputs" / "genome.fa").write_text(
    ">chr1\n" + chr1 + "\n>chrY\n" + chry + "\n",
    encoding="ascii",
)

expected_y = []
expected_noy = []
with (out / "inputs" / "mixed_R1.fastq").open("wt", encoding="ascii") as r1, \
     (out / "inputs" / "mixed_R2.fastq").open("wt", encoding="ascii") as r2:
    for i in range(20):
        is_y = i % 2 == 1
        chrom = chry if is_y else chr1
        name = f"{'Y' if is_y else 'noY'}_{i:02d}"
        start1 = 100 + i * 120
        start2 = start1 + 70
        seq1 = chrom[start1:start1 + 50]
        seq2 = rc(chrom[start2:start2 + 50])
        qchar = chr(ord("F") + (i % 20))
        r1.write(f"@{name}/1 1:N:0:ACGT\n{seq1}\n+\n{qchar * len(seq1)}\n")
        r2.write(f"@{name}/2 2:N:0:ACGT\n{seq2}\n+\n{qchar * len(seq2)}\n")
        target = expected_y if is_y else expected_noy
        target.extend([name, name])

(out / "expected_y_names.txt").write_text("\n".join(expected_y) + "\n", encoding="ascii")
(out / "expected_noy_names.txt").write_text("\n".join(expected_noy) + "\n", encoding="ascii")
PY

CBQ="$OUT_ROOT/inputs/mixed.cbq"
"$ENCODER_BIN" \
    --readFilesIn "$OUT_ROOT/inputs/mixed_R1.fastq" "$OUT_ROOT/inputs/mixed_R2.fastq" \
    --outFile "$CBQ" \
    --compressionLevel 0 \
    --blockSize 512 \
    > "$OUT_ROOT/logs/encode.stdout" 2> "$OUT_ROOT/logs/encode.stderr"

"$STAR_BIN" \
    --runMode genomeGenerate \
    --genomeDir "$OUT_ROOT/genome" \
    --genomeFastaFiles "$OUT_ROOT/inputs/genome.fa" \
    --runThreadN 1 \
    --genomeSAindexNbases 3 \
    --genomeChrBinNbits 4 \
    > "$OUT_ROOT/logs/genome_generate.stdout" \
    2> "$OUT_ROOT/logs/genome_generate.stderr"

run_case() {
    local case_name="$1"
    shift
    local star_dir="$OUT_ROOT/star_${case_name}"
    local dump_dir="$OUT_ROOT/dumps/${case_name}"
    mkdir -p "$star_dir" "$dump_dir"

    "$STAR_BIN" \
        --runMode alignReads \
        --genomeDir "$OUT_ROOT/genome" \
        --runThreadN 2 \
        --readNameSeparator / \
        "$@" \
        --emitYNoY yes \
        --emitYNoYFormat cbq \
        --emitYNoYCbqCompressionLevel 0 \
        --emitYNoYCbqBlockSize 512 \
        --outSAMtype None \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 0 \
        --alignIntronMax 1 \
        --alignMatesGapMax 1000 \
        --alignEndsType EndToEnd \
        --limitIObufferSize 205000 500000 \
        --outFileNamePrefix "$star_dir/" \
        > "$OUT_ROOT/logs/star_${case_name}.stdout" \
        2> "$OUT_ROOT/logs/star_${case_name}.stderr"

    local y_cbq="$star_dir/y_separated/mixed_Y.cbq"
    local noy_cbq="$star_dir/y_separated/mixed_noY.cbq"
    [[ -s "$y_cbq" ]] || { echo "ERROR: missing Y CBQ output: $y_cbq" >&2; exit 1; }
    [[ -s "$noy_cbq" ]] || { echo "ERROR: missing noY CBQ output: $noy_cbq" >&2; exit 1; }

    "$CBQ_BIN" --readNameSeparator / --mateCount 2 --readFilesIn "$y_cbq" --dump-fastq \
        > "$dump_dir/y.fastq"
    "$CBQ_BIN" --readNameSeparator / --mateCount 2 --readFilesIn "$noy_cbq" --dump-fastq \
        > "$dump_dir/noy.fastq"

    awk 'NR % 4 == 1 { name=$1; sub(/^@/, "", name); print name }' "$dump_dir/y.fastq" \
        > "$dump_dir/y.names"
    awk 'NR % 4 == 1 { name=$1; sub(/^@/, "", name); print name }' "$dump_dir/noy.fastq" \
        > "$dump_dir/noy.names"

    cmp -s "$OUT_ROOT/expected_y_names.txt" "$dump_dir/y.names"
    cmp -s "$OUT_ROOT/expected_noy_names.txt" "$dump_dir/noy.names"

    if grep -q '^noY_' "$dump_dir/y.names"; then
        echo "ERROR: noY read present in Y CBQ output for $case_name" >&2
        exit 1
    fi
    if grep -q '^Y_' "$dump_dir/noy.names"; then
        echo "ERROR: Y read present in noY CBQ output for $case_name" >&2
        exit 1
    fi
}

run_case fastq_input \
    --readFilesType Fastx \
    --readFilesIn "$OUT_ROOT/inputs/mixed_R1.fastq" "$OUT_ROOT/inputs/mixed_R2.fastq"

run_case cbq_input \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ"

echo "PASS: ordered CBQ Y/noY smoke output at $OUT_ROOT"
