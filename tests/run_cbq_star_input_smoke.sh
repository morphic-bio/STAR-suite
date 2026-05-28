#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${STAR_BIN:-$ROOT_DIR/core/legacy/source/STAR}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_star_input_smoke}"

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

if [[ ! -x "$STAR_BIN" ]]; then
    make -C "$ROOT_DIR/core/legacy/source" STAR
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{inputs,genome,fastx,fastx_b,fastx_multi,binseq,binseq_manifest,binseq_multi,binseq_multi_manifest,fastx_se,binseq_se,batch_binseq,batch_binseq_manifest}

python3 - "$OUT_ROOT" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(91231)
alphabet = "ACGT"
genome = "".join(rng.choice(alphabet) for _ in range(1200))

def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

reads_a = [
    ("read001", genome[101:151], rc(genome[211:261]), "I"),
    ("read002", genome[403:453], rc(genome[552:602]), "H"),
    ("read003", genome[735:785], rc(genome[894:944]), "J"),
]
reads_b = [
    ("read101", genome[167:217], rc(genome[318:368]), "I"),
    ("read102", genome[626:676], rc(genome[809:859]), "H"),
    ("read103", genome[944:994], rc(genome[1043:1093]), "J"),
]

(out / "inputs" / "genome.fa").write_text(">chrSynthetic\n" + genome + "\n", encoding="ascii")

def write_pe(prefix, reads):
    with (out / "inputs" / f"{prefix}_R1_001.fastq").open("wt", encoding="ascii") as r1, \
         (out / "inputs" / f"{prefix}_R2_001.fastq").open("wt", encoding="ascii") as r2:
        for name, seq1, seq2, qchar in reads:
            r1.write(f"@{name} 1:N:0:ACGT\n{seq1}\n+\n{qchar * len(seq1)}\n")
            r2.write(f"@{name} 2:N:0:ACGT\n{seq2}\n+\n{qchar * len(seq2)}\n")

write_pe("sampleA_S1_L001", reads_a)
write_pe("sampleB_S2_L001", reads_b)

with (out / "inputs" / "reads_SE.fastq").open("wt", encoding="ascii") as se:
    for name, seq1, _seq2, qchar in reads_a:
        se.write(f"@{name} 1:N:0:ACGT\n{seq1}\n+\n{qchar * len(seq1)}\n")
PY

FASTQ_A_R1="$OUT_ROOT/inputs/sampleA_S1_L001_R1_001.fastq"
FASTQ_A_R2="$OUT_ROOT/inputs/sampleA_S1_L001_R2_001.fastq"
FASTQ_B_R1="$OUT_ROOT/inputs/sampleB_S2_L001_R1_001.fastq"
FASTQ_B_R2="$OUT_ROOT/inputs/sampleB_S2_L001_R2_001.fastq"
CBQ="$OUT_ROOT/inputs/sampleA_S1_L001.cbq"
CBQ_B="$OUT_ROOT/inputs/sampleB_S2_L001.cbq"
SE_CBQ="$OUT_ROOT/inputs/reads_single.cbq"

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

if ! encode_pair_cbq "$FASTQ_A_R1" "$FASTQ_A_R2" "$CBQ"; then
    echo "ERROR: failed to encode paired FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

if ! encode_pair_cbq "$FASTQ_B_R1" "$FASTQ_B_R2" "$CBQ_B"; then
    echo "ERROR: failed to encode second paired FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

if ! encode_single_cbq "$OUT_ROOT/inputs/reads_SE.fastq" "$SE_CBQ"; then
    echo "ERROR: failed to encode single-end FASTQ to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_single.stderr" >&2
    exit 1
fi

"$STAR_BIN" \
    --runMode genomeGenerate \
    --genomeDir "$OUT_ROOT/genome" \
    --genomeFastaFiles "$OUT_ROOT/inputs/genome.fa" \
    --runThreadN 1 \
    --genomeSAindexNbases 3 \
    --genomeChrBinNbits 4 \
    > "$OUT_ROOT/genome_generate.stdout" \
    2> "$OUT_ROOT/genome_generate.stderr"

COMMON_ARGS=(
    --runMode alignReads
    --genomeDir "$OUT_ROOT/genome"
    --runThreadN 1
    --readNameSeparator /
    --outSAMtype SAM
    --outFilterMultimapNmax 1
    --outFilterMismatchNmax 0
    --alignIntronMax 1
    --alignMatesGapMax 1000
)
BASE_ATTRS=(--outSAMattributes NH HI AS nM)
RG_ATTRS=(--outSAMattributes NH HI AS nM RG)

"$STAR_BIN" "${COMMON_ARGS[@]}" "${BASE_ATTRS[@]}" \
    --readFilesType Fastx \
    --readFilesIn "$FASTQ_A_R1" "$FASTQ_A_R2" \
    --outFileNamePrefix "$OUT_ROOT/fastx/" \
    > "$OUT_ROOT/fastx/star.stdout" \
    2> "$OUT_ROOT/fastx/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" "${BASE_ATTRS[@]}" \
    --readFilesType Fastx \
    --readFilesIn "$FASTQ_B_R1" "$FASTQ_B_R2" \
    --outFileNamePrefix "$OUT_ROOT/fastx_b/" \
    > "$OUT_ROOT/fastx_b/star.stdout" \
    2> "$OUT_ROOT/fastx_b/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" "${BASE_ATTRS[@]}" \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ" \
    --outFileNamePrefix "$OUT_ROOT/binseq/" \
    > "$OUT_ROOT/binseq/star.stdout" \
    2> "$OUT_ROOT/binseq/star.stderr"

printf '%s\t-\tID:lane1\n' "$CBQ" > "$OUT_ROOT/inputs/manifest.tsv"
"$STAR_BIN" "${COMMON_ARGS[@]}" "${BASE_ATTRS[@]}" \
    --readFilesType Binseq PE \
    --readFilesManifest "$OUT_ROOT/inputs/manifest.tsv" \
    --outFileNamePrefix "$OUT_ROOT/binseq_manifest/" \
    > "$OUT_ROOT/binseq_manifest/star.stdout" \
    2> "$OUT_ROOT/binseq_manifest/star.stderr"

RG_LINES=(--outSAMattrRGline ID:sampleA , ID:sampleB)

"$STAR_BIN" "${COMMON_ARGS[@]}" "${RG_ATTRS[@]}" "${RG_LINES[@]}" \
    --readFilesType Fastx \
    --readFilesIn "$FASTQ_A_R1,$FASTQ_B_R1" "$FASTQ_A_R2,$FASTQ_B_R2" \
    --outFileNamePrefix "$OUT_ROOT/fastx_multi/" \
    > "$OUT_ROOT/fastx_multi/star.stdout" \
    2> "$OUT_ROOT/fastx_multi/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" "${RG_ATTRS[@]}" "${RG_LINES[@]}" \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ,$CBQ_B" \
    --outFileNamePrefix "$OUT_ROOT/binseq_multi/" \
    > "$OUT_ROOT/binseq_multi/star.stdout" \
    2> "$OUT_ROOT/binseq_multi/star.stderr"

printf '%s\t-\tID:sampleA\n%s\t-\tID:sampleB\n' "$CBQ" "$CBQ_B" > "$OUT_ROOT/inputs/manifest_multi.tsv"
"$STAR_BIN" "${COMMON_ARGS[@]}" "${RG_ATTRS[@]}" \
    --readFilesType Binseq PE \
    --readFilesManifest "$OUT_ROOT/inputs/manifest_multi.tsv" \
    --outFileNamePrefix "$OUT_ROOT/binseq_multi_manifest/" \
    > "$OUT_ROOT/binseq_multi_manifest/star.stdout" \
    2> "$OUT_ROOT/binseq_multi_manifest/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" "${BASE_ATTRS[@]}" \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ,$CBQ_B" \
    --batchMode 1 \
    --outFileNamePrefix "$OUT_ROOT/batch_binseq/" \
    > "$OUT_ROOT/batch_binseq/star.stdout" \
    2> "$OUT_ROOT/batch_binseq/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" "${RG_ATTRS[@]}" \
    --readFilesType Binseq PE \
    --readFilesManifest "$OUT_ROOT/inputs/manifest_multi.tsv" \
    --batchMode 1 \
    --outFileNamePrefix "$OUT_ROOT/batch_binseq_manifest/" \
    > "$OUT_ROOT/batch_binseq_manifest/star.stdout" \
    2> "$OUT_ROOT/batch_binseq_manifest/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" "${BASE_ATTRS[@]}" \
    --readFilesType Fastx \
    --readFilesIn "$OUT_ROOT/inputs/reads_SE.fastq" \
    --outFileNamePrefix "$OUT_ROOT/fastx_se/" \
    > "$OUT_ROOT/fastx_se/star.stdout" \
    2> "$OUT_ROOT/fastx_se/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" "${BASE_ATTRS[@]}" \
    --readFilesType Binseq SE \
    --readFilesIn "$SE_CBQ" \
    --outFileNamePrefix "$OUT_ROOT/binseq_se/" \
    > "$OUT_ROOT/binseq_se/star.stdout" \
    2> "$OUT_ROOT/binseq_se/star.stderr"

awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/fastx/Aligned.out.sam" > "$OUT_ROOT/fastx/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/fastx_b/Aligned.out.sam" > "$OUT_ROOT/fastx_b/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/binseq/Aligned.out.sam" > "$OUT_ROOT/binseq/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/binseq_manifest/Aligned.out.sam" > "$OUT_ROOT/binseq_manifest/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/fastx_multi/Aligned.out.sam" > "$OUT_ROOT/fastx_multi/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/binseq_multi/Aligned.out.sam" > "$OUT_ROOT/binseq_multi/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/binseq_multi_manifest/Aligned.out.sam" > "$OUT_ROOT/binseq_multi_manifest/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/fastx_se/Aligned.out.sam" > "$OUT_ROOT/fastx_se/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/binseq_se/Aligned.out.sam" > "$OUT_ROOT/binseq_se/body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/batch_binseq/alignments/sampleA/sampleA_Aligned.out.sam" > "$OUT_ROOT/batch_binseq/sampleA.body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/batch_binseq/alignments/sampleB/sampleB_Aligned.out.sam" > "$OUT_ROOT/batch_binseq/sampleB.body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/batch_binseq_manifest/alignments/sampleA/sampleA_Aligned.out.sam" > "$OUT_ROOT/batch_binseq_manifest/sampleA.body.sam"
awk 'substr($0, 1, 1) != "@" { print }' "$OUT_ROOT/batch_binseq_manifest/alignments/sampleB/sampleB_Aligned.out.sam" > "$OUT_ROOT/batch_binseq_manifest/sampleB.body.sam"

cmp -s "$OUT_ROOT/fastx/body.sam" "$OUT_ROOT/binseq/body.sam"
cmp -s "$OUT_ROOT/fastx/body.sam" "$OUT_ROOT/binseq_manifest/body.sam"
cmp -s "$OUT_ROOT/fastx_multi/body.sam" "$OUT_ROOT/binseq_multi/body.sam"
cmp -s "$OUT_ROOT/fastx_multi/body.sam" "$OUT_ROOT/binseq_multi_manifest/body.sam"
cmp -s "$OUT_ROOT/fastx/body.sam" "$OUT_ROOT/batch_binseq/sampleA.body.sam"
cmp -s "$OUT_ROOT/fastx_b/body.sam" "$OUT_ROOT/batch_binseq/sampleB.body.sam"
cmp -s "$OUT_ROOT/fastx_se/body.sam" "$OUT_ROOT/binseq_se/body.sam"

grep -q $'\tRG:Z:sampleA' "$OUT_ROOT/batch_binseq_manifest/sampleA.body.sam"
grep -q $'\tRG:Z:sampleB' "$OUT_ROOT/batch_binseq_manifest/sampleB.body.sam"
if grep -q $'\tRG:Z:sampleB' "$OUT_ROOT/batch_binseq_manifest/sampleA.body.sam"; then
    echo "ERROR: sampleA batch manifest output contains sampleB read group" >&2
    exit 1
fi
if grep -q $'\tRG:Z:sampleA' "$OUT_ROOT/batch_binseq_manifest/sampleB.body.sam"; then
    echo "ERROR: sampleB batch manifest output contains sampleA read group" >&2
    exit 1
fi

echo "PASS: STAR CBQ input smoke completed at $OUT_ROOT"
