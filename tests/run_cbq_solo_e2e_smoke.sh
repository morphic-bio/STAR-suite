#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${STAR_BIN:-$ROOT_DIR/core/legacy/source/STAR}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_cbq_solo_e2e_smoke}"

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
mkdir -p "$OUT_ROOT"/{inputs,genome,fastq,cbq,cbq_level0,cbq_manifest,cbq_range,cbq_range_bam,cbq_range_cap}

python3 - "$OUT_ROOT" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(57421)
alphabet = "ACGT"
genome = "".join(rng.choice(alphabet) for _ in range(1200))
gene_start = 101
gene_end = 900
reads = [
    ("read001", 140, "ACGTACGTACGTACGT", "ACGTACGTACGT", "I"),
    ("read002", 330, "TGCATGCATGCATGCA", "TGCATGCATGCA", "H"),
    ("read003", 610, "ACGTACGTACGTACGT", "GATCGATCGATC", "J"),
]

(out / "inputs" / "genome.fa").write_text(">chrSynthetic\n" + genome + "\n", encoding="ascii")
(out / "inputs" / "genes.gtf").write_text(
    "\n".join(
        [
            f'chrSynthetic\tSTARsuite\tgene\t{gene_start}\t{gene_end}\t.\t+\t.\tgene_id "GENE1"; gene_name "GENE1";',
            f'chrSynthetic\tSTARsuite\ttranscript\t{gene_start}\t{gene_end}\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";',
            f'chrSynthetic\tSTARsuite\texon\t{gene_start}\t{gene_end}\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; gene_name "GENE1";',
        ]
    )
    + "\n",
    encoding="ascii",
)
(out / "inputs" / "whitelist.txt").write_text(
    "ACGTACGTACGTACGT\nTGCATGCATGCATGCA\n",
    encoding="ascii",
)

with (out / "inputs" / "barcodes_R1.fastq").open("wt", encoding="ascii") as r1, \
     (out / "inputs" / "cdna_R2.fastq").open("wt", encoding="ascii") as r2:
    for name, start, barcode, umi, qchar in reads:
        barcode_read = barcode + umi
        cdna_read = genome[start:start + 60]
        r1.write(f"@{name} 1:N:0:ACGT\n{barcode_read}\n+\n{qchar * len(barcode_read)}\n")
        r2.write(f"@{name} 2:N:0:ACGT\n{cdna_read}\n+\n{qchar * len(cdna_read)}\n")
PY

BARCODE_FASTQ="$OUT_ROOT/inputs/barcodes_R1.fastq"
CDNA_FASTQ="$OUT_ROOT/inputs/cdna_R2.fastq"
CBQ="$OUT_ROOT/inputs/solo_cdna_barcode.cbq"
CBQ_LEVEL0="$OUT_ROOT/inputs/solo_cdna_barcode_level0.cbq"

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

# STARsolo expects cDNA as mate 1 and CB/UMI as mate 2 in --readFilesIn.
# Encode the CBQ in that same mate order.
if ! encode_pair_cbq "$CDNA_FASTQ" "$BARCODE_FASTQ" "$CBQ"; then
    echo "ERROR: failed to encode STARsolo FASTQ pair to CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

if ! encode_pair_cbq "$CDNA_FASTQ" "$BARCODE_FASTQ" "$CBQ_LEVEL0" 0; then
    echo "ERROR: failed to encode STARsolo FASTQ pair to level-0 CBQ with bqtools" >&2
    echo "See $OUT_ROOT/encode_pair.stderr" >&2
    exit 1
fi

"$STAR_BIN" \
    --runMode genomeGenerate \
    --genomeDir "$OUT_ROOT/genome" \
    --genomeFastaFiles "$OUT_ROOT/inputs/genome.fa" \
    --sjdbGTFfile "$OUT_ROOT/inputs/genes.gtf" \
    --sjdbOverhang 59 \
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
    --outSAMtype None
    --outFilterMultimapNmax 1
    --outFilterMismatchNmax 0
    --alignIntronMax 1
    --alignMatesGapMax 1000
    --clipAdapterType CellRanger4
    --soloType CB_UMI_Simple
    --soloCBstart 1
    --soloCBlen 16
    --soloUMIstart 17
    --soloUMIlen 12
    --soloBarcodeReadLength 0
    --soloCBwhitelist "$OUT_ROOT/inputs/whitelist.txt"
    --soloCBmatchWLtype Exact
    --soloUMIdedup Exact
    --soloCellFilter None
    --soloFeatures Gene
)

RANGE_ARGS=("${COMMON_ARGS[@]}")
for i in "${!RANGE_ARGS[@]}"; do
    if [[ "${RANGE_ARGS[$i]}" == "--runThreadN" ]]; then
        RANGE_ARGS[$((i + 1))]=4
        break
    fi
done

BAM_RANGE_ARGS=()
skip_next=0
for i in "${!RANGE_ARGS[@]}"; do
    if (( skip_next )); then
        skip_next=0
        continue
    fi
    if [[ "${RANGE_ARGS[$i]}" == "--outSAMtype" ]]; then
        BAM_RANGE_ARGS+=(--outSAMtype BAM SortedByCoordinate)
        skip_next=1
    else
        BAM_RANGE_ARGS+=("${RANGE_ARGS[$i]}")
    fi
done

"$STAR_BIN" "${COMMON_ARGS[@]}" \
    --readFilesType Fastx \
    --readFilesIn "$CDNA_FASTQ" "$BARCODE_FASTQ" \
    --outFileNamePrefix "$OUT_ROOT/fastq/" \
    > "$OUT_ROOT/fastq/star.stdout" \
    2> "$OUT_ROOT/fastq/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ" \
    --outFileNamePrefix "$OUT_ROOT/cbq/" \
    > "$OUT_ROOT/cbq/star.stdout" \
    2> "$OUT_ROOT/cbq/star.stderr"

"$STAR_BIN" "${COMMON_ARGS[@]}" \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ_LEVEL0" \
    --outFileNamePrefix "$OUT_ROOT/cbq_level0/" \
    > "$OUT_ROOT/cbq_level0/star.stdout" \
    2> "$OUT_ROOT/cbq_level0/star.stderr"

"$STAR_BIN" "${RANGE_ARGS[@]}" \
    --readFilesCbqRangeMode range \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ" \
    --outFileNamePrefix "$OUT_ROOT/cbq_range/" \
    > "$OUT_ROOT/cbq_range/star.stdout" \
    2> "$OUT_ROOT/cbq_range/star.stderr"

"$STAR_BIN" "${BAM_RANGE_ARGS[@]}" \
    --readFilesCbqRangeMode range \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ" \
    --outFileNamePrefix "$OUT_ROOT/cbq_range_bam/" \
    > "$OUT_ROOT/cbq_range_bam/star.stdout" \
    2> "$OUT_ROOT/cbq_range_bam/star.stderr"

"$STAR_BIN" "${RANGE_ARGS[@]}" \
    --readFilesCbqRangeMode range \
    --readMapNumber 2 \
    --readFilesType Binseq PE \
    --readFilesIn "$CBQ" \
    --outFileNamePrefix "$OUT_ROOT/cbq_range_cap/" \
    > "$OUT_ROOT/cbq_range_cap/star.stdout" \
    2> "$OUT_ROOT/cbq_range_cap/star.stderr"

printf '%s\t-\tID:solo_cbq\n' "$CBQ" > "$OUT_ROOT/inputs/manifest.tsv"
"$STAR_BIN" "${COMMON_ARGS[@]}" \
    --readFilesType Binseq PE \
    --readFilesManifest "$OUT_ROOT/inputs/manifest.tsv" \
    --outFileNamePrefix "$OUT_ROOT/cbq_manifest/" \
    > "$OUT_ROOT/cbq_manifest/star.stdout" \
    2> "$OUT_ROOT/cbq_manifest/star.stderr"

compare_solo_raw() {
    local observed_dir="$1"
    local rel
    for rel in barcodes.tsv features.tsv matrix.mtx; do
        cmp -s "$OUT_ROOT/fastq/Solo.out/Gene/raw/$rel" "$observed_dir/Solo.out/Gene/raw/$rel" || {
            echo "ERROR: STARsolo CBQ parity failed for $rel in $observed_dir" >&2
            diff -u "$OUT_ROOT/fastq/Solo.out/Gene/raw/$rel" "$observed_dir/Solo.out/Gene/raw/$rel" >&2 || true
            exit 1
        }
    done
}

require_nonzero_matrix() {
    python3 - "$1" <<'PY'
import sys

matrix_path = sys.argv[1]
with open(matrix_path, "rt", encoding="ascii") as handle:
    for line in handle:
        if line.startswith("%"):
            continue
        fields = line.split()
        if len(fields) != 3:
            raise SystemExit(f"invalid MatrixMarket size line in {matrix_path}: {line.rstrip()}")
        if int(fields[2]) <= 0:
            raise SystemExit(f"empty STARsolo matrix: {matrix_path}")
        break
    else:
        raise SystemExit(f"missing MatrixMarket size line: {matrix_path}")
PY
}

compare_solo_raw "$OUT_ROOT/cbq"
compare_solo_raw "$OUT_ROOT/cbq_level0"
compare_solo_raw "$OUT_ROOT/cbq_manifest"
compare_solo_raw "$OUT_ROOT/cbq_range"
require_nonzero_matrix "$OUT_ROOT/fastq/Solo.out/Gene/raw/matrix.mtx"
grep -q "CBQ indexed range reader: active" "$OUT_ROOT/cbq_range/Log.out"
grep -q "CBQ indexed range reader: active" "$OUT_ROOT/cbq_range_bam/Log.out"
grep -q "CBQ indexed range reader: active" "$OUT_ROOT/cbq_range_cap/Log.out"
[[ -s "$OUT_ROOT/cbq_range_bam/Aligned.sortedByCoord.out.bam" ]]
awk -F'|' '/Number of input reads/ { gsub(/[ \t]/, "", $2); if ($2 == "2") found=1 } END { exit(found ? 0 : 1) }' \
    "$OUT_ROOT/cbq_range_cap/Log.final.out"

avg_input_read_length() {
    awk -F'|' '/Average input read length/ { gsub(/[ \t]/, "", $2); print $2 }' "$1/Log.final.out"
}

compare_avg_input_read_length() {
    local observed_dir="$1"
    local expected observed
    expected="$(avg_input_read_length "$OUT_ROOT/fastq")"
    observed="$(avg_input_read_length "$observed_dir")"
    if [[ -z "$expected" || -z "$observed" || "$expected" != "$observed" ]]; then
        echo "ERROR: STARsolo CBQ average input read length mismatch for $observed_dir: expected $expected observed $observed" >&2
        exit 1
    fi
}

compare_avg_input_read_length "$OUT_ROOT/cbq"
compare_avg_input_read_length "$OUT_ROOT/cbq_level0"
compare_avg_input_read_length "$OUT_ROOT/cbq_manifest"
compare_avg_input_read_length "$OUT_ROOT/cbq_range"

echo "PASS: STARsolo CBQ E2E smoke completed at $OUT_ROOT"
