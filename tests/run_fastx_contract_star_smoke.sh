#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${STAR_BIN:-$ROOT_DIR/core/legacy/source/STAR}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_fastx_contract_star_smoke}"

if [[ ! -x "$STAR_BIN" ]]; then
    echo "ERROR: STAR binary not found or not executable: $STAR_BIN" >&2
    echo "Build it with: make -C core/legacy/source STAR" >&2
    exit 1
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{genome,reads,runs}

cat > "$OUT_ROOT/genome/genome.fa" <<'FASTA'
>chr1
ACGTTGCAACGATCGTACGTTAGCTAGGCTAACCGTTAACCGGTTAGCTAGCTAGGATCCGATCGTACCTAGGCTAACCGTTAACCGGTTAGCTAGCATCGATCGTACGATCGTAGCTAGCTAGGCTA
FASTA

"$STAR_BIN" \
    --runMode genomeGenerate \
    --runThreadN 1 \
    --genomeDir "$OUT_ROOT/genome" \
    --genomeFastaFiles "$OUT_ROOT/genome/genome.fa" \
    --genomeSAindexNbases 2 \
    --genomeChrBinNbits 4 \
    > "$OUT_ROOT/genome_generate.stdout" 2> "$OUT_ROOT/genome_generate.stderr"

seq1="ACGTTGCAACGATCGTACGTTAGCTAGGCTAACCGTTAACCGGTTA"
seq2="GCTAGCTAGGATCCGATCGTACCTAGGCTAACCGTTAACCGGTTAG"
qual1="$(printf '%*s' "${#seq1}" '' | tr ' ' I)"
qual2="$(printf '%*s' "${#seq2}" '' | tr ' ' I)"

{
    printf '@read001 1:N:0:1\n'
    printf '%s\n+\n%s\n' "$seq1" "$qual1"
} > "$OUT_ROOT/reads/lane1.fastq"

{
    printf '@read002 1:N:0:1\n'
    printf '%s\n+\n%s\n' "$seq2" "$qual2"
} > "$OUT_ROOT/reads/lane2.fastq"

gzip -c "$OUT_ROOT/reads/lane1.fastq" > "$OUT_ROOT/reads/lane1.fastq.gz"
printf '%s\t-\tID:lane1\n' "$OUT_ROOT/reads/lane1.fastq" > "$OUT_ROOT/reads/manifest.tsv"
cat > "$OUT_ROOT/reads/not_fastq.fa" <<'FASTA'
>read001
ACGTTGCAACGATCGTACGTTAGCTAGGCTAACCGTTAACCGGTTA
FASTA

run_star_case() {
    local name="$1"
    local expected_reads="$2"
    shift 2

    local out_dir="$OUT_ROOT/runs/$name"
    mkdir -p "$out_dir"

    "$STAR_BIN" \
        --runThreadN 1 \
        --genomeDir "$OUT_ROOT/genome" \
        --outFileNamePrefix "$out_dir/" \
        --outSAMtype None \
        "$@" \
        > "$out_dir/stdout" 2> "$out_dir/stderr"

    local observed_reads
    observed_reads="$(awk -F'|' '/Number of input reads/ {gsub(/[ \t]/, "", $2); print $2}' "$out_dir/Log.final.out")"
    if [[ "$observed_reads" != "$expected_reads" ]]; then
        echo "ERROR: $name expected $expected_reads input reads, observed ${observed_reads:-missing}" >&2
        echo "See $out_dir/Log.final.out" >&2
        exit 1
    fi

    if ! grep -F "Fastx input path: direct STAR chunk reader" "$out_dir/Log.out" >/dev/null; then
        echo "ERROR: $name did not use the direct STAR Fastx input path" >&2
        echo "See $out_dir/Log.out" >&2
        exit 1
    fi
}

run_star_case plain_fastq 1 \
    --readFilesIn "$OUT_ROOT/reads/lane1.fastq"

run_star_case gzip_internal 1 \
    --readFilesIn "$OUT_ROOT/reads/lane1.fastq.gz"

run_star_case gzip_command 1 \
    --readFilesIn "$OUT_ROOT/reads/lane1.fastq.gz" \
    --readFilesCommand zcat

run_star_case comma_lanes 2 \
    --readFilesIn "$OUT_ROOT/reads/lane1.fastq,$OUT_ROOT/reads/lane2.fastq"

run_star_case manifest 1 \
    --readFilesManifest "$OUT_ROOT/reads/manifest.tsv"

gate_dir="$OUT_ROOT/runs/y_gate_fasta"
mkdir -p "$gate_dir"
set +e
"$STAR_BIN" \
    --runThreadN 1 \
    --genomeDir "$OUT_ROOT/genome" \
    --outFileNamePrefix "$gate_dir/" \
    --outSAMtype None \
    --readFilesIn "$OUT_ROOT/reads/not_fastq.fa" \
    --emitYNoYFastq yes \
    > "$gate_dir/stdout" 2> "$gate_dir/stderr"
gate_status=$?
set -e

if [[ "$gate_status" -eq 0 ]]; then
    echo "ERROR: --emitYNoYFastq accepted a non-FASTQ input path" >&2
    exit 1
fi

if ! grep -R "requires FASTQ input" "$gate_dir" >/dev/null 2>&1; then
    echo "ERROR: --emitYNoYFastq non-FASTQ failure did not include the expected gate message" >&2
    exit 1
fi

echo "PASS: Fastx STAR contract smoke completed at $OUT_ROOT"
