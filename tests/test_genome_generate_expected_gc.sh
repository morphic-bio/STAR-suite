#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
OUTDIR="${OUTDIR:-$(mktemp -d /tmp/star_suite_genome_expected_gc_XXXXXX)}"

[[ -x "${STAR_BIN}" ]] || {
    echo "ERROR: STAR binary not executable: ${STAR_BIN}" >&2
    exit 2
}

mkdir -p "${OUTDIR}/ref" "${OUTDIR}/index"

cat > "${OUTDIR}/ref/genome.fa" <<'FASTA'
>chr1
AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT
FASTA

cat > "${OUTDIR}/ref/genes.gtf" <<'GTF'
chr1	synthetic	gene	1	120	.	+	.	gene_id "GENE1"; gene_name "Gene1";
chr1	synthetic	transcript	1	120	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "Gene1";
chr1	synthetic	exon	1	120	.	+	.	gene_id "GENE1"; transcript_id "TX1"; gene_name "Gene1";
GTF

run_index() {
    local label="$1"
    (
        cd "${OUTDIR}"
        "${STAR_BIN}" \
            --runMode genomeGenerate \
            --runThreadN 1 \
            --genomeDir "${OUTDIR}/index" \
            --genomeFastaFiles "${OUTDIR}/ref/genome.fa" \
            --sjdbGTFfile "${OUTDIR}/ref/genes.gtf" \
            --sjdbOverhang 20 \
            --genomeSAindexNbases 2 \
            --genomeChrBinNbits 4 \
            --genomeGenerateTranscriptome Yes \
            > "${OUTDIR}/${label}.stdout" 2> "${OUTDIR}/${label}.stderr"
    )
}

run_index first

EXPECTED_GC="${OUTDIR}/index/expected_gc.tsv"
[[ -s "${EXPECTED_GC}" ]] || {
    echo "ERROR: expected GC sidecar was not created: ${EXPECTED_GC}" >&2
    exit 1
}

rows="$(wc -l < "${EXPECTED_GC}")"
[[ "${rows}" == "101" ]] || {
    echo "ERROR: expected 101 GC bins, found ${rows}" >&2
    exit 1
}

awk '
    { sum += $2 }
    END {
        if (sum < 0.999 || sum > 1.001) {
            printf("ERROR: expected_gc.tsv probabilities sum to %.9f\n", sum) > "/dev/stderr";
            exit 1;
        }
    }
' "${EXPECTED_GC}"

grep -q "Successfully computed expected GC distribution" "${OUTDIR}/index/Log.out" || {
    echo "ERROR: genomeGenerate log did not report expected-GC creation" >&2
    exit 1
}

before_sha="$(sha256sum "${EXPECTED_GC}" | awk '{print $1}')"
run_index second
after_sha="$(sha256sum "${EXPECTED_GC}" | awk '{print $1}')"

[[ "${before_sha}" == "${after_sha}" ]] || {
    echo "ERROR: expected GC sidecar changed on second genomeGenerate run" >&2
    exit 1
}

grep -q "Expected GC distribution already exists; reusing" "${OUTDIR}/index/Log.out" || {
    echo "ERROR: second genomeGenerate did not report expected-GC reuse" >&2
    exit 1
}

echo "PASS: genomeGenerate expected-GC sidecar created and reused at ${OUTDIR}"
