#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
OUTDIR="${BAMSORT_BYTE_PARITY_OUTDIR:-${SCRIPT_DIR}/bamsort_samtools_byte_parity_output_$(date +%Y%m%d_%H%M%S)}"
THREADS="${BAMSORT_BYTE_PARITY_THREADS:-1}"
SORT_RAM="${BAMSORT_BYTE_PARITY_SORT_RAM:-1000}"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || fail "Missing required command: $1"
}

require_file() {
  [[ -f "$1" ]] || fail "Missing required file: $1"
}

compare_decompressed_bam() {
  local label="$1"
  local expected="$2"
  local observed="$3"
  local case_dir="$4"

  require_file "${expected}"
  require_file "${observed}"
  samtools quickcheck "${expected}" "${observed}"

  if ! cmp -s <(gzip -dc "${expected}") <(gzip -dc "${observed}"); then
    samtools view --no-PG -H "${expected}" > "${case_dir}/${label}.samtools.header.sam" 2>/dev/null || true
    samtools view --no-PG -H "${observed}" > "${case_dir}/${label}.star.header.sam" 2>/dev/null || true
    samtools view "${expected}" > "${case_dir}/${label}.samtools.body.sam" 2>/dev/null || true
    samtools view "${observed}" > "${case_dir}/${label}.star.body.sam" 2>/dev/null || true
    fail "${label}: decompressed BAM byte stream differs; see ${case_dir}"
  fi
}

run_case() {
  local mode="$1"
  shift
  local case_dir="${OUTDIR}/${mode}"
  local star_dir="${case_dir}/star"
  local samtools_bam="${case_dir}/samtools.sorted.bam"
  mkdir -p "${star_dir}"

  echo ">>> ${mode}: STAR unsorted+sorted run"
  "${STAR_BIN}" \
    --runMode alignReads \
    --runThreadN "${THREADS}" \
    --genomeDir "${OUTDIR}/ref/star_index" \
    --readFilesIn "$@" \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --outBAMsortMethod samtools \
    --limitBAMsortRAM "${SORT_RAM}" \
    --outSAMattributes NH HI AS nM NM \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 1 \
    --outFileNamePrefix "${star_dir}/" \
    --outTmpDir "${case_dir}/star_tmp" \
    > "${case_dir}/star.log" 2>&1

  require_file "${star_dir}/Aligned.out.bam"
  require_file "${star_dir}/Aligned.sortedByCoord.out.bam"

  echo ">>> ${mode}: samtools sort comparator"
  samtools sort --no-PG -@ 1 \
    -o "${samtools_bam}" \
    -T "${case_dir}/samtools_tmp" \
    "${star_dir}/Aligned.out.bam"

  compare_decompressed_bam "${mode}" "${samtools_bam}" "${star_dir}/Aligned.sortedByCoord.out.bam" "${case_dir}"

  local records
  records="$(samtools view -c "${star_dir}/Aligned.sortedByCoord.out.bam")"
  local coord_ties
  coord_ties="$(samtools view "${star_dir}/Aligned.sortedByCoord.out.bam" \
    | awk '{key=$3 ":" $4; n[key]++} END {ties=0; for (key in n) if (n[key] > 1) ties++; print ties+0}')"
  [[ "${coord_ties}" -gt 0 ]] || fail "${mode}: synthetic fixture did not produce coordinate ties"
  printf '%s\tPASS\t%s\t%s\n' "${mode}" "${records}" "${coord_ties}" >> "${OUTDIR}/summary.tsv"
  echo "PASS: ${mode} decompressed BAM byte stream matches samtools sort (${records} records, ${coord_ties} tied coordinates)"
}

require_cmd python3
require_cmd samtools
require_cmd gzip
require_cmd cmp
require_file "${STAR_BIN}"

SAMTOOLS_SORT_HELP="$(samtools sort --help 2>&1 || true)"
grep -q -- '--no-PG' <<< "${SAMTOOLS_SORT_HELP}" \
  || fail "samtools sort does not support --no-PG; cannot compare headers byte-for-byte"

rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}/ref" "${OUTDIR}/fastq"
printf 'mode\tstatus\trecords\tcoord_ties\n' > "${OUTDIR}/summary.tsv"

python3 - "${OUTDIR}" <<'PY'
from pathlib import Path
import random
import sys

root = Path(sys.argv[1])
random.seed(1)
bases = "ACGT"
seq = "".join(random.choice(bases) for _ in range(3000))
chr1 = seq[:1500]
chr2 = seq[1500:]
comp = str.maketrans("ACGT", "TGCA")

def revcomp(s):
    return s.translate(comp)[::-1]

(root / "ref" / "genome.fa").write_text(f">chr1\n{chr1}\n>chr2\n{chr2}\n")
qual = "I" * 60

with (root / "fastq" / "se.fastq").open("w") as handle:
    read_no = 1
    for start in [10, 70, 130, 190, 250, 310, 1570, 1630, 1690, 1750]:
        chrom = chr1 if start < 1500 else chr2
        local = start if start < 1500 else start - 1500
        read = chrom[local:local + 60]
        handle.write(f"@se{read_no:03d}\n{read}\n+\n{qual}\n")
        read_no += 1
    # Same-position forward/reverse ties exercise the samtools strand tie-break.
    for start in [430, 1810]:
        chrom = chr1 if start < 1500 else chr2
        local = start if start < 1500 else start - 1500
        read = chrom[local:local + 60]
        handle.write(f"@se{read_no:03d}\n{read}\n+\n{qual}\n")
        read_no += 1
        handle.write(f"@se{read_no:03d}\n{revcomp(read)}\n+\n{qual}\n")
        read_no += 1

with (root / "fastq" / "pe_R1.fastq").open("w") as r1, (root / "fastq" / "pe_R2.fastq").open("w") as r2:
    # Fragment starts are chosen so an R1 from one fragment shares coordinate
    # with an R2 from a prior fragment, exercising PE coordinate ties.
    for i, start in enumerate([20, 70, 120, 180, 240, 310, 390, 1570, 1640, 1710, 1780, 1850], 1):
        chrom = chr1 if start < 1500 else chr2
        local = start if start < 1500 else start - 1500
        mate1 = chrom[local:local + 60]
        mate2 = revcomp(chrom[local + 140:local + 200])
        r1.write(f"@pe{i:03d}\n{mate1}\n+\n{qual}\n")
        r2.write(f"@pe{i:03d}\n{mate2}\n+\n{qual}\n")
PY

echo "=== BAM Sort Samtools Byte Parity ==="
echo "STAR: ${STAR_BIN}"
echo "samtools: $(samtools --version | head -1)"
echo "Output: ${OUTDIR}"
echo "Threads: ${THREADS}"
echo "Sorter RAM: ${SORT_RAM}"
echo

"${STAR_BIN}" \
  --runMode genomeGenerate \
  --runThreadN 1 \
  --genomeDir "${OUTDIR}/ref/star_index" \
  --genomeFastaFiles "${OUTDIR}/ref/genome.fa" \
  --genomeSAindexNbases 5 \
  > "${OUTDIR}/genomeGenerate.log" 2>&1

run_case "SE" "${OUTDIR}/fastq/se.fastq"
run_case "PE" "${OUTDIR}/fastq/pe_R1.fastq" "${OUTDIR}/fastq/pe_R2.fastq"

echo
echo "PASS: SE and PE STAR-suite samtools-sorter BAMs are byte-identical to samtools sort after BGZF decompression"
echo "Summary: ${OUTDIR}/summary.tsv"
