#!/usr/bin/env bash
# =============================================================================
# STAR Suite vs nf-core/rnaseq star_salmon — bulk RNA-seq comparison harness
# -----------------------------------------------------------------------------
# Runs, on the SAME real human reads and the SAME reference, both:
#
#   (A) the reference chain  : Trim Galore -> STAR (2.7.x) -> Salmon
#                              (the tools nf-core/rnaseq --aligner star_salmon uses)
#   (B) STAR Suite integrated: one binary doing trim -> align -> quant
#       (B') STAR Suite BAM -> Salmon : exact star_salmon-equivalent quant
#
# It EMITS every intermediate (trimmed reads, genome BAM, transcriptome BAM,
# quant.sf) and COMPARES everything (wall time + peak RSS, read counts, BAM
# stats, per-transcript and per-gene quant concordance), then writes a Markdown
# report.
#
# Modes:
#   --mode chr22 : real human reads aligned to a chr22-only reference (fast).
#   --mode full  : full GRCh38 reference (the production number; slower).
#
# Nothing here is bespoke to our machines beyond the reference/binary paths at
# the top, which are CLI-overridable. Run `./run_compare.sh --help`.
# =============================================================================
set -euo pipefail

# ---- defaults (all overridable) ---------------------------------------------
MODE="chr22"
THREADS="$(nproc 2>/dev/null || echo 8)"
SALMON_THREADS="${SALMON_THREADS:-}"  # Default is set after --mode parsing.
SPOTS="2000000"                 # reads (pairs) to fetch for the fast path
ACCESSION="SRR4422207"          # public human PE bulk (GSE88509 / GSM2344101)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR=""                       # default set after arg parse
READS_DIR=""                    # provide pre-downloaded R1/R2 to skip SRA

STAR_SUITE_BIN="${STAR_SUITE_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}"
STAR_UPSTREAM_BIN="${STAR_UPSTREAM_BIN:-/usr/local/bin/STAR}"
SALMON_BIN="${SALMON_BIN:-salmon}"
TRIMGALORE_BIN="${TRIMGALORE_BIN:-trim_galore}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
GFFREAD_BIN="${GFFREAD_BIN:-gffread}"

# self-consistent chr-prefixed GRCh38 (cellranger-style) used for chr22 build
SRC_GENOME="${SRC_GENOME:-/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-cellranger/fasta/genome.fa}"
SRC_GTF="${SRC_GTF:-/storage/autoindex_110_44/refdata-gex-GRCh38-autoindex11044-cellranger/genes/genes.gtf}"
# full-mode prebuilt STAR index + transcriptome
FULL_STAR_INDEX="${FULL_STAR_INDEX:-/storage/autoindex_110_44/bulk_index}"
FULL_TRANSCRIPTOME="${FULL_TRANSCRIPTOME:-/storage/autoindex_110_44/bulk_index/transcriptome.fa}"

usage() { sed -n '2,33p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0; }
log()   { printf '\n\033[1;34m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }
timed() { # timed <label> <logfile> -- cmd... : capture wall + Max RSS
  local label="$1" tlog="$2"; shift 2
  log "RUN ${label}: $*"
  /usr/bin/time -v -o "${tlog}" "$@"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --salmon-threads) SALMON_THREADS="$2"; shift 2;;
    --spots) SPOTS="$2"; shift 2;;
    --accession) ACCESSION="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --reads-dir) READS_DIR="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
[[ -z "${OUTDIR}" ]] && OUTDIR="${SCRIPT_DIR}/runs/${MODE}_$(date +%Y%m%d_%H%M%S)"
if [[ -z "${SALMON_THREADS}" ]]; then
  if [[ "${MODE}" == "chr22" ]]; then
    SALMON_THREADS=1  # Salmon's online FLD is thread-order sensitive on chr22.
  else
    SALMON_THREADS="${THREADS}"
  fi
fi
mkdir -p "${OUTDIR}"
REF_DIR="${OUTDIR}/reference"; READS_OUT="${OUTDIR}/reads"
REFCHAIN="${OUTDIR}/A_reference_chain"; SUITE="${OUTDIR}/B_starsuite"
SUITE_SALMON="${OUTDIR}/B2_starsuite_bam_salmon"; CMP="${OUTDIR}/compare"
mkdir -p "${REF_DIR}" "${READS_OUT}" "${REFCHAIN}/trim" "${REFCHAIN}/star" \
         "${REFCHAIN}/salmon" "${SUITE}" "${SUITE_SALMON}" "${CMP}"

log "mode=${MODE} threads=${THREADS} salmon_threads=${SALMON_THREADS} outdir=${OUTDIR}"
echo "${MODE}"   > "${OUTDIR}/.mode"

# ---- version guard --------------------------------------------------------
# The reference baseline must be the SAME upstream STAR that STAR Suite is built
# on, otherwise an alignment difference could be the upstream version bump rather
# than STAR Suite's changes. STAR Suite reports its base via --upstream-version,
# and nf-core/rnaseq pins exactly that version (star=2.7.11b), so matching it is
# doubly correct.
SUITE_BASE="$("${STAR_SUITE_BIN}" --upstream-version 2>/dev/null | tr -d '[:space:]')"
UP_VER="$("${STAR_UPSTREAM_BIN}" --version 2>/dev/null | tr -d '[:space:]')"
log "STAR Suite base=${SUITE_BASE:-?}; reference baseline STAR=${UP_VER:-?}; nf-core/rnaseq pins star=${SUITE_BASE:-2.7.11b}"
if [[ -n "${SUITE_BASE}" && "${UP_VER}" != "${SUITE_BASE}" ]]; then
  echo "ERROR: reference baseline STAR (${STAR_UPSTREAM_BIN}) is '${UP_VER}', but STAR Suite is built on '${SUITE_BASE}'." >&2
  echo "       For a clean comparison use upstream STAR ${SUITE_BASE} (also what nf-core/rnaseq pins)." >&2
  echo "       e.g. STAR_UPSTREAM_BIN=/path/to/STAR-${SUITE_BASE} ./run_compare.sh ..." >&2
  exit 3
fi

# =============================================================================
# 1. Reference
# =============================================================================
if [[ "${MODE}" == "chr22" ]]; then
  GENOME="${REF_DIR}/chr22.fa"; GTF="${REF_DIR}/chr22.gtf"
  TRANSCRIPTOME="${REF_DIR}/transcriptome.fa"; STAR_INDEX="${REF_DIR}/star_index"
  if [[ ! -s "${STAR_INDEX}/SAindex" ]]; then
    log "Building chr22 reference (genome + gtf + transcriptome + STAR index)"
    "${SAMTOOLS_BIN}" faidx "${SRC_GENOME}" chr22 > "${GENOME}"
    "${SAMTOOLS_BIN}" faidx "${GENOME}"
    awk -F'\t' 'BEGIN{OFS="\t"} $1=="chr22"' "${SRC_GTF}" > "${GTF}"
    "${GFFREAD_BIN}" -w "${TRANSCRIPTOME}" -g "${GENOME}" "${GTF}"
    mkdir -p "${STAR_INDEX}"
    "${STAR_UPSTREAM_BIN}" --runMode genomeGenerate --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX}" --genomeFastaFiles "${GENOME}" \
      --sjdbGTFfile "${GTF}" --sjdbOverhang 100 --genomeSAindexNbases 11 \
      --outTmpDir "${STAR_INDEX}/_tmp"
  else
    log "Reusing chr22 reference at ${REF_DIR}"
  fi
  # STAR Suite TranscriptVB reads the transcriptome FASTA from inside the genome dir
  ln -sf "$(readlink -f "${TRANSCRIPTOME}")" "${STAR_INDEX}/transcriptome.fa"
else
  GTF="${SRC_GTF}"; TRANSCRIPTOME="${FULL_TRANSCRIPTOME}"; STAR_INDEX="${FULL_STAR_INDEX}"
  log "Full mode: STAR index=${STAR_INDEX} transcriptome=${TRANSCRIPTOME}"
fi
# transcript -> gene map for gene-level aggregation
python3 - "${GTF}" "${REF_DIR}/tx2gene.tsv" <<'PY'
import re, sys
gtf, out = sys.argv[1], sys.argv[2]
seen = {}
with open(gtf) as fh, open(out, "w") as o:
    for line in fh:
        if line.startswith("#") or "\ttranscript\t" not in line: continue
        t = re.search(r'transcript_id "([^"]+)"', line)
        g = re.search(r'gene_id "([^"]+)"', line)
        if t and g and t.group(1) not in seen:
            seen[t.group(1)] = 1
            o.write(f"{t.group(1)}\t{g.group(1)}\n")
print(f"tx2gene: {len(seen)} transcripts -> {out}")
PY

# =============================================================================
# 2. Reads (public human PE fixture; deterministic first N spots)
# =============================================================================
if [[ -n "${READS_DIR}" ]]; then
  R1="$(ls "${READS_DIR}"/*_1*.f*q.gz | head -1)"
  R2="$(ls "${READS_DIR}"/*_2*.f*q.gz | head -1)"
  log "Using provided reads: ${R1} ; ${R2}"
else
  R1="${READS_OUT}/${ACCESSION}_1.fastq.gz"; R2="${READS_OUT}/${ACCESSION}_2.fastq.gz"
  if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
    log "Fetching ${ACCESSION} first ${SPOTS} spots (fastq-dump)"
    ( cd "${READS_OUT}" && fastq-dump --split-files --gzip -X "${SPOTS}" "${ACCESSION}" )
  else
    log "Reusing reads at ${READS_OUT}"
  fi
fi

# =============================================================================
# 2b. chr22 read fixture — keep only reads that map to chr22 (dense parity signal)
#     Built once and cached under fixtures/chr22/, then reused. This is what makes
#     the chr22 parity check fast and dense (like nf-core's own chr22 test data),
#     rather than feeding whole-transcriptome reads to a chr22 reference.
# =============================================================================
if [[ "${MODE}" == "chr22" ]]; then
  FIX_DIR="${SCRIPT_DIR}/fixtures/chr22"; mkdir -p "${FIX_DIR}"
  FR1="${FIX_DIR}/${ACCESSION}_chr22_1.fastq.gz"; FR2="${FIX_DIR}/${ACCESSION}_chr22_2.fastq.gz"
  if [[ ! -s "${FR1}" || ! -s "${FR2}" ]]; then
    log "Building chr22 read fixture (align -> collect chr22 read IDs -> subset original FASTQs)"
    FX="${OUTDIR}/_fixture"; mkdir -p "${FX}"
    "${STAR_UPSTREAM_BIN}" --runMode alignReads --runThreadN "${THREADS}" \
      --genomeDir "${STAR_INDEX}" --readFilesIn "${R1}" "${R2}" --readFilesCommand zcat \
      --outSAMtype BAM Unsorted --outFileNamePrefix "${FX}/" --outTmpDir "${FX}/_tmp"
    "${SAMTOOLS_BIN}" view -@ "${THREADS}" -F 0x4 "${FX}/Aligned.out.bam" \
      | cut -f1 | sort -u > "${FX}/chr22_ids.txt"
    log "chr22-mapping reads: $(wc -l < "${FX}/chr22_ids.txt") unique IDs"
    python3 "${SCRIPT_DIR}/subset_fastq_by_id.py" "${FX}/chr22_ids.txt" "${R1}" "${FR1}"
    python3 "${SCRIPT_DIR}/subset_fastq_by_id.py" "${FX}/chr22_ids.txt" "${R2}" "${FR2}"
  else
    log "Reusing cached chr22 read fixture at ${FIX_DIR}"
  fi
  R1="${FR1}"; R2="${FR2}"
  log "chr22 fixture reads: $(zcat "${R1}" | wc -l | awk '{print $1/4}') pairs"
fi

# =============================================================================
# 3. (A) Reference chain: Trim Galore -> STAR -> Salmon
# =============================================================================
timed "A.trim_galore" "${REFCHAIN}/trim/time.log" \
  "${TRIMGALORE_BIN}" --paired --gzip --cores "${THREADS}" \
  --output_dir "${REFCHAIN}/trim" "${R1}" "${R2}"
TV1="$(ls "${REFCHAIN}"/trim/*_val_1.fq.gz | head -1)"
TV2="$(ls "${REFCHAIN}"/trim/*_val_2.fq.gz | head -1)"

timed "A.STAR_upstream" "${REFCHAIN}/star/time.log" \
  "${STAR_UPSTREAM_BIN}" --runMode alignReads --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX}" --readFilesIn "${TV1}" "${TV2}" --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
  --outFileNamePrefix "${REFCHAIN}/star/" --outTmpDir "${REFCHAIN}/star/_tmp"

timed "A.salmon" "${REFCHAIN}/salmon/time.log" \
  "${SALMON_BIN}" quant -t "${TRANSCRIPTOME}" -l A -p "${SALMON_THREADS}" \
  -a "${REFCHAIN}/star/Aligned.toTranscriptome.out.bam" -o "${REFCHAIN}/salmon"

# =============================================================================
# 4. (B) STAR Suite integrated: one binary, trim -> align -> quant
# =============================================================================
timed "B.starsuite_integrated" "${SUITE}/time.log" \
  "${STAR_SUITE_BIN}" --runMode alignReads --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX}" --readFilesIn "${R1}" "${R2}" \
  --trimCutadapt Yes \
  --outSAMtype BAM SortedByCoordinate --outBAMsortMethod samtools \
  --quantMode TranscriptomeSAM TranscriptVB --quantVBLibType A \
  --transcriptomeFasta "${TRANSCRIPTOME}" \
  --outFileNamePrefix "${SUITE}/" --outTmpDir "${SUITE}/_tmp"

# =============================================================================
# 5. (B') STAR Suite BAM -> Salmon : exact star_salmon-equivalent quant
# =============================================================================
timed "B2.starsuite_bam_salmon" "${SUITE_SALMON}/time.log" \
  "${SALMON_BIN}" quant -t "${TRANSCRIPTOME}" -l A -p "${SALMON_THREADS}" \
  -a "${SUITE}/Aligned.toTranscriptome.out.bam" -o "${SUITE_SALMON}"

# =============================================================================
# 6. Emit BAM stats (everything is kept on disk under ${OUTDIR})
# =============================================================================
for tag in "A_reference_chain/star" "B_starsuite"; do
  bam="${OUTDIR}/${tag}/Aligned.sortedByCoord.out.bam"
  [[ -s "${bam}" ]] || continue
  "${SAMTOOLS_BIN}" flagstat "${bam}" > "${OUTDIR}/${tag}/flagstat.txt"
  "${SAMTOOLS_BIN}" index "${bam}" 2>/dev/null || true
done

# =============================================================================
# 7. Compare + report
# =============================================================================
log "Comparing outputs -> ${CMP}/report.md"
python3 "${SCRIPT_DIR}/compare_outputs.py" \
  --outdir "${OUTDIR}" --mode "${MODE}" --threads "${THREADS}" \
  --salmon-threads "${SALMON_THREADS}" \
  --tx2gene "${REF_DIR}/tx2gene.tsv"

log "DONE. Report: ${CMP}/report.md  (all intermediates kept under ${OUTDIR})"
cat "${CMP}/report.md"
