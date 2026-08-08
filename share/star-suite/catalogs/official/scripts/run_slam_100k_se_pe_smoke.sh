#!/usr/bin/env bash
set -euo pipefail

# Reproducible 100K STAR-SLAM smoke for comparing R1-only SE vs R1/R2 PE.
# Dataset and reference locations are deliberately supplied by the caller.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"

FASTQ_DIR="${FASTQ_DIR:-}"
OUT_BASE="${OUT_BASE:-/tmp/slam_pe_100k_smoke_${STAMP}}"
STAR_BIN="${STAR_BIN:-STAR}"
STAR_INDEX="${STAR_INDEX:-}"
SNP_MASK="${SNP_MASK:--}"
THREADS="${THREADS:-16}"
READS="${READS:-100000}"
MIN_READ_COUNT="${MIN_READ_COUNT:-100}"
MIN_NTR_PEARSON="${MIN_NTR_PEARSON:-0.90}"
MAX_ABS_DELTA="${MAX_ABS_DELTA:-1.0}"
SLAM_CB_FORMAT="${SLAM_CB_FORMAT:-star}"
SLAM_MIN_CALLABLE_LENGTH="${SLAM_MIN_CALLABLE_LENGTH:-30}"
SLAM_STRANDNESS="${SLAM_STRANDNESS:-Sense}"

# Fixed trims from the noSU 100K PE calibration.
SE_TRIM5P="${SE_TRIM5P:-8}"
SE_TRIM3P="${SE_TRIM3P:-12}"
PE_TRIM5P_M1="${PE_TRIM5P_M1:-8}"
PE_TRIM3P_M1="${PE_TRIM3P_M1:-13}"
PE_TRIM5P_M2="${PE_TRIM5P_M2:-19}"
PE_TRIM3P_M2="${PE_TRIM3P_M2:-14}"

DRY_RUN=0
RESUME=0
COMPARE_ONLY=0
declare -a REQUESTED_SAMPLES=()

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Runs a 100K downsampled SLAM smoke for each selected sample in both modes:
R1-only SE and R1/R2 PE. It writes manifests, STAR commands, compare logs, and
compare/smoke_metrics.tsv under the output root.

The input directory, genome directory, and sample names have no public
dataset-specific defaults. Supply them explicitly. Output defaults to:
  ${OUT_BASE}

Selection:
  --sample NAME              Add one sample. May be repeated.
  --samples CSV              Add comma-separated sample names.

Paths:
  --fastq-dir DIR            Input FASTQ directory (required)
  --outdir DIR               Output root (default: ${OUT_BASE})
  --star-bin PATH            STAR binary (default: ${STAR_BIN})
  --genome-dir DIR           STAR genomeDir (required)
  --snp-mask PATH            SLAM SNP mask BED.GZ; use '-' to disable

Runtime:
  --threads N                STAR threads per run (default: ${THREADS})
  --reads N                  Reads to keep per sample/mate (default: ${READS})
  --resume                   Reuse existing downsampled FASTQs and completed runs.
  --compare-only             Recompute compare logs/metrics from an existing outdir.
  --dry-run                  Write commands/manifests but do not execute STAR.

Gates:
  --min-read-count N         ReadCount floor for NTR Pearson (default: ${MIN_READ_COUNT})
  --min-ntr-pearson X        Treatment-sample NTR Pearson gate (default: ${MIN_NTR_PEARSON})
  --max-abs-delta X          compare_slam_summary |delta NTR| gate (default: ${MAX_ABS_DELTA})
  --slam-min-callable N      Minimum callable bases after trimming (default: ${SLAM_MIN_CALLABLE_LENGTH})

Examples:
  $(basename "$0") --fastq-dir ./fastq --genome-dir ./index --sample control --sample treatment
  $(basename "$0") --outdir /tmp/slam_pe_100k_smoke_rerun --compare-only
EOF
}

log() {
  printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*"
}

die() {
  log "ERROR: $*" >&2
  exit 1
}

trap 'die "Command failed at line ${LINENO}: ${BASH_COMMAND}"' ERR

quote_cmd() {
  printf '%q ' "$@"
}

require_file() {
  [[ -f "$1" ]] || die "Missing file: $1"
}

require_dir() {
  [[ -d "$1" ]] || die "Missing directory: $1"
}

require_exec() {
  [[ -x "$1" ]] || die "Missing executable: $1"
}

csv_to_samples() {
  tr ',' '\n' <<< "$1" | sed '/^$/d'
}

is_nosu_sample() {
  local sample_lc
  sample_lc="$(tr '[:upper:]' '[:lower:]' <<< "$1")"
  [[ "$sample_lc" == *nosu* || "$sample_lc" == *no4su* || "$sample_lc" == *no-su* ]]
}

find_mate_fastq() {
  local sample="$1"
  local mate="$2"
  local direct
  local sample_id
  local -a candidates=()

  direct="${FASTQ_DIR}/${sample}_${mate}_001.fastq.gz"
  if [[ -f "$direct" ]]; then
    echo "$direct"
    return 0
  fi

  sample_id="$(sed -nE 's/^.*_(S[0-9]+)$/\1/p' <<< "$sample")"
  if [[ -n "$sample_id" ]]; then
    mapfile -t candidates < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*_${sample_id}_${mate}_001.fastq.gz" | sort)
    if [[ "${#candidates[@]}" -eq 1 ]]; then
      echo "${candidates[0]}"
      return 0
    fi
  fi

  die "Could not find ${mate} FASTQ for ${sample} under ${FASTQ_DIR}"
}

count_fastq_reads() {
  python3 - "$1" <<'PY'
import gzip
import sys

path = sys.argv[1]
lines = 0
with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
    for lines, _ in enumerate(handle, 1):
        pass
if lines % 4:
    raise SystemExit(f"{path}: FASTQ line count is not divisible by 4: {lines}")
print(lines // 4)
PY
}

downsample_fastq() {
  local src="$1"
  local dst="$2"
  local reads="$3"
  local observed

  if [[ "$RESUME" == "1" && -s "$dst" ]]; then
    observed="$(count_fastq_reads "$dst")"
    if [[ "$observed" == "$reads" ]]; then
      log "Reusing downsampled FASTQ: ${dst}"
      return 0
    fi
  fi

  mkdir -p "$(dirname "$dst")"
  python3 - "$src" "$dst" "$reads" <<'PY'
import gzip
import sys

src, dst, reads_s = sys.argv[1:4]
reads = int(reads_s)
limit = reads * 4
written = 0
with gzip.open(src, "rt", encoding="utf-8", errors="replace") as inp, \
     gzip.open(dst, "wt", encoding="utf-8") as out:
    for line in inp:
        if written >= limit:
            break
        out.write(line)
        written += 1
if written != limit:
    raise SystemExit(f"{src}: requested {reads} reads but found only {written // 4}")
PY

  observed="$(count_fastq_reads "$dst")"
  [[ "$observed" == "$reads" ]] || die "${dst}: expected ${reads} reads, observed ${observed}"
}

write_run_config() {
  mkdir -p "${OUT_BASE}/manifests" "${OUT_BASE}/compare"
  cat > "${OUT_BASE}/manifests/run_config.env" <<EOF
FASTQ_DIR=${FASTQ_DIR}
OUT_BASE=${OUT_BASE}
STAR_BIN=${STAR_BIN}
STAR_INDEX=${STAR_INDEX}
SNP_MASK=${SNP_MASK}
THREADS=${THREADS}
READS=${READS}
MIN_READ_COUNT=${MIN_READ_COUNT}
MIN_NTR_PEARSON=${MIN_NTR_PEARSON}
MAX_ABS_DELTA=${MAX_ABS_DELTA}
SLAM_CB_FORMAT=${SLAM_CB_FORMAT}
SLAM_MIN_CALLABLE_LENGTH=${SLAM_MIN_CALLABLE_LENGTH}
SLAM_STRANDNESS=${SLAM_STRANDNESS}
SE_TRIM5P=${SE_TRIM5P}
SE_TRIM3P=${SE_TRIM3P}
PE_TRIM5P_M1=${PE_TRIM5P_M1}
PE_TRIM3P_M1=${PE_TRIM3P_M1}
PE_TRIM5P_M2=${PE_TRIM5P_M2}
PE_TRIM3P_M2=${PE_TRIM3P_M2}
EOF
  printf 'sample\tsource_r1\tsource_r2\tdownsample_r1\tdownsample_r2\n' \
    > "${OUT_BASE}/manifests/samples.tsv"
  printf 'sample\tmode\tcommand\n' > "${OUT_BASE}/manifests/commands.tsv"
}

append_command_manifest() {
  local sample="$1"
  local mode="$2"
  shift 2
  printf '%s\t%s\t' "$sample" "$mode" >> "${OUT_BASE}/manifests/commands.tsv"
  quote_cmd "$@" >> "${OUT_BASE}/manifests/commands.tsv"
  printf '\n' >> "${OUT_BASE}/manifests/commands.tsv"
}

write_command_script() {
  local path="$1"
  shift
  {
    printf '#!/usr/bin/env bash\nset -euo pipefail\n\n'
    quote_cmd "$@"
    printf '\n'
  } > "$path"
  chmod +x "$path"
}

build_star_command() {
  local sample="$1"
  local mode="$2"
  local r1="$3"
  local r2="$4"
  local run_dir="$5"
  local -n out_cmd="$6"

  out_cmd=(
    "$STAR_BIN"
    --runThreadN "$THREADS"
    --genomeDir "$STAR_INDEX"
    --readFilesIn "$r1"
  )
  if [[ "$mode" == "pe" ]]; then
    out_cmd+=("$r2")
  fi
  out_cmd+=(
    --readFilesCommand zcat
    --outFileNamePrefix "${run_dir}/star_"
    --outTmpDir "${run_dir}/tmp"
    --outSAMtype None
    --quantMode TranscriptVB
    --quantTranscriptomeSAMoutput BanSingleEnd
    --quantVBgcBias 1
    --quantVBgenes 1
    --quantVBgenesMode Tximport
    --quantVBLibType A
    --slamQuantMode 1
    --slamGrandSlamOut 1
    --slamCbOut 1
    --slamCbFormat "$SLAM_CB_FORMAT"
    --slamCbOutFile "${run_dir}/SlamQuant.cB.${SLAM_CB_FORMAT}.tsv"
    --slamStrandness "$SLAM_STRANDNESS"
    --slamOutFile "${run_dir}/SlamQuant.out"
    --slamMinCallableLength "$SLAM_MIN_CALLABLE_LENGTH"
  )
  if [[ -n "$SNP_MASK" && "$SNP_MASK" != "-" ]]; then
    out_cmd+=(--slamSnpMaskIn "$SNP_MASK")
  fi
  if [[ "$mode" == "pe" ]]; then
    out_cmd+=(
      --slamCompatTrim5pMate1 "$PE_TRIM5P_M1"
      --slamCompatTrim3pMate1 "$PE_TRIM3P_M1"
      --slamCompatTrim5pMate2 "$PE_TRIM5P_M2"
      --slamCompatTrim3pMate2 "$PE_TRIM3P_M2"
    )
  else
    out_cmd+=(
      --slamCompatTrim5p "$SE_TRIM5P"
      --slamCompatTrim3p "$SE_TRIM3P"
    )
  fi
}

run_star_mode() {
  local sample="$1"
  local mode="$2"
  local r1="$3"
  local r2="$4"
  local run_dir="${OUT_BASE}/runs/${sample}/${mode}"
  local -a cmd

  mkdir -p "$run_dir"
  build_star_command "$sample" "$mode" "$r1" "$r2" "$run_dir" cmd
  write_command_script "${run_dir}/RUN_COMMAND.sh" "${cmd[@]}"
  append_command_manifest "$sample" "$mode" "${cmd[@]}"

  if [[ "$DRY_RUN" == "1" ]]; then
    log "DRY-RUN ${sample} ${mode}"
    quote_cmd "${cmd[@]}"
    printf '\n'
    return 0
  fi

  if [[ "$RESUME" == "1" && -s "${run_dir}/SlamQuant.out" && -s "${run_dir}/star_Log.final.out" ]]; then
    log "Reusing completed run: ${sample} ${mode}"
    return 0
  fi

  rm -rf "${run_dir}/tmp"
  log "Running ${sample} ${mode}"
  /usr/bin/time -v -o "${run_dir}/star.time.log" \
    "${cmd[@]}" \
    > "${run_dir}/stdout.log" \
    2> "${run_dir}/stderr.log"
}

compare_sample() {
  local sample="$1"
  local se="${OUT_BASE}/runs/${sample}/se/SlamQuant.out"
  local pe="${OUT_BASE}/runs/${sample}/pe/SlamQuant.out"
  local log_file="${OUT_BASE}/compare/${sample}.se_vs_pe.read${MIN_READ_COUNT}.log"
  local err_file="${OUT_BASE}/compare/${sample}.se_vs_pe.read${MIN_READ_COUNT}.stderr"
  local -a flags=(
    --allow-gene-set-mismatch
    --min-read-count "$MIN_READ_COUNT"
    --max-abs-delta "$MAX_ABS_DELTA"
  )

  require_file "$se"
  require_file "$pe"
  if ! is_nosu_sample "$sample"; then
    flags+=(--min-ntr-pearson "$MIN_NTR_PEARSON")
  fi

  log "Comparing ${sample} SE vs PE"
  if python3 "${REPO_ROOT}/scripts/compare_slam_summary.py" "$se" "$pe" "${flags[@]}" \
      > "$log_file" 2> "$err_file"; then
    return 0
  fi
  {
    printf '\n=== stderr ===\n'
    cat "$err_file"
  } >> "$log_file"
  return 1
}

write_metrics() {
  python3 - "$OUT_BASE" "$MIN_READ_COUNT" "${REQUESTED_SAMPLES[@]}" <<'PY'
import csv
import math
import sys
from pathlib import Path

root = Path(sys.argv[1])
min_read_count = float(sys.argv[2])
samples = sys.argv[3:]

def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0.0 or vy == 0.0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)

def ranks(vals):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and vals[order[j]] == vals[order[i]]:
            j += 1
        avg = (i + j - 1) / 2.0 + 1.0
        for k in range(i, j):
            out[order[k]] = avg
        i = j
    return out

def read_slam(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            row["Gene"]: {
                "ReadCount": float(row["ReadCount"]),
                "Conversions": float(row["Conversions"]),
                "Coverage": float(row["Coverage"]),
                "NTR": float(row["NTR"]),
            }
            for row in reader
        }

def read_sf(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["Name"]: float(row["NumReads"]) for row in reader}

def safe_round(x):
    if isinstance(x, float):
        if math.isnan(x):
            return "nan"
        return f"{x:.6f}"
    return str(x)

out = root / "compare" / "smoke_metrics.tsv"
fields = [
    "sample",
    "slam_shared_genes",
    "slam_se_genes",
    "slam_pe_genes",
    "slam_both_readcount_ge_min",
    "slam_ntr_pearson_ge_min",
    "slam_ntr_spearman_ge_min",
    "slam_readcount_pearson_shared",
    "slam_coverage_pearson_shared",
    "slam_conversions_pearson_shared",
    "slam_se_readcount_total",
    "slam_pe_readcount_total",
    "slam_se_coverage_total",
    "slam_pe_coverage_total",
    "slam_se_conversions_total",
    "slam_pe_conversions_total",
    "tximport_gene_numreads_pearson_all",
    "tximport_gene_numreads_pearson_ge1_both",
    "tximport_gene_ge1_both",
    "tximport_gene_se_total",
    "tximport_gene_pe_total",
]
with open(out, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(fields)
    for sample in samples:
        se_dir = root / "runs" / sample / "se"
        pe_dir = root / "runs" / sample / "pe"
        se = read_slam(se_dir / "SlamQuant.out")
        pe = read_slam(pe_dir / "SlamQuant.out")
        shared = sorted(set(se) & set(pe))
        high = [g for g in shared if se[g]["ReadCount"] >= min_read_count and pe[g]["ReadCount"] >= min_read_count]
        ntr_se = [se[g]["NTR"] for g in high]
        ntr_pe = [pe[g]["NTR"] for g in high]

        tx_se = read_sf(se_dir / "star_quant.genes.tximport.sf")
        tx_pe = read_sf(pe_dir / "star_quant.genes.tximport.sf")
        tx_shared = sorted(set(tx_se) & set(tx_pe))
        tx_ge1 = [g for g in tx_shared if tx_se[g] >= 1.0 and tx_pe[g] >= 1.0]

        row = [
            sample,
            len(shared),
            len(se),
            len(pe),
            len(high),
            pearson(ntr_se, ntr_pe),
            pearson(ranks(ntr_se), ranks(ntr_pe)),
            pearson([se[g]["ReadCount"] for g in shared], [pe[g]["ReadCount"] for g in shared]),
            pearson([se[g]["Coverage"] for g in shared], [pe[g]["Coverage"] for g in shared]),
            pearson([se[g]["Conversions"] for g in shared], [pe[g]["Conversions"] for g in shared]),
            sum(v["ReadCount"] for v in se.values()),
            sum(v["ReadCount"] for v in pe.values()),
            sum(v["Coverage"] for v in se.values()),
            sum(v["Coverage"] for v in pe.values()),
            sum(v["Conversions"] for v in se.values()),
            sum(v["Conversions"] for v in pe.values()),
            pearson([tx_se[g] for g in tx_shared], [tx_pe[g] for g in tx_shared]),
            pearson([tx_se[g] for g in tx_ge1], [tx_pe[g] for g in tx_ge1]),
            len(tx_ge1),
            sum(tx_se.values()),
            sum(tx_pe.values()),
        ]
        writer.writerow([safe_round(v) for v in row])
PY
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample) REQUESTED_SAMPLES+=("$2"); shift 2 ;;
    --samples)
      while IFS= read -r sample; do REQUESTED_SAMPLES+=("$sample"); done < <(csv_to_samples "$2")
      shift 2
      ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --outdir|--out-base) OUT_BASE="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; shift 2 ;;
    --genome-dir) STAR_INDEX="$2"; shift 2 ;;
    --snp-mask) SNP_MASK="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --reads) READS="$2"; shift 2 ;;
    --min-read-count) MIN_READ_COUNT="$2"; shift 2 ;;
    --min-ntr-pearson) MIN_NTR_PEARSON="$2"; shift 2 ;;
    --max-abs-delta) MAX_ABS_DELTA="$2"; shift 2 ;;
    --slam-min-callable) SLAM_MIN_CALLABLE_LENGTH="$2"; shift 2 ;;
    --slam-cb-format) SLAM_CB_FORMAT="$2"; shift 2 ;;
    --resume) RESUME=1; shift ;;
    --compare-only) COMPARE_ONLY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

((${#REQUESTED_SAMPLES[@]} > 0)) || die "at least one --sample or --samples value is required"

[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || die "--threads must be a positive integer"
[[ "$READS" =~ ^[0-9]+$ && "$READS" -gt 0 ]] || die "--reads must be a positive integer"
[[ "$MIN_READ_COUNT" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "--min-read-count must be numeric"
[[ "$SLAM_MIN_CALLABLE_LENGTH" =~ ^[0-9]+$ ]] || die "--slam-min-callable must be a non-negative integer"
[[ "$SLAM_CB_FORMAT" == "star" || "$SLAM_CB_FORMAT" == "ezbakr" ]] || die "--slam-cb-format must be star or ezbakr"

if [[ "$COMPARE_ONLY" == "0" ]]; then
  require_dir "$FASTQ_DIR"
  require_exec "$STAR_BIN"
  require_dir "$STAR_INDEX"
  if [[ -n "$SNP_MASK" && "$SNP_MASK" != "-" ]]; then
    require_file "$SNP_MASK"
  fi
fi

write_run_config

log "SLAM 100K SE/PE smoke"
log "Output root: ${OUT_BASE}"
log "Samples: ${REQUESTED_SAMPLES[*]}"

if [[ "$COMPARE_ONLY" == "0" ]]; then
  mkdir -p "${OUT_BASE}/fastq"
  for sample in "${REQUESTED_SAMPLES[@]}"; do
    src_r1="$(find_mate_fastq "$sample" R1)"
    src_r2="$(find_mate_fastq "$sample" R2)"
    ds_r1="${OUT_BASE}/fastq/${sample}_100k_R1_001.fastq.gz"
    ds_r2="${OUT_BASE}/fastq/${sample}_100k_R2_001.fastq.gz"

    if [[ "$DRY_RUN" == "0" ]]; then
      downsample_fastq "$src_r1" "$ds_r1" "$READS"
      downsample_fastq "$src_r2" "$ds_r2" "$READS"
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$src_r1" "$src_r2" "$ds_r1" "$ds_r2" \
      >> "${OUT_BASE}/manifests/samples.tsv"

    run_star_mode "$sample" se "$ds_r1" ""
    run_star_mode "$sample" pe "$ds_r1" "$ds_r2"
  done
fi

if [[ "$DRY_RUN" == "0" ]]; then
  for sample in "${REQUESTED_SAMPLES[@]}"; do
    compare_sample "$sample"
  done
  write_metrics
  log "Metrics: ${OUT_BASE}/compare/smoke_metrics.tsv"
fi

log "Done"
