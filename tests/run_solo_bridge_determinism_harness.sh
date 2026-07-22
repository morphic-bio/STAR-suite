#!/usr/bin/env bash
# Run two non-Flex Solo bridge FASTQ runs with stage-level determinism
# instrumentation and report where the outputs first diverge.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
ROOT="${A375_ROOT:-/storage/A375}"
FASTQ_ROOT="${A375_FASTQ_ROOT:-${ROOT}/fastqs/1k_CRISPR_5p_gemx_fastqs}"
GEX_DIR="${A375_GEX_DIR:-${FASTQ_ROOT}/gex}"
WHITELIST="${A375_WHITELIST:-${ROOT}/3M-5pgex-jan-2023.txt}"
GENOME_DIR="${A375_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
OUT_ROOT="${OUT_ROOT:-/mnt/pikachu/solo_bridge_determinism_$(date -u +%Y%m%dT%H%M%SZ)}"
THREADS_A="${THREADS_A:-32}"
THREADS_B="${THREADS_B:-${THREADS_A}}"
READ_MAP_NUMBER="${READ_MAP_NUMBER:-0}"
BRIDGE_DETERMINISM_CB="${BRIDGE_DETERMINISM_CB:-0}"
RUN_STAR="${RUN_STAR:-1}"
FORCE="${FORCE:-0}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --out-root DIR          Output root (default: ${OUT_ROOT})
  --threads-a N           Threads for run_a (default: ${THREADS_A})
  --threads-b N           Threads for run_b (default: ${THREADS_B})
  --read-map-number N     Limit STAR to N reads; 0 means full input (default: ${READ_MAP_NUMBER})
  --per-cb                Write and compare per-cell-barcode bridge digests
  --no-per-cb             Disable per-cell-barcode bridge digests
  --skip-run              Only compare existing run_a/run_b outputs
  --force                 Remove existing run_a/run_b before running
  --star-bin PATH         STAR binary (default: ${STAR_BIN})
  -h, --help              Show this help

Env overrides:
  A375_ROOT, A375_FASTQ_ROOT, A375_GEX_DIR, A375_WHITELIST, A375_GENOME_DIR,
  OUT_ROOT, THREADS_A, THREADS_B, READ_MAP_NUMBER, BRIDGE_DETERMINISM_CB,
  RUN_STAR, FORCE
EOF
}

log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }
die() { echo "FAIL: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-root)        OUT_ROOT="$2"; shift 2 ;;
    --threads-a)       THREADS_A="$2"; shift 2 ;;
    --threads-b)       THREADS_B="$2"; shift 2 ;;
    --read-map-number) READ_MAP_NUMBER="$2"; shift 2 ;;
    --per-cb)          BRIDGE_DETERMINISM_CB=1; shift ;;
    --no-per-cb)       BRIDGE_DETERMINISM_CB=0; shift ;;
    --skip-run)        RUN_STAR=0; shift ;;
    --force)           FORCE=1; shift ;;
    --star-bin)        STAR_BIN="$2"; shift 2 ;;
    -h|--help)         usage; exit 0 ;;
    *)                 die "Unknown argument: $1" ;;
  esac
done

[[ -x "${STAR_BIN}" ]] || die "STAR binary not executable: ${STAR_BIN}"
[[ -d "${GEX_DIR}" ]] || die "GEX FASTQ directory missing: ${GEX_DIR}"
[[ -f "${WHITELIST}" ]] || die "Whitelist missing: ${WHITELIST}"
[[ -d "${GENOME_DIR}" ]] || die "Genome directory missing: ${GENOME_DIR}"

R1_FILES=$(find -L "${GEX_DIR}" -maxdepth 1 -type f -name '*R1*.fastq.gz' | sort | paste -sd, -)
R2_FILES=$(find -L "${GEX_DIR}" -maxdepth 1 -type f -name '*R2*.fastq.gz' | sort | paste -sd, -)
[[ -n "${R1_FILES}" ]] || die "No R1 FASTQs under ${GEX_DIR}"
[[ -n "${R2_FILES}" ]] || die "No R2 FASTQs under ${GEX_DIR}"

prepare_outdir() {
  local out="$1"
  if [[ "${RUN_STAR}" -eq 0 ]]; then
    [[ -d "${out}" ]] || die "Missing existing output dir for --skip-run: ${out}"
    return
  fi
  if [[ -e "${out}" ]]; then
    if [[ "${FORCE}" -eq 1 ]]; then
      rm -rf "${out}"
    else
      die "Output dir already exists: ${out} (use --force or set OUT_ROOT to a fresh path)"
    fi
  fi
  mkdir -p "${out}"
  rm -rf "${out}/tmp"
}

build_star_cmd() {
  local out="$1"
  local threads="$2"
  STAR_CMD=(
    "${STAR_BIN}"
    --runThreadN "${threads}"
    --genomeDir "${GENOME_DIR}"
    --readFilesIn "${R2_FILES}" "${R1_FILES}"
    --readFilesCommand zcat
    --outFileNamePrefix "${out}/"
    --outTmpDir "${out}/tmp"
    --clipAdapterType CellRanger4
    --clip3pPolyG yes
    --alignEndsType Local
    --chimSegmentMin 1000000
    --soloType CB_UMI_Simple
    --soloCBstart 1 --soloCBlen 16
    --soloUMIstart 17 --soloUMIlen 12
    --soloBarcodeReadLength 0
    --soloCBwhitelist "${WHITELIST}"
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloMultiMappers Unique
    --soloCellFilter EmptyDrops_CR
    --soloCbUbRequireTogether no
    --soloStrand Unstranded
    --soloFeatures GeneFull
    --soloCrGexFeature genefull
    --soloCrMultimapRescue yes
    --dynamicThreadInterface 1
    --dynamicThreadConstMapPermits "${threads}"
    --dynamicThreadTelemetry 1
    --outSAMtype None
    --soloInlineHashMode yes
  )
  if [[ "${READ_MAP_NUMBER}" != "0" ]]; then
    STAR_CMD+=(--readMapNumber "${READ_MAP_NUMBER}")
  fi
}

write_run_command() {
  local out="$1"
  {
    echo '#!/usr/bin/env bash'
    echo 'set -euo pipefail'
    echo 'export STAR_SOLO_NONFLEX_HASH_BRIDGE=1'
    printf 'export STAR_SOLO_BRIDGE_SOURCE_DETERMINISM_OUT=%q\n' "${out}/bridge_source_determinism.tsv"
    printf 'export STAR_SOLO_BRIDGE_DETERMINISM_OUT=%q\n' "${out}/bridge_determinism.tsv"
    if [[ "${BRIDGE_DETERMINISM_CB}" -eq 1 ]]; then
      echo 'export STAR_SOLO_BRIDGE_DETERMINISM_CB=1'
    fi
    printf '%q ' "${STAR_CMD[@]}"
    echo
  } > "${out}/RUN_COMMAND.sh"
  chmod +x "${out}/RUN_COMMAND.sh"
}

run_one() {
  local label="$1"
  local out="$2"
  local threads="$3"
  prepare_outdir "${out}"
  build_star_cmd "${out}" "${threads}"
  write_run_command "${out}"
  if [[ "${RUN_STAR}" -eq 0 ]]; then
    return
  fi

  log "Starting ${label}: threads=${threads}, out=${out}"
  local start
  start=$(date +%s)
  (
    export STAR_SOLO_NONFLEX_HASH_BRIDGE=1
    export STAR_SOLO_BRIDGE_SOURCE_DETERMINISM_OUT="${out}/bridge_source_determinism.tsv"
    export STAR_SOLO_BRIDGE_DETERMINISM_OUT="${out}/bridge_determinism.tsv"
    if [[ "${BRIDGE_DETERMINISM_CB}" -eq 1 ]]; then
      export STAR_SOLO_BRIDGE_DETERMINISM_CB=1
    else
      unset STAR_SOLO_BRIDGE_DETERMINISM_CB 2>/dev/null || true
    fi
    "${STAR_CMD[@]}"
  ) > "${out}/star.log" 2>&1
  local elapsed=$(( $(date +%s) - start ))
  printf 'label=%s\nthreads=%s\nwall_seconds=%s\noutdir=%s\n' \
    "${label}" "${threads}" "${elapsed}" "${out}" > "${out}/HARNESS_RUN_SUMMARY.txt"
  [[ -f "${out}/bridge_determinism.tsv" ]] || die "Missing bridge determinism digest for ${label}: ${out}"
  log "Finished ${label}: ${elapsed}s"
}

compare_outputs() {
  local a="$1"
  local b="$2"
  local report="${OUT_ROOT}/comparison_summary.tsv"
  local source_diff="${OUT_ROOT}/bridge_source_determinism.diff"
  local digest_diff="${OUT_ROOT}/bridge_determinism.diff"
  local cb_diff="${OUT_ROOT}/bridge_determinism.cb.diff"

  set +e
  if [[ -f "${a}/bridge_source_determinism.tsv" && -f "${b}/bridge_source_determinism.tsv" ]]; then
    diff -u "${a}/bridge_source_determinism.tsv" "${b}/bridge_source_determinism.tsv" > "${source_diff}"
  fi
  diff -u "${a}/bridge_determinism.tsv" "${b}/bridge_determinism.tsv" > "${digest_diff}"
  local digest_status=$?
  if [[ "${BRIDGE_DETERMINISM_CB}" -eq 1 && -f "${a}/bridge_determinism.tsv.cb.tsv" && -f "${b}/bridge_determinism.tsv.cb.tsv" ]]; then
    diff -u "${a}/bridge_determinism.tsv.cb.tsv" "${b}/bridge_determinism.tsv.cb.tsv" > "${cb_diff}"
  fi
  set -e

  python3 - "$a" "$b" "$report" "${BRIDGE_DETERMINISM_CB}" <<'PY'
import gzip
import hashlib
import sys
from pathlib import Path

run_a = Path(sys.argv[1])
run_b = Path(sys.argv[2])
report = Path(sys.argv[3])
per_cb = sys.argv[4] == "1"
source_stages = [
    "thread_hash_pre_resolve",
    "thread_hash_logical_pre_resolve",
    "pending_ambig_umi_gene_pre_resolve",
    "pending_ambig_candidates_pre_resolve",
    "pending_ambig_context_representative_pre_resolve",
    "pending_ambig_context_pin_quals_pre_resolve",
    "pending_ambig_context_evidence_pre_resolve",
    "pending_ambig_context_accounting_counts_pre_resolve",
    "pending_ambig_context_accounting_flags_pre_resolve",
    "cb_read_count_pre_resolve",
    "resolved_ambig_hash_post_resolve",
    "post_accounting_pin_reads",
    "post_accounting_feature_stats",
]
stages = ["pre_cr_exact", "post_cr", "resolved_umi_gene", "final_gene_count"]
mask = (1 << 64) - 1

def open_text(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")

def open_bytes(path):
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return open(path, "rb")

def splitmix64(x):
    x = (x + 0x9e3779b97f4a7c15) & mask
    x = ((x ^ (x >> 30)) * 0xbf58476d1ce4e5b9) & mask
    x = ((x ^ (x >> 27)) * 0x94d049bb133111eb) & mask
    return (x ^ (x >> 31)) & mask

def resolve(mex_dir, name):
    for candidate in (mex_dir / name, mex_dir / f"{name}.gz"):
        if candidate.exists():
            return candidate
    return None

def read_stage_digest(path):
    out = {}
    with open(path, "rt", encoding="utf-8") as handle:
        header = None
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            out[parts[0]] = tuple(int(x) for x in parts[1:])
    return out

def read_cb_digest(path):
    out = {}
    icb = {}
    with open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            key = (parts[idx["stage"]], int(parts[idx["wl_cb"]]), parts[idx["barcode"]])
            out[key] = tuple(int(parts[idx[name]]) for name in ("records", "total_count", "hash_sum", "hash_xor", "hash_sum2"))
            icb[key] = int(parts[idx["iCB"]])
    return out, icb

def file_digest(mex_dir, name):
    path = resolve(mex_dir, name)
    if path is None:
        return ("missing", "", "0")
    sha = hashlib.sha256()
    lines = 0
    with open_bytes(path) as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(chunk)
            lines += chunk.count(b"\n")
    return (str(path.relative_to(mex_dir)), sha.hexdigest(), str(lines))

def matrix_digest(mex_dir):
    path = resolve(mex_dir, "matrix.mtx")
    if path is None:
        return {"status": "missing"}
    rows = cols = nnz_header = None
    records = 0
    value_sum = 0
    hash_sum = 0
    hash_xor = 0
    hash_sum2 = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.split()
            if rows is None:
                rows, cols, nnz_header = map(int, parts[:3])
                continue
            if len(parts) < 3:
                continue
            r, c, v = int(parts[0]), int(parts[1]), int(parts[2])
            records += 1
            value_sum += v
            h = splitmix64(r | (c << 32)) ^ splitmix64(v | 0x517cc1b727220a95)
            h = splitmix64(h)
            hash_sum = (hash_sum + h) & mask
            hash_xor ^= h
            hash_sum2 = (hash_sum2 + ((h * h) & mask) + 0x9e3779b97f4a7c15) & mask
    return {
        "status": "ok",
        "path": str(path.relative_to(mex_dir)),
        "rows": rows,
        "cols": cols,
        "nnz_header": nnz_header,
        "records": records,
        "value_sum": value_sum,
        "hash_sum": hash_sum,
        "hash_xor": hash_xor,
        "hash_sum2": hash_sum2,
    }

def mex_summary(run):
    out = {}
    for layer in ("raw", "filtered"):
        mex = run / "Solo.out" / "GeneFull" / layer
        out[layer] = {
            "matrix": matrix_digest(mex),
            "barcodes": file_digest(mex, "barcodes.tsv"),
            "features": file_digest(mex, "features.tsv"),
        }
    return out

stage_a = read_stage_digest(run_a / "bridge_determinism.tsv")
stage_b = read_stage_digest(run_b / "bridge_determinism.tsv")
first_stage = "none"
for stage in stages:
    if stage_a.get(stage) != stage_b.get(stage):
        first_stage = stage
        break

lines = []
lines.append("section\tname\trun_a\trun_b\tmatch")
source_path_a = run_a / "bridge_source_determinism.tsv"
source_path_b = run_b / "bridge_source_determinism.tsv"
first_source_stage = "missing"
if source_path_a.exists() and source_path_b.exists():
    source_a = read_stage_digest(source_path_a)
    source_b = read_stage_digest(source_path_b)
    first_source_stage = "none"
    for stage in source_stages:
        if source_a.get(stage) != source_b.get(stage):
            first_source_stage = stage
            break
    for stage in source_stages:
        va = source_a.get(stage)
        vb = source_b.get(stage)
        lines.append(f"bridge_source_stage\t{stage}\t{va}\t{vb}\t{va == vb}")
lines.append(f"bridge_source_stage\tfirst_different_source_stage\t{first_source_stage}\t{first_source_stage}\t{first_source_stage == 'none'}")
for stage in stages:
    va = stage_a.get(stage)
    vb = stage_b.get(stage)
    lines.append(f"bridge_stage\t{stage}\t{va}\t{vb}\t{va == vb}")
lines.append(f"bridge_stage\tfirst_different_stage\t{first_stage}\t{first_stage}\t{first_stage == 'none'}")

if per_cb:
    path_a = run_a / "bridge_determinism.tsv.cb.tsv"
    path_b = run_b / "bridge_determinism.tsv.cb.tsv"
    if path_a.exists() and path_b.exists():
        cb_a, icb_a = read_cb_digest(path_a)
        cb_b, icb_b = read_cb_digest(path_b)
        for stage in stages:
            keys = sorted({k for k in cb_a if k[0] == stage} | {k for k in cb_b if k[0] == stage})
            diffs = [k for k in keys if cb_a.get(k) != cb_b.get(k)]
            missing_a = sum(1 for k in diffs if k not in cb_a)
            missing_b = sum(1 for k in diffs if k not in cb_b)
            changed = len(diffs) - missing_a - missing_b
            lines.append(f"bridge_per_cb\t{stage}_diff_count\t{len(diffs)}\t{len(diffs)}\t{len(diffs) == 0}")
            lines.append(f"bridge_per_cb\t{stage}_changed_present_in_both\t{changed}\t{changed}\t{changed == 0}")
            lines.append(f"bridge_per_cb\t{stage}_only_run_a\t{missing_b}\t{missing_b}\t{missing_b == 0}")
            lines.append(f"bridge_per_cb\t{stage}_only_run_b\t{missing_a}\t{missing_a}\t{missing_a == 0}")
            for key in diffs[:20]:
                label = (
                    f"{stage}:wl={key[1]}:bc={key[2]}"
                    f":iCB_a={icb_a.get(key, 'NA')}:iCB_b={icb_b.get(key, 'NA')}"
                )
                lines.append(f"bridge_per_cb\t{label}\t{cb_a.get(key)}\t{cb_b.get(key)}\tFalse")
    else:
        lines.append("bridge_per_cb\tstatus\tmissing_per_cb_files\tmissing_per_cb_files\tFalse")

mex_a = mex_summary(run_a)
mex_b = mex_summary(run_b)
for layer in ("raw", "filtered"):
    ma = mex_a[layer]["matrix"]
    mb = mex_b[layer]["matrix"]
    lines.append(f"mex_matrix\tGeneFull/{layer}\t{ma}\t{mb}\t{ma == mb}")
    for name in ("barcodes", "features"):
        va = mex_a[layer][name]
        vb = mex_b[layer][name]
        lines.append(f"mex_file\tGeneFull/{layer}/{name}\t{va}\t{vb}\t{va == vb}")

report.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"first_different_source_stage={first_source_stage}")
print(f"first_different_bridge_stage={first_stage}")
for line in lines:
    if (line.startswith("mex_matrix\t")
            or line.startswith("bridge_stage\tfirst_different_stage")
            or line.startswith("bridge_source_stage\tfirst_different_source_stage")):
        print(line)
PY

  if [[ -f "${source_diff}" ]]; then
    if [[ -s "${source_diff}" ]]; then
      log "Bridge source digests differ; diff: ${source_diff}"
    else
      log "Bridge source digests match exactly"
    fi
  fi
  if [[ "${digest_status}" -eq 0 ]]; then
    log "Bridge stage digests match exactly"
  else
    log "Bridge stage digests differ; diff: ${digest_diff}"
  fi
  log "Comparison report: ${report}"
}

mkdir -p "${OUT_ROOT}"
RUN_A="${OUT_ROOT}/run_a"
RUN_B="${OUT_ROOT}/run_b"

log "Harness output root: ${OUT_ROOT}"
log "Input: ${GEX_DIR}"
log "Threads: run_a=${THREADS_A}, run_b=${THREADS_B}; read_map_number=${READ_MAP_NUMBER}; per_cb=${BRIDGE_DETERMINISM_CB}"

run_one "run_a" "${RUN_A}" "${THREADS_A}"
run_one "run_b" "${RUN_B}" "${THREADS_B}"
compare_outputs "${RUN_A}" "${RUN_B}"
