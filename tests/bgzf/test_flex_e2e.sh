#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
STAR_BIN="${STAR_BIN:-${ROOT_DIR}/core/legacy/source/STAR}"
CASE="${BGZF_E2E_CASE:-all}"
THREADS="${BGZF_E2E_THREADS:-4}"
READ_LIMIT="${BGZF_E2E_READ_LIMIT:-2000}"
WORKDIR="${BGZF_E2E_OUT_ROOT:-/tmp/star_suite_bgzf_flex_e2e}"
TINY_FASTQ_DIR="${TINY_FASTQ_DIR:-/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_fastq}"
TINY_REF_ROOT="${TINY_REF_ROOT:-/home/lhhung/cellranger-9.0.1/external/cellranger_tiny_ref}"

die() { echo "FAIL: $*" >&2; exit 1; }
[[ -x "${STAR_BIN}" ]] || die "STAR binary is absent: ${STAR_BIN}"
[[ -d "${TINY_FASTQ_DIR}" ]] || die "tiny FASTQ fixture is absent: ${TINY_FASTQ_DIR}"
[[ -d "${TINY_REF_ROOT}" ]] || die "tiny reference fixture is absent: ${TINY_REF_ROOT}"

rm -rf "${WORKDIR}"
mkdir -p "${WORKDIR}"
python3 "${ROOT_DIR}/scripts/codespaces/generate_tiny_flex_demo.py" \
    --tiny-fastq-dir "${TINY_FASTQ_DIR}" \
    --tiny-ref-root "${TINY_REF_ROOT}" \
    --outdir "${WORKDIR}/assets_base" \
    --read-limit "${READ_LIMIT}"

gtf="${TINY_REF_ROOT}/genes/genes.gtf"
[[ -f "${gtf}" ]] || gtf="${TINY_REF_ROOT}/genes/genes.gtf.gz"
"${ROOT_DIR}/flex/scripts/build_filtered_reference.sh" \
    --probe-set "${WORKDIR}/assets_base/probe_set.csv" \
    --base-fasta "${TINY_REF_ROOT}/fasta/genome.fa" \
    --base-gtf "${gtf}" \
    --work-dir "${WORKDIR}/refwork" --skip-filter --quiet
"${ROOT_DIR}/flex/scripts/make_filtered_star_index.sh" \
    --filtered-reference "${WORKDIR}/refwork/filtered_reference" \
    --output-dir "${WORKDIR}/star_index" --threads "${THREADS}" \
    --sa-index-n-bases 11 --star-bin "${STAR_BIN}"

base_fastq="${WORKDIR}/assets_base/gex/tiny_flex"
for layout in plain blocked mixed mixed_control; do
    mkdir -p "${WORKDIR}/${layout}_fastq"
done
for mate in R1 R2; do
    src="${base_fastq}/tinyflex_S1_L001_${mate}_001.fastq.gz"
    cp "${src}" "${WORKDIR}/plain_fastq/tinyflex_S1_L001_${mate}_001.fastq.gz"
    "${ROOT_DIR}/tools/make_bgzf_fixture.sh" --block-bytes 1021 "${src}" \
        "${WORKDIR}/blocked_fastq/tinyflex_S1_L001_${mate}_001.fastq.gz"
    cp "${WORKDIR}/blocked_fastq/tinyflex_S1_L001_${mate}_001.fastq.gz" \
        "${WORKDIR}/mixed_fastq/tinyflex_S1_L001_${mate}_001.fastq.gz"
    cp "${src}" "${WORKDIR}/mixed_fastq/tinyflex_S1_L002_${mate}_001.fastq.gz"
    cp "${src}" "${WORKDIR}/mixed_control_fastq/tinyflex_S1_L001_${mate}_001.fastq.gz"
    cp "${src}" "${WORKDIR}/mixed_control_fastq/tinyflex_S1_L002_${mate}_001.fastq.gz"
done

make_config() {
    local layout="$1"
    python3 - "${WORKDIR}/assets_base/config.csv" "${WORKDIR}/${layout}.config.csv" \
        "${base_fastq}" "${WORKDIR}/${layout}_fastq" <<'PY'
import sys
from pathlib import Path
source, dest, old, new = sys.argv[1:]
text = Path(source).read_text(encoding="utf-8")
if old not in text:
    raise SystemExit(f"FASTQ path {old} was absent from {source}")
Path(dest).write_text(text.replace(old, new), encoding="utf-8")
PY
}
for layout in plain blocked mixed mixed_control; do make_config "${layout}"; done

# The BGZF adapter is intentionally wired only into the fully-fused Flex
# consumers. Build a tiny H0/H1 cache and make every comparison run exercise
# that path rather than the standard align-everything path.
fused_inputs="${WORKDIR}/fused_inputs"
mkdir -p "${fused_inputs}" "${WORKDIR}/hash_cache_run"
python3 "${ROOT_DIR}/scripts/flex_compat/render_flex_inputs_from_cr_config.py" \
    --config "${WORKDIR}/plain.config.csv" \
    --sample-whitelist-out "${fused_inputs}/sample_whitelist.tsv" \
    --sample-probes-out "${fused_inputs}/sample_probes.tsv" \
    --probe-list-out "${fused_inputs}/probe_list.txt" \
    --probe-catalog "${WORKDIR}/assets_base/sample_probe_catalog.tsv" \
    --emit-env > "${fused_inputs}/inputs.env"
"${STAR_BIN}" --runMode hashCacheGenerate --runThreadN "${THREADS}" \
    --genomeDir "${WORKDIR}/star_index" \
    --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 10 --soloBarcodeReadLength 0 \
    --soloCBwhitelist "${WORKDIR}/assets_base/whitelist.txt" \
    --flex yes --soloFlexExpectedCellsPerTag 3000 --soloFeatures Gene \
    --soloProbeList "${fused_inputs}/probe_list.txt" \
    --soloSampleWhitelist "${fused_inputs}/sample_whitelist.tsv" \
    --soloSampleProbes "${fused_inputs}/sample_probes.tsv" \
    --soloSampleProbeOffset 68 \
    --hashCacheOutput "${fused_inputs}/hash_cache.bin" \
    --hashCacheTiers H0,H1 --outSAMtype None \
    --outFileNamePrefix "${WORKDIR}/hash_cache_run/"

run_layout() {
    local layout="$1"
    local mode="$2"
    local wrapper="${WORKDIR}/STAR-${mode}"
    python3 - "${wrapper}" "${STAR_BIN}" "${mode}" <<'PY'
import os
import shlex
import sys
from pathlib import Path
wrapper, star, mode = sys.argv[1:]
Path(wrapper).write_text(
    "#!/usr/bin/env bash\nexec " + shlex.quote(star) +
    " --readFilesBgzfMode " + shlex.quote(mode) +
    " --soloHashScreenFile " + shlex.quote(str(Path(wrapper).parent / "fused_inputs/hash_cache.bin")) +
    " --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 1 \"$@\"\n",
    encoding="utf-8")
os.chmod(wrapper, 0o755)
PY
    STAR_BIN="${wrapper}" "${ROOT_DIR}/scripts/run_flex_cr_config.sh" \
        --cr-config "${WORKDIR}/${layout}.config.csv" \
        --genome-dir "${WORKDIR}/star_index" \
        --cb-whitelist "${WORKDIR}/assets_base/whitelist.txt" \
        --solo-cb-start 1 --solo-cb-len 16 --solo-umi-start 17 --solo-umi-len 10 \
        --sample-probe-catalog "${WORKDIR}/assets_base/sample_probe_catalog.tsv" \
        --sample-probe-offset 68 --out-samtype none \
        --out-base "${WORKDIR}/runs" --run-id "${layout}" --threads "${THREADS}"
}

canonical_manifest() {
    local run_dir="$1"
    local output="$2"
    python3 - "${run_dir}" "${output}" <<'PY'
import hashlib
import sys
from pathlib import Path
root = Path(sys.argv[1])
out = Path(sys.argv[2])
selected = []
for base_name in ("Solo.out", "per_sample"):
    base = root / base_name
    if not base.exists():
        continue
    for path in base.rglob("*"):
        if path.is_file() and path.name in {
            "matrix.mtx", "matrix.mtx.gz", "barcodes.tsv", "barcodes.tsv.gz",
            "features.tsv", "features.tsv.gz", "Barcodes.stats", "Summary.csv",
        }:
            selected.append(path)
if not selected:
    raise SystemExit(f"no canonical Flex outputs found under {root}")
lines = []
for path in sorted(selected):
    if path.suffix == ".gz":
        import gzip
        with gzip.open(path, "rb") as handle:
            payload = handle.read()
    else:
        payload = path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    lines.append(f"{path.relative_to(root)}\t{digest}")
out.write_text("\n".join(lines) + "\n", encoding="utf-8")
PY
}

if [[ "${CASE}" == T4 || "${CASE}" == all ]]; then
    run_layout plain off
    run_layout blocked range
    grep -F "BGZF parallel range readers: active" "${WORKDIR}/runs/blocked/Log.out" >/dev/null \
        || die "T4 range-reader activation was not logged"
    canonical_manifest "${WORKDIR}/runs/plain" "${WORKDIR}/plain.manifest"
    canonical_manifest "${WORKDIR}/runs/blocked" "${WORKDIR}/blocked.manifest"
    diff -u "${WORKDIR}/plain.manifest" "${WORKDIR}/blocked.manifest"
fi

if [[ "${CASE}" == T6 || "${CASE}" == all ]]; then
    run_layout mixed_control off
    run_layout mixed auto
    grep -F "BGZF input lane 0: range" "${WORKDIR}/runs/mixed/Log.out" >/dev/null \
        || die "T6 BGZF lane did not use a range reader"
    grep -F "BGZF input lane 1: zlib" "${WORKDIR}/runs/mixed/Log.out" >/dev/null \
        || die "T6 plain-gzip lane did not use the zlib stream"
    canonical_manifest "${WORKDIR}/runs/mixed_control" "${WORKDIR}/mixed_control.manifest"
    canonical_manifest "${WORKDIR}/runs/mixed" "${WORKDIR}/mixed.manifest"
    diff -u "${WORKDIR}/mixed_control.manifest" "${WORKDIR}/mixed.manifest"
fi

echo "PASS: Flex BGZF ${CASE} end-to-end fixture (${WORKDIR})"
