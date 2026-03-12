#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_ucsf_full_compat_forward_rescue_guides.sh"
HELPER="${REPO_ROOT}/scripts/ucsf_parity/render_star_inputs_from_cr_config.py"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

mkdir -p "${tmpdir}/gex/sample1" "${tmpdir}/guides/sample1" "${tmpdir}/out"

make_fastq() {
  local path="$1"
  python3 - <<'PY' "${path}"
import gzip, sys
with gzip.open(sys.argv[1], "wt") as handle:
    handle.write("@r1\nACGT\n+\nFFFF\n")
PY
}

make_fastq "${tmpdir}/gex/sample1/sample1_S1_L001_R1_001.fastq.gz"
make_fastq "${tmpdir}/gex/sample1/sample1_S1_L001_R2_001.fastq.gz"
make_fastq "${tmpdir}/gex/sample1/sample1_S1_L002_R1_001.fastq.gz"
make_fastq "${tmpdir}/gex/sample1/sample1_S1_L002_R2_001.fastq.gz"
make_fastq "${tmpdir}/guides/sample1/sample1_S1_L001_R1_001.fastq.gz"
make_fastq "${tmpdir}/guides/sample1/sample1_S1_L001_R2_001.fastq.gz"

config="${tmpdir}/config.csv"
cat > "${config}" <<EOF
[gene-expression]
reference,/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A
create-bam,false
chemistry,auto

[feature]
reference,${tmpdir}/feature_ref.csv

[libraries]
fastq_id,fastqs,feature_types
sample1,${tmpdir}/gex/sample1,Gene Expression
sample1,${tmpdir}/guides/sample1,CRISPR Guide Capture
EOF

cat > "${tmpdir}/feature_ref.csv" <<EOF
id,name,sequence,feature_type
g1,G1,ACGT,CRISPR Guide Capture
EOF

helper_env="${tmpdir}/helper.env"
python3 "${HELPER}" --config "${config}" --pf-multi-out "${tmpdir}/pf_multi.csv" --emit-env > "${helper_env}"

grep -F "CR_CONFIG=${config}" "${helper_env}" >/dev/null
grep -F "CR_FEATURE_REF=${tmpdir}/feature_ref.csv" "${helper_env}" >/dev/null
grep -F "GEX_R1=${tmpdir}/gex/sample1/sample1_S1_L001_R1_001.fastq.gz,${tmpdir}/gex/sample1/sample1_S1_L002_R1_001.fastq.gz" "${helper_env}" >/dev/null
grep -F "GUIDE_R2=${tmpdir}/guides/sample1/sample1_S1_L001_R2_001.fastq.gz" "${helper_env}" >/dev/null
grep -F "${tmpdir}/gex/sample1,sample1,Gene Expression,Gene Expression" "${tmpdir}/pf_multi.csv" >/dev/null
grep -F "${tmpdir}/guides/sample1,sample1,CRISPR Guide Capture,CRISPR Guide Capture" "${tmpdir}/pf_multi.csv" >/dev/null

STAR_BIN=/bin/true \
  "${RUN_SCRIPT}" \
  --cr-config "${config}" \
  --out-base "${tmpdir}/out" \
  --run-id smoke \
  --dry-run

manifest="${tmpdir}/out/smoke/RUN_MANIFEST.txt"
command="${tmpdir}/out/smoke/RUN_COMMAND.sh"
generated_pf="${tmpdir}/out/smoke/pf_multi_config.from_cr.csv"

[[ -f "${manifest}" ]]
[[ -f "${command}" ]]
[[ -f "${generated_pf}" ]]
grep -F "cr_config=${config}" "${manifest}" >/dev/null
grep -F "cr_gene_expression_chemistry=auto" "${manifest}" >/dev/null
grep -F -- "--pfMultiConfig ${generated_pf}" "${command}" >/dev/null
grep -F -- "${tmpdir}/gex/sample1/sample1_S1_L001_R2_001.fastq.gz" "${command}" >/dev/null
grep -F -- "${tmpdir}/guides/sample1/sample1_S1_L001_R2_001.fastq.gz" "${command}" >/dev/null

echo "PASS"
