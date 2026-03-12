#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_ucsf_full_compat_forward_rescue_guides.sh"
HELPER="${REPO_ROOT}/scripts/ucsf_parity/render_star_inputs_from_cr_config.py"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

mkdir -p "${tmpdir}/gex/sample1" "${tmpdir}/guides/sample1" "${tmpdir}/custom/sample1" "${tmpdir}/out"

make_fastq() {
  local path="$1"
  python3 - <<'PY' "${path}"
import gzip, sys
with gzip.open(sys.argv[1], 'wt') as handle:
    handle.write('@r1\nACGT\n+\nFFFF\n')
PY
}

for dir in "${tmpdir}/gex/sample1" "${tmpdir}/guides/sample1" "${tmpdir}/custom/sample1"; do
  make_fastq "${dir}/sample1_S1_L001_R1_001.fastq.gz"
  make_fastq "${dir}/sample1_S1_L001_R2_001.fastq.gz"
done

cat > "${tmpdir}/guide_ref.csv" <<EOF
id,name,sequence,feature_type
g1,G1,ACGT,CRISPR Guide Capture
EOF

cat > "${tmpdir}/custom_ref.csv" <<EOF
id,name,sequence,feature_type
c1,C1,TGCA,Custom
EOF

config="${tmpdir}/config.csv"
cat > "${config}" <<EOF
[gene-expression]
reference,/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A
create-bam,false
chemistry,auto

[star]
sample-id,DE_30KO

[feature]
reference,${tmpdir}/guide_ref.csv

[libraries]
fastq_id,fastqs,feature_types,star_chemistry,star_feature_ref,star_library_id
mRNA_DE_30KO,${tmpdir}/gex/sample1,Gene Expression,, ,gex_main
PolyIII_DE_30KO,${tmpdir}/guides/sample1,CRISPR Guide Capture,NXT,,guide_main
LARRY_DE_30KO,${tmpdir}/custom/sample1,Custom,TRU,${tmpdir}/custom_ref.csv,custom_main
EOF

helper_env="${tmpdir}/helper.env"
python3 "${HELPER}" --config "${config}" --pf-multi-out "${tmpdir}/pf_multi.csv" --emit-env > "${helper_env}"

grep -F "SAMPLE_ID=DE_30KO" "${helper_env}" >/dev/null
grep -F "CR_FEATURE_REF=${tmpdir}/guide_ref.csv" "${helper_env}" >/dev/null
grep -F "FEATURE_FASTQ_DIRS=${tmpdir}/guides/sample1,${tmpdir}/custom/sample1" "${helper_env}" >/dev/null
grep -F "GUIDE_FASTQ_DIRS=${tmpdir}/guides/sample1" "${helper_env}" >/dev/null
grep -F "FEATURE_R2=${tmpdir}/guides/sample1/sample1_S1_L001_R2_001.fastq.gz,${tmpdir}/custom/sample1/sample1_S1_L001_R2_001.fastq.gz" "${helper_env}" >/dev/null
grep -F "${tmpdir}/custom/sample1,DE_30KO,Custom,Custom,TRU,${tmpdir}/custom_ref.csv,custom_main" "${tmpdir}/pf_multi.csv" >/dev/null
grep -F "${tmpdir}/guides/sample1,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,,guide_main" "${tmpdir}/pf_multi.csv" >/dev/null

python3 "${HELPER}" --config "${config}" --sample-id OVERRIDE_SAMPLE --pf-multi-out "${tmpdir}/pf_multi.override.csv" --emit-env > "${tmpdir}/helper.override.env"
grep -F "SAMPLE_ID=OVERRIDE_SAMPLE" "${tmpdir}/helper.override.env" >/dev/null
grep -F "${tmpdir}/custom/sample1,OVERRIDE_SAMPLE,Custom,Custom,TRU,${tmpdir}/custom_ref.csv,custom_main" "${tmpdir}/pf_multi.override.csv" >/dev/null

STAR_BIN=/bin/true \
  "${RUN_SCRIPT}" \
  --cr-config "${config}" \
  --out-base "${tmpdir}/out" \
  --run-id smoke_multi \
  --dry-run

manifest="${tmpdir}/out/smoke_multi/RUN_MANIFEST.txt"
command="${tmpdir}/out/smoke_multi/RUN_COMMAND.sh"

grep -F "cr_config=${config}" "${manifest}" >/dev/null
grep -F -- "${tmpdir}/guides/sample1/sample1_S1_L001_R2_001.fastq.gz" "${command}" >/dev/null
grep -F -- "${tmpdir}/custom/sample1/sample1_S1_L001_R2_001.fastq.gz" "${command}" >/dev/null
grep -F -- "--pfMultiConfig ${tmpdir}/out/smoke_multi/pf_multi_config.from_cr.csv" "${command}" >/dev/null

echo "PASS"
