#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
HELPER="${REPO_ROOT}/scripts/flex_compat/render_flex_inputs_from_cr_config.py"
RUN_SCRIPT="${REPO_ROOT}/scripts/run_flex_cr_config.sh"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT
mkdir -p "${tmpdir}/gex/sample1" "${tmpdir}/out"

make_fastq() {
  local path="$1"
  python3 - <<'PY' "${path}"
import gzip, sys
with gzip.open(sys.argv[1], 'wt') as handle:
    handle.write('@r1\nACGT\n+\nFFFF\n')
PY
}

make_fastq "${tmpdir}/gex/sample1/sample1_S1_L001_R1_001.fastq.gz"
make_fastq "${tmpdir}/gex/sample1/sample1_S1_L001_R2_001.fastq.gz"
make_fastq "${tmpdir}/gex/sample1/sample1_S1_L002_R1_001.fastq.gz"
make_fastq "${tmpdir}/gex/sample1/sample1_S1_L002_R2_001.fastq.gz"

cat > "${tmpdir}/probe_set.csv" <<EOF
#probe_set_file_format=2.0
gene_id,probe_seq,probe_id,included,region
ENSG000001,AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA,ENSG000001|GENE1|a,TRUE,spliced
ENSG000002,CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC,ENSG000002|GENE2|b,TRUE,spliced
ENSG000001,GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG,ENSG000001|GENE1|c,TRUE,unspliced
EOF

cat > "${tmpdir}/sample_probe_catalog.tsv" <<EOF
AAAACCCC	AAAACCCC	BC004
AAAACCCG	AAAACCCC	BC004
TTTTGGGG	TTTTGGGG	BC006
TTTTGGGA	TTTTGGGG	BC006
EOF

config="${tmpdir}/config.csv"
cat > "${config}" <<EOF
[gene-expression]
reference,/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A
probe-set,${tmpdir}/probe_set.csv
create-bam,false

[star]
sample-probe-catalog,${tmpdir}/sample_probe_catalog.tsv
sample-probe-offset,71

[libraries]
fastq_id,fastqs,feature_types
sample1,${tmpdir}/gex/sample1,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
BC004,BC004,Sample 4
BC006,BC006,Sample 6
EOF

python3 "${HELPER}" \
  --config "${config}" \
  --sample-whitelist-out "${tmpdir}/sample_whitelist.tsv" \
  --sample-probes-out "${tmpdir}/sample_probes.tsv" \
  --probe-list-out "${tmpdir}/probe_list.txt" \
  --emit-env > "${tmpdir}/helper.env"

grep -F "FLEX_SAMPLE_WHITELIST=${tmpdir}/sample_whitelist.tsv" "${tmpdir}/helper.env" >/dev/null
grep -F "FLEX_PROBE_LIST=${tmpdir}/probe_list.txt" "${tmpdir}/helper.env" >/dev/null
grep -F "FLEX_SAMPLE_PROBE_CATALOG=${tmpdir}/sample_probe_catalog.tsv" "${tmpdir}/helper.env" >/dev/null
grep -F "FLEX_SAMPLE_PROBE_OFFSET=71" "${tmpdir}/helper.env" >/dev/null
grep -F $'BC004	AAAACCCC' "${tmpdir}/sample_whitelist.tsv" >/dev/null
grep -F $'BC006	TTTTGGGG' "${tmpdir}/sample_whitelist.tsv" >/dev/null
grep -F $'AAAACCCG	AAAACCCC	BC004' "${tmpdir}/sample_probes.tsv" >/dev/null
grep -F 'ENSG000001' "${tmpdir}/probe_list.txt" >/dev/null
grep -F 'ENSG000002' "${tmpdir}/probe_list.txt" >/dev/null

STAR_BIN=/bin/true \
  "${RUN_SCRIPT}" \
  --cr-config "${config}" \
  --out-base "${tmpdir}/out" \
  --run-id smoke_flex \
  --dry-run

manifest="${tmpdir}/out/smoke_flex/RUN_MANIFEST.txt"
command="${tmpdir}/out/smoke_flex/RUN_COMMAND.sh"

grep -F "cr_config=${config}" "${manifest}" >/dev/null
grep -F "cr_gene_expression_probe_set=${tmpdir}/probe_set.csv" "${manifest}" >/dev/null
grep -F "sample_probe_catalog=${tmpdir}/sample_probe_catalog.tsv" "${manifest}" >/dev/null
grep -F "sample_probe_offset=71" "${manifest}" >/dev/null
grep -F -- "--soloSampleWhitelist ${tmpdir}/out/smoke_flex/sample_whitelist.from_cr.tsv" "${command}" >/dev/null
grep -F -- "--soloSampleProbes ${tmpdir}/out/smoke_flex/sample_probes.from_cr.tsv" "${command}" >/dev/null
grep -F -- "--soloProbeList ${tmpdir}/out/smoke_flex/probe_list.from_cr.txt" "${command}" >/dev/null
grep -F -- "--soloSampleProbeOffset 71" "${command}" >/dev/null

echo "PASS"
