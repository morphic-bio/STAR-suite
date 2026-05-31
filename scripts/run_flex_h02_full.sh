#!/usr/bin/env bash
# Full-scale H2 cache build: all H0 probes → H2 variants → STAR-Flex → merged binary.
# Uses --genomeLoad LoadAndKeep to amortize the 23s genome load across all shards.
#
# Env:
#   THREADS        — STAR threads (default 8)
#   STAR_BIN       — STAR binary (default: repo core/legacy/source/STAR)
#   OUT_DIR        — output root (default: /storage/.../flex_h02_full_<stamp>)
#   NUM_SHARDS     — override auto-computed shard count
#   MAX_READS_PER_SHARD — reads per shard (default 600000, must be <= CB whitelist)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
THREADS="${THREADS:-8}"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"

STAMP="$(date -u +%Y%m%d_%H%M%S)"
OUT_DIR="${OUT_DIR:-/storage/downsampled_100K/SC2300771/results/flex_h02_full_${STAMP}}"
GENOME_DIR="${GENOME_DIR:-/storage/flex_filtered_reference/star_index}"
PROBE_LIST="${PROBE_LIST:-/storage/flex_filtered_reference/filtered_reference/probe_list.txt}"
BINARY_CACHE="${BINARY_CACHE:-/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin}"
CB_WHITELIST="${CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${SAMPLE_WHITELIST:-/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv}"
SAMPLE_PROBES="${SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
MAX_READS_PER_SHARD="${MAX_READS_PER_SHARD:-600000}"

# Compute shard count from cache if not overridden
if [[ -z "${NUM_SHARDS:-}" ]]; then
    NUM_SHARDS=$(python3 -c "
from scripts.flex_h01_pilot import load_binary_hash_cache
import math, sys
cache = load_binary_hash_cache('${BINARY_CACHE}')
n_h0 = len(set((r['seq_lo'], r['seq_hi']) for r in cache if r['cache_class'] == 0))
total = n_h0 * 11025
ns = max(1, math.ceil(total / ${MAX_READS_PER_SHARD}))
print(ns)
" 2>/dev/null)
fi

mkdir -p "${OUT_DIR}"

echo "================================================================"
echo "H2 full build: ${OUT_DIR}"
echo "STAR: ${STAR_BIN}  threads=${THREADS}"
echo "Shards: ${NUM_SHARDS}  max_reads/shard: ${MAX_READS_PER_SHARD}"
echo "================================================================"

T0=$(date +%s)

# ── Load genome into shared memory ──────────────────────────────────────────
echo "[$(date +%H:%M:%S)] Loading genome into shared memory..."
"${STAR_BIN}" --genomeDir "${GENOME_DIR}" --genomeLoad LoadAndKeep \
    --outFileNamePrefix "${OUT_DIR}/_genome_load/" 2>/dev/null || true

STAR_COMMON_ARGS=(
    --runThreadN "${THREADS}"
    --genomeDir "${GENOME_DIR}"
    --genomeLoad LoadAndKeep
    --soloType CB_UMI_Simple
    --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0
    --soloCBwhitelist "${CB_WHITELIST}"
    --flex yes
    --soloRunFlexFilter no
    --soloFlexExpectedCellsPerTag 3000
    --soloSampleWhitelist "${SAMPLE_WHITELIST}"
    --soloSampleProbes "${SAMPLE_PROBES}"
    --soloSampleProbeOffset 68
    --soloSampleSearchNearby no
    --soloProbeList "${PROBE_LIST}"
    --limitIObufferSize 50000000 50000000
    --outSJtype None
    --outBAMcompression 0
    --soloMultiMappers Rescue
    --alignIntronMax 500000
    --outFilterMismatchNmax 6
    --outFilterMismatchNoverReadLmax 1.0
    --outFilterMatchNmin 25
    --outSAMunmapped None
    --outFilterMatchNminOverLread 0
    --outFilterMultimapNmax 10000
    --outFilterMultimapScoreRange 4
    --outSAMmultNmax 10000
    --winAnchorMultimapNmax 200
    --outSAMprimaryFlag AllBestScore
    --outFilterScoreMin 0
    --outFilterScoreMinOverLread 0
    --outSAMattributes NH HI AS nM NM GX GN ZG
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    --soloUMIfiltering MultiGeneUMI_CR
    --soloUMIdedup 1MM_CR
    --soloCellFilter None
    --clipAdapterType CellRanger4
    --soloFeatures Gene
    --alignEndsType Local
    --soloStrand Unstranded
    --chimSegmentMin 1000000
    --soloKeysCompat cr
    --outSAMtype BAM Unsorted
)

echo "[$(date +%H:%M:%S)] Genome loaded. Starting ${NUM_SHARDS} shards..."

FAIL_COUNT=0

for (( S=0; S<NUM_SHARDS; S++ )); do
    SHARD_DIR="${OUT_DIR}/shard_$(printf '%04d' $S)"
    mkdir -p "${SHARD_DIR}/synthetic" "${SHARD_DIR}/alignment"

    # Generate synthetic FASTQ for this shard
    python3 "${SCRIPT_DIR}/flex_h01_pilot.py" h2-make-synth-fastq \
        --binary-cache "${BINARY_CACHE}" \
        --probe-list "${PROBE_LIST}" \
        --cb-whitelist "${CB_WHITELIST}" \
        --sample-whitelist "${SAMPLE_WHITELIST}" \
        --probe-offset 0 --sample-offset 68 --read-length 90 --fill-base A \
        --num-shards "${NUM_SHARDS}" \
        --shard-id "${S}" \
        --max-reads-per-shard "${MAX_READS_PER_SHARD}" \
        --r1-fastq "${SHARD_DIR}/synthetic/h2_R1.fastq.gz" \
        --r2-fastq "${SHARD_DIR}/synthetic/h2_R2.fastq.gz" \
        --manifest-tsv "${SHARD_DIR}/synthetic/h2_manifest.tsv" \
        2>/dev/null || { echo "  SKIP shard ${S} (no variants)"; continue; }

    # STAR alignment (genome already in shared memory)
    "${STAR_BIN}" "${STAR_COMMON_ARGS[@]}" \
        --outFileNamePrefix "${SHARD_DIR}/alignment/" \
        --readFilesIn "${SHARD_DIR}/synthetic/h2_R2.fastq.gz" "${SHARD_DIR}/synthetic/h2_R1.fastq.gz" \
        2>/dev/null

    # Build per-shard sequence cache
    python3 "${SCRIPT_DIR}/flex_h01_pilot.py" h2-build-cache-from-mex \
        --manifest-tsv "${SHARD_DIR}/synthetic/h2_manifest.tsv" \
        --probe-list "${PROBE_LIST}" \
        --barcodes-tsv "${SHARD_DIR}/alignment/Solo.out/Gene/raw/barcodes.tsv" \
        --features-tsv "${SHARD_DIR}/alignment/Solo.out/Gene/raw/features.tsv" \
        --matrix-mtx "${SHARD_DIR}/alignment/Solo.out/Gene/raw/matrix.mtx" \
        --qname-cache-tsv "${SHARD_DIR}/h2_qname_cache.tsv" \
        --sequence-cache-tsv "${SHARD_DIR}/h2_sequence_cache.tsv" \
        2>/dev/null

    if (( (S + 1) % 50 == 0 || S + 1 == NUM_SHARDS )); then
        ELAPSED=$(( $(date +%s) - T0 ))
        echo "[$(date +%H:%M:%S)] Shard $((S+1))/${NUM_SHARDS} done  (${ELAPSED}s elapsed)"
    fi

    # Clean up large intermediate files to save disk
    rm -f "${SHARD_DIR}/alignment/Aligned.out.bam"
    rm -f "${SHARD_DIR}/synthetic/h2_R1.fastq.gz" "${SHARD_DIR}/synthetic/h2_R2.fastq.gz"
done

# ── Unload genome from shared memory ───────────────────────────────────────
echo "[$(date +%H:%M:%S)] Unloading genome from shared memory..."
"${STAR_BIN}" --genomeDir "${GENOME_DIR}" --genomeLoad Remove \
    --outFileNamePrefix "${OUT_DIR}/_genome_unload/" 2>/dev/null || true

# ── Merge per-shard sequence caches ────────────────────────────────────────
echo "[$(date +%H:%M:%S)] Merging ${NUM_SHARDS} shard TSVs..."

python3 -c "
import csv, sys, os
from collections import defaultdict

out_dir = '${OUT_DIR}'
num_shards = ${NUM_SHARDS}
NEG_PROBE_AMBIG = 1
CACHE_CLASS_H2_KEEP = 3

by_seq = defaultdict(lambda: {'keep_genes': set(), 'has_deny': False,
                               'source_genes': set(), 'n_recs': 0, 'n_keep': 0, 'n_deny': 0})

for s in range(num_shards):
    tsv = os.path.join(out_dir, f'shard_{s:04d}', 'h2_sequence_cache.tsv')
    if not os.path.exists(tsv):
        continue
    with open(tsv, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq = row['sequence']
            cc = int(row['cache_class'])
            gene = int(row['resolved_gene_idx15'])
            e = by_seq[seq]
            e['source_genes'].update(row['source_gene_ids'].split(','))
            e['n_recs'] += int(row['n_records'])
            e['n_keep'] += int(row['n_keep'])
            e['n_deny'] += int(row['n_deny'])
            if cc == CACHE_CLASS_H2_KEEP and gene > 0:
                e['keep_genes'].add(gene)
            if cc == 2 or int(row['n_deny']) > 0:
                e['has_deny'] = True

rows_out = []
for seq, e in by_seq.items():
    if e['has_deny'] or len(e['keep_genes']) != 1:
        cc = 2; gene = 0; neg = NEG_PROBE_AMBIG
    else:
        cc = CACHE_CLASS_H2_KEEP; gene = next(iter(e['keep_genes'])); neg = 0
    rows_out.append({
        'sequence': seq, 'cache_class': cc, 'negative_code': neg,
        'resolved_gene_idx15': gene,
        'source_gene_ids': ','.join(sorted(e['source_genes'])),
        'n_records': e['n_recs'], 'n_keep': e['n_keep'], 'n_deny': e['n_deny'],
    })

rows_out.sort(key=lambda r: (r['cache_class'], r['resolved_gene_idx15'], r['sequence']))

merged_tsv = os.path.join(out_dir, 'h2_sequence_cache_merged.tsv')
with open(merged_tsv, 'w', newline='') as f:
    w = csv.DictWriter(f, delimiter='\t',
        fieldnames=['sequence','cache_class','negative_code','resolved_gene_idx15',
                    'source_gene_ids','n_records','n_keep','n_deny'])
    w.writeheader()
    w.writerows(rows_out)

n_keep = sum(1 for r in rows_out if r['cache_class'] == CACHE_CLASS_H2_KEEP)
n_deny = sum(1 for r in rows_out if r['cache_class'] == 2)
print(f'Merged: {len(rows_out)} sequences ({n_keep} KEEP, {n_deny} DENY)')
"

# ── Write final binary cache ───────────────────────────────────────────────
python3 "${SCRIPT_DIR}/flex_h01_pilot.py" h2-write-binary-cache \
    --sequence-cache-tsv "${OUT_DIR}/h2_sequence_cache_merged.tsv" \
    --output-bin "${OUT_DIR}/h2_keep_only.bin" \
    --sample-idx 0 \
    --keep-only

ELAPSED=$(( $(date +%s) - T0 ))
echo "================================================================"
echo "Done in ${ELAPSED}s."
echo "Merged TSV: ${OUT_DIR}/h2_sequence_cache_merged.tsv"
echo "H2 binary:  ${OUT_DIR}/h2_keep_only.bin"
echo "================================================================"
