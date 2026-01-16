# EC Filter Test Suite

This directory contains tests for the EC (Equivalence Class) Filter Module, which replicates Salmon's alignment filtering and EC construction pipeline.

## Overview

The EC Filter Module implements:
- **Alignment Filtering**: Salmon-compatible filtering with minAlnProb gating
- **EC Construction**: Equivalence class building with range factorization
- **Extended Pruning**: Optional non-Salmon pruning modes (disabled by default for parity)

## Test Types

### 1. Unit Tests

Test individual components with synthetic fixtures.

**Run:**
```bash
cd tools/ec_filter_test
make test
python3 generate_fixtures.py > fixtures.json
# Process fixtures with ec_filter_cli (after CLI implementation)
python3 verify_filtering.py results.json
```

### 2. Parity Tests

Compare EC output against Salmon's `--dumpEqWeights` output to verify >=99% parity.

#### Synthetic Data Parity Test

Uses synthetic reads aligned with minimap2.

**Run:**
```bash
cd tools/ec_filter_test
./generate_parity_bam.sh  # Generate BAM file
./run_salmon_parity.sh /tmp/ec_parity_test
cat /tmp/ec_parity_test/parity_report.txt
```

**Data Location:**
- Fixtures: `test/fixtures/salmon_eq/`
- Transcriptome: `test/fixtures/salmon_eq/chr22_trans.fa`
- Reads: `test/fixtures/salmon_eq/synthetic_reads{1,2}.fq`
- BAM: `test/fixtures/salmon_eq/alignments.bam` (generated)

#### nf-core/test-datasets Parity Test

Uses real RNA-seq data from nf-core/test-datasets (rnaseq branch).

**Prerequisites:**
- nf-core/test-datasets repository cloned at `/mnt/pikachu/test-datasets`
- STAR binary available
- Salmon binary available
- samtools available

**Setup (one-time):**
```bash
# Clone nf-core/test-datasets if not already present
git clone https://github.com/nf-core/test-datasets.git /mnt/pikachu/test-datasets

# The script will automatically checkout the rnaseq branch worktree
```

**Run:**
```bash
cd tools/ec_filter_test
./run_nfcore_parity.sh [output_dir] [threads]
```

**Example:**
```bash
./run_nfcore_parity.sh /tmp/nfcore_parity_test 8
```

**Data Paths:**
- Test datasets root: `/mnt/pikachu/test-datasets`
- rnaseq branch worktree: `/mnt/pikachu/test-datasets-rnaseq`
- FASTQ files: `data/rna-seq/srr6357070/SRR6357070_{1,2}.fastq.gz`
- Reference transcriptome: `reference/transcriptome.fasta.gz`
- Reference GTF: `reference/genes.gtf.gz`

**What the script does:**
1. Checks out nf-core/test-datasets rnaseq branch (creates worktree if needed)
2. Builds STAR index from transcriptome + GTF
3. Aligns reads with STAR in transcriptome mode (`--quantMode TranscriptomeSAM`)
4. Runs Salmon quant in alignment mode with `--dumpEqWeights`
5. Runs our EC filter CLI on the same BAM + transcriptome
6. Compares EC files using `compare_ecs.py`
7. Generates parity report

**Output:**
- STAR index: `{output_dir}/star_idx_nfcore/`
- STAR alignment: `{output_dir}/star_align/Aligned.toTranscriptome.out.bam`
- Salmon output: `{output_dir}/salmon_out/`
- Our EC output: `{output_dir}/our_out/eq_classes.txt`
- Parity report: `{output_dir}/parity_report.txt`

**Expected Result:** >=99% parity on EC labels, weights, and counts.

## Build

### Build Library
```bash
cd source/libem
make clean
make -j4
```

### Build CLI
```bash
cd tools/ec_filter_test
make
```

## Parameters

All tests use Salmon-compatible defaults:

- `score_exp = 1.0`
- `min_aln_prob = 1e-5`
- `min_score_fraction = 0.65`
- `range_factorization_bins = 4`
- `hard_filter = false` (soft filter default)

**Important**: Extended pruning (local/global) is **disabled** by default for parity testing. Use `--no-local-pruning --no-global-pruning` flags.

## Comparison Tool

`compare_ecs.py` compares EC files and reports:
- Matching labels, weights, and counts
- Label-only mismatches
- Weight mismatches (with tolerance)
- Count mismatches
- Parity percentage

**Usage:**
```bash
python3 compare_ecs.py \
    --salmon salmon_eq_classes.txt \
    --ours our_eq_classes.txt \
    --tolerance 1e-6 \
    --report parity_report.txt
```

## CLI Status

**Current Status**: CLI skeleton exists (`ec_filter_cli.cpp`) but requires implementation of:
- BAM file reading (htslib/samtools integration)
- Transcriptome loading (FASTA parser)
- Alignment-to-RawAlignment conversion
- EC output in Salmon format

See `.cursor/plans/ec_filter_handoff.md` for implementation details.

## Troubleshooting

### nf-core/test-datasets not found
```bash
git clone https://github.com/nf-core/test-datasets.git /mnt/pikachu/test-datasets
```

### FASTQ files not found
The script will list available FASTQ files if the expected ones are missing. You can modify the script to use a different pair.

### STAR index build fails
- Ensure transcriptome and GTF files are valid
- Check disk space
- Verify STAR binary is working: `STAR --version`

### Salmon not found
Set `SALMON` environment variable:
```bash
export SALMON=/path/to/salmon
./run_nfcore_parity.sh
```

### CLI not implemented
If `ec_filter_cli` is still a skeleton, the script will skip comparison and exit gracefully. This is expected until CLI integration is complete.

## References

- Salmon source: `SalmonMappingUtils.hpp`, `SalmonQuantify.cpp`
- Plan document: `.cursor/plans/ec_gating_cleanup_module_489701da.plan.md`
- Handoff document: `.cursor/plans/ec_filter_handoff.md`
