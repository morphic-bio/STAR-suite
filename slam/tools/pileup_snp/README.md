# pileup_snp: External htslib-based SNP mask caller

A fast external tool for calling SNPs from BAM files using htslib pileup, designed to match GEDI/STAR SNP detection semantics for parity testing and eventual integration into STAR.

## Building

```bash
cd tools/pileup_snp
make
```

Requires:
- htslib (vendored in `source/htslib/`)
- Reference FASTA with `.fai` index (run `samtools faidx ref.fa`)
- BAM file with `.bai` index (run `samtools index bam.bam`)

## Usage

```bash
./pileup_snp \
  --bam input.bam \
  --bed candidate_loci.bed \
  --ref reference.fa \
  --output output_mask.bed.gz \
  [--debug-tsv debug.tsv] \
  [--kMode conv|any] \
  [--includeSecondary 0|1] \
  [--includeSupplementary 0|1] \
  [--minMapQ INT] \
  [--minBaseQ INT] \
  [--pErr FLOAT] \
  [--pval FLOAT] \
  [--minTcRatio FLOAT] \
  [--minCov INT] \
  [--minAlt INT]
```

### Parameters

**Inputs:**
- `--bam`: Input BAM file (must be indexed)
- `--bed`: Candidate loci BED file (0-based start, exclusive end)
- `--ref`: Reference FASTA (must be indexed with `samtools faidx`)

**Outputs:**
- `--output`: Output BED file (bgzip-compressed, tabix-indexed)
- `--debug-tsv`: Optional TSV with per-site counts for debugging

**Counting modes:**
- `--kMode`: `conv` (T→C/A→G only) or `any` (any mismatch) [default: conv]
- `--includeSecondary`: Include secondary alignments (FLAG 256) [default: 0]
- `--includeSupplementary`: Include supplementary alignments [default: 0]
- `--minMapQ`: Minimum MAPQ [default: 0]
- `--minBaseQ`: Minimum base quality (Phred) [default: 0]

**Calling params (GEDI/STAR compat):**
- `--pErr`: Error rate for binomial test [default: 0.001]
- `--pval`: P-value threshold [default: 0.001]
- `--minTcRatio`: Minimum conversion ratio (k/n) [default: 0.3]
- `--minCov`: Minimum coverage [default: 6]
- `--minAlt`: Minimum alternative count [default: 1]

## Output Format

**BED output** (bgzip+tabix):
```
#chrom	start	end	ref	alt	n	k	f	pval
chr1	12345	12346	T	C	20	8	0.400000	1.2e-5
```

**Debug TSV** (if `--debug-tsv` provided):
```
chrom	pos0	ref	n	k_any	k_conv	f_any	f_conv	pval_any	pval_conv	called_any	called_conv
chr1	12345	T	20	8	8	0.400000	0.400000	1.2e-5	1.2e-5	1	1
```

## GEDI Compatibility Mode

To match GEDI defaults:
```bash
./pileup_snp \
  --bam wt.bam \
  --bed candidates.bed \
  --ref ref.fa \
  --output gedi_compat.bed.gz \
  --kMode any \
  --includeSecondary 1 \
  --minMapQ 0 \
  --minBaseQ 0 \
  --pErr 0.001 \
  --pval 0.001 \
  --minTcRatio 0.3 \
  --minCov 6 \
  --minAlt 1
```

## STAR Compatibility Mode

To match STAR `--slamSnpMaskCompat gedi` defaults:
```bash
./pileup_snp \
  --bam wt.bam \
  --bed candidates.bed \
  --ref ref.fa \
  --output star_compat.bed.gz \
  --kMode conv \
  --includeSecondary 0 \
  --minMapQ 0 \
  --minBaseQ 0 \
  --pErr 0.001 \
  --pval 0.001 \
  --minTcRatio 0.3 \
  --minCov 6 \
  --minAlt 1
```

## Testing

See `test_parity.sh` for parity tests against GEDI and STAR outputs.
