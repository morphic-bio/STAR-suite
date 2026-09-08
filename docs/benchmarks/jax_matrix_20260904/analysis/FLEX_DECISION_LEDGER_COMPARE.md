# Flex decision-ledger versus BAM audit

`flex_decision_ledger_compare` joins STAR's optional fused-Flex decision
ledger to one or more Cell Ranger BAMs by read name. It reports read-level
assignment agreement, cache route, rejection reason, BAM feature status, and
pre-collapse molecule balance. A STAR BAM is optional because the decision
ledger is authoritative for cache-only and no-align runs.

Build it with:

```bash
make -C docs/benchmarks/jax_matrix_20260904/analysis \
  flex_decision_ledger_compare
```

Run STAR with `STAR_FLEX_DECISION_LEDGER=/path/decision_ledger.tsv`, then use:

```bash
docs/benchmarks/jax_matrix_20260904/analysis/flex_decision_ledger_compare \
  --star-ledger /path/decision_ledger.tsv \
  --cr-bam sample-a=/path/sample-a.bam \
  --cr-bam unassigned=/path/unassigned.bam \
  --probe-set /path/probe-set.csv \
  --probe-list /path/probe-list.txt \
  --tag-map /path/sample-whitelist.tsv \
  --reference /path/genome.fa \
  --out-dir /path/audit
```

Add `--star-bam /path/star.bam` when alignment geometry is part of the audit.
The tool invokes `samtools calmd` read-only and accepts `--samtools` and
`--samtools-threads` overrides. The output directory contains the complete
`read_ledger.tsv` plus category, route, gene, CIGAR, MAPQ, feature-status, and
molecule-balance summaries.

For Bayesian-resolved cell barcodes, the tool reconstructs the molecule key
using the final ledger row's tag token when the resolution event does not
repeat the tag. This is required to avoid treating one biological molecule as
separate tag-0 and tag-N records.
