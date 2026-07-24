# Runbook: in-memory Visium HD GEX `1MM_CR` resolution in STAR

Date: 2026-07-24

Branch: `feature/visium-hd-gex-in-memory-1mm-cr-20260724`

Base: STAR Suite `master` at `d956a898108e8db99be4ce0ee09da0ed982e875c`

Design handoff:
`docs/HANDOFF_VISIUM_HD_GEX_IN_MEMORY_STAR_1MM_CR_20260724.md`

## Objective

Implement the production Visium HD 3-prime GEX molecule path inside one STAR
invocation. STAR reads each paired FASTQ record once, maps only R2, decodes raw
R1 on the mapping worker, combines the spatial candidates with the modern
post-rescue GeneFull result, and retains compact evidence through molecule
resolution and matrix materialization.

The production graph is:

```text
paired FASTQ input
  -> raw R1 H0-H2 candidate decode + raw 9-nt UMI
  -> R2 alignment + annotated CR-compatible GeneFull rescue
  -> compact combined evidence
  -> read cliques + frozen spatial posterior
  -> candidate-specific 1MM_CR
  -> GEX-only MultiGeneUMI_CR
  -> strict / soft_expected / hard / gated_hard
  -> 2 / 8 / 16 micrometer matrices
```

No normal production stage requires a BAM, GeneFull sidecar, FIFO, read-name
join, TSV evidence ledger, standalone resolver, or standalone materializer.

## Scope

This implementation supports exactly:

```text
--soloUMIdedup 1MM_CR
--soloUMIfiltering MultiGeneUMI_CR
--soloMultiMappers Unique
```

It produces the four assignment products `strict`, `soft_expected`, `hard`,
and `gated_hard` at 2, 8, and 16 micrometers. Exact-UMI sensitivity is not
computed or stored.

The existing `--soloSpatialFeatureSidecar` and
`--soloSpatialR1FastqTap` interfaces remain temporary diagnostic oracles. They
are not part of the integrated execution graph and are intended for eventual
retirement.

## Integrated mode contract

Working option names are:

```text
--soloSpatialGexIntegrated yes|no
--soloSpatialBarcodeContract DIR
--soloSpatialBc1Oligos FILE
--soloSpatialBc2Oligos FILE
--soloSpatialAssignmentProducts strict,soft_expected,hard,gated_hard|all
--soloSpatialBinSizes 2,8,16
--soloSpatialExpectedReads N
--soloSpatialExpectedCandidates N
--soloSpatialMemoryFraction F
--soloSpatialOverflowPolicy Fail|Spill
```

The mode defaults to `no`. With it absent, STAR allocates no decoder index,
spatial accumulator, overflow buffer, or matrix state and does not enter any
spatial branch.

Enabled mode fails closed unless all of the following hold:

- paired FASTQ input is ordered R2 then raw R1;
- only R2 is mapped;
- `GeneFull` is the sole Solo feature and CR GEX feature;
- CR-compatible genomic multimapper rescue is enabled;
- annotated score-first retained-GTF rescue evidence is selected;
- `--soloMultiMappers Unique` is selected;
- UMI policy is exactly `1MM_CR` plus `MultiGeneUMI_CR`;
- no cell filtering precedes raw spatial matrix materialization;
- two-pass, restart, batch, SLAM, and transcript-VB replay modes are disabled;
- `GX`, `GN`, `UR`, `UB`, `CB`, and `CR` are not requested as computational
  evidence.

`--outSAMtype None` is the production default. A diagnostic BAM without the
forbidden compatibility tags may be allowed only after BAM-on/off matrix
identity is covered by tests.

## Read-lifetime integration

STAR's FASTQ input retains both physical ends in `Read0` and `Qual0` even when
`readNmates=1` excludes raw R1 from alignment. Decode raw R1 at the beginning
of the mapping worker's `ReadAlign::oneReadLoaded()` lifetime, before any
operation can modify its sequence or qualities. Do not decode in the input
callback, because that callback runs under STAR's input mutex and would
serialize decoder CPU.

Keep the decoded result on the same `ReadAlign` object. In
`ReadAlign::outputAlignments()`, after annotated CR rescue and modern GeneFull
annotation, append the combined result to that worker's accumulator. No
ordinal join is necessary because R1 and R2 remain in the same read object.

The optional legacy feature sidecar may be written from the same GeneFull hook.
It must not affect the integrated accumulator or any final result.

## Decoder contract

Factor the accepted direct tiered H0-H2 decoder from companion commit
`5e9af58f9b86e9d3b95f8612a4a77e40f9a3ed86` into a C++11 STAR module.

Use one immutable shared oligo/index object and one scratch object per mapping
worker. For each read retain:

- raw 9-nt UMI packed in two bits per base;
- every coordinate in the best edit-distance envelope;
- BC1 and BC2 edit distances and observed lengths;
- full-barcode start and minimum tier;
- the accepted Phred-derived alignment log likelihood;
- compact audit flags needed for parity counters.

The barcode contract is a complete row-major universe. Represent candidates by
the compact coordinate index, not a string or hash-map key. Validate the
contract and oligo dimensions before mapping.

Accumulate PCR-inclusive exact-H0 counts for both oligo axes while decoding.
Freeze the merged run-level counts after input closes. Candidate prior support
is computed as:

```text
(bc1_h0[col] + 1) * (bc2_h0[row] + 1) - 1
```

with checked overflow. Do not reread R1 or emit/reload an H0 table.

## Compact storage

Use contiguous thread-local pools with fixed-width records. The logical
minimum is:

```cpp
struct SpatialReadEvidence {
    uint32_t geneIndex;
    uint32_t rawUmi;
    uint32_t candidateBegin;
    uint16_t candidateCount;
    uint16_t flags;
};

struct SpatialCandidateEvidence {
    uint32_t coordinateIndex;
    uint32_t auditBits;
    double logSequenceLikelihood;
};
```

The exact ABI may change only if size assertions and memory-model tests are
updated. Do not store read names, feature strings, coordinate strings, member
lists, or an STL node per read/candidate. A compact source ordinal may be kept
where needed for canonical clique ordering and audit hashes.

## Memory model

Integrated mode requires declared upper bounds for reads and candidate rows.
Use checked arithmetic to calculate the peak of every phase, including:

1. decoder indices, per-thread scratch, read evidence, and candidate pool;
2. grouping indices plus read and candidate evidence;
3. clique support plus candidate-specific UMI correction state;
4. corrected support plus flat cross-gene reconciliation state;
5. final molecule and matrix-entry state.

Compare the maximum phase estimate with the minimum applicable host and cgroup
memory limit multiplied by `--soloSpatialMemoryFraction`. Log the complete
model before reading the first FASTQ record. Runtime counters must fail before
exceeding their declared capacities.

Full-slide planning uses the observed upper bounds:

```text
474,131,092 paired reads
529,580,381 candidate rows
386,530,812 unique-gene reads
351,886,421 joined gene-plus-coordinate reads
```

## Lazy overflow spool

The overflow mechanism is a bounded-memory safety path, not the legacy feature
sidecar.

`--soloSpatialOverflowPolicy Fail` stops before exceeding the configured
budget. `Spill` permits lazy compact binary overflow under STAR's temporary
directory. No overflow file is created unless a stage reaches its high-water
mark.

Overflow requirements:

- use fixed-width or length-delimited binary records with magic, schema,
  source revision, contract digests, record counts, and checksums;
- write immutable sorted runs through temporary names and commit them
  atomically;
- sort by the canonical key for the current stage;
- merge in-memory and spilled runs with deterministic k-way grouping;
- support read evidence, clique support, corrected support, and final matrix
  entries without converting to text;
- reject truncated, corrupt, wrong-schema, wrong-reference, or incomplete
  runs;
- remove committed overflow runs after final matrices and summaries are
  atomically committed;
- report triggered stages, records, bytes, merge time, and peak live memory.

Spill partitioning must not affect scientific output. A forced-spill run with a
small test budget must produce the same declared molecule and matrix hashes as
the all-memory run.

During initial validation `Fail` remains the default. `Spill` may become the
production default only after forced-spill parity and a reviewed full-slide
benchmark.

## In-memory molecule phases

### A. Read cliques and posterior

Merge worker pools with deterministic counts and prefix sums. Group by
`(gene, raw UMI)` and schedule that pair—not an entire gene—as the parallel
unit. Sort reads and candidates by canonical numeric keys before intersection
or floating-point accumulation.

Reproduce the accepted global-intersection clique rule. Pairwise chained
overlap must not manufacture a clique. Apply sequence likelihood and the
frozen H0 prior exactly once per read clique.

Write candidate support directly to a flat table keyed by
`(gene, coordinate, raw UMI)`. Keep all four policy values together.

### B. Candidate-specific `1MM_CR`

Group support by `(gene, coordinate)`. Reuse STAR's existing CR UMI correction
kernel. Do not correct per read, across coordinates, or before spatial
posterior resolution. Preserve the accepted soft expected correction and
occupancy semantics separately from integer count products.

### C. GEX-only `MultiGeneUMI_CR`

Group corrected provisional support by `(coordinate, corrected UMI, policy)`.
Use `multi_gene_umi_cr::resolve` and the direct bridge's flat `MgRow` approach.
Preserve corrected-count ties, original-at-corrected-UMI dominance, and all
rejection counters.

### D. Direct materialization

Feed final numeric molecules directly into a reusable matrix writer. Traverse
the final array once and aggregate 2-micrometer coordinates to 8 and 16
micrometers by integer transforms. Preserve the accepted feature and barcode
axis ordering.

Write only final MEX products, axes, summaries, logs, and explicitly requested
diagnostics. Do not write ordinary evidence, clique, support, reconciliation,
or exact-UMI artifacts.

## Implementation sequence and commits

1. Runbook and design review.
2. CLI, compact structures, memory model, and default-off tests.
3. In-process decoder and exact decoder-parity tests.
4. Combined read-lifetime accumulator and accounting tests.
5. Clique/posterior plus shared `1MM_CR` kernel.
6. Flat `MultiGeneUMI_CR`, materialization, and forced-spill support.
7. Regression, determinism, and 100K validation report.

Keep shared-core commits narrow and preserve history with a non-squash merge
only after user review.

## Frozen 100K gate

Use only the declared fixture sources:

```text
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1
/mnt/pikachu/star-spatial/references/refdata-gex-GRCh38-2024-A_STAR-2.7.11a/star
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc1_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc2_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/barcode_contract/
```

Create a fresh output root. Never use a previous candidate, sidecar, joined
evidence, clique, support, resolver, or matrix as computational input. The
accepted annotated bridge may be read only as a parity oracle after the
integrated output is complete.

Accepted counts are:

| Quantity | Value |
|---|---:|
| Paired reads | 100,000 |
| Unique-gene reads | 81,012 |
| Joined eligible reads | 73,892 |
| Read cliques | 68,212 |
| Strict `1MM_CR` molecules | 55,437 |
| Soft `1MM_CR` mass | 68,171.458563 |
| Hard `1MM_CR` molecules | 68,172 |
| Gated-hard `1MM_CR` molecules | 59,895 |

Required 100K executions are:

1. all-memory, sidecar off;
2. all-memory, optional diagnostic sidecar on;
3. forced-spill with a deliberately small spatial memory budget;
4. one-thread and multithread determinism runs.

All four products and all 2/8/16 matrices must match the accepted molecule
sets and axes. Sidecar on/off and all-memory/forced-spill must have identical
declared hashes.

## Regression gates

At minimum run after a mandatory clean rebuild:

```text
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
make -C core/legacy/source test_SpatialFeatureSidecar
make -C core/legacy/source test_SpatialR1FastqTap
make -C core/legacy/source test_MultiGeneUmiCr
make -C core/legacy/source test_CrMultimapRescuePolicy
tests/run_scrna_sidecar_off_golden.sh
tests/run_molecule_first_native_smoke.sh
```

Add focused tests for decoder parity, complete accounting, global clique
intersection, candidate-specific UMI correction, cross-gene ties/dominance,
soft conservation, hierarchy conservation, corrupt overflow rejection,
forced-spill identity, worker utilization, and default-off isolation.

## Stop conditions

Stop on any candidate, likelihood, H0, GeneFull, clique, posterior, UMI,
cross-gene, policy-mass, scale, thread-count, sidecar-on/off, spill/no-spill,
or default-off mismatch. Diagnose the smallest fixture. Do not change the
scientific policy to make an integration mismatch disappear.

Do not launch the 474-million-read full slide and do not merge to `master`
until the user reviews the complete 100K parity, timing, memory, and overflow
report.

