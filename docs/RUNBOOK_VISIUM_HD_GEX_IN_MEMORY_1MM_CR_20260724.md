# Runbook: in-memory Visium HD GEX 1MM-CR inside STAR

Date: 2026-07-24

Branch: `feature/visium-hd-gex-in-memory-1mm-cr-20260724`

Base: STAR Suite `master` at `d956a898108e8db99be4ce0ee09da0ed982e875c`
plus handoff commit `7bfd1df`.

Status: implementation and validation in progress. Do not merge before the
100K parity and performance review.

## Objective

Implement the complete production Visium HD 3-prime GEX molecule path inside
one STAR invocation:

```text
paired FASTQ input, read once
  raw R1 -> all best-envelope spatial candidates + raw 9-nt UMI + H0 counts
  R2     -> alignment -> annotated score-first CR rescue -> GeneFull evidence
                               |
                    combined compact read evidence
                               |
              clique/posterior spatial resolution
                               |
                  candidate-specific 1MM_CR
                               |
                   GEX-only MultiGeneUMI_CR
                               |
       strict / soft_expected / hard / gated_hard
                               |
                       2 / 8 / 16 um MEX
```

The normal execution graph is in memory. It must not require the historical
GeneFull feature sidecar, raw-R1 FIFO, ordinal join, normalized evidence TSV,
feature shard, external text sort, standalone resolver, or standalone
materializer.

## Scientific scope

This tranche implements only:

```text
--soloUMIdedup 1MM_CR
--soloUMIfiltering MultiGeneUMI_CR
--soloMultiMappers Unique
```

It produces all four policies from the same candidate universe:

- `strict`;
- `soft_expected`;
- `hard`;
- `gated_hard`.

Each product is materialized at 2, 8, and 16 micrometers. Exact-UMI products
are deliberately absent. Xenium comparisons remain downstream work.

Do not change the frozen R1 decoder, posterior, GeneFull rescue, UMI
correction, or assignment policies to resolve a parity discrepancy.

## Mode and isolation contract

Working options:

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

The option is explicit and default-off. Supporting options use the same
`soloSpatial` namespace and provide:

- the complete row-major barcode/coordinate contract;
- BC1 and BC2 full oligos;
- requested products and scales;
- declared maximum read and candidate capacities;
- the allowed fraction of available host/cgroup memory;
- overflow policy `Fail` or `Spill`.

The enabled mode fails closed unless the command uses:

- paired FASTQ input ordered as R2 then raw R1;
- R2-only mapping;
- `GeneFull` for `soloFeatures` and `soloCrGexFeature`;
- CR-compatible genomic multimapper rescue with annotated score-first retained-GTF evidence;
- `--soloMultiMappers Unique`;
- exactly `--soloUMIdedup 1MM_CR`;
- exactly `--soloUMIfiltering MultiGeneUMI_CR`;
- no cell filtering before spatial matrix materialization;
- no two-pass, restart, batch, SLAM, or transcript-VB replay mode.

Production examples use `--outSAMtype None`. Optional BAM output may be tested
separately, but `GX`, `GN`, `UR`, `UB`, `CB`, and `CR` are never accepted as
integrated-mode computational evidence.

With the option absent, no decoder index, candidate pool, spatial accumulator,
memory reservation, overflow file, or matrix writer is created. Existing
scRNA-seq, Solo, Flex, perturb-seq, SLAM, and paired-end behavior must remain
unchanged.

## Paired-read ownership

STAR's FASTQ loader retains both physical ends because `readNends == 2`, while
the sidecar-compatible spatial recipe sets `readNmates == 1` so only R2 is
mapped. `ReadAlign::oneRead()` therefore loads raw R1 into `Read0[1]` and
`Qual0[1]` even though alignment consumes only `Read0[0]`.

Decode R1 at the beginning of the `ReadAlign::oneReadLoaded()` lifetime,
before any operation can mutate the raw sequence or qualities. This runs on
mapping workers, outside the serialized FASTQ input mutex. Retain the compact
decode result on that `ReadAlign` object until `outputAlignments()` completes
the modern GeneFull and CR-rescue decision. Append combined evidence once to
the owning worker's spatial accumulator.

Do not decode inside `fastxAppendRecord()`: that callback executes while STAR
owns the input lock and would serialize the decoder.

## Decoder integration

Factor the accepted direct-tiered H0-H2 primitives from
`visium-hd-processing` commit
`5e9af58f9b86e9d3b95f8612a4a77e40f9a3ed86` into a C++11-compatible STAR
module. The module has:

- one immutable run-level oligo/index object;
- one scratch object per mapping worker;
- a packed 9-nt raw UMI;
- every coordinate in the best edit-distance envelope;
- exact tier, edit, observed-length, frame/start, and Phred likelihood fields;
- exact-H0 BC1 and BC2 counters accumulated per worker.

The barcode contract is validated as a complete row-major universe. The hot
coordinate key is its compact integer index, not a barcode sequence or
row/column string.

The run-level H0 arrays are merged deterministically after input closes. The
posterior uses the frozen coordinate prior:

```text
(bc1_h0[col] + 1) * (bc2_h0[row] + 1) - 1
```

FASTQs are not reopened and an H0 TSV is not emitted in production.

## Compact storage contract

Representative logical records are:

```cpp
struct SpatialReadEvidence {
    uint32_t gene;
    uint32_t rawUmi;
    uint32_t candidateBegin;
    uint16_t candidateCount;
    uint16_t flags;
    uint64_t sourceOrdinal;
};

struct SpatialCandidateEvidence {
    uint32_t coordinate;
    uint32_t auditBits;
    double logSequenceLikelihood;
};
```

The exact packing may change after `static_assert` and capacity tests, but the
numeric likelihood representation must remain exact. Production storage must
not contain read-name strings, coordinate strings, semicolon member lists, an
STL node per read, or a hash map per candidate. Source ordinal is retained
only where canonical scheduling or audit accounting requires it.

Thread-local append pools are finalized into contiguous arrays with checked
count and prefix-sum arithmetic. Every capacity multiplication uses overflow-
checked 64-bit or wider arithmetic.

## Memory preflight and ownership

The enabled command declares maximum reads and candidate rows. Preflight
computes the worst live peak for every logical phase:

1. decoder/read/candidate accumulation;
2. clique grouping and policy support;
3. candidate-specific UMI correction;
4. cross-gene reconciliation;
5. final molecule and matrix aggregation.

The limit is the smaller of physical RAM and the active cgroup memory limit,
multiplied by the configured safe fraction after reserving the STAR genome and
runtime allowance. If `Fail` is selected and the model does not fit, STAR
exits before reading FASTQs.

Buffers are released or reused in this order:

1. decoder scratch after the input pass;
2. read/candidate evidence after clique support is complete;
3. raw support after UMI correction;
4. provisional gene support after `MultiGeneUMI_CR`;
5. final molecule entries after matrices and summaries commit.

The full-slide sizing gate uses at least these observed bounds:

- 474,131,092 paired reads;
- 386,530,812 unique-gene reads;
- 422,945,305 candidate-bearing reads;
- 529,580,381 candidate rows;
- 351,886,421 joined gene-plus-coordinate reads.

## Lazy overflow spool

The historical GeneFull sidecar is not the overflow mechanism. It lacks raw
R1 candidates and cannot resume molecule processing.

Integrated mode instead supports:

```text
--soloSpatialOverflowPolicy Fail|Spill
--soloSpatialMemoryFraction FLOAT
```

`Fail` remains the default until forced-spill parity and full-scale review are
complete. With `Spill`, each flat accumulator has a deterministic high-water
mark. On reaching it, STAR:

1. sorts the completed segment by that stage's canonical numeric key;
2. writes a compact, versioned binary run below `outTmpDir`;
3. records schema, key, element size/count, source revision, checksums, and a
   completion marker;
4. releases the segment memory;
5. later performs a deterministic k-way merge over memory and completed runs.

The generic run mechanism may be used at these boundaries:

- combined read evidence keyed by `(gene, raw UMI, source ordinal)`;
- candidate/policy support keyed by `(gene, coordinate, raw UMI)`;
- provisional cross-gene support keyed by
  `(coordinate, corrected UMI, policy, gene)`;
- final matrix entries keyed by `(policy, scale, coordinate, gene)`.

It never writes TSV, read names, BAM tags, feature strings, or legacy Solo
spool records. A partial, corrupt, wrong-version, wrong-key, or wrong-source
run is fatal. Completed temporary runs are deleted only after final matrices
and summaries commit atomically. The run summary reports spill stages, run
counts, bytes, merge time, and peak/live memory.

No automatic transition to the old sidecar/FIFO/standalone graph is permitted.

## Molecule phases

### A. Read cliques and posterior

Group by `(gene, raw UMI)` using deterministic histogram/prefix-sum or radix
partitioning. The parallel work unit is one such group, never an entire gene.
Within each group, reproduce the accepted global-intersection clique rule and
frozen posterior. Sort reads and candidates by canonical numeric keys before
intersection or floating-point accumulation.

Populate one flat support table keyed by `(gene, coordinate, raw UMI)` with
the four policy values held together.

### B. Candidate-specific 1MM_CR

Group support by `(gene, coordinate)`. Use bounded UMI arrays and the same
`umiArrayCorrect_CR` policy kernel used by ordinary Solo. Correction never
crosses coordinates and is never performed independently per read.

Integer policy support follows STAR count and directional-dominance rules.
Soft support preserves the accepted weighted correction and occupancy
contract exactly, including its declared numeric tolerance.

### C. GEX-only MultiGeneUMI_CR

Group corrected provisional support by
`(coordinate, corrected UMI, policy)`. Call
`multi_gene_umi_cr::resolve` using the direct bridge's flat `MgRow` pattern.
Preserve corrected-count ties, original-at-corrected-UMI dominance, rejection
rules, and GEX Unique-read semantics.

### D. Direct materialization

Factor the numeric axis, aggregation, atomic output, and summary logic from
the optimized molecule-first materializer into a library interface. Consume
final numeric molecules directly. Build all selected policies and 2/8/16
micrometer parents in one traversal. Preserve accepted feature ordering,
decimal coordinate ordering, MEX headers, numeric precision, and mass
conservation.

## Legacy sidecar retirement

`--soloSpatialFeatureSidecar` and `--soloSpatialR1FastqTap` remain unchanged
only as temporary diagnostic oracles. They are not enabled by integrated mode,
not used for overflow, and absent from production recipes. Documentation must
label them diagnostic and scheduled for retirement.

An explicitly requested feature sidecar may coexist with integrated mode for
the 100K gene-evidence audit. Sidecar on/off must produce identical integrated
molecule and matrix hashes. The primary accepted integrated run has the
sidecar disabled and creates no FIFO, digest, join, or evidence artifact.

## Implementation sequence and commits

Keep changes reviewable:

1. runbook, CLI contract, compact records, checked memory arithmetic, and
   fixture-free isolation tests;
2. C++11 decoder library and exact synthetic/companion parity tests;
3. `ReadAlign` decode lifetime and thread-local combined accumulator;
4. deterministic clique/posterior and candidate-specific `1MM_CR`;
5. flat `MultiGeneUMI_CR` and direct matrix materialization;
6. lazy overflow runs and forced-spill parity tests;
7. source-only 100K validation and performance report.

Do not merge shared core changes until the user reviews step 7. Preserve the
DAG with `git merge --no-ff` if integration is later approved.

## Fixture-free tests

Add focused tests for:

- integrated mode absent: zero spatial allocations and unchanged normal CLI;
- invalid mixed recipes rejected before input opens;
- checked arithmetic overflow and host/cgroup memory limits;
- decoder candidate, audit, UMI, likelihood, and H0 behavior;
- coordinate contract completeness and row-major mapping;
- global intersection versus invalid chained pairwise overlap;
- posterior/prior applied once per clique;
- candidate-specific correction and no cross-coordinate correction;
- `MultiGeneUMI_CR` ties and original-UMI dominance;
- deterministic soft accumulation;
- hierarchy and policy mass conservation;
- empty, one-run, multi-run, truncated, corrupt, and wrong-schema spill cases;
- byte-identical all-memory and forced-spill results;
- sidecar on/off and optional BAM on/off integrated parity.

Retain at least these existing gates:

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

Run the relevant standard Solo, direct bridge, Flex, SLAM, and binary test
matrix gates before accepting the branch.

## Frozen 100K validation

Primary inputs are the frozen ovarian fixture, declared GRCh38 index, full
oligos, and barcode contract listed in the handoff. Do not use earlier STAR,
decoder, sidecar, resolver, or matrix artifacts as computational input.

Run in fresh output roots:

1. integrated all-memory, sidecar off;
2. integrated forced-spill with a deliberately small spatial budget;
3. integrated diagnostic sidecar on;
4. one-thread and configured multithread determinism runs as required.

Compare with the accepted annotated bridge oracle only after the new source-
only run finishes. Required accepted values are:

| Quantity | Accepted value |
|---|---:|
| Paired reads | 100,000 |
| Unique-gene reads | 81,012 |
| Joined eligible reads | 73,892 |
| Read cliques | 68,212 |
| Strict 1MM-CR molecules | 55,437 |
| Soft 1MM-CR mass | 68,171.458563 |
| Hard 1MM-CR molecules | 68,172 |
| Gated-hard 1MM-CR molecules | 59,895 |

All four final molecule sets and all 12 matrices must match the accepted
annotated policy. One-thread, multithread, all-memory, forced-spill, and
sidecar-enabled runs must have identical declared hashes. Report per-stage
wall/CPU time, peak RSS, live-byte high-water marks, spill bytes, and merge
time.

## Full-slide gate

Do not launch the 474,131,092-read full slide until the user reviews the 100K
result. If approved, use a fresh output root and serialize the benchmark. The
run must pass the preflight model, complete without text evidence artifacts,
keep post-read workers utilized, and report mapping, post-read, spill, total
wall time, and peak RSS separately.

The performance objective remains greater than 90% elimination of rejected
bridge intermediate bytes and at least 4x faster post-mapping processing. A
spill run is an emergency capacity path, not permission to recreate the text
bridge.

## Mandatory stop conditions

Stop on any:

- paired-read, ordinal, lane, or raw-UMI mismatch;
- decoder candidate, audit, likelihood, or H0 mismatch;
- coordinate outside the complete contract;
- GeneFull feature-axis, strand, or rescue mismatch;
- clique, posterior, `1MM_CR`, or `MultiGeneUMI_CR` mismatch;
- raw, policy, or hierarchy mass non-conservation;
- thread, memory/spill, sidecar, or BAM nondeterminism;
- default-off regression;
- checked memory-model failure;
- incomplete or corrupt overflow run.

Diagnose the smallest fixture. Do not relax the policy, tolerance, or oracle to
make an integration mismatch disappear.

## Review deliverables

Report:

- source commits and clean-built STAR SHA-256;
- final CLI and checked memory formula;
- decoder/collapse/materializer factoring map;
- unit, regression, isolation, spill, and determinism results;
- 100K molecule, matrix, and summary hashes;
- per-stage timing and peak/live memory;
- spill trigger, bytes, runs, and merge time;
- confirmation that exact-UMI and text bridge outputs are absent;
- unresolved correctness or performance caveats.

Do not commit binaries or generated outputs. Record new artifact roots in
`tests/ARTIFACTS.md`.
