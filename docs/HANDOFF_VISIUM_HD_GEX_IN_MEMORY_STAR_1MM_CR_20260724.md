# Handoff: in-memory Visium HD GEX molecule resolution inside STAR

Date: 2026-07-24

Status: implementation handoff; production design is approved, implementation
and validation are not complete.

Base: fetched `origin/master`
`d956a898108e8db99be4ce0ee09da0ed982e875c`.

Suggested implementation branch:
`feature/visium-hd-gex-in-memory-1mm-cr-20260724`.

## Decision

Implement the complete Visium HD 3-prime GEX production molecule path inside
STAR. STAR already owns the paired read, the R2 alignment and modern GeneFull
evidence, and the production `1MM_CR`/`MultiGeneUMI_CR` algorithms. It must
also decode the paired raw R1 spatial candidates and retain the combined
evidence in compact memory until molecule resolution and matrix
materialization finish.

The production path is:

```text
one STAR paired-read pass
  R1 -> spatial coordinate candidates + raw 9-nt UMI
  R2 -> alignment -> CR rescue -> modern GeneFull evidence
                    |
             combined read evidence
                    |
        spatial read cliques and posterior
                    |
          candidate-specific 1MM_CR
                    |
      GEX-only MultiGeneUMI_CR reconciliation
                    |
     final (gene, corrected UMI, coordinate, policy, weight)
                    |
      strict / soft / hard / gated-hard matrices
                    |
                 2 / 8 / 16 um
```

All hot-path evidence remains in memory after it is read. Do not build a new
sidecar converter, standalone join, feature text sharder, external sort,
disk-hash reconciliation, or provisional-support ledger.

The existing sidecar and standalone tools become optional parity and
diagnostic oracles. They are not the production execution graph.

## Explicit scope

This tranche implements only the production UMI policy:

```text
--soloUMIdedup 1MM_CR
--soloUMIfiltering MultiGeneUMI_CR
```

Materialize these four assignment products from the same candidate universe:

- `strict`;
- `soft_expected`;
- `hard`;
- `gated_hard`.

Materialize each product at 2, 8, and 16 micrometers.

Exact-UMI sensitivity is out of scope. Do not compute it, store it, add a flag
for it, or include it in a performance estimate.

The adjacent-Xenium comparisons remain downstream analysis. They do not
belong in STAR and may run concurrently after the conserved 16-micrometer
matrices are complete.

## Why the previous full-slide graph is rejected

The stopped full-slide attempt used tested components as process boundaries:

```text
STAR binary sidecar
  -> external R1 join TSV
  -> feature text sort
  -> standalone resolver
  -> provisional-support TSV
  -> disk hash reconciliation
  -> standalone materializer
```

Measured costs from
`/mnt/pikachu/star-spatial/gex_sidecar_tests/20260724_ovarian_full_gex_annotated_v1`
were:

| Stage | Wall time | Large artifact |
|---|---:|---:|
| STAR alignment + sidecar | 24:48 | 3.79 GB binary |
| external R1/feature join | 50:07 | 50.99 GB TSV |
| feature text grouping | 53:10 | second 50.99 GB TSV |
| standalone resolver before termination | 2:01:02 | 273.13 GB support TSV plus 52.30 GB clique TSV |
| prior measured disk-hash reconciliation | 1:35:11 | text/hash bridge |
| prior measured 1MM-CR materialization | 1:11:06 | final matrices |

The run was deliberately terminated at 2026-07-24T15:41:53Z and is recorded
as `pipeline failed_exit_143`. It did not begin reconciliation or
materialization. Preserve it as a failed bridge audit; do not reuse or delete
it as part of this implementation.

None of the large text artifacts is a scientific requirement. The molecule
operations are required; their serialization between separate programs is
not.

## Read these before changing code

Read completely:

- `AGENTS.md` and the host's `AGENTS.local.md`;
- `docs/RUNBOOK_VISIUM_HD_GEX_FUSED_READ_PROCESSING_20260724.md`;
- `docs/VALIDATION_VISIUM_HD_GEX_FUSED_READ_PROCESSING_20260724.md`;
- `docs/RUNBOOK_VISIUM_HD_GEX_FEATURE_SIDECAR_20260722.md`;
- `docs/SPATIAL_FEATURE_SIDECAR.md`;
- `docs/HANDOFF_SOLO_BRIDGE_TUPLE_AMBIG_REDESIGN_NEXT_AGENT_20260325.md`;
- `docs/RUNBOOK_SOLO_BRIDGE_AMBIGUOUS_CB_DESIGN_20260325.md`;
- the relevant direct-bridge sections of
  `docs/HANDOFF_SOLO_OPTIMIZATION_20260324.md`.

Also inspect the validated companion R1 decoder and molecule reference at
`visium-hd-processing` commit
`5e9af58f9b86e9d3b95f8612a4a77e40f9a3ed86`:

- `native/hd_r1_anchored_decode.cpp`;
- `native/hd_candidate_preserving_reference.cpp`;
- `tests/test_hd_r1_native_tiered_decode.py`;
- `tests/test_molecule_first_paper_validation.py`.

Do not assume a STAR defect or change the scientific policy before consulting
these sources and their tests.

## Existing code to reuse

Do not reimplement algorithms that are already tested.

### STAR paired-read and GeneFull path

- `core/legacy/source/ReadAlignChunk_processChunks.cpp` owns the validated
  paired-input hook used by the current R1 tap.
- `core/legacy/source/ReadAlign_outputAlignments.cpp` owns the current modern
  post-rescue GeneFull sidecar hook.
- `core/legacy/source/SpatialFeatureSidecar.{h,cpp}` defines the validated
  modern GeneFull evidence record and remains an optional audit sink.
- `core/legacy/source/SpatialR1FastqTap.{h,cpp}` proves that STAR owns the
  correctly paired raw R1 record. The FIFO is not used in production after
  this work.

Validated fused-read implementation commit: `459028a`.

### Existing in-memory Solo collapse

Use the data-structure and lifetime approach already proven in:

- `core/legacy/source/SoloFeature_collapseUMI_fromBridgeHash.cpp`;
- `core/legacy/source/SoloFeature_collapseUMIall.cpp`, especially
  `SoloFeature::umiArrayCorrect_CR`;
- `core/legacy/source/SoloFeature_countCBgeneUMI.cpp`;
- `core/legacy/source/SoloReadFeature_record_base.cpp`;
- `core/legacy/source/MultiGeneUmiCr.{h,cpp}`.

The direct bridge already demonstrates the intended production style:
thread-local aggregation, flat/CSR grouping, bounded per-group UMI arrays,
flat `MgRow` cross-gene reconciliation, explicit buffer release, and no
per-read replay when BAM output is disabled. Generalize that machinery to
spatial candidate coordinates rather than constructing a second framework.

### Spatial molecule semantics

Factor reusable algorithms, not the command-line TSV interface, from:

- `flex/tools/molecule_first_resolver/molecule_first_resolver.cpp`;
- the molecule-first resolver library used by that executable;
- `flex/tools/molecule_first_resolver/molecule_first_materialize.cpp`.

Historical optimization commits `aa52c71`, `f5b34c6`, `5905208`, `ad8a20c`,
and `420dd55` are useful for semantics and regression fixtures. Their disk
sharding, TSV, and external hash interfaces are not the target architecture.

## Required STAR-owned read path

### 1. Decode R1 inside the paired read

Move or factor the stable decoder primitives from
`native/hd_r1_anchored_decode.cpp` into a STAR Suite module with a small
library interface. Do not invoke the standalone executable and do not write a
FIFO.

For every paired record, while STAR still owns both mates:

1. map only R2;
2. decode raw R1 into every best-envelope spatial coordinate candidate;
3. retain the raw 9-nt UMI from R1;
4. retain the current Phred-derived sequence likelihood and edit/frame audit
   fields required by the frozen posterior;
5. obtain modern post-rescue GeneFull evidence from the current R2 policy;
6. append the combined evidence once to a thread-local spatial accumulator.

Prefer combining the R1 decode result with GeneFull evidence in the same
`ReadAlign` lifetime. If the two hooks cannot meet in one object without an
invasive change, use an ordinal-indexed segmented in-memory store. Do not use
a file merely to join evidence already present in the same STAR process.

The compatibility sidecar may still be written when
`--soloSpatialFeatureSidecar` is explicitly requested. The integrated mode
must not require it.

### 2. Use compact, flat storage

Use packed indices and contiguous pools. Do not allocate an STL object or hash
map per read or candidate. A representative logical layout is:

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
    // Preserve the frozen numeric likelihood representation exactly.
    LikelihoodType logSequenceLikelihood;
    uint32_t auditBits;
};
```

This is illustrative, not a mandated ABI. Preserve all fields needed to
reproduce the accepted 100K posterior and audit counters. Do not store read
names or semicolon-separated member lists in production memory. Stable source
ordinals may be retained compactly only where they are required for
deterministic IDs or an explicitly enabled audit.

Candidate coordinates should use the compact coordinate index from the
validated barcode contract. Row/column strings are materialization metadata,
not hot-path keys.

### 3. Build the H0 prior without rereading R1

Accumulate the PCR-inclusive exact-H0 observation counts while decoding R1.
Merge thread-local counts deterministically after input closes, freeze the H0
prior, and apply it during clique/posterior resolution.

The read group being resolved must not train itself beyond the already frozen
run-level rule. Preserve the existing prior semantics and counts exactly. Do
not reopen FASTQs and do not emit/reload `h0_read_prior.tsv` in production.

## In-memory molecule stages

These are logical phases over memory owned by one STAR invocation. They are
not external commands and must not spill ordinary full-slide evidence to
disk.

### Phase A: read cliques and posterior

Group by `(gene, raw UMI)` using a parallel histogram/prefix-sum or the
equivalent flat-array technique already used by the direct bridge. Within each
group, construct the existing global-intersection candidate-overlap read
cliques and frozen posterior.

The parallel work unit is `(gene, raw UMI)`, not a whole gene. This prevents
an abundant gene from becoming the single-threaded multi-minute tail observed
in the rejected full-slide resolver.

Populate candidate support directly into an in-memory flat table keyed by:

```text
(gene, coordinate, raw UMI)
```

Store the four policy values together in one value record. Do not expand one
logical support into separate text rows by policy or UMI mode.

### Phase B: candidate-specific 1MM_CR

Group the support table by `(gene, coordinate)`. For each group, build the
bounded UMI array and call the existing STAR `umiArrayCorrect_CR` semantics.

The order is mandatory:

```text
spatial clique/posterior combination
  -> candidate-specific 1MM_CR
  -> cross-gene reconciliation
```

Do not correct the same raw UMI independently per read. Do not correct across
coordinates. Do not change the count, dominance, or Hamming-one rules.

### Phase C: GEX-only MultiGeneUMI_CR

After candidate-specific correction, group provisional gene support in memory
by:

```text
(coordinate, corrected UMI, policy)
```

Use `multi_gene_umi_cr::resolve` and the flat `MgRow` approach from the direct
bridge. Preserve corrected-count ties, original-at-corrected-UMI dominance,
and rejection rules exactly.

Canonical `--soloMultiMappers Unique` behavior remains:

- a read with one post-rescue GeneFull gene carries that provisional gene;
- a read still assigned to several genes is excluded, though counters may
  record it;
- different uniquely assigned reads that collapse to the same corrected UMI
  are reconciled here.

### Phase D: materialize before releasing final molecules

Once a final gene is selected, a GEX gene is an ordinary `feature_id`. Feed
the final in-memory molecule records directly to the existing optimized Flex
materialization logic.

Materialize all requested policies and all three scales in one traversal of
the final molecule array. Parent coordinates are deterministic integer
transforms of the 2-micrometer coordinate.

Only final MEX/H5 products, feature/barcode axes, compact summaries, logs, and
optional explicitly requested diagnostic binaries go to disk. Production
must not emit:

- normalized evidence TSV;
- feature-sorted evidence TSV;
- read-clique/member TSV;
- provisional-support TSV;
- hash-reconciliation partitions;
- exact-UMI products.

## Memory ownership and concurrency

The full ovarian fixture has 474,131,092 paired reads. The accepted annotated
sidecar run observed 386,530,812 unique-gene reads, 422,945,305 reads with
spatial candidates, 529,580,381 candidate rows, and 351,886,421 joined
gene-plus-coordinate reads.

Before a full run, add a checked memory model using these observed upper
bounds. Flat arrays and reusable buffers are required; general-purpose node
maps for hundreds of millions of elements are not acceptable.

Required lifetime plan:

1. thread-local paired-read evidence is finalized into contiguous global
   arrays;
2. release alignment-only spatial scratch before molecule grouping;
3. reuse or release read-evidence storage after clique support is built;
4. reuse or release raw support after `1MM_CR` correction;
5. reuse or release provisional gene support after `MultiGeneUMI_CR`;
6. release final molecules after matrices and summaries are committed.

Use all configured post-read worker threads. Build deterministic buckets with
parallel counts and prefix sums; never depend on unordered-map traversal
order. Within a group, order evidence by its canonical compact key before any
floating-point accumulation.

If the preflight memory model exceeds the configured safe fraction of host
RAM, fail before reading the full dataset. Do not silently fall back to the
rejected text bridge. Any optional emergency spill design requires separate
review and is outside this handoff.

## Mode and isolation contract

Add one explicit default-off Visium HD GEX integrated mode. Choose final CLI
names consistently with current STAR conventions, but require inputs for:

- the slide barcode/coordinate contract;
- BC1 and BC2 full oligos or their frozen compiled index;
- the frozen R1 decoder configuration;
- the requested assignment products and spatial scales.

The mode must require the current fail-closed scientific configuration:

- paired FASTQ input ordered as R2 then raw R1;
- R2-only mapping;
- `GeneFull` for `soloFeatures` and `soloCrGexFeature`;
- CR-compatible genomic multimapper rescue;
- current annotated score-first retained-GTF evidence policy;
- `--soloMultiMappers Unique`;
- `--soloUMIdedup 1MM_CR`;
- `--soloUMIfiltering MultiGeneUMI_CR`;
- no cell filtering before raw spatial matrix materialization;
- no SAM/BAM carrier requirement.

Do not use or emit legacy `GX`, `GN`, `UR`, `UB`, `CB`, or Space Ranger tags
as computational evidence. Optional BAM output must not alter molecule
results.

With the integrated option absent, ordinary scRNA-seq, Flex, perturb-seq, and
paired-end alignment behavior must remain byte-identical.

## Implementation sequence

Keep commits reviewable and do not merge until the 100K gate is complete.

1. **Design and memory contract**
   - add the final CLI contract, compact structures, capacity arithmetic, and
     fixture-free tests;
   - no production behavior yet.
2. **In-process R1 decoder**
   - factor the validated companion decoder into a STAR module;
   - prove per-read candidate and likelihood parity on synthetic and 100K
     fixtures;
   - retain the old sidecar/FIFO path as the oracle.
3. **Combined in-memory evidence**
   - join R1 candidates and GeneFull inside the paired-read lifetime;
   - add exact accounting and thread-determinism tests.
4. **In-memory clique/posterior and 1MM_CR**
   - factor molecule-first algorithms into reusable library code;
   - use direct-bridge flat grouping and STAR's existing correction routine.
5. **In-memory MultiGeneUMI_CR and materialization**
   - reuse the flat `MgRow` reconciliation;
   - call the materializer directly for four policies and three scales.
6. **100K conformance and performance**
   - compare against the frozen accepted bridge outputs;
   - fix the smallest fixture on any mismatch.
7. **Fresh full-slide benchmark**
   - only after all parity and memory gates pass;
   - use a new run ID and output root.

Do not register this experimental mode in publication recipes until its 100K
parity and full-slide performance are reviewed.

## Frozen 100K validation fixture

Inputs:

```text
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/fastqs
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/checksums.sha256
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/summary.json
/mnt/pikachu/star-spatial/references/refdata-gex-GRCh38-2024-A_STAR-2.7.11a/star
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc1_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc2_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/barcode_contract/
```

Current annotated-policy 100K bridge oracle:

```text
/mnt/pikachu/star-spatial/gex_sidecar_tests/20260724_ovarian_100k_annotated_v2
```

Thread-one repeat:

```text
/mnt/pikachu/star-spatial/gex_sidecar_tests/20260724_ovarian_100k_annotated_thread1_v1
```

The accepted bridge reported:

- 100,000 paired reads;
- 81,012 unique-gene reads;
- 73,892 joined eligible reads;
- 68,212 read cliques;
- 55,437 strict 1MM-CR molecules;
- 68,171.458563 soft 1MM-CR mass;
- 68,172 hard 1MM-CR molecules;
- 59,895 gated-hard 1MM-CR molecules.

The accepted 16-thread and one-thread bridge artifacts were byte-identical.
Use the bridge only as a parity oracle; do not reuse any generated evidence as
input to the integrated STAR run.

## Required tests

### Decoder and evidence

- every frozen 100K R1 candidate set, coordinate, edit/frame field, likelihood,
  and raw UMI matches the native decoder oracle;
- every modern GeneFull status/gene matches the optional binary sidecar from
  the same STAR invocation;
- paired read count, lane order, and names reconcile without reopening FASTQs;
- no-candidate, no-gene, multi-gene, and unmapped accounting sums to the raw
  universe;
- H0 counts and frozen priors match the accepted oracle.

### Molecule semantics

- global candidate intersection is required for a read clique;
- chained pairwise overlap cannot manufacture a clique;
- clique output is invariant to input and thread scheduling;
- sequence likelihood and H0 prior are applied exactly once per read clique;
- candidate-specific `1MM_CR` matches STAR's existing implementation;
- UMI correction never crosses coordinates;
- `MultiGeneUMI_CR` corrected-count ties and original-UMI dominance match the
  direct bridge;
- strict/shared support remains separate from the ambiguous increment;
- soft mass is not renormalized after occupancy deduplication;
- all policy mass and 2/8/16 hierarchy mass reconcile.

### Determinism and isolation

- 1, 8, 16, and the intended full-slide thread count produce identical
  integer molecule/member sets and canonical hashes;
- soft accumulation is deterministic under the declared numeric contract;
- high-count synthetic features divide across `(gene, raw UMI)` work units
  and cannot force a whole feature onto one worker;
- sidecar on/off gives identical integrated matrices;
- optional BAM on/off gives identical integrated matrices;
- existing scRNA direct-bridge, standard Solo, Flex, and sidecar-off golden
  tests remain unchanged.

Run at least:

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

Add focused integrated-mode tests rather than relying only on the legacy
standalone wrapper.

## Acceptance gates

The 100K integrated run is accepted only if:

- it reads each FASTQ record once in STAR;
- no production sidecar, FIFO, join, TSV, sort, disk-hash partition, or
  standalone resolver/materializer is required;
- all four 1MM-CR products match the accepted annotated bridge molecule sets
  and matrices;
- all accounting, posterior, UMI, cross-gene, policy, and scale conservation
  checks pass;
- one-thread and multithread results have identical declared hashes;
- default-off regressions pass after a mandatory clean build;
- peak memory and post-read timing are reported from the integrated process.

The first full-slide integrated run is accepted only if:

- the checked memory model passes before launch;
- no hot-path spill or text evidence artifact appears;
- the run completes on the full 474,131,092-read universe;
- post-read workers remain utilized on high-count features;
- total wall time, post-read wall time, and peak RSS are reported separately;
- all final policy/scale conservation and determinism gates pass;
- the four downstream Xenium comparisons can consume the final 16-micrometer
  matrices without conversion.

Performance objective: eliminate more than 90% of the rejected bridge's
intermediate bytes and at least 4x of its post-mapping bridge wall time. The
agent should profile and optimize the integrated data structures before
accepting a slower result; it must not restore the text bridge to satisfy the
functional gate.

## Mandatory stop conditions

Stop and diagnose the smallest fixture on any:

- paired-read or ordinal mismatch;
- reference, GeneFull feature-axis, strand, or rescue-policy mismatch;
- coordinate outside the barcode contract;
- candidate-set or H0-prior parity failure;
- `1MM_CR` or `MultiGeneUMI_CR` parity failure;
- posterior or raw molecule mass non-conservation;
- 2/8/16 hierarchy mismatch;
- thread-count nondeterminism;
- default-off regression;
- preflight memory-model failure.

Do not change the decoder, posterior, GeneFull rescue, UMI, or assignment
policy to make an integration mismatch disappear.

## Handoff deliverables

At review, report:

- exact source commits and clean-built STAR binary SHA-256;
- final integrated-mode CLI and memory model;
- code reuse/factoring map;
- unit, regression, 100K parity, and determinism results;
- per-stage in-process timing and peak/live memory;
- final 100K molecule and matrix hashes;
- fresh full-slide result ID only if every 100K gate passes;
- all failures and unresolved performance caveats;
- confirmation that exact UMI and all full-scale text bridges are absent.

Do not merge to `master` until the user reviews the 100K parity and performance
results. Never modify or build in the dirty `/mnt/pikachu/STAR-suite`
checkout; use a clean worktree and clean rebuild.
