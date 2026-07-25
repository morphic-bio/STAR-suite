# Runbook: in-memory Visium HD GEX 1MM-CR inside STAR

Date: 2026-07-24

Branch: `feature/visium-hd-gex-in-memory-1mm-cr-20260724`

Base: STAR Suite `master` at `d956a898108e8db99be4ce0ee09da0ed982e875c`
plus handoff commit `7bfd1df`.

Status: implementation and source-only 100K validation complete. Do not merge
or launch the full slide before the 100K review.

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

Do not change the frozen all-ACGT R1 fast path, posterior, GeneFull rescue,
UMI correction, or assignment policies to resolve a parity discrepancy. The
bounded non-ACGT fallback specified below is an intentional compatibility
extension, not permission to retune the accepted all-ACGT decoder.

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
--soloSpatialSpillHighWaterCandidates N
```

The option is explicit and default-off. Supporting options use the same
`soloSpatial` namespace and provide:

- the complete row-major barcode/coordinate contract;
- BC1 and BC2 full oligos;
- requested products and scales;
- declared maximum read and candidate capacities;
- the allowed fraction of available host/cgroup memory;
- overflow policy `Fail` or `Spill`.

The advanced spill high-water is a per-worker candidate count. Zero derives a
limit from the checked memory budget; a nonzero value is accepted only with
`Spill` and exists for deterministic capacity testing and explicit tuning.

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

Integrated mode requires `--outSAMtype None`. BAM comparisons belong to a
separate non-integrated diagnostic run; `GX`, `GN`, `UR`, `UB`, `CB`, and `CR`
are never accepted as integrated-mode computational evidence.

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

### Language and ambiguous-base contract

The decoder and its STAR-facing interface compile as C++11 with the rest of
the canonical STAR core. Packed hits, spans, candidates, edit metrics, and
scratch caches are performance records, not C++17 data structures. Give
records explicit C++11 constructors where brace initialization would
otherwise require newer aggregate rules; do not add an isolated
`-std=c++17` translation-unit exception.

Spatial-barcode `N` bases are recoverable uncertainty, not an automatic read
rejection. Classify the bounded barcode decode envelope once per raw R1:

- all A/C/G/T selects the unchanged packed H0-H2 fast path;
- one or more `N` bases selects bounded alignment/DP over the original bases
  and qualities;
- any other unsupported base is invalid input for the spatial decoder.

The hot candidate and molecule records are not widened for this rare case. A
small transient mask may record `N` positions and count for dispatch,
accounting, and audit, following the existing packed-key-plus-`nMask`
convention. It must not substitute `N` with A/C/G/T and then call a matching
substitution an exact base: every `N` contributes one unit of uncertain
mismatch evidence. DP retains or rejects coordinates under the same edit
budget, best-envelope, uniqueness, and ambiguity rules as ordinary spatial
barcode evidence.

The raw 9-nt UMI has a deliberately different contract. Any `N` or other
non-ACGT base sets `rawUmiValid=false`. Barcode decoding may still complete
for read accounting, but that read cannot enter clique formation, UMI
correction, or molecule output. No fifth UMI alphabet symbol is stored.

In the frozen source FASTQs, 768 raw 9-nt UMIs contain `N`. The broader
possible spatial-decode envelope contains 12 reads and 14 `N` bases. The
matched retained Space Ranger 100K BAM is a read-only oracle for the overlap:
eight reads have `N` in emitted CR, six have CB and two do not; 768 reads have
`N` in UR and none has UB. These observations are validation expectations
only and are never computational inputs.

Report at least barcode reads with `N`, barcode `N` bases, DP-recovered,
ambiguous, and unassigned reads, plus UMI reads with `N`. These are scalar
counters, not per-candidate fields.

## Compact storage contract

Representative logical records are:

```cpp
struct SpatialReadEvidence {
    uint32_t geneIndex;
    uint32_t rawUmi;
    uint32_t candidateBegin;
    uint32_t sourceOrdinal;
    uint16_t candidateCount;
    uint16_t flags;
};

struct SpatialCandidateEvidence {
    uint32_t coordinateIndex;
    uint32_t auditBits;
    double logSequenceLikelihood;
};
```

The implemented records are 20 and 16 bytes respectively and are protected by
`static_assert`. The declared read and candidate capacities, including source
ordinals, must fit `uint32_t`. Production storage must
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

The limit is the smaller of host `MemAvailable` and remaining active cgroup
memory, multiplied by the configured safe fraction. Pipeline construction
occurs after the genome/transcriptome is resident, so those allocations are
already reflected in the available-memory observation. If `Fail` is selected
and the model does not fit, STAR exits before reading FASTQs.

The correction estimate deliberately includes conservative bounds for the
current flat-vector and soft-occupancy map workspaces: 224 bytes per declared
candidate, 96 bytes per declared read, and 8 MiB per worker. This is an
admission guard, not a promise that later phases have already been converted
to bounded external-memory algorithms.

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

> **Updated 2026-07-25:** the accumulation-only boundary described below is
> retained as implementation history. Production `Spill` now covers the full
> downstream correction, reconciliation, and materialization path. See
> `docs/RUNBOOK_VISIUM_HD_GEX_DOWNSTREAM_SPOOL_20260725.md`.

The historical GeneFull sidecar is not the overflow mechanism. It lacks raw
R1 candidates and cannot resume molecule processing.

Integrated mode instead supports:

```text
--soloSpatialOverflowPolicy Fail|Spill
--soloSpatialMemoryFraction FLOAT
```

`Fail` remains the default. The implemented `Spill` boundary is the combined
read/candidate accumulator only. Each worker has a deterministic candidate
high-water mark; on reaching it, STAR:

1. sorts the completed segment by `(gene, raw UMI, source ordinal)`;
2. writes a compact, versioned binary run below `outTmpDir`;
3. records schema, key, element size/count, source revision, checksums, and a
   completion marker;
4. releases the segment memory;
5. later performs a deterministic k-way merge over memory and completed runs.

The original implementation stopped at this boundary. The 2026-07-25
extension adds versioned coordinate-contribution and matrix-entry runs, with
bounded per-shard correction and reconciliation. The diagnostic sidecar is
still never an implicit overflow fallback.

It never writes TSV, read names, BAM tags, feature strings, or legacy Solo
spool records. A partial, corrupt, wrong-version, wrong-key, or wrong-source
run is fatal. Completed temporary runs are deleted only after final matrices
and summaries commit atomically. The run summary reports run count, bytes,
merge time, selected high-water mark, and peak resident evidence.

No automatic transition to the old sidecar/FIFO/standalone graph is permitted.

## Molecule phases

### A. Read cliques and posterior

Group by `(gene, raw UMI)` in canonical numeric order. The current finalizer
processes these groups deterministically; when this phase is parallelized for
full-scale throughput, the work unit must be one group, never an entire gene.
Within each group, reproduce the accepted global-intersection clique rule and
frozen posterior. Sort reads and candidates by canonical numeric keys before
intersection or floating-point accumulation.

Retain one shared flat clique/candidate table and derive all four policy
products from that identical candidate universe.

### B. Candidate-specific 1MM_CR

Group support by `(gene, coordinate)`. The flat resolver reproduces the
accepted `1MM_CR` count and directional-dominance contract and is tested
against the shared `MultiGeneUmiCr` policy helper; it does not replay the
legacy barcode-owned Solo arrays. Correction never crosses coordinates and is
never performed independently per read.

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
micrometer parents from the shared numeric molecule set. Preserve accepted feature ordering,
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
- one- and two-`N` barcode recovery through DP, including `N` in BC1, BC2,
  and both halves;
- `N`-bearing barcode ambiguity and no-match cases remain ambiguous or
  unassigned under the normal rules;
- an `N`-bearing UMI remains invalid and produces no molecule even when its
  spatial barcode resolves;
- all-ACGT packed-path output is unchanged when the rare fallback is enabled;
- coordinate contract completeness and row-major mapping;
- global intersection versus invalid chained pairwise overlap;
- posterior/prior applied once per clique;
- candidate-specific correction and no cross-coordinate correction;
- `MultiGeneUMI_CR` ties and original-UMI dominance;
- deterministic soft accumulation;
- hierarchy and policy mass conservation;
- empty, one-run, multi-run, truncated, corrupt, and wrong-schema spill cases;
- byte-identical all-memory and forced-spill results;
- sidecar on/off integrated parity.

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
2. integrated forced-spill with an explicit small candidate high-water mark;
3. integrated diagnostic sidecar on;
4. one-thread and configured multithread determinism runs as required.

Compare with the accepted annotated bridge oracle only after the new source-
only run finishes. The pre-N-policy bridge values are retained as an ACGT and
gene-policy control:

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

The N-aware result is not required to reproduce molecule totals from a bridge
that admitted N-bearing raw UMIs. It must preserve the ACGT decoder and gene
policy invariants, retain spatial-N candidates under the ordinary
uniqueness/ambiguity rules, and exclude every N-bearing raw UMI before clique
formation. One-thread, multithread, all-memory, forced-spill, and
sidecar-enabled runs must have identical declared hashes within the same
N-aware policy. Report per-stage wall/CPU time, peak RSS, live-byte high-water
marks, spill bytes, and merge time.

### Completed C++11/N-policy validation

The final audit clean-built STAR SHA-256 is
`9773a1a95e1ec19ce0d7acbce82f316eb7836051e562e94272df2808afa61d7a`.
The 16-thread acceptance runs used the preceding source-identical clean build
(`96bc1eb1529b3a56bdc383e8cd89bbc7c8e5365c582510c20898ecd3b1e2e998`),
and the accepted one-thread run used another source-identical clean build
(`d90c92fcbf2c835408e34c3cb65ef586e90217e56982beef5d9b57c86f9a898d`).
These binary hash differences are build provenance/compilation time, not
source changes. The one-thread output is byte-identical to the 16-thread
outputs; the final audit build passed the focused and default-off regression
gates listed below.

The reproducible runner is
`scripts/run_visium_hd_gex_integrated_100k.sh`. It refuses an existing output
root or an active STAR process, verifies the fixture checksum manifest,
records STAR and index hashes, owns the child process group, and writes
`BENCHMARK_SUMMARY.txt` and `RUN_COMPLETE` only after successful STAR and
integrated completion markers. Examples:

```bash
scripts/run_visium_hd_gex_integrated_100k.sh /fresh/all_memory

SPATIAL_OVERFLOW_POLICY=Spill \
SPATIAL_SPILL_HIGH_WATER_CANDIDATES=1000 \
scripts/run_visium_hd_gex_integrated_100k.sh /fresh/forced_spill

SPATIAL_DIAGNOSTIC_SIDECAR_PREFIX=/fresh/sidecar/gene_evidence \
scripts/run_visium_hd_gex_integrated_100k.sh /fresh/sidecar

STAR_THREADS=1 \
scripts/run_visium_hd_gex_integrated_100k.sh /fresh/one_thread
```

The final source-only roots are:

- `20260724_ovarian_100k_integrated_accept_allmem_1802`;
- `20260724_ovarian_100k_integrated_accept_spill_1804`;
- `20260724_ovarian_100k_integrated_accept_sidecar_1806`;
- `20260724_ovarian_100k_integrated_accept_1t_retry2_1816`.

Their accepted accounting is:

| Quantity | N-aware value |
|---|---:|
| Paired reads | 100,000 |
| Reads with spatial candidates | 89,159 |
| Unique-gene reads | 81,012 |
| Joined eligible reads | 73,326 |
| Candidate rows | 91,832 |
| Exact-H0 reads | 37,953 |
| Read cliques | 67,665 |
| Strict 1MM-CR molecules | 54,996 |
| Soft 1MM-CR mass | 67,636.488767 |
| Hard 1MM-CR molecules | 67,637 |
| Gated-hard 1MM-CR molecules | 59,421 |

The bounded decode envelope contained 12 reads and 14 bases with N: eight
were uniquely decoded, one retained an ambiguous candidate set, and three had
no candidate. The raw UMI contained N in 768 reads; all 768 were invalidated.
No unsupported spatial base was observed. The matched read-only Space Ranger
100K BAM independently contains eight CR-N reads: six have CB and two do not;
it contains 768 UR-N reads and none has UB. STAR retains candidates for the
same six CR-N reads (five unique plus one ambiguous) and has no candidates for
the same two. The four additional STAR DP dispatches contain N in the broader
possible decode envelope outside Space Ranger's emitted CR span.

All 36 policy MEX component files are byte-identical across one and 16
threads, all-memory and forced spill, and sidecar off/on. Their relative-path
hash-manifest SHA-256 is
`581657e47d7b5b32a803d565c8d2c4694d9509e5c909846813748de2cf80b368`.
The `summary.tsv` SHA-256 is
`cf9b0930dd4918a854461428c3abf1ce50c0cc8c6c975e4d91cdc8ee968c21c7`.
STAR elapsed time was 21.85 seconds at one thread and 18.99 seconds for the
primary 16-thread all-memory run, with peak RSS 15,126,120 and 16,823,548 KiB
respectively. Index/fixture hashing is outside the timed STAR child. These
small-fixture timings are determinism checks, not throughput benchmarks.

The forced-spill run completed in fresh root
`20260724_ovarian_100k_integrated_accept_spill_1804` with a test-only high-water mark
of 1,000 candidates per thread. It wrote 92 validated runs totaling 2,944,664
bytes, observed peak resident evidence of 1,581 reads and 1,980 candidates,
and merged in 0.01434 seconds. Its STAR wall time was 19.12 seconds and peak
RSS was 16,820,276 KiB. The completed spill directory and partial files were
removed after matrix/summary commit.

The diagnostic-sidecar run completed in 19.34 seconds with peak RSS
16,827,340 KiB. Its 100,000-record sidecar completed with binary SHA-256
`41a2b2d20521bac76a30934be7489ea0e6713a01dd05c12ad0f3540d0ad136d0`,
while every integrated MEX component remained byte-identical to sidecar-off.
The primary all-memory run emitted no BAM, sidecar, FIFO, normalized evidence,
text join, standalone resolver output, or exact-UMI product.

The final clean C++11 build and focused tests passed:

- `test_SpatialGex` and `test_SpatialGexSpill`;
- `test_SpatialFeatureSidecar` and `test_SpatialR1FastqTap`;
- `test_MultiGeneUmiCr`, including 20,000 legacy/bridge parity cases;
- `test_CrMultimapRescuePolicy` and the R1-tap/integrated-mode fail-closed CLI guards;
- normal scRNA sidecar-off golden MEX parity;
- the standalone molecule-first native smoke after building its local tools.

Two late one-thread launch roots are non-computational diagnostics: the first
exited 127 while a concurrent clean rebuild temporarily removed STAR, and the
second stopped during a concurrent edit that briefly left the runner
incomplete. Neither opened STAR input or contributed validation data.

### Scaling-only 2M run

A balanced scaling fixture was created from the full ovarian source FASTQs at
`/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_2m_scaling_v1`.
It contains the first 500,000 paired records from each of four lanes: exactly
2,000,000 pairs total, with paired identifiers and fixed 43-nt R1/75-nt R2
lengths verified. `checksums.sha256` and `summary.json` describe this fixture.
It is a throughput fixture, not a count oracle.

The serialized all-memory scaling matrix is at
`/mnt/pikachu/star-spatial/gex_in_memory_scaling/20260724_ovarian_2m_clean_d90c92fc`.
Every accepted point used STAR SHA-256
`d90c92fcbf2c835408e34c3cb65ef586e90217e56982beef5d9b57c86f9a898d`,
external `zcat`, expected capacities of 2,000,000 reads and 2,300,000
candidates, and fresh output roots. All wrapper and integrated completion
markers are present. No oracle comparison was performed.

| Threads | STAR wall (s) | Speedup | Reads/s | Mean CPU | Peak RSS (GiB) | STAR mapping M reads/h | Mutex wait (thread-s) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 134.72 | 1.00x | 14,846 | 101% | 14.76 | 60.50 | 0.00 |
| 8 | 36.73 | 3.67x | 54,451 | 374% | 15.53 | 327.27 | 3.10 |
| 16 | 32.98 | 4.08x | 60,643 | 497% | 16.39 | 400.00 | 15.65 |
| 32 | 31.27 | 4.31x | 63,959 | 575% | 18.13 | 450.00 | 69.00 |

Initialization before mapping remained about 15--16 seconds. Eight threads is
the practical knee for this 2M fixture: 16 threads improve end-to-end time by
only 10%, and 32 threads by another 5%. Mean CPU occupancy reaches only about
5.75 cores at 32 requested threads, while accumulated mutex wait rises from
3.10 to 69.00 thread-seconds and aggregate map-chunk CPU time rises from
111.04 to 152.23 thread-seconds. This points to shared pipeline/input
contention and fixed initialization/finalization work rather than exhausted
alignment compute. A larger fixture is still required to measure full-slide
steady-state scaling.

As a non-oracle execution sanity check, all four runs decoded 2,000,000 reads
and reported identical scalar accounting. Of these, 829 barcode envelopes
contained N; 435 were uniquely recovered by bounded DP, 130 remained
ambiguous, and 264 were unassigned. The rare string-keyed non-ACGT memo path
therefore handled only 0.041% of reads and is not a visible scaling driver.
All 36 MEX matrix/axis components are byte-identical across 1, 8, 16, and 32
threads; their relative-path hash-manifest SHA-256 is
`044d4e8760c0b7c00c0d9560c58b5e4555f47e204131b40e063edd045a334091`.
The policy `summary.tsv` SHA-256 is
`c213cb90e679e3ac3281c3a28afab614632469df2dfc27ea31abd1cbc1a06e64`.
The run summaries intentionally differ in thread-dependent memory and
resident-evidence fields.

An optional four-thread add-on did not launch STAR: a concurrent clean rebuild
from another session removed the binary between wrapper preflight and child
exec. Its exit-127 diagnostic is retained as
`failed_threads_4_concurrent_rebuild_1810` and is excluded from the matrix.

### Pre-merge regression audit

The serialized audit root is
`/mnt/pikachu/star-spatial/gex_in_memory_tests/20260724_merge_gate_1820`.
Diagnostic reproductions are retained alongside the accepted results; a
failed exploratory row is not an acceptance result. The final clean
Chromap-enabled build passed and has STAR SHA-256
`1cd667014caf3250ffba8ff7335bfc4a4693b7203677091626723591d1ab2ff5`.

The focused spatial units, spill units, feature-sidecar and R1-tap units,
`MultiGeneUmiCr` 20,000-case parity test, CR multimap policy test, integrated
CLI guards, normal-scRNA sidecar-off golden, and Solo smoke all pass. The full
CR-compatible CRISPR regression initially exposed an unrelated shared-helper
regression: valid legacy Solo rows with zero corrected support had been
rejected before their original-UMI support could participate in the dominance
veto. The helper now preserves that legacy contract, has explicit accept,
veto, and all-zero unit cases, and the full CRISPR regression passes.

The missing derived Flex cache
`sequence_cache_bc004_full16.bin` is no longer a hidden prerequisite.
`scripts/restamp_flex_hash_cache.py` checksum-validates the retained v1 source
and reconstructs the historical cache byte-for-byte (SHA-256
`9adbe006723fbedc58a3d21f1bd1b33ccf8794275fb9a9f402894f989358f51c`)
by sample-stamping H0 while leaving H1/deny rows global. Exact classifier
correctness is now a separate registered gate:
`tests/run_flex_hash_screen_replay_regression.sh` passes 15 synthetic cases and
reproduces all 800,000 retained decisions with zero flat-vs-truth,
tiered-vs-truth, or flat-vs-tiered mismatches. The downstream 100K integration
gate compares keyed MEX outputs under tight bounds rather than incorrectly
conflating aggregation with classifier exactness. On the reconstructed-cache
BAM-path audit, pooled mismatch and total-count deltas are 0.0320% and 0.0237%;
the worst sample is BC004 at 0.1814% and 0.1557%, with barcode-set difference
0.0872% and maximum coordinate delta one.

The external `bqtools` oracle was originally installed under `/tmp`, so it
disappeared on host reboot even though native C++ CBQ support remained built
and tested. `scripts/provision_bqtools_oracle.sh` now installs pinned v0.5.6 in
the persistent user prefix with builder/binary provenance. All seven local CBQ
conversion, reader, ordered-encoder, STAR, STARsolo, process_features, and
Chromap cases pass, and the registered production contract tier is 12/12
runnable and passing with zero skips. The CB/UB quick regression passes for
sorted and unsorted BAM (identical 13,106 cells and 140,925 UMIs), the
synthetic GeneFull GEX smoke passes 11/11, the Flex CR-config 100K smoke passes,
and the corrected-index SLAM fixture passes its registered parity thresholds
(reported NTR correlation 0.998865).

## Full-slide gate

Do not launch the 474,131,092-read full slide until the user reviews the 100K
result. If approved, use a fresh output root and serialize the benchmark. The
run must pass the preflight model, complete without text evidence artifacts,
keep post-read workers utilized, and report mapping, post-read, spill, total
wall time, and peak RSS separately.

The current conservative full-slide correction estimate is 164,679,461,088
bytes (153.37 GiB) at 64 workers. It exceeds the 0.8 memory budget observed on
this host after the index is resident, so the current implementation is
expected to fail closed before reading full-slide FASTQs. Approval to test the
full slide is not approval to weaken that guard. Before full-scale execution,
either reduce and revalidate the later flat/map workspaces or implement
additional compact spill boundaries for correction/reconciliation and matrix
aggregation. The diagnostic feature sidecar is not an allowed workaround.

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

## Interrupted 100K attempt: failure boundary

The first integrated all-memory attempt used fresh output root
`/mnt/pikachu/star-spatial/gex_in_memory_tests/20260724_100k_all_memory_v1`.
It was manually stopped after the runtime was clearly abnormal. This partial
root is diagnostic evidence only and must never be resumed or used as a
computational input.

The run completed argument checking and its 206,080,896-byte spatial memory
estimate. It then selected STAR's private-genome mode, read the 3,216,071,051-
byte `Genome` file, and stalled while reading the 8,738,781,846-byte `SA`
file. The log ends at `Loading SA ...`; it contains no finished-genome,
mapping-start, FASTQ-consumption, decoder, materializer, or STAR fatal-error
record. Therefore this attempt did not exercise the new spatial path.

This was not a change in STAR's index reader. The accepted sidecar run used
the same sparse-D3 index path and private-genome mode, loading `Genome`, `SA`,
and `SAindex` in about four seconds. The index files predate both runs, and
the current branch has no change from the accepted STAR commit in
`Genome_genomeLoad.cpp`, `streamFuns.cpp`, or `fstreamReadBig()`.

The abnormal attempt overlapped broad recursive `rg` and `find` jobs on
`/mnt/pikachu`. During diagnosis, two of those `rg` jobs alone sustained about
150 and 180 MB/s, `md0` was 100% busy with roughly 230--240 ms read latency,
and STAR received only about 6--12 MB/s. The failure was therefore benchmark
orchestration under saturated storage, before the spatial code, rather than a
measured spatial implementation runtime.

The first interactive interrupt also stopped only the foreground control
path; its STAR child remained alive and continued reading until its exact PID
was terminated during diagnosis. Future wrappers must own a process group,
forward `INT`/`TERM` to the complete group, wait for child exit, and write the
completion summary only after that exit. A partial `Log.out` is never a
completion signal.

A separate session subsequently started
`/mnt/pikachu/star-spatial/gex_in_memory_tests/20260724_ovarian_100k_integrated_v2`
while this diagnosis was active. It was also stopped during private index
loading so that it would not compete for the saturated array. That root is
partial diagnostic evidence, not a resumable run.

Before creating the next fresh attempt:

1. confirm no STAR or competing benchmark process is active;
2. reject the launch while broad scans or sustained I/O saturate `md0`;
3. record index file identities and the selected `genomeLoad` mode;
4. use a wrapper with process-group cancellation and an authoritative
   completion marker; and
5. never resume either interrupted output root.

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
