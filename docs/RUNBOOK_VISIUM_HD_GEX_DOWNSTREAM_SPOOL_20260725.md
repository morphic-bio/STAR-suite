# Runbook: bounded downstream spool for integrated Visium HD GEX

Date: 2026-07-25
Branch: `feature/visium-hd-gex-downstream-spool-20260725`
Base: `26e6d05a48537761be330fd926b4bef9b68b3e8e`

## Objective

Extend `--soloSpatialOverflowPolicy Spill` from accumulation-only overflow to
an end-to-end bounded-memory implementation of integrated Visium HD GEX. The
spooled path must form read cliques, perform candidate-specific `1MM_CR` UMI
correction, apply `MultiGeneUMI_CR`, and materialize strict, soft-expected,
hard, and gated-hard matrices without retaining slide-sized downstream vectors.

The forced-spool 100K result must be byte-identical to the accepted in-memory
result for all 36 MEX components. The all-memory `Fail` path remains the
reference implementation and must not change when the spool is not selected.

This is an internal binary overflow implementation. It is not the diagnostic
feature sidecar, BAM replay, the historical Solo binary spool, or a resumable
user-facing intermediate format.

## Motivation and measured boundary

The complete ovarian slide declares 474,131,092 paired reads and 529,580,381
candidate rows. At the merged implementation's conservative capacity model:

| Phase | Estimated peak |
|---|---:|
| Read/candidate accumulation | 16.91 GiB |
| Clique formation | 39.98 GiB |
| Candidate-specific correction | 152.99 GiB |
| Multigene reconciliation | 29.91 GiB |
| Matrix materialization | 56.52 GiB |

After the reference index was resident, the host offered an 83.76-GiB
configured spatial budget. The mandatory gate therefore exited 102 before
mapping. The guard was correct: the existing spill merge reconstructs all
cliques and candidates in memory before correction.

The downstream spool should target no more than 24--32 GiB of spatial state so
the assay is practical on 64--128-GiB systems, subject to the index and other
STAR allocations.

## Frozen correctness contracts

Do not change these policies to obtain a smaller memory footprint:

- Gene evidence is captured after modern GeneFull annotation and
  CR-compatible multimap rescue.
- Raw R1 is decoded once and is never aligned or sent through Solo barcode or
  UMI correction.
- Spatial barcode `N` bases use the bounded DP fallback; raw UMIs containing
  `N` remain invalid.
- Read cliques are formed globally by `(gene, raw UMI)` using the current
  deterministic candidate-intersection rule.
- H0 priors are applied once while normalizing each clique posterior.
- UMI correction is candidate-specific and scoped to `(coordinate, gene)`.
- `MultiGeneUMI_CR` is scoped to `(coordinate, corrected UMI)` and retains
  zero-corrected-support rows as original-support veto evidence.
- Soft expected occupancy and original-UMI dominance retain the current
  deterministic floating-point semantics.
- The feature axis is canonical gene-ID order. Spatial axes retain the current
  row-major coordinate contract and 2/8/16-um barcode naming.
- Cell calling remains downstream of integer post-collapse matrices; soft
  expected matrices are never cell-called.

## Streaming design

```text
mapping workers
    -> compact read/candidate runs sorted by (gene, raw UMI, ordinal)
    -> k-way merge and one complete read clique at a time
    -> posterior + strict/hard/gated/soft contribution records
    -> coordinate-partitioned binary spool
    -> bounded sort/group for one coordinate shard
    -> per-(coordinate,gene) UMI correction
    -> per-(coordinate,corrected UMI) MultiGeneUMI_CR
    -> final matrix-entry runs for each product and scale
    -> ordered k-way assembly of atomic MEX outputs
```

Everything after a clique posterior is computed is local to its 2-um
coordinate. Contributions may therefore be partitioned without changing the
policy. Partition by a deterministic hash of the 16-um parent coordinate so
all 2-um children of an 8- or 16-um bin remain in the same shard. This permits
local scale aggregation while distributing tissue-dense regions.

### Phase 1: stream cliques

Refactor clique construction to emit a completed clique to a sink rather than
requiring global `vector<Clique>` and `vector<CliqueCandidate>` ownership.

The existing in-memory sink preserves the `Fail` path. The spool sink records,
for each clique candidate:

- gene and raw UMI;
- coordinate;
- deterministic clique/candidate ordinal;
- read member count;
- posterior;
- strict, hard, and gated-hard selection bits.

Soft correction needs both summed posterior support and the ordered product of
absence probabilities. Retain enough ordering information to reproduce the
current operation order exactly; do not replace products with log-space math
unless both paths are deliberately changed and the accepted MEX hash remains
unchanged.

### Phase 2: coordinate-local correction

For each shard, use bounded sorted runs keyed first by
`(coordinate, gene, raw UMI, clique ordinal, candidate ordinal)`. Aggregate
duplicate raw support and apply the current directional UMI correction
independently for strict, soft, hard, and gated-hard support.

The largest irreducible correction group is one `(coordinate, gene)` raw-UMI
universe. The packed 9-base UMI namespace bounds this far below a slide-sized
vector. Check group size and arithmetic explicitly and fail closed on schema,
capacity, or memory-budget violations.

### Phase 3: coordinate-local multigene reconciliation

Emit corrected records ordered by
`(product, coordinate, corrected UMI, gene)`. Aggregate corrected and original
support per gene, then call the shared `MultiGeneUmiCr` resolver. For soft
expected support, apply the same winner, tie, and original-dominance policy as
the current resolver.

Emit final molecules or already aggregated matrix entries; do not accumulate a
slide-sized `vector<FinalMolecule>`.

### Phase 4: streaming materialization

Keep only the compact feature-used vector and coordinate bitset globally.
Build the canonical axes after clique streaming completes. Produce bounded,
sorted matrix-entry runs per product and scale. Record `nnz` and mass while
merging, write the MatrixMarket body to a temporary file, then atomically
assemble the header and body.

All selected scales must conserve mass. A matrix, summary, or completion marker
is committed only after every contributing run validates successfully.

## Binary-spool safety contract

Every run must contain:

- fixed magic and schema version;
- explicit record type, record size, and sort key;
- source-revision hash and policy version;
- shard, run, and declared record counts;
- payload byte count and checksum;
- completion trailer written before atomic rename.

Readers reject partial, truncated, corrupt, wrong-schema, wrong-record-size,
wrong-key, wrong-policy, wrong-source, over-capacity, or trailing-byte files.
Spools are run-owned under STAR's temporary directory. Successful completion
removes them only after all MEX and summary files commit. Failed roots are
diagnostic evidence and are never accepted as resumable computational input.

## Admission and configuration

`--soloSpatialOverflowPolicy Fail` retains the existing all-memory admission
model and behavior.

`--soloSpatialOverflowPolicy Spill` selects the end-to-end spool path. Its
admission model must include:

- fixed decoder, H0, feature, coordinate-bitset, and merge-cursor state;
- accumulation high-water memory;
- clique-group high-water memory;
- contribution sort buffers and one correction group;
- reconciliation group and matrix-run buffers;
- configured fraction of memory available after the genome is resident;
- a conservative disk requirement for simultaneously live spool generations.

The production implementation uses 256 coordinate shards and admits at most
the checked 8-GiB downstream workspace plus fixed decoder/thread state. Each
shard is independently checked at a conservative 1,024 resident bytes per
contribution before it is loaded. Internal test-only shard/buffer overrides
must remain clearly labeled and deterministic. The production path must refuse
an unsuitable temporary filesystem rather than falling back to a diagnostic
sidecar or weakening capacity declarations.

Before mapping, `Spill` also requires available space below `outTmpDir` for a
conservative all-product/all-scale bound of 256 bytes per declared read, 64
bytes per declared candidate, plus 1 GiB of fixed run/header margin. This bound
covers simultaneously live read, contribution, matrix-body, and MEX files; a
disk-capacity failure is fatal rather than deferred until a partial write.

## Implementation order

1. Add fixed-width contribution/run schemas and corruption tests.
2. Factor clique emission behind in-memory and spool sinks.
3. Implement coordinate partitioning and bounded run sorting.
4. Implement integer and soft correction over grouped contribution streams.
5. Implement multigene reconciliation and final-entry runs.
6. Implement atomic streaming MEX assembly at 2/8/16 um.
7. Replace the accumulation-only spill admission guard with a checked
   end-to-end spool model and disk admission.
8. Add spool metrics to `run_summary.tsv`: generations, shards, runs, bytes,
   sort/merge wall time, largest group, and resident high-water marks.

## Unit and focused tests

At minimum cover:

- empty, singleton, multi-run, many-shard, and forced-small-buffer spools;
- truncated, corrupt, wrong-schema, wrong-key, wrong-record-size,
  wrong-source, incomplete, and trailing-byte runs;
- clique groups split across every input run and mapping thread;
- candidate sets crossing downstream shards;
- strict, hard, and gated selection from the same clique;
- candidate-specific correction with no cross-coordinate correction;
- integer directional UMI ties and chains;
- soft support ordering, ties, absence products, and original dominance;
- zero-corrected-support multigene veto rows;
- 2/8/16-um aggregation across shard and run boundaries;
- all-ACGT and spatial-`N` DP paths plus invalid-UMI accounting;
- one versus multiple mapping threads and multiple shard counts;
- cleanup only after successful atomic output commit.

## Frozen 100K acceptance

Use source files only:

```text
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/fastqs
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/checksums.sha256
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/summary.json
/mnt/pikachu/star-spatial/references/refdata-gex-GRCh38-2024-A_STAR-2.7.11a/star
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc1_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc2_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/barcode_contract
```

Use fresh output roots. Ignore previous STAR-derived computational artifacts;
Space Ranger artifacts remain allowed for sanity checks.

Run:

1. all-memory `Fail` reference;
2. production `Spill` with automatic memory/shard settings;
3. forced-small-buffer/many-run spool;
4. one-thread forced spool;
5. diagnostic sidecar on with spool, verifying matrix parity.

Accepted molecule accounting:

| Quantity | Value |
|---|---:|
| Read cliques | 67,665 |
| Strict molecules | 54,996 |
| Soft expected mass | 67,636.4887669825 |
| Hard molecules | 67,637 |
| Gated-hard molecules | 59,421 |

All 36 MEX components must reproduce:

```text
581657e47d7b5b32a803d565c8d2c4694d9509e5c909846813748de2cf80b368
```

Integer biological accounting must be identical. All MEX components and the
policy `summary.tsv` must be byte-identical. The top-level diagnostic
`soft_expected_mass` may differ by at most `1e-9` solely because coordinate
shards change the order of floating-point addition; this tolerance does not
permit any matrix-entry or policy-summary difference.

## Implementation and 100K acceptance result

The end-to-end implementation passed the frozen acceptance on 2026-07-25.
Fresh output roots were:

- all-memory reference:
  `/mnt/pikachu/star-spatial/gex_in_memory_tests/20260725_ovarian_100k_downstream_spool_fail_v1`;
- production 256-shard spool:
  `/mnt/pikachu/star-spatial/gex_in_memory_tests/20260725_ovarian_100k_downstream_spool_spill_v3`;
- one-thread, 100-candidate high-water run (922 accumulation runs):
  `/mnt/pikachu/star-spatial/gex_in_memory_tests/20260725_ovarian_100k_downstream_spool_1t_manyrun_v1`;
- spool plus diagnostic sidecar:
  `/mnt/pikachu/star-spatial/gex_in_memory_tests/20260725_ovarian_100k_downstream_spool_sidecar_v1`.

All four runs reproduced the accepted 36-component MEX manifest above and the
same policy-summary SHA-256:

```text
cf9b0930dd4918a854461428c3abf1ce50c0cc8c6c975e4d91cdc8ee968c21c7
```

The production spool emitted 84,731 contribution records in 256 runs
(2,738,016 bytes), with 432 records in its largest shard. It emitted 3,072
matrix runs totaling 12,739,744 bytes. No spool or `.partial` output remained
after successful commit. The top-level soft diagnostic was
`67636.488766982337`, versus `67636.488766982526` in the all-memory run; the
absolute `1.89e-10` difference is floating-point addition order only, as shown
by byte-identical matrices and policy summaries.

The focused unit targets and all runnable regression/smoke gates listed below
passed. `run_cbub_regression_test.sh` was an environmental skip because its
configured `/storage/SC2300771_filtered_2M/sample_whitelist.tsv` fixture was
absent; no substitute fixture was found in the configured storage roots. The
aggregate `make -C core/legacy/source test` target is not runnable on the base
revision because its first pre-existing rule references the absent
`core/legacy/test/PackedReadInfo_test.cpp` and
`core/legacy/source/PackedReadInfo.cpp`; the explicit unit targets below are
the authoritative gates for this change.

## Regression and smoke gates

After a mandatory clean rebuild, run at least:

```text
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
make -C core/legacy/source test_SpatialGex
make -C core/legacy/source test_SpatialGexSpill
make -C core/legacy/source test_SpatialGexDownstreamSpool
make -C core/legacy/source test_SpatialFeatureSidecar
make -C core/legacy/source test_SpatialR1FastqTap
make -C core/legacy/source test_MultiGeneUmiCr
make -C core/legacy/source test_CrMultimapRescuePolicy
tests/run_spatial_r1_tap_guard.sh
tests/run_scrna_sidecar_off_golden.sh
tests/run_molecule_first_native_smoke.sh
tests/run_solo_smoke.sh
tests/run_gex_binary_smoke.sh
tests/run_cbub_regression_test.sh
tests/test_cr_compat_crispr_calling.sh
tests/run_slam_parity_smoke.sh
tests/run_flex_hash_screen_replay_regression.sh
```

Retain the exact Flex classifier replay gate. Treat reconstructed historical
cache E2E drift as a separately reported provenance diagnostic, not as a
substitute for the checksum-pinned 800K decision replay.

## Full-slide readiness gate

Do not launch the complete ovarian slide until the forced-spool 100K matrix,
determinism, corruption, memory-model, and regression gates pass. Before the
full slide, preflight both available memory and temporary-disk capacity using
the declared 474,131,092 reads and 529,580,381 candidates.

The full-slide run must use a fresh root, authoritative wrapper completion
markers, and serialized benchmark execution. After successful MEX commit,
rerun the existing adjacent-Xenium generators and comparison suite.

## Stop conditions

Stop on any:

- change to accepted 100K MEX bytes or biological scalar accounting;
- nondeterminism by thread count, shard count, or sort-buffer size;
- cross-coordinate UMI correction or incomplete multigene groups;
- soft mass/order drift;
- unsupported spool corruption or partial-file acceptance;
- memory or disk use exceeding the admitted budget;
- incomplete atomic output or premature spool cleanup;
- normal scRNA, Flex, CR-compatible, CB/UB, SLAM, or sidecar-off regression.

Do not bypass the guard, reduce declared full-slide capacities, commit binary
spools/MEX outputs, or use previous run artifacts as computational input.
