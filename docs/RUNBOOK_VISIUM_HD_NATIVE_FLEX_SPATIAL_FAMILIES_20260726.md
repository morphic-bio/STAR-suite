# Runbook: native Visium HD Flex with spatial barcode families

Date: 2026-07-26

Status: implemented and release-validated on the feature branch; pending local
v1.6.0 release integration

Branch: `feature/visium-hd-native-flex-spatial-families-20260726`

Base: STAR Suite v1.6.0
`ddc75a47de372504cd0ae1f887a6b7eee73f3803`

Governing production handoff:

```text
/mnt/pikachu/visium-hd-processing-recipes/docs/
  HANDOFF_NATIVE_FLEX_CRC_CODEX_20260726.md
```

## Objective

Add an explicit native Visium HD Flex mode in which STAR processes R1 and R2
as one paired-read transaction:

```text
R1 -> raw UMI + spatial barcode candidate family
R2 -> Flex H0/H1 or alignment-fallback feature
                         |
                         v
       feature + raw UMI + spatial candidate family
                         |
                         v
        delayed candidate-specific UMI resolution
                         |
                         v
       strict / soft expected / hard / gated hard
```

The new mode is Flex with a different barcode resolver. It is not a new
expression pipeline:

- ordinary Flex resolves R1 to a single corrected cell barcode;
- Visium HD Flex preserves the bounded family of valid spatial coordinates;
- feature assignment remains the existing Flex H0/H1 plus alignment-fallback
  policy;
- UMI correction remains delayed until the candidate-specific molecule phase;
- the accepted v1.6.0 spatial clique, posterior, bounded-spool, and
  materialization engine remains the single downstream implementation.

No BAM, GX/UR tags, text ledger, external sort, companion resolver, or
companion materializer is part of the native route.

## Audit conclusion

The required implementation is narrow. The following capabilities already
exist at the v1.6.0 base:

1. `ReadAlign::oneRead()` reads R2 and raw R1 together and assigns one global
   input ordinal.
2. `FlexHashScreenCache::classifyReadH0H1Offset0()` implements the frozen
   no-sample CRC H0/H1 decision.
3. `flexResolveGeneIdx15_inlineResolver()` and
   `resolveGeneFromCandidates()` implement the normal Flex probe/genomic
   alignment fallback and return one 1-based probe-list feature index.
4. `SpatialR1Decoder` returns the raw nine-base UMI and a bounded spatial
   coordinate candidate family, including the non-ACGT barcode DP fallback.
5. `spatial_gex::Pipeline` already performs:
   - globally candidate-consistent read-clique formation;
   - spatial posterior normalization;
   - candidate-specific `1MM_CR`;
   - strict, soft expected-count, hard, and gated-hard policies;
   - 2, 8, and 16 micrometer matrix materialization;
   - deterministic all-memory and bounded-spool execution.

The ordinary Flex inline hash cannot be used as the spatial molecule store.
Its key contains one corrected CB index and therefore collapses the spatial
candidate family too early. Native spatial Flex must branch immediately before
that single-CB insertion and send the same Flex feature result to the existing
spatial family pipeline.

The existing internal `spatial_gex` name does not imply GEX-specific molecule
semantics. Renaming the mature engine is outside this tranche. The new
user-facing output directory will be `SpatialFlex.out`.

## Implementation record

The implementation remained inside the audited boundary:

- `--soloSpatialFlexIntegrated yes` is default-off and mutually exclusive with
  integrated spatial GEX.
- The ordinary paired reader decodes raw R1 into a candidate family and uses
  the same global source ordinal for the terminal R2 feature decision.
- H0 and H1 cache keeps append the zero-based probe feature directly; denies
  consume the pair as feature-unassigned; misses continue to alignment.
- Alignment misses invoke `flexResolveGeneIdx15_inlineResolver()`, the same
  probe/genomic candidate builder and modern resolver used by ordinary Flex.
- Native spatial mode bypasses `SoloRead::record()` and the ordinary
  single-cell-barcode insertion path.
- `ProbeListIndex` now preserves the ordered probe feature axis used for final
  MEX features.
- The existing spatial engine performs delayed candidate-specific `1MM_CR`,
  clique formation, all four assignment policies, all three bin sizes, and
  bounded spilling.
- Run summaries identify the native Flex schema, feature-axis and cache hashes,
  route counters, feature/spatial cross-tabs, barcode-`N` fallback outcomes,
  and invalid-UMI accounting.
- Explicit `--soloRunFlexFilter no` and `--soloMapqMode off` values are no
  longer overwritten by the Flex omnibus defaults. Omitted values retain
  ordinary Flex production defaults.

No BAM, GX/UR bridge, per-read text ledger, external resolver, or external
materializer participates in the native route.

## Frozen semantic boundary

### Early read-pair completion

Each input ordinal has exactly one R1 decode and exactly one terminal R2
feature outcome:

- H0 keep;
- H1 keep;
- hash deny;
- hash miss followed by alignment-resolved;
- hash miss followed by alignment-unresolved.

Completing the read-pair evidence record does not complete a molecule.

### Delayed molecule completion

The adapter must retain:

- the zero-based Flex probe feature index;
- the raw packed UMI;
- every valid spatial coordinate candidate;
- the per-candidate sequence likelihood;
- the original input ordinal.

It must not produce a corrected UMI, choose a hard coordinate, or collapse
reads during mapping. Candidate-specific UMI resolution and all hard/soft
policies run only in the existing downstream molecule phase.

### Barcode and UMI base handling

- A non-ACGT spatial barcode window enters the existing bounded DP fallback.
  `N` bases are uncertain mismatch evidence, not an automatic rejection.
- A raw UMI containing `N` or another unsupported base remains invalid and
  cannot create a molecule.
- Both outcomes remain visible in accounting.

## Implementation surface

### 1. Explicit opt-in mode

Add:

```text
--soloSpatialFlexIntegrated yes|no
```

Default: `no`.

The mode must be mutually exclusive with:

```text
--soloSpatialGexIntegrated yes
```

The first implementation is intentionally limited to the validated ordinary
paired reader:

```text
--flex yes
--flexPipeline no
--soloInlineHashMode yes
--soloSkipProcessing yes
--soloRunFlexFilter no
--soloFeatures Gene
--outSAMtype None
```

The fused Flex packet pipeline is not required by the frozen CRC method.
Supporting it later requires carrying the spatial candidate family through its
reader/triage/alignment queues and will be a separate extension.

The mode must fail before reading FASTQs unless all of the following are true:

- `--runMode alignReads`;
- two paired FASTQ ends ordered as R2 then raw R1;
- `--flex yes`;
- `--flexPipeline no`;
- a readable H0/H1 cache;
- a readable probe gene list;
- no sample whitelist or sample-tag requirement;
- a complete spatial barcode contract and both oligo files;
- `--soloUMIdedup 1MM_CR`;
- no BAM output and no GX/GN/UR/UB/CB/CR output attributes;
- no two-pass, restart, batch, SLAM, or TranscriptVB mode;
- positive declared read and candidate capacities;
- a valid all-memory or spill memory admission result.

Ordinary Flex, ordinary STARsolo, integrated spatial GEX, and bulk STAR remain
unchanged when the flag is absent.

### 2. Coupled R1/R2 state

Extend the current mapping-thread spatial state so it records the source
ordinal at R1 decode time. Require the same ordinal when R2 reaches its
terminal feature outcome.

The state needs two explicit terminal operations:

```text
append assigned feature and consume current R1
record unassigned feature and consume current R1
```

Starting a new read while an earlier R1 remains unconsumed, completing an
ordinal twice, or completing a different ordinal is a fatal consistency error.

This closes a current accounting gap in which an unassigned GEX read is
implicitly overwritten by the next decode. The change must preserve existing
GEX matrices byte-for-byte.

### 3. Flex H0/H1 route

For native spatial Flex, run the no-sample classifier directly on raw R2
before alignment:

- `Keep`, class 0: record H0 and append `geneIdx15 - 1`;
- `Keep`, class 1: record H1 and append `geneIdx15 - 1`;
- `Deny`: record hash-denied and consume the R1 state without a feature;
- `Pass`: record hash miss and continue to STAR alignment.

Do not call the ordinary single-barcode Flex hash insertion for this mode.
A valid feature decision remains valid even when R1 has no spatial candidates
or has an invalid UMI; those dimensions are accounted independently.

### 4. Flex alignment fallback

After STAR has populated the alignment MAPQ, score, mismatch, CIGAR, chromosome
name, and Gene annotation surfaces, invoke the same candidate builder and
`resolveGeneFromCandidates()` policy used by ordinary Flex.

The shared resolver must remain one implementation. If necessary, extract a
small barcode-independent helper from
`flex/SoloReadFeature_record_flex.cpp`; the existing ordinary Flex wrapper and
the new spatial call site must both call it.

- resolved `geneIdx15 > 0`: append `geneIdx15 - 1`;
- unresolved `geneIdx15 == 0`: record alignment-unresolved and consume R1
  without a feature.

Do not call ordinary `SoloRead::record()` for the native spatial product. That
call performs the single-CB insertion that this mode replaces.

### 5. Probe feature axis

Extend `ProbeListIndex` to retain the ordered gene IDs in addition to its
gene-ID-to-1-based-index map.

Use that same loaded object for:

- H0/H1 index interpretation;
- probe chromosome lookup;
- transcriptome-gene-to-probe mapping;
- alignment-fallback validation;
- final `features.tsv` ordering.

Every accepted `geneIdx15` must be within the ordered axis. The conversion to
the spatial engine is exactly:

```text
spatial_feature_index = geneIdx15 - 1
```

The CRC v2 probe axis has 18,085 entries. Record the path and SHA-256 of the
probe list and cache in the STAR log and run summary. The present cache schema
does not contain a probe-axis digest, so the command preflight must continue
to pin both files independently.

### 6. Accounting

Add a native Flex feature section to `run_summary.tsv` with at least:

```text
feature_hash_h0
feature_hash_h1
feature_hash_deny
feature_hash_miss
feature_alignment_resolved
feature_alignment_unresolved
feature_assigned_total
feature_unassigned_total
```

Also report the cross-tab:

```text
feature assigned   x spatial candidates present
feature assigned   x spatial candidates absent
feature unassigned x spatial candidates present
feature unassigned x spatial candidates absent
```

and retain the existing valid/invalid UMI accounting.

Required conservation:

```text
H0 + H1 + deny + miss = input read pairs
alignment_resolved + alignment_unresolved = miss
feature_assigned = H0 + H1 + alignment_resolved
feature_assigned + feature_unassigned = input read pairs
four feature/spatial cross-tab cells = input read pairs
```

### 7. Output and diagnostics

Write primary matrices under:

```text
<out-prefix>/SpatialFlex.out/
```

The output contract remains the existing four products and selected scales.
The production route emits no per-read text sidecar.

For deterministic validation, prefer canonical phase hashes or compact audit
digests for:

- eligible paired evidence;
- read-clique member sets;
- corrected support;
- each policy/scale MEX component.

An optional small-fixture audit may be added if a discrepancy cannot be
localized from counters and hashes. It must default off and must not become an
input to resolution or materialization.

## Source locations expected to change

The likely minimal STAR Suite surface is:

```text
core/legacy/source/Parameters.{h,cpp}
core/legacy/source/parametersDefault
core/legacy/source/STAR.cpp
core/legacy/source/ReadAlign_oneRead.cpp
core/legacy/source/ReadAlign_outputAlignments.cpp
core/legacy/source/SpatialGex.{h,cpp}
core/legacy/source/SpatialGexPipeline.cpp
core/legacy/source/ProbeListIndex.{h,cpp}
core/legacy/source/SoloReadFeature_record_shared.h
core/legacy/source/flex/SoloReadFeature_record_flex.cpp
core/legacy/test/test_spatial_gex.cpp
```

A small `SpatialFlex` adapter source is acceptable if it keeps feature routing
and accounting out of the mature downstream numeric kernels. Duplicating the
spatial resolver, Flex gene resolver, UMI correction, or MEX materializer is
not acceptable.

## Unit and component tests

Perform the mandatory clean rebuild before interpreting any result:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
```

Add deterministic tests for:

1. default-off parsing and invalid yes/no values;
2. mutual exclusion with integrated GEX;
3. fail-closed Flex recipe validation;
4. source ordinal pairing, double completion, missing completion, and
   cross-thread isolation;
5. probe-axis order, one-based/zero-based conversion, and bounds;
6. H0, H1, deny, and miss routing;
7. alignment-resolver parity between ordinary Flex and spatial Flex;
8. feature/spatial cross-tab conservation;
9. invalid UMI accounting without molecule creation;
10. spatial barcode `N` fallback and two-`N` parity cases;
11. candidate-family clique construction;
12. candidate-specific `1MM_CR`;
13. strict, soft expected, hard, and gated-hard products;
14. all-memory versus forced-spill byte parity;
15. one-thread versus multithread determinism.

Run the existing native components:

```bash
make -C core/legacy/source test_SpatialGex
make -C core/legacy/source test_SpatialGexSpill
make -C core/legacy/source test_SpatialGexDownstreamSpool
make -C core/legacy/source test_MultiGeneUmiCr
make -C core/legacy/source test_CrMultimapRescuePolicy
make -C flex/tools/molecule_first_resolver clean test
bash tests/run_molecule_first_native_smoke.sh
bash tests/run_molecule_first_reference_conformance.sh
```

## Regression and smoke gates

At minimum, run:

```bash
bash tests/run_flex_hash_alignment_fallback_smoke.sh
bash tests/run_flex_hash_screen_replay_regression.sh
bash tests/run_flex_skip_processing_inline_smoke.sh
bash tests/run_flex_skip_processing_parameter_test.sh
bash tests/run_flex_cr_config_smoke.sh
bash tests/flex_smoke/run_flex_smoke.sh
bash tests/run_solo_smoke.sh
bash tests/run_cbub_regression_test.sh
```

Run the repository parameter-generation regression and the full applicable
core regression suite. An unavailable external fixture is an explicit skip,
not a pass.

The default-off gate compares selected normal Flex, STARsolo, and integrated
spatial GEX products with the v1.6.0 baseline. There must be no intentional
scientific change in those modes.

## CRC 100K gate

### Immutable inputs

Fixture:

```text
/mnt/pikachu/star-spatial/10x/visium_hd_human_crc/downsample_100k
```

It contains exactly 100,000 pairs across lanes L006-L008:

```text
L006: 33,334
L007: 33,333
L008: 33,333
```

Use only R2 and raw R1 FASTQs. Recheck
`fastq_checksums.md5` before the first run.

Reference and feature inputs:

```text
STAR index:
  /mnt/pikachu/autoindex_98_32/visium_hd_human_crc_v2_flex_index
probe axis:
  /mnt/pikachu/autoindex_98_32/visium_hd_human_crc_v2_flex_index/
    probe_gene_list.txt
H0/H1 cache:
  /mnt/pikachu/autoindex_98_32/visium_hd_human_crc_v2_flex_index/
    flex_h01_sequence_cache.bin
barcode contract:
  /mnt/pikachu/star-spatial/10x/visium_hd_human_crc/hd_barcode_contract/
    visium_hd_v1_readspace_suffix_prefix_20260504
```

The oligo inputs must be the same validated BC1/BC2 files selected by the
archived CRC preflight. Re-audit their exact paths and hashes before rendering
the command.

### Frozen settings

```text
R1 UMI: offset 0, length 9
R1 BC1: offset 11, length 14
R1 BC2: offset 25, length 14
R2 probe: offset 0, length 50
feature route: H0/H1 plus alignment fallback
sample tags: none
barcode candidates: frozen spatial decoder through H2
UMI mode: candidate-specific 1MM_CR
temperature: 1
prior alpha: 1
prior beta: 1
spatial lambda: 0
gate posterior: 0.95
gate margin: 0.90
products: strict,soft_expected,hard,gated_hard
scales: 2,8,16
```

### Execution matrix

Use a fresh output directory for every case:

1. one thread, all memory;
2. production thread count, all memory;
3. production thread count, forced spill;
4. independent repeat of the accepted production configuration.

Before each run, render the complete argv and run:

```bash
python3 /mnt/pikachu/visium-hd-processing-recipes/scripts/guard_full_slide_execution.py \
  --scale fixture \
  --route star_integrated_spatial \
  --commands <rendered-command-record>
```

The fixture guard must be tightened to require the integrated route and an
approved integrated switch, as the full-slide guard already does.

The rendered record must not contain:

```text
sort_hd_candidate_feature_ledger.py
normalize_hd_molecule_first_ledger.py
molecule_first_resolver
molecule_first_materialize
```

### Diagnostic oracle

The following directory is read-only diagnostic evidence:

```text
/mnt/pikachu/star-spatial/runs/
  20260725_human_crc_flex_100k_v160_release_gate_v1
```

Do not execute its bridge or use its derived artifacts as native inputs.

Expected feature accounting:

| Quantity | Archived value |
| --- | ---: |
| Input pairs | 100,000 |
| H0 | 90,832 |
| H1 | 1,843 |
| Hash deny | 966 |
| Hash miss | 6,359 |
| Alignment resolved | 4,130 |
| Alignment unresolved | 2,229 |
| Feature assigned | 96,805 |

Expected spatial/molecule diagnostics:

| Quantity | Archived value |
| --- | ---: |
| Raw spatial candidate rows | 113,399 |
| Joined candidate rows | 109,852 |
| Joined reads | 91,283 |
| Read cliques | 79,644 |
| Strict molecules | 66,845 |
| Soft expected mass | 79,164.873884433939 |
| Hard molecules | 79,144 |
| Gated-hard molecules | 71,262 |

These values identify the first divergent cohort; they are not targets to
force. Any difference must be traced to one of:

- H0/H1 decision;
- hash-miss alignment resolution;
- probe-axis translation;
- spatial candidate family;
- UMI eligibility;
- clique membership;
- candidate-specific correction;
- policy materialization.

### Completed CRC 100K result

The accepted native configuration explicitly uses the ordinary Flex production
policy:

```text
--soloMapqMode genomic
```

Fresh output roots:

```text
/mnt/pikachu/star-spatial/native_flex_tests/
  20260726_crc_100k_native_spatial_flex_release_1t_memory_v1
  20260726_crc_100k_native_spatial_flex_release_32t_memory_v1
  20260726_crc_100k_native_spatial_flex_release_32t_forced_spill_v1
  20260726_crc_100k_native_spatial_flex_release_32t_memory_repeat_v1
```

Accepted 32-thread all-memory accounting:

| Quantity | Native value |
| --- | ---: |
| Input pairs | 100,000 |
| R1 pairs with spatial candidates | 94,266 |
| H0 | 90,832 |
| H1 | 1,843 |
| Hash deny | 966 |
| Hash miss | 6,359 |
| Alignment resolved | 971 |
| Alignment unresolved | 5,388 |
| Feature assigned | 93,646 |
| Joined candidate rows | 106,118 |
| Joined reads | 88,241 |
| Read cliques | 76,945 |
| Strict molecules | 64,617 |
| Soft expected mass | 76,495.315452584386 |
| Hard molecules | 76,476 |
| Gated-hard molecules | 68,892 |

Conservation passed exactly:

```text
90832 + 1843 + 966 + 6359 = 100000
971 + 5388 = 6359
90832 + 1843 + 971 = 93646
93646 + 6354 = 100000
```

The input identities were recorded as:

```text
probe axis SHA-256:
  0dc8d7fd2fc0cb8cf54e6707bc5450d853fc890a39fa6626836e332644eae2b8
H0/H1 cache SHA-256:
  7351167d5ad3c0e65fb5f9f5e0f4f4b5885dcc70f6abea3822523ce7ebfc5a7f
probe features:
  18,085
```

All four runs produced the same canonical 36-component MEX manifest:

```text
ce1c010321cb7249d732222a25f9ec311359a8cf07e570e76672d18d27f79296
```

The forced-spill case exercised 1,666 accumulation runs, 256 downstream
contribution runs, and 3,072 matrix runs. Its diagnostic soft-mass accumulator
differed from all-memory by `2.0e-13`, while every emitted soft and integer MEX
component remained byte-identical.

The difference from the archived 4,130 alignment assignments is intentional
and localized. The archive accepted the single primary BAM `GX` annotation
after mapping; it did not run the current shared Flex probe/genomic resolver.
The native H0/H1/deny/miss surface matches the archive exactly. A controlled
`--soloMapqMode off` trace raised alignment resolution to 3,875, demonstrating
that most of the gap is the ordinary Flex genomic MAPQ policy; the remaining
255 reads are current probe/genomic conflict and rescue decisions. Forcing
4,130 would reintroduce the retired alignment-tag compatibility semantics.

### Acceptance

Require:

- exactly 100,000 paired ordinals consumed once;
- exact feature-accounting conservation;
- no BAM or external bridge in the primary ancestry;
- identical feature axis and declared input hashes;
- identical clique/member-set digest across thread and spill modes;
- byte-identical integer MEX components across thread, spill, and repeat runs;
- byte-identical soft MEX components under the existing deterministic
  reduction contract;
- identical declared scientific hashes in the independent repeat;
- an explanation for every material difference from the archived diagnostic
  surface.

Stop after the 100K report. A full CRC or CODEX launch requires separate
authorization.

## Completed regression record

The following release gates passed:

- clean Chromap-enabled STAR build;
- full `core/legacy/source` unit target, including the new ordered-probe-axis
  test and all spatial, spill, downstream-spool, `MultiGeneUMI_CR`, and
  multimap-rescue tests;
- molecule-first resolver unit/CLI tests and frozen Python-reference
  conformance;
- 800,000-read H0/H1 hash-screen replay with zero mismatches;
- Flex hash/alignment fallback smoke;
- Flex skip-processing inline and parameter-policy smokes;
- ordinary Flex CR-config 100K smoke;
- Solo smoke;
- full standard parameter/regression suite: 18/18 passed, including GEX,
  CR-multi tiers, Flex inline, Y-chromosome, SLAM parity, and quick tests.

The CB/UB full fixture script reported an explicit skip because
`/storage/SC2300771_filtered_2M/sample_whitelist.tsv` is absent.

Because terminal-read accounting is shared with integrated GEX, a fresh
ovarian GEX 100K regression was also run:

```text
/mnt/pikachu/star-spatial/gex_in_memory_tests/
  20260726_ovarian_100k_integrated_native_flex_regression_v1
```

All 36 GEX MEX components were byte-identical to the accepted
`26e6d05a` oracle, with canonical manifest:

```text
a6e497bf8e0aecb450e123602fee824968cccb8e3365bb5fb4f55a65a8df5f9b
```

## Recipe and guard follow-up

After the STAR entrypoint and 100K command are stable, update the private
`visium-hd-processing-recipes` repository in its own feature branch:

- add the canonical native Flex recipe and workflow schema;
- include a `COMPOSITION` header;
- add dry-run and tiny executed smoke profiles;
- require `route=star_integrated_spatial`;
- require `--soloSpatialFlexIntegrated yes` at fixture and full scale;
- reject every retired bridge token;
- add the recipe to the catalog only after deterministic rendering tests pass.

Do not place host datasets, run outputs, generated binaries, or archived
oracle files in either repository.

## Mandatory stop conditions

Stop immediately for:

- paired-name, paired-count, or source-ordinal mismatch;
- more or fewer than one feature outcome per read;
- probe-list/cache/reference hash drift;
- a feature index outside the probe axis;
- any normal Flex single-CB insertion in the native spatial primary;
- premature corrected-UMI or hard-coordinate selection;
- loss of feature-unassigned, spatial-unassigned, or invalid-UMI accounting;
- posterior or mass non-conservation;
- thread, spill, or repeat nondeterminism;
- a BAM/tag or retired-bridge dependency;
- default-off Flex, STARsolo, or GEX regression;
- reuse of an existing or partial output directory.

Diagnose on the smallest deterministic fixture and do not proceed to full
scale.
