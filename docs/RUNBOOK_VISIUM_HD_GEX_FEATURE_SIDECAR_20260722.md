# Runbook: modern GeneFull sidecar for Visium HD 3-prime GEX

Date: 2026-07-22  
Branch: `feature/visium-hd-gex-sidecar-20260722`  
Base commit: `a996107e271c013e39f9398151deda0017da35d6`  
Companion processing source: `visium-hd-processing` commit `7b581c8`  
Fixture: human ovarian Visium HD 3-prime, deterministic 100K paired-read prefix

## Objective

Extend the working molecule-first Visium HD Flex path to 3-prime GEX without
using BAM `GX`/`UR` records or the legacy Solo replay spool.

STAR will emit a compact per-read GeneFull evidence sidecar after the current
CR-compatible genomic multimapper rescue and GeneFull annotation logic, but
before barcode correction and UMI collapse. The raw-R1 spatial decoder will
join that evidence by read ordinal to its coordinate candidates and raw UMI.
The existing molecule-first resolver will perform spatial candidate resolution
and candidate-specific 1MM-CR UMI correction. A GEX-only reconciliation step
will then apply `MultiGeneUMI_CR` per spatial candidate and corrected UMI.
After the final gene is selected, the existing Flex molecule and materializer
path is reused unchanged.

The target operation order is:

```text
R2 alignment
  -> CR-compatible genomic multimapper rescue
  -> per-read modern GeneFull annotation
  -> fixed-width GeneFull evidence sidecar

raw R1
  -> all supported HD coordinate candidates + likelihoods
  -> raw UMI
  -> join GeneFull evidence by read ordinal
  -> read cliques and spatial posterior
  -> candidate-specific 1MM-CR UMI correction
  -> candidate-specific MultiGeneUMI_CR
  -> final (gene, corrected UMI, coordinate, policy, weight)
  -> existing Flex resolver/materializer output path
```

## Clean source boundary (mandatory)

For this fixture, ignore every previous STAR, decoder, resolver, materializer,
and comparison output. Prior outputs are neither caches nor baselines.

Permitted primary inputs are only:

- the eight FASTQs under the frozen 100K fixture;
- that fixture's `checksums.sha256` and `summary.json`;
- the declared STAR GRCh38 index;
- the two declared full-oligo lists;
- the declared spatial barcode/coordinate contract;
- committed source files from the two development worktrees at recorded
  commits.

Hard prohibitions:

- Do not inspect, copy, link, or consume any existing file below
  `/mnt/pikachu/star-spatial/runs/`.
- Do not reuse any old STAR output, ordered SAM, BAM, GX/UR ledger, candidate
  shard, H0 prior, resolved directory, matrix, summary, or completion marker.
- Do not use existing binaries under either development worktree. Rebuild from
  committed source and record executable SHA-256 values.
- Do not use generated files under the companion processing worktree's
  `build/` directory as inputs.
- Do not use Space Ranger assignments, corrected barcodes, corrected UMIs,
  molecule definitions, or cell/bin calls as computational evidence.

The extant Space Ranger run is permitted only as a final, read-only sanity
comparison after every source-only invariant passes. Its path and checksum
must be recorded before use. It must not supply candidates, priors, genes,
UMIs, filtering decisions, expected counts, or tuning targets. If its location
is under the prohibited `/mnt/pikachu/star-spatial/runs/` tree, stop and obtain
an explicit exception or an allowed copy before reading it.

## Frozen source fixture

Root:

```text
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1
```

FASTQ directory:

```text
/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1/fastqs
```

The lane order is `L001`, `L002`, `L003`, `L004`. Each lane contains exactly
25,000 paired reads. The total is exactly 100,000 pairs. R1 is 43 nt and R2 is
75 nt. The frozen summary states that paired IDs match and within-lane order is
preserved.

Before every fixture run:

```bash
FIXTURE=/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1
(cd "${FIXTURE}" && sha256sum -c checksums.sha256)
```

All nine entries passed on branch creation. Do not regenerate or edit the
fixture manifest.

## Reference and barcode contract

STAR index:

```text
/mnt/pikachu/star-spatial/references/refdata-gex-GRCh38-2024-A_STAR-2.7.11a/star
```

The recorded index has `versionGenome 2.7.4a`, `sjdbOverhang 100`, and was
built from the GRCh38-2024-A FASTA/GTF. The STAR Suite source binary must be
rebuilt cleanly; the index itself is the declared immutable input.

Full oligos:

```text
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc1_full_oligos.txt
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc2_full_oligos.txt
```

Coordinate contract:

```text
/storage/star-spatial/runs/cleanroom_hd_mouse_brain/barcode_contract/
```

Required files:

```text
barcode_coords.tsv
raw28_barcode_coords.tsv
bc1_whitelist.txt
bc2_whitelist.txt
summary.json
```

Frozen contract facts:

- 11,222,500 coordinate records;
- 3,168 distinct BC1 whitelist segments;
- 3,350 BC2 whitelist segments;
- BC1 is R1 offset 11, length 14;
- BC2 is R1 offset 25, length 14;
- `row2` is BC2/Y and `col2` is BC1/X;
- 10,612,800 distinct raw-28 sequence keys;
- 596,300 raw-28 keys have multiple valid coordinates;
- BC1 has 178 colliding segments and a maximum collision multiplicity of 3.

The candidate decoder must therefore map a sequence key to a vector of
coordinates. Any one-sequence-to-one-coordinate map is invalid for this work.

## Historical 3-prime path that must not be reused

The companion processing source currently contains an older prototype:

```text
native/hd_3prime_ordered_candidate_feature_shard.cpp
scripts/run_hd_3prime_candidate_feature_sharding.py
```

That prototype:

- creates a synthetic feature-carrier R1 FASTQ and whitelist;
- asks STAR for ordered SAM output;
- reads `GX` and `UR` from the SAM stream;
- uses `GeneFull_ExonOverIntron` rather than this runbook's modern GeneFull
  contract;
- couples feature extraction to an alignment output format.

It may be read to understand the raw-candidate shard schema, but it must not be
executed, copied into the new driver, or treated as a correctness reference.
No carrier FASTQ, ordered SAM, BAM, `GX`, or `UR` is part of the new path.

The reusable source-side pieces are limited to behavior that is independently
revalidated from the frozen inputs, especially:

- `native/hd_r1_anchored_decode.cpp` for candidate-preserving raw-R1 decode;
- the H0 read-prior aggregation source;
- deterministic candidate shard sorting/merging;
- the native molecule-first resolver and optimized materializer.

## Exact STAR interception point

The current read flow in `ReadAlign::outputAlignments()` is:

1. CR-compatible genomic multimapper rescue mutates the read's alignment set;
2. `ReadAlign::alignedAnnotation()` creates the modern GeneFull read evidence
   and applies the CR exonic-over-intronic filter;
3. `SoloReadBarcode::getCBandUMI()` parses and corrects barcode/UMI fields;
4. Solo records and later collapses barcode/gene/UMI state.

The sidecar hook belongs after step 2 and before step 3. It must execute exactly
once per input read, including unmapped and rejected reads, and only in the
final mapping pass. It must not call `SoloReadFeature::record()`,
`outputReadCB_base()`, the inline bridge hash, the binary Solo spool, or BAM tag
code.

The sidecar record is post-rescue per-read GeneFull evidence. It is not the
final molecule-level gene call. The latter does not exist until
candidate-specific UMI correction and `MultiGeneUMI_CR` complete.

## Proposed sidecar v1

Working CLI name:

```text
--soloSpatialFeatureSidecar <output-prefix>
```

The option is explicit and default-off. The enabled path must require the
modern GEX policy used in this runbook and must fail closed for unsupported
feature/multimapper combinations.

Files:

```text
<prefix>.bin
<prefix>.features.tsv
<prefix>.read_name_digests.tsv
<prefix>.summary.json
```

### Fixed-width read record

Use one little-endian fixed-width slot per global STAR read ordinal:

```c
struct SpatialFeatureRecordV1 {
    uint32_t gene_index;   // UINT32_MAX unless uniquely GeneFull eligible
    uint16_t status_flags;
    uint16_t reserved;     // zero in v1
};
```

The read ordinal is implicit from the record offset. This permits concurrent
`pwrite()` without a shared append lock and allows the R1 adapter to scan the
sidecar sequentially without sorting or hashing read names. A 100-million-read
run uses approximately 800 MB for records, far smaller than an alignment BAM.

Required status information:

- record present;
- mapped versus unmapped/filtered;
- exactly one eligible GeneFull gene;
- no GeneFull gene;
- unresolved multi-gene read rejected by `--soloMultiMappers Unique`;
- genomic multimapper whose GeneFull union is one gene;
- CR exonic rescue;
- CR intronic-fallback rescue;
- GeneFull overlap class needed for diagnostics.

Status values must be explicit and mutually validated. Zero-filled or missing
slots are invalid, not aliases for gene index zero.

### Header and completion contract

The binary header must contain at least:

- magic and schema version;
- byte order and record size;
- total read count;
- feature count;
- STAR Suite version and source commit;
- GeneFull feature type and strand policy;
- CR multimapper-rescue and intronic-fallback settings;
- input-manifest digest;
- gene dictionary digest;
- completion state.

Write to a temporary path. On successful finalization, verify exactly one
present record for every input ordinal, write the final header, flush/fsync,
and atomically rename. An interrupted file must never advertise completion.

`features.tsv` maps the compact `gene_index` to the canonical gene ID and name.
Its checksum is stored in both the binary header and summary. Internal gene
indices must never be interpreted against a different reference.

The summary records outcome counts, rescue counts, first/last ordinal, input
lane boundaries, parameter values, source commit, binary SHA-256, and failure
status.

### Ordering validation

The sidecar join relies on the same global lane/read order in STAR and the raw
R1 decoder. The driver must:

- pass `L001` through `L004` in that exact order to both processes;
- assign a global ordinal across lane boundaries;
- validate the frozen FASTQ checksums before launch;
- calculate a per-lane or blockwise normalized-read-name digest independently
  in both streams;
- reject the join on any count, boundary, or digest mismatch.

Do not store every read name merely to make the join work. Names may appear in
an opt-in diagnostic trace, but the production join key is the ordinal.

## Modern GeneFull policy

The first 100K implementation must use:

```text
GeneFull
--soloCrMultimapRescue yes
--soloUMIfiltering MultiGeneUMI_CR
--soloUMIdedup 1MM_CR
--soloMultiMappers Unique
--soloStrand Forward
--soloCellFilter None
--outSAMtype None
```

The final rendered command must not contain `GX`, `GN`, `UR`, `UB`, a BAM/SAM
output type, `Gene`, `GeneFull_ExonOverIntron`, a carrier whitelist, or the
experimental Solo binary-spool environment variables.

The sidecar mode may enable the annotation structures needed for GeneFull, but
its output must not depend on a successful Solo barcode match. If existing
parameter initialization cannot express annotation-without-barcode-collapse,
add a narrowly scoped sidecar mode rather than fabricating a compatible
barcode read.

Per-read policy:

- a post-rescue GeneFull set of size one records that gene as eligible;
- a set of size zero records `no_gene`;
- a set larger than one records `multi_gene_rejected` for canonical Unique
  mode and does not invent a winner;
- rescue provenance is orthogonal to the GeneFull outcome.

## R1 candidate join

The raw-R1 decoder owns:

- raw 9-nt UMI extraction;
- BC1/BC2 observed sequence and qualities;
- every supported coordinate candidate;
- sequence log likelihoods;
- exact H0 read counts used by the frozen read prior.

The STAR sidecar owns only post-rescue GeneFull evidence. The joiner must never
read a corrected barcode/UMI or alignment tag.

For each raw-R1 candidate-bearing read ordinal with one eligible gene, emit the
existing normalized evidence fields:

```text
read_id
feature_id
raw_umi
candidate
log_sequence_likelihood
exact_read_count
```

Candidate rows for one read must preserve the complete coordinate set. The
join may discard no candidate because another coordinate has the same raw-28
sequence. Rejected/no-gene sidecar states contribute to explicit accounting
but do not become feature rows.

For scale, prefer a compact internal candidate stream plus a separate exact
prior table over repeating large strings and prior counts in every temporary
row. Any compact v2 input is private to execution; the existing normalized TSV
route remains available for conformance and small fixtures.

## GEX-specific post-UMI reconciliation

The current molecule-first resolver already provides:

- read-clique construction by provisional feature and raw UMI;
- clique-level spatial posterior and strict/hard/gated calls;
- candidate-specific exact or 1MM-CR UMI correction;
- policy molecule rows consumed by the optimized materializer.

GEX adds one stage after candidate-specific UMI correction. For each spatial
candidate and corrected UMI, collect support across provisional genes and
apply STAR's current `MultiGeneUMI_CR` winner/tie/original-UMI checks. Factor
that policy into a shared helper used by both the normal direct bridge and the
molecule-first GEX path. Do not copy or approximate the legacy implementation.

For the integer strict, hard, and gated-hard policies, the result must exactly
match the shared integer helper. Soft expected counts use candidate-weighted
support and must be labeled as molecule-first expected-count semantics rather
than byte-identical standard Solo MEX semantics. Its rule and mass behavior
must be frozen in tests before a publication run.

After reconciliation, output the same molecule shape as Flex:

```text
(feature_id, corrected_umi, candidate_coordinate, product, weight)
```

From that point onward, reuse the existing optimized resolver/materializer.
Do not add a second GEX-specific matrix writer.

## Implementation tranches

### Tranche 0: runbook and clean-room preflight

- Commit this runbook before implementation.
- Record both repository commits and dirty-state checks.
- Verify all frozen FASTQ checksums.
- Verify reference and contract paths.
- Verify the output root is new and empty.
- Confirm no process or command reads the prohibited run tree.

### Tranche 1: sidecar writer and reader

- Add the explicit default-off parameter and fail-closed validation.
- Add the fixed-width writer at the post-annotation/pre-barcode hook.
- Add atomic completion and feature dictionary/summary outputs.
- Add an independent reader/dumper used only by tests and the spatial adapter.
- Keep all standard STAR/Solo outputs unchanged when the option is absent.

### Tranche 2: ordinal R1 join

- Add global read ordinal to the raw-R1 candidate stream.
- Validate lane boundaries and read-name digests.
- Replace the ordered-SAM/GX/UR adapter with the sidecar join.
- Preserve every raw-R1 coordinate candidate.
- Build H0 priors only from the current raw R1 inputs.

### Tranche 3: candidate-specific MultiGeneUMI_CR

- Factor the current production rule into a shared tested helper.
- Apply it after candidate-specific 1MM-CR correction for each policy.
- Emit standard molecule rows for the existing materializer.
- Preserve the Flex behavior and default-off GEX behavior exactly.

### Tranche 4: 100K source-only integration

- Clean-build STAR Suite and all companion tools from recorded sources.
- Run only the frozen 100K FASTQs into the fresh output root.
- Materialize strict, soft expected, hard, and gated-hard matrices at 2, 8,
  and 16 micrometers.
- Complete all internal gates before inspecting Space Ranger.

## Required tests

### Sidecar unit and contract tests

- little-endian encode/decode round trip;
- exact header and record-size checks;
- gene index zero is distinct from a missing record;
- truncated, incomplete, duplicate-write, and missing-slot files fail closed;
- unique, no-gene, multi-gene, unmapped, same-gene multimapper, exonic rescue,
  and intronic-fallback records have correct status;
- one versus many mapping threads produce byte-identical finalized sidecars;
- two-pass/restart modes cannot duplicate records;
- default-off synthetic Solo outputs remain byte-identical.

### Candidate join tests

- ordinal and per-block read-name digests match;
- a lane-order permutation is rejected;
- candidate-less and feature-less read accounting reconciles;
- all candidates for colliding raw-28 keys survive;
- raw UMI comes from R1, never `UR` or `UB`;
- no `GX`, `CB`, `UB`, Space Ranger field, or cell-call field is accepted;
- exact-read priors are PCR-inclusive and derived only from this invocation.

### GEX policy tests

- post-rescue unique-gene evidence matches the in-memory GeneFull set at the
  sidecar hook;
- unresolved multi-gene reads remain rejected in Unique mode;
- candidate-specific 1MM-CR correction precedes cross-gene reconciliation;
- shared `MultiGeneUMI_CR` tuples match the normal direct bridge exactly;
- equal gene support rejects the molecule;
- original-UMI dominance checks match STAR's current rule;
- Flex resolver outputs do not change.

### 100K integration gates

- fixture checksum verification passes;
- exactly 100,000 sidecar records finalize;
- exactly 25,000 R1/R2 pairs reconcile in each lane;
- R1 and STAR blockwise read-name digests match;
- all observed candidates are in the declared 11,222,500-coordinate universe;
- every candidate set is nonempty and contains no duplicate coordinate;
- no unsupported coordinate appears after resolution;
- strict/hard/gated integer products are deterministic across repeated runs;
- soft expected mass is finite, bounded, and conserved within the frozen
  tolerance;
- every policy conserves 2-micrometer mass at 8 and 16 micrometers;
- feature/barcode axes, dimensions, and NNZ counts are valid for all twelve
  matrices;
- the wrapper writes its completion summary only after every child exits zero.

Space Ranger comparisons are descriptive sanity checks only. Report overlap,
mass, correlations, and obvious axis/reference mismatches without changing a
threshold or policy to improve agreement.

## Build and execution hygiene

Before testing any STAR behavior, perform the mandatory clean rebuild:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
```

Build companion native tools from their recorded source commit into a fresh
`bin/` directory under the test root. Record command lines and SHA-256 values.

Do not run another benchmark concurrently. Use wrapper-written completion
artifacts as authoritative; do not infer completion from `pgrep`, `Log.out`,
or partial matrices.

## Fresh output root

Created empty on 2026-07-22:

```text
/mnt/pikachu/star-spatial/gex_sidecar_tests/20260722_ovarian_100k_v1
```

This path is outside the prohibited historical run tree. Do not reuse it after
a failed or partial integration run. Preserve the failed directory for audit
and create a new versioned sibling such as `_v2` for the next attempt.

Planned layout:

```text
source_manifest/
bin/
star/
  GeneFull.read_evidence.bin
  GeneFull.read_evidence.features.tsv
  GeneFull.read_evidence.summary.json
r1_candidates/
joined_shards/
resolved/
policy_mex_1mm_cr/
logs/
commands.json
RUN_COMPLETE.json
```

Large binaries, sidecars, candidate shards, matrices, FASTQ derivatives, and
logs remain untracked. Only small source manifests and final validation
summaries may be proposed for version control after review.

## Execution ledger

| Gate | Command or artifact | Status |
|---|---|---|
| STAR Suite feature branch | `feature/visium-hd-gex-sidecar-20260722` from `a996107` | complete |
| Companion processing source | branch `feature/visium-hd-gex-sidecar-20260722` from `7b581c8` | ordinal source change implemented; no output reused |
| Frozen fixture checksum | `(cd "$FIXTURE" && sha256sum -c checksums.sha256)` | passed, 9/9 |
| Fixture shape | four lanes x 25,000 pairs; R1 43 nt; R2 75 nt | verified from frozen summary |
| Barcode contract | 11,222,500 coordinates; 596,300 ambiguous raw-28 keys | verified from contract summary |
| Fresh output root | `/mnt/pikachu/star-spatial/gex_sidecar_tests/20260722_ovarian_100k_v1` | created empty |
| Runbook committed before implementation | commit `68fb69e` | complete |
| Clean STAR build | command above | passed after implementation; warnings only |
| Sidecar unit tests | `make -C core/legacy/source test_SpatialFeatureSidecar` | passed before clean-build gate |
| Candidate join tests | `test_spatial_feature_sidecar_join` including lane permutation/raw-UMI rejection | passed before clean-build gate |
| Shared MultiGeneUMI_CR tests | core helper plus GEX resolver correction/tie/original-dominance cases | passed before clean-build gate |
| Source-only 100K run | fresh output root | not run |
| Space Ranger sanity check | comparator path/checksum to record after internal gates | not run |

## Commit plan

Keep the work reviewable and preserve source-policy provenance:

1. runbook and artifact-location declaration;
2. sidecar parameter, binary contract, writer/reader, and unit tests;
3. raw-R1 ordinal join and source-only wrapper tests;
4. shared `MultiGeneUMI_CR` factor plus GEX resolver integration;
5. 100K validation summary and final documentation.

Never add `Co-Authored-By` trailers. Do not commit generated binaries or test
outputs.
