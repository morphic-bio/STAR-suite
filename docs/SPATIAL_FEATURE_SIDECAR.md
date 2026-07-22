# Spatial GeneFull feature sidecar

`--soloSpatialFeatureSidecar PREFIX` is a default-off, annotation-only output
for molecule-first Visium HD 3-prime GEX. It records modern post-rescue
GeneFull evidence by global input-read ordinal before any spatial barcode or
UMI correction.

This is not a BAM-tag export. `GX` is an alignment-level annotation and is not
the final modern GeneFull/CR-compatible molecule decision. `UR` is the raw UMI,
not a corrected UMI. The sidecar recipe therefore forbids `GX`, `GN`, `UR`,
`UB`, `CB`, and `CR`, and requires `--outSAMtype None`. Raw UMI and coordinate
candidates come independently from raw R1.

## Required STAR policy

The option currently fails closed unless the command uses:

```text
--soloType None
--soloFeatures GeneFull
--soloCrGexFeature GeneFull
--soloCrMultimapRescue yes
--soloUMIdedup 1MM_CR
--soloUMIfiltering MultiGeneUMI_CR
--soloMultiMappers Unique
--soloStrand Forward
--soloCellFilter None
--outSAMtype None
```

Pass R2 first and raw R1 second in `--readFilesIn`. Only R2 is aligned. The
sidecar hook runs once after CR genomic-multimapper rescue and modern GeneFull
annotation, but before Solo barcode/UMI parsing.

## Output contract

The prefix produces:

```text
PREFIX.bin
PREFIX.features.tsv
PREFIX.read_name_digests.tsv
PREFIX.summary.json
```

The binary begins with a 1,024-byte little-endian header and then one 8-byte
record per global input ordinal:

```c
uint32_t gene_index;   // UINT32_MAX unless exactly one GeneFull gene is eligible
uint16_t status_flags;
uint16_t reserved;     // zero in schema v1
```

The header carries schema/record sizes, completion state, read and feature
counts, policy/provenance text, the input-manifest digest, and the canonical
gene-dictionary digest. A completed reader requires an exact file size, a
present valid record at every ordinal, and `records_written == total_reads`.
Writers use ordinal `pwrite` slots and commit the binary by atomic rename only
after a full scan succeeds. Duplicate writes and sparse/missing slots fail
finalization.

`features.tsv` maps the compact index to canonical reference gene ID/name.
`read_name_digests.tsv` contains SHA-256 digests over normalized-name hashes in
16,384-read blocks, reset at each input lane. The raw-R1 joiner recomputes those
digests from paired FASTQs and rejects lane permutations, count differences,
and name mismatches.

Status flags distinguish mapped/unmapped, unique/no-gene/multi-gene outcomes,
same-gene genomic multimappers, exonic versus intronic-fallback CR rescue, and
the GeneFull overlap class. Rescue provenance is orthogonal to gene
eligibility. Zero-filled slots are invalid and gene index zero is a valid gene.

## Spatial adapter and GEX reconciliation

The standalone tools are in `flex/tools/molecule_first_resolver/`:

- `spatial_feature_sidecar_dump` validates and dumps the fixed records.
- `spatial_feature_sidecar_join` joins sidecar ordinals to the current raw-R1
  candidate stream, validates the coordinate contract and block digests, and
  emits the normalized molecule-first evidence schema. It accepts no BAM/tag
  fields and verifies the 9-nt UMI directly against R1.
- `molecule_first_resolver --gex-multigene-umi-cr` performs candidate-specific
  exact or `1MM_CR` UMI correction and then calls the same shared
  `MultiGeneUMI_CR` helper used by ordinary Solo collapse.

Integer strict/hard/gated-hard products use read counts and the production
winner/tie/original-UMI rule. Soft output is explicitly
candidate-weighted expected-occupancy semantics: it reconciles provisional
gene occupancies after weighted UMI correction and reports input/output mass in
`summary.tsv`. It is not represented as byte-identical standard Solo MEX
semantics.

The frozen clean-room implementation and validation procedure is
`docs/RUNBOOK_VISIUM_HD_GEX_FEATURE_SIDECAR_20260722.md`; its executable driver
is `scripts/run_visium_hd_gex_sidecar_100k.py`. The producer fork/join schedule,
thread-budget contract, and serial/concurrent parity procedure are documented in
`docs/RUNBOOK_VISIUM_HD_GEX_CONCURRENT_PRODUCERS_20260722.md`.
