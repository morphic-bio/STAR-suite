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
--soloCrMultimapRescueEvidence compatibility|annotated
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

`compatibility` preserves the vendor-oriented exon-first rescue contract.
`annotated` is the alignment-aware STAR policy: alignments without a retained
GTF gene are ignored, only the highest-STAR-score annotated alignments enter
the feature-uniqueness decision, and every tied alignment must support the
same single countable gene. A higher-scoring intronic GeneFull alignment can
therefore win over a lower-scoring paralogous exon, and vice versa. Equal-score
different-gene evidence remains unresolved. Exonic/intronic status selects a
representative only after gene uniqueness has been established.

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
16,384-read blocks, reset at each input lane. The contracts joiner recomputes
those digests from paired FASTQs and rejects lane permutations, count
differences, and name mismatches.

For fused input, add `--soloSpatialR1FastqTap FIFO`. STAR remains the only
reader of the paired source FASTQs: it maps the first end (R2) and writes the
already-paired second end (raw R1) to an existing FIFO in canonical input
order. The optimized decoder reads that FIFO and emits its complete decode and
candidate tables. `spatial_feature_sidecar_join --decode-reads FILE` validates
the decoder's one-row-per-read stream against the STAR block digests without
reopening R1 or R2. The tap requires the full fail-closed sidecar recipe, FASTQ
input, and an active FIFO reader; it is absent by default.

Status flags distinguish mapped/unmapped, unique/no-gene/multi-gene outcomes,
same-gene genomic multimappers, exonic versus intronic-fallback CR rescue, and
the GeneFull overlap class. Annotated-mode rejection is explicit for
conflicting best-score genes and a best alignment that itself overlaps
multiple countable genes. Rescue provenance is orthogonal to gene eligibility.
Zero-filled slots are invalid and gene index zero is a valid gene. The legacy
best-score-NA status bit remains reserved in schema v1 but is not emitted by
the annotated policy.

## Spatial adapter and GEX reconciliation

The standalone tools are in `flex/tools/molecule_first_resolver/`:

- `spatial_feature_sidecar_dump` validates and dumps the fixed records.
- `spatial_feature_sidecar_join` joins sidecar ordinals to the current raw-R1
  candidate stream, validates the coordinate contract and block digests, and
  emits the normalized molecule-first evidence schema. Contracts mode verifies
  the 9-nt UMI directly against R1. Fused mode accepts `--decode-reads` instead
  of FASTQ paths and verifies the STAR-owned read order without a source rescan.
  Neither mode accepts BAM/tag fields.
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
The fused single-input-stream implementation and 100K parity procedure are in
`docs/RUNBOOK_VISIUM_HD_GEX_FUSED_READ_PROCESSING_20260724.md`.
