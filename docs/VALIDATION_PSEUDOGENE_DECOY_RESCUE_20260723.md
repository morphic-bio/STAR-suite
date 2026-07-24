# Pseudogene-decoy rescue validation, 2026-07-23

> **Withdrawn 2026-07-24.** The `decoy` policy described below treated every
> annotation-free genomic alignment as `NA` evidence. Full-slide validation
> showed that this was much broader than the motivating pseudogene case. The
> option has been removed. The replacement `annotated` policy ignores
> annotation-free alignments, ranks retained GeneFull evidence by STAR score,
> and resolves equal-best different-gene evidence conservatively. See
> `VALIDATION_ANNOTATED_SCORE_RESCUE_20260724.md`.

## Outcome

STAR can now keep an unannotated or non-countable locus as `NA` alignment
evidence without making it a countable feature. The new explicit policy is:

```text
--soloCrMultimapRescueEvidence decoy
```

Only maximum-STAR-score alignments enter the feature-uniqueness decision. All
tied best alignments must support the same single countable gene. A best-score
`NA` alignment vetoes feature uniqueness; a lower-score `NA` alignment does
not erase a better annotated alignment. After uniqueness is established, an
exonic alignment is preferred as the representative, followed by an intronic
alignment when intronic rescue is enabled. `NA` is never emitted as a gene.

The existing vendor-oriented behavior remains the default and is selected
explicitly with:

```text
--soloCrMultimapRescueEvidence compatibility
```

This avoids silently changing existing CR-parity workflows. The compatibility
implementation also corrects an alignment-order bug: a later intronic
candidate can no longer overwrite the index of the single exonic winner.

## Code and binary identity

- Base: `2fc05770d1b0528abf3709da66199b708ca17afb`
- Implementation commit: `61c29379aa95d67aeb6d21b11aae6ecaf5505aac`
- Branch: `fix/pseudogene-exon-rescue-20260723`
- STAR-suite version: `1.5.0`
- Clean release binary:
  `core/legacy/source/STAR`
- Binary SHA-256:
  `6f88215a889c6479cb9754574d9261bd8acfb237079cb0a040f78ff1a15a145f`

The binary's embedded source identity is the implementation commit with an
empty `diff files` field.

## External validation fixture

The fixture contains 200 paired reads: 100 ordinary unique MECOM reads and
100 reads from the RPL22/pseudogene family identified during the ovarian
Xenium audit.

| Input | SHA-256 |
|---|---|
| `sample_R1.fastq` | `d86e4b31d8e98b333a57a114f0a2f0efbcd7fd1191d5c76625ecbe9eb3cb0bb8` |
| `sample_R2.fastq` | `311979b9ec3e8c673df8e8f1b4dd510835677bce4b70bfae80ffe7b50d553b1c` |
| `sample.tsv` | `5c267fa65918d91004eaf3f63300860cefbf3aee9436e9f071daf0e56f01b128` |

Fixture directory:

```text
/mnt/pikachu/star-spatial/runs/20260723_ovarian_xenium_sr_first_v1/mecom_audit_v1/read_sample
```

The index is:

```text
/mnt/pikachu/star-spatial/references/refdata-gex-GRCh38-2024-A_STAR-2.7.11a/star
```

Selected index hashes are:

| Index file | SHA-256 |
|---|---|
| `Genome` | `9d97df99c23bc0d84323ece5bac1c72200343eab459a767786c14eff7ecd9d4c` |
| `genomeParameters.txt` | `c0a34434da7a29257036109a7593968b9b8c7faf3538f8c877c43de2c6f7663c` |
| `geneInfo.tab` | `40f6ac3112cdb0c36ccef1543eb40c075904f02efc2aae33644f1d90424a6e08` |
| `exonGeTrInfo.tab` | `5eb3694cb41c0d9159cbea7f46fdd2e812701f3385dbd3b08e121c4d9cfd2445` |
| `exonInfo.tab` | `2d4e2260991cd97219776fe89130a262dac70e87a6b61ea428e21d13223289bd` |
| `transcriptInfo.tab` | `575664d55782dc3209837b8783daf3fe000eea89ab5f531883d2aacbd91d1666` |

## Replay command

The two policies used the same command, changing only `MODE` and the fresh
output directory:

```bash
core/legacy/source/STAR \
  --runThreadN 4 \
  --genomeDir "$INDEX" \
  --readFilesIn "$R2" "$R1" \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloType None \
  --soloFeatures GeneFull \
  --soloCrGexFeature GeneFull \
  --soloCrMultimapRescue yes \
  --soloCrMultimapRescueIntronic auto \
  --soloCrMultimapRescueEvidence "$MODE" \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloMultiMappers Unique \
  --soloStrand Forward \
  --soloCellFilter None \
  --outSAMtype None \
  --soloSpatialFeatureSidecar "$OUT/features" \
  --outFileNamePrefix "$OUT/"
```

The checked-in `spatial_feature_sidecar_dump` source was used to validate and
summarize every fixed-width record.

## Results

| Fixture cohort | Reads | Compatibility | Decoy |
|---|---:|---:|---:|
| Ordinary unique MECOM | 100 | 100 MECOM | 100 MECOM |
| RPL22/pseudogene family | 100 | 100 RPL22 | 13 RPL22, 87 unassigned |
| Total unique feature calls | 200 | 200 | 113 |

For the 100 RPL22/pseudogene-family reads, decoy mode found:

- 13 annotated RPL22 alignments whose STAR score was one point higher than
  the retained decoy alternative; these remain RPL22.
- 84 reads with conflicting countable genes among equal best-score
  alignments; these remain unassigned.
- 3 reads with `NA` present in the best-score evidence set; these remain
  unassigned. `NA` has deterministic diagnostic precedence if another
  ambiguity is present in the same best-score set.

There were no per-alignment multi-gene winners and no intronic-only rescues.
The 100 ordinary MECOM calls were unchanged. This is the intended correction:
removed annotation changes feature eligibility but does not remove the locus
from genomic uniqueness or score comparison.

Final output directories:

```text
/mnt/pikachu/star-spatial/runs/20260723_pseudogene_decoy_validation/compatibility_commit_v1
/mnt/pikachu/star-spatial/runs/20260723_pseudogene_decoy_validation/decoy_commit_v1
/mnt/pikachu/star-spatial/runs/20260723_pseudogene_decoy_validation/decoy_commit_repeat_v1
```

## Determinism

The two fresh decoy runs are byte-identical for every declared sidecar
artifact:

| Artifact | SHA-256 |
|---|---|
| `features.bin` | `39edc0444309d040c5a9c4e4351e2ffe033adf7fcebf1aee029d536a96e00487` |
| `features.features.tsv` | `f60d9c932ae7d5a933af84a2d23e23570fd14c81afce4c3c1f4d7ebcca78f0b3` |
| `features.read_name_digests.tsv` | `0619cdf27f8984f0de5cf412ab1648ab4992ea0f14853781adf60acd04326980` |
| `features.summary.json` | `5fe9c7f5cc581fc60d85aa59bbd16776f3e7e87998586f837ab79a2dd1e8adab` |

The compatibility binary and summary hashes are respectively
`78d1c4d1d0508df5f868e8de6b8a07f03ed448a60c33b9a59eff51113254f2c5`
and
`9c4c76c22af13a763d66738e32cef4d1c8d7d4895cd95eea2b051f654f233043`.

## Regression gates

Passed:

- Clean release build: `make clean && make -j8 STAR`.
- `make test_CrMultimapRescuePolicy`.
- `make test_SpatialFeatureSidecar`.
- `make test_MultiGeneUmiCr`, including 20,000 legacy/bridge parity cases.
- Invalid `--soloCrMultimapRescueEvidence` value rejected with parameter exit
  code 102.
- Real 200-read compatibility and decoy runs.
- Fresh decoy replay with identical declared hashes.
- `tests/run_scrna_sidecar_off_golden.sh`; normal sidecar-off scRNA count
  artifacts remain byte-identical to golden commit
  `a996107e271c013e39f9398151deda0017da35d6`.

Repository-level test limitations, unrelated to this change:

- `make test` stops at the existing `test_packed_readinfo` target because
  `../test/PackedReadInfo_test.cpp` and source-tree `PackedReadInfo.cpp` do not
  exist at the master base.
- `make compat-unit-tests` currently fails to link because its target omits
  required existing HTSlib, chromap-contract, libscrna, probe-filter, and
  related libraries.
- `tests/run_cr_compat_integration_smoke.sh` completed a clean ASAN STAR build
  but could not execute because the referenced
  `tests/solo_smoke/fastq/R2.fastq` fixture is absent from the checkout.

No test failure implicated the new rescue policy, sidecar record validation,
normal scRNA output, or the real pseudogene-decoy replay.

## Master integration validation

The feature branch was merged without conflict onto fetched master
`2fc05770d1b0528abf3709da66199b708ca17afb` in a separate clean worktree.
The tested merge commit is
`3dfec3c4cd3f184d94791d0d64eb3cf4aaa63c3f`.

From that merge commit:

- A clean release build passed. Its binary SHA-256 is
  `9263645ef24dd002fd8f6f1ce135a91ad43a75ffe23affb033d552755eee939b`,
  and its embedded source identity has an empty `diff files` field.
- The policy, sidecar, and 20,000-case `MultiGeneUMI_CR` parity tests passed.
- The normal sidecar-off scRNA golden test passed.
- The real decoy replay at
  `/mnt/pikachu/star-spatial/runs/20260723_pseudogene_decoy_validation/decoy_merge_v1`
  reproduced 100 MECOM, 13 RPL22, and 87 unassigned records with the same
  3 best-score-NA and 84 conflicting-gene diagnostics.
- The checked-in dump of all 200 fixed-width records was byte-identical as
  text to the implementation-commit replay. The merge sidecar binary differs
  only in declared build/source provenance and has SHA-256
  `a6cbbf879c4a2d68b56978bd9036421dc362bcc3e09f19e98cc99c84db11d4fb`.
