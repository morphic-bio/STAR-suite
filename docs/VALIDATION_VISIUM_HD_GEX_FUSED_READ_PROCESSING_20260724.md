# Validation: fused paired-read processing for Visium HD 3-prime GEX

Date: 2026-07-24

Implementation commit: `459028af4116ed124ed12ca1bee0a596bcc8d7a1`

Companion decoder commit: `5e9af58f9b86e9d3b95f8612a4a77e40f9a3ed86`

Runbook: `docs/RUNBOOK_VISIUM_HD_GEX_FUSED_READ_PROCESSING_20260724.md`

## Result

The fused mode passed the 100K ovarian Visium HD 3-prime parity gate. STAR was
the sole reader of the paired source FASTQs, mapped R2 only, and streamed the
already-paired raw R1 records through the bounded FIFO tap. The decoder used
that FIFO as its only read input. The fused join used the complete decoder row
stream and STAR sidecar digests, with zero R1/R2 FASTQ arguments.

All 49 scheduling-independent files compared between the fresh contracts and
fused runs had identical byte sizes and SHA-256 digests:

- decoder rows, complete candidate rows, oligo statistics, and H0 prior;
- STAR sidecar binary, feature dictionary, and read-name digest table;
- normalized molecule evidence and resolver read cliques;
- strict, soft-expected, hard, and gated-hard resolver products; and
- all 36 feature, barcode, and matrix files for four policies at 2, 8, and 16
  micrometers.

No mismatch was observed. Command manifests also confirm that the fused FIFO
was removed after producer completion and that the wrapper wrote
`RUN_COMPLETE.json` only after both repeated resolver and materializer gates.

## Fresh source-only runs

Contracts oracle:

```text
/mnt/pikachu/star-spatial/gex_fused_tests/20260724_ovarian_100k_contracts_v1
```

Fused run:

```text
/mnt/pikachu/star-spatial/gex_fused_tests/20260724_ovarian_100k_fused_v1
```

Both wrappers recorded the implementation and companion commits above. They
validated the frozen fixture checksums and consumed no earlier STAR, decoder,
candidate, sidecar, resolver, or matrix artifact.

Completion-marker SHA-256 values are:

```text
79a18cd18f76b6d1a3e4d8b5846c3cc4b8635e13aab941deb862a1117f1c5e12  contracts/RUN_COMPLETE.json
cde6325afa0d06c25033ccad956204c2f1dab4b73acda44ac5b7fe434a8694be  fused/RUN_COMPLETE.json
```

## Count reconciliation

Both modes reported exactly:

| Quantity | Count |
|---|---:|
| Input read pairs | 100,000 |
| Sidecar records | 100,000 |
| Reads with coordinate candidates | 89,156 |
| Candidate rows | 111,744 |
| Gene-eligible reads | 85,953 |
| Emitted gene-plus-coordinate reads | 78,453 |
| Emitted normalized evidence rows | 98,269 |
| No-gene reads | 4,579 |
| Multi-gene-rejected reads | 2,971 |
| Unmapped or filtered reads | 6,497 |

The fused completion invariants are all true, including
`single_star_input_stream`, `join_reopened_fastqs=false`,
`lane_order_and_name_digests_match`, `raw_umi_from_r1`,
`gex_multigene_umi_cr_enabled`, repeated integer-product determinism, and mass
conservation at all three scales.

## Command and scheduling audit

The contracts join reopened four R1 and four R2 files. In fused mode:

- the decoder had one `--r1-fastq` value, the output-local FIFO;
- STAR had the eight fixture FASTQs plus `--soloSpatialR1FastqTap FIFO`;
- the join had `--decode-reads` and no R1/R2 FASTQ argument; and
- the FIFO did not remain after successful completion.

The contracts producer phase was serial with 16 threads reused by each branch.
The fused producer phase respected the total 16-thread cap as 8 decoder plus 8
STAR threads. Producer wall times were 59.31 seconds for contracts and 13.92
seconds for fused. The latter ran after the contracts run and therefore used a
warm genome-index cache. These numbers validate scheduling and do not establish
a performance claim; representative full-slide warm/cold benchmarks remain
separate work.

## Default-off and fixture-free gates

The following passed before the 100K runs:

```text
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
make -C core/legacy/source test_SpatialR1FastqTap
make -C flex/tools/molecule_first_resolver all test
python3 tests/test_visium_hd_gex_sidecar_concurrency.py
tests/run_molecule_first_native_smoke.sh
tests/run_spatial_r1_tap_guard.sh
tests/run_scrna_sidecar_off_golden.sh
```

The normal scRNA test omitted all spatial options and reproduced the frozen
GeneFull MEX byte-for-byte. The tap-only guard failed before STAR opened normal
read input. The FIFO unit test verified exact record order/content and rejection
of a regular-file target. The join tests verified byte-identical contracts and
fused evidence and rejected reordered and short decode streams. The native
materializer smoke also verified that `--products soft_expected` emits only the
selected real-valued policy directory.

## Read-only Space Ranger sanity check

After internal parity passed, the fused matrices were compared read-only with
the published ovarian Space Ranger 4.0.1 raw binned matrices. All 12
policy/scale comparisons completed and every observed feature ID existed in the
Space Ranger feature axis.

For the hard policy, gene-total raw Pearson was `0.986420`, common-detected-gene
Spearman was `0.807402`, and 93 of the top 100 genes overlapped. Fractions of
our occupied bins present in Space Ranger were `0.997937`, `0.999907`, and
`0.999979` at 2, 8, and 16 micrometers. Raw spatial-bin Pearson rose from
`0.066684` to `0.216267` and `0.386860` across those scales, as expected for a
100K prefix compared with the full 474-million-read run.

The complete descriptive report is:

```text
/mnt/pikachu/star-spatial/gex_fused_tests/20260724_ovarian_100k_fused_v1_space_ranger_sanity.json
```

Its SHA-256 is
`b19e44c8d6721bddba78ea8d38aeb7e2a029a680fce98e9f4dd6770f23e56c68`.
Space Ranger evidence was not used for resolution, tuning, or parity decisions.
