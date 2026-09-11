# LARRY search corrections — 2026-09-11

Worktree: `/mnt/pikachu/STAR-suite-larry-regression-20260911`, branch
`fix-larry-feature-regression-20260911`. These are development results, not a
replacement for the full MSK paper benchmark. No Cell Ranger/cyto source was inspected.

## Finding

The large feature hash is populated and queried correctly. All 245,979 reference
sequences resolve to themselves. On the first 1,000 instrumented learning calls,
both the preserved current implementation and a clean v1.4.3 control resolved 861
reads and performed 644,214,716 cached-offset accesses. The expensive work surrounds
the cache: learning iterates over every feature after a lookup has identified the
only possible winner. A broad fallback also scans every feature on reads that the
complete Hamming prehash can prove have no match.

The v1.4.3 commit is `e4292f1c15449a60e62de2ff142f3de54be67546` (the annotated tag
object is different). Both diagnostic runs were deliberately stopped after tracing;
neither is a completed timing control. The historical April paper log says
`cc84382-dirty`; it is not proof that the tagged v1.4.3 binary produced that run.
The matched learning configuration reproduces the repeated work in v1.4.3 too.

## Changes

1. For a single uniform anchor group covering the full reference, visit only the
   features that cached offset lookups can return. Preserve original feature and
   offset ordering, adaptive distance ceilings, first positions, ambiguity rules,
   and learned histograms. Other reference layouts keep the original loop.
2. Before the broad single-thread scan, use a complete cumulative H0/H1/H2 hash as
   a negative filter. Possible hits, including ambiguous payloads, still run the
   original scan. This guard is restricted to uniform, byte-aligned feature lengths
   and decodes the same packed windows, including existing padded-tail semantics.

Learning remains 100,000 reads in the production test. Barcode correction, input
selection, feature acceptance, and post-learning policy are unchanged. Earlier
zero-N and deterministic learning corrections are committed separately in `e2bb279`.

## Input and validation

Artifacts root: `/home/lhhung/pf_larry_regression_20260911/`.
The `larry_subset_200k` fixture is the first 25,000 paired records from each of eight
MSK ES LARRY lanes: 200,000 pairs total, 16 ordinary-gzip FASTQs, full feature reference
and whitelist. It is a lane-balanced prefix, not a random sample. Manifest verification
checks input/output identities, paired names, counts, gzip CRCs, and SHA256 digests.
Auto detection selects gzip on all eight lanes and starts zero BGZF inflater workers.

`larry_fix/parity2` compares preserved pre-change and corrected libraries on 27 cases:
suffix/prefix/dual anchors, mixed lengths, 64/128-bit keys, small/large references,
H0/H1/H2, allowed/forbidden Ns, repeated/multiple anchors, broad fallback, short tails,
and learning cutover. All per-read assignments, distances, positions, learned histograms
and modes are exactly equal. This includes 1,600 actual LARRY reads with the full reference.
`larry_fix/asan2` instruments the changed translation unit and driver with AddressSanitizer;
all 27 cases pass without sanitizer errors and retain the same ledgers.

The isolated real-read driver (including reference/prehash initialization) takes
28.49 s before and 7.78 s after both changes. This is a correctness fixture, not a
full-library performance claim.

## Completed PF-only subset runs

No genome index, GEX mapping, BAM, sidecar, or EmptyDrops is included.
Whole-process time includes initialization and cleanup; all rows exited successfully.

| Learning reads | Implementation | Wall | Deduplicated UMIs | Raw barcodes |
|---|---|---:|---:|---:|
| 1,000 (diagnostic control) | Preserved pre-change | 28.76 s | 159,027 | 35,626 |
| 1,000 (diagnostic control) | Uniform-anchor change only | 26.68 s | 159,027 | 35,626 |
| 1,000 (diagnostic control) | Both changes | 16.57 s | 159,027 | 35,626 |
| 100,000 (normal policy) | Both changes | 279.12 s | 159,035 | 35,626 |

The matched 1,000-learning comparison has zero differing sparse matrix entries after
barcode reordering: 42.4% lower wall time. It is not the production learning setting.
The 100,000-learning run assigns 171,279 feature counts, has 26,442 unmatched reads,
and peaks at 1,731,056 KiB RSS. There is no completed matching 100,000-learning baseline,
so do not assign a production speedup from these rows. The eight additional UMIs across
thresholds reflect a different learning configuration and are not a same-policy parity test.
The earlier `after_100000` attempt was stopped for a stack diagnosis; its elapsed time
must not be presented as a completed benchmark.

The normal-learning run still takes 4m39s and averages only 1.31 CPU cores. Serialized
learning and remaining positive-hit broad scans remain possible follow-up costs;
no additional scheduling or search-policy changes were bundled into this correction.
The full 50.68M-read LARRY/P02 paper benchmark has not been rerun.

## Build and reproduction

Clean commands: `make -C core/features/process_features clean`,
`make -C core/legacy/source clean`, `make -C core/legacy/source -j8 STAR`.
Frozen binary: `larry_fix/build2/STAR`, SHA256
`84ed32dc8130a325f9292c4dce46427078dd70013d9f3d604637a5e6ba7b4fa4`.
Frozen assignment source SHA256:
`a26ce9f5d627149d6d9b41d8d15cb5082c4f44734554f15aee1d4c1e00baf4b7`.
Build logs, source patches, libraries, commands, environments, hashes and exit statuses
are retained with each run. `tests/run_bootstrap_anchor_regression.py` is the checked-in
differential driver; its baseline directory must contain the preserved PF and libscrna
archives. Reuse recorded controls when validating a later binary.

Full A375 integration passed in `hierarchy/a375_larry_fix_fixed4`: 146.22 s versus
145.47 s for the preserved flat-hash control. All six raw/filtered/combined MEX matrices
and all three guide-call CSVs are exactly equal. Comparison: `larry_fix/a375_comparison.json`.
The benchmark uses native BGZF through auto detection, so both compression paths have
completed output validation.

## Ordered follow-up

After committing these corrections, inspect
`/mnt/pikachu/STAR-suite/docs/AUDIT_ALLOCATION_PATTERNS_20260911.md` and implement its
simple memory fixes with independent validation. Preserve this LARRY binary and evidence
so subsequent allocation changes can be compared without conflating the causes.
