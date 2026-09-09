# Full 320K STAR all-feature caller integration — 2026-09-09

The fix is implemented in `/mnt/pikachu/STAR-suite` and validated with one full CBQ and one full BGZF cloud run. All eight samples improve against the matched 100K-simulation active-only control. Pooled cell Jaccard rises from 0.954180 to **0.973961**; CR-missed calls fall from 7,813 to **331**. The new deprecated counts are independently quantified by STAR from reads. No CR counts enter the STAR model or export.

## Cell concordance

All identities are full CB16+TAG8; only a terminal `-1` is removed. Every result below is after whole-library occupancy. CBQ and BGZF have exactly identical pre-occupancy and final callsets. CR results contain 325,410 cells across eight samples / sixteen tags.

| Sample | Active Jaccard | New STAR Jaccard | CR cells | STAR cells | Shared | CR-absent | CR-missed | Saved cyto Jaccard |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Glioblastoma_BC1-2 | 0.972255 | 0.974652 | 33,342 | 34,122 | 33,299 | 823 | 43 | 0.922388 |
| Colorectal_BC3-4 | 0.934325 | 0.970922 | 52,560 | 54,059 | 52,523 | 1,536 | 37 | 0.875132 |
| LungCancer2_BC5-6 | 0.970723 | 0.984363 | 69,330 | 70,393 | 69,311 | 1,082 | 19 | 0.940730 |
| BreastCancer1_BC7-8 | 0.961588 | 0.967407 | 15,134 | 15,591 | 15,108 | 483 | 26 | 0.696350 |
| LNReactive_BC9-10 | 0.862962 | 0.957777 | 31,668 | 32,876 | 31,576 | 1,300 | 92 | 0.744279 |
| Kidney_BC11-12 | 0.973544 | 0.976322 | 33,687 | 34,423 | 33,647 | 776 | 40 | 0.951631 |
| SkinMelanoma_BC13-14 | 0.969374 | 0.972206 | 45,139 | 46,326 | 45,088 | 1,238 | 51 | 0.924822 |
| Endo_BC15-16 | 0.970232 | 0.974930 | 44,550 | 45,649 | 44,527 | 1,122 | 23 | 0.883761 |
| POOLED | 0.954180 | 0.973961 | 325,410 | 333,439 | 325,079 | 8,360 | 331 | 0.882754 |

Precision is 97.4928% versus 97.7116% in the control; recall is 99.8983% versus 97.5990%. CR-absent calls increase 7,438 → 8,360. “CR-absent” describes reference disagreement, not established biological false positives.

Against the active control, the new calls gain 7,603 CR identities and lose 121; they gain 2,038 CR-absent identities and lose 1,116. Exact per-sample identity counts are in `identity_changes.json`. The prior automatic hybrid had Jaccard 0.973943; `comparison.tsv` also records exact native-STAR versus hybrid callset overlap. Similar pooled scores do not imply identical callsets.

The active and hybrid controls are the completed all-sample 100K-simulation replays in `../full320k_all_samples_deprecated_20260909/`. This avoids confounding the cell comparison with the earlier full-read benchmark’s 10K simulation depth. Cyto numbers reuse the completed full-320K cyto 0.4.7 callsets; cyto was not rerun. Tag suffixes are converted through the actual whitelist before comparison.

## Count concordance and audit

Expression correlations use shared called cells and the 18,129 included export genes. Gene correlations require ≥20 total UMIs and detection in ≥1% of shared cells in both matrices, using the canonical STAR parity helper. Correlations on the previous CR-only gene mask are also retained in the JSON/TSV. Cell correlation compares per-cell total UMIs; gene correlation compares per-gene total UMIs. Count-weighted Jaccard is sum(min(STAR, CR)) / sum(max(STAR, CR)) across cell–gene coordinates in shared cells.

| Sample | Cell UMI Pearson | Cell UMI Spearman | Gene UMI Pearson | Gene UMI Spearman | Count-weighted Jaccard |
| --- | --- | --- | --- | --- | --- |
| Glioblastoma_BC1-2 | 0.999992973 | 0.999978144 | 0.999998103 | 0.999984101 | 0.994601 |
| Colorectal_BC3-4 | 0.999981647 | 0.999995224 | 0.999999838 | 0.999996232 | 0.993281 |
| LungCancer2_BC5-6 | 0.999978424 | 0.999996653 | 0.999999805 | 0.999999106 | 0.993097 |
| BreastCancer1_BC7-8 | 0.999995579 | 0.999973439 | 0.999998954 | 0.999989153 | 0.995599 |
| LNReactive_BC9-10 | 0.999994372 | 0.999984671 | 0.999998339 | 0.999994154 | 0.995014 |
| Kidney_BC11-12 | 0.999995190 | 0.999987782 | 0.999999688 | 0.999993002 | 0.995426 |
| SkinMelanoma_BC13-14 | 0.999987603 | 0.999991351 | 0.999999690 | 0.999995029 | 0.994687 |
| Endo_BC15-16 | 0.999987746 | 0.999996076 | 0.999998682 | 0.999995828 | 0.994324 |
| POOLED | 0.999980266 | 0.999993486 | 0.999999663 | 0.999999740 | — |

Full CBQ/BGZF raw matrices have 7,924,268 barcodes, 19,068 features, 725,480,265 nonzero coordinates, and 1,310,330,528 UMIs. **Every coordinate is identical between formats**. All eight filtered exports exactly equal the included-feature subsets of the new raw matrix.

The original 18,526 active features have 628 changed coordinates versus the previous STAR full raw matrix; absolute UMI difference 639, with totals 1,283,927,962 → 1,283,927,323. Barcode-axis growth includes cells detected only through deprecated features. The read classifier, barcode correction and UMI-deduplication C++ code is unchanged; `source_changes_vs_previous_benchmark.json` records the source audit.

Independent deprecated-feature count audit against CR raw output HDF5:

| Sample | STAR deprecated UMIs | CR deprecated UMIs | STAR / CR | Gene-total Pearson | Different coordinates |
| --- | --- | --- | --- | --- | --- |
| Glioblastoma_BC1-2 | 571,917 | 568,369 | 1.006242 | 0.999969848 | 3,645 |
| Colorectal_BC3-4 | 9,397,591 | 9,345,235 | 1.005602 | 0.999999889 | 33,944 |
| LungCancer2_BC5-6 | 9,480,932 | 9,424,963 | 1.005938 | 0.999999988 | 36,516 |
| BreastCancer1_BC7-8 | 350,406 | 348,700 | 1.004892 | 0.999984443 | 1,733 |
| LNReactive_BC9-10 | 1,747,151 | 1,737,124 | 1.005772 | 0.999999961 | 6,357 |
| Kidney_BC11-12 | 782,461 | 778,419 | 1.005193 | 0.999986143 | 4,143 |
| SkinMelanoma_BC13-14 | 1,609,958 | 1,601,126 | 1.005516 | 0.999985042 | 9,300 |
| Endo_BC15-16 | 2,462,789 | 2,448,356 | 1.005895 | 0.999981737 | 15,074 |

The deprecated audit compares the full raw barcode union and 542 separate deprecated feature IDs. It measures count differences; no read-level BAM/sidecar comparison was performed in these count-only benchmarks.

## Caller stages and occupancy

| Sample | Primary | Candidates | Tail tested | Pre-occupancy | Removed | Final |
| --- | --- | --- | --- | --- | --- | --- |
| Glioblastoma_BC1-2 | 33,157 | 37,025 | 3,868 | 34,631 | 509 | 34,122 |
| Colorectal_BC3-4 | 7,652 | 71,871 | 64,219 | 54,843 | 784 | 54,059 |
| LungCancer2_BC5-6 | 5,130 | 76,855 | 71,725 | 71,062 | 669 | 70,393 |
| BreastCancer1_BC7-8 | 13,104 | 17,603 | 4,499 | 15,966 | 375 | 15,591 |
| LNReactive_BC9-10 | 5,797 | 44,041 | 38,244 | 33,284 | 408 | 32,876 |
| Kidney_BC11-12 | 27,874 | 37,670 | 9,796 | 34,979 | 556 | 34,423 |
| SkinMelanoma_BC13-14 | 24,533 | 59,197 | 34,664 | 46,959 | 633 | 46,326 |
| Endo_BC15-16 | 10,742 | 58,051 | 47,309 | 46,235 | 586 | 45,649 |

Occupancy uses all 337,959 pre-calls across sixteen tags: 129,844 occupied GEMs, mean observed tags 2.602807985, fitted λ=2.356086655, ordinary-Poisson 0.999 cutoff **8**. It rejects 474 GEMs and removes 4,520 calls; **0 removed calls match CR**. An independent zero-truncated Poisson calculation reproduces the production final sets exactly. The matched active control had λ=2.296315606, the same cutoff 8, and 3,915 removals.

Minimum included-export UMIs among final calls: 114. Across all eight samples there are 0 zero-UMI, 0 one-UMI, and 0 below-100-UMI exports. The 500-UMI floor applies to the model counts; export counts can be lower after excluding modeling-only genes.

## Implementation and validation

- Main STAR now integrates fixed ambient ranks, observed-tag occupancy after all samples, and an export-only `--soloFlexFilteredGeneList` allowlist. The existing primary-floor fix remains enforced.
- Tag-aware STAR and standalone defaults use 100,000 EmptyDrops simulations, BH FDR 0.01, and fixed ambient ranks 45,000–90,000 per tag. The legacy caller retains its 10K default.
- `flex/scripts/extend_model_probe_cache.py` extends a matching active cache from the public 2024-A probe CSV. It preserves old gene indices and adds 623 deprecated probe parents / 542 feature IDs. No new-key collisions occurred for this panel. The model has 19,068 features; filtered export has 18,129.
- The all-feature cache must be paired with its gene-ID list and `--soloRemoveDeprecated No`. Invocation and cache collision rules are documented in `docs/FLEX_MODEL_FEATURES.md`. Reusing an old active-only cache will not add the missing counts.
- A fresh isolated Chromap-enabled STAR build passed seven shared caller tests, the grouped fixture, exact internal/standalone full-call and candidate-ledger parity, and the model-cache unit fixture.
- One CBQ and one BGZF read fixture produce exact raw counts and eight intended cells across four real tags, with zero calls in twelve unused tags containing deprecated-only 1–3-UMI barcodes. Large grouped-caller fixture coverage also includes 90,500 such low-count barcodes. The full JAX deprecated-feature read benchmark was not repeated.
- Initial build validation needed a missing standalone helper after STAR compilation/tests had already passed. Initial CBQ read-fixture validation incorrectly assumed empty samples had no output directories; they contained gDNA metrics only. Both validation repairs preserve the original outputs/failure records. No successful STAR or caller execution was repeated.
- CR source was not inspected. Only public probe data, CR output matrices and saved reference calls were used. The excluded earlier source-exposed investigation was not read.

## Full benchmark and provenance

Cloud instance `i-06de289faa5d78117`, us-west-2. One execution per input format, serialized, cold caches, 48 threads, full caller diagnostics, 7,303,142,230 read pairs per arm. `--flexNoAlign 1` uses an empty genome directory; no STAR genome-index load, BAM or decision sidecar. The small read fixture additionally traces file opens to confirm no index files are accessed.

| Input | Wall seconds | Wall time | Peak RSS GiB | Pairs/s |
| --- | --- | --- | --- | --- |
| CBQ | 1217.150 | 20m 17.15s | 138.852 | 6,000,199 |
| BGZF | 1452.071 | 24m 12.07s | 139.944 | 5,029,466 |

These timings include the 100K caller and diagnostics. The earlier active-only full-read benchmark used 10K simulations, so its timings are not an isolated estimate of deprecated-feature overhead.

- Local artifacts: `/mnt/pikachu/star_suite_paper/analysis/full320k_star_deprecated_20260909`.
- Cloud artifacts: `/scratch/full320k_star_deprecated_20260909_v1`.
- S3 root: `s3://star-suite-320k-benchmark-alt-171440768238-us-west-2-20260904/analysis-tools/full320k_star_deprecated_20260909_v1/`.
- Binary SHA-256: `1b7fc3f961fbccae2ef3b84e67bdaec9c20b487b1d162974a764d7657bc622ec`. Local binary: `source/core/legacy/source/STAR`.
- The same tested binary is installed at `/mnt/pikachu/STAR-suite/core/legacy/source/STAR` (reports version 1.8.4). All 1,588 compiled source/configuration files match the main checkout. The previous executable is preserved at `before/main_STAR.binary`; see `installed_binary.json` and `final_source_verification.json`.
- Model-cache SHA-256: `c18e89e6ea6c6c5dc0c1375425af868512b8fb6df9a1722b4bd1bdd6f64a410d`.
- Verified compact result archive SHA-256: `16c612747689b44ebad11d46eb3e4a7da98bd356a63b95d58651c38a25141a8f`.
- `implementation.patch` records this task’s source changes over the inherited dirty checkout; source/build/bundle manifests and the uploaded source snapshot preserve exact provenance. Existing unrelated edits were preserved. No commit or merge was made.
- `BENCHMARK_COMPLETE.json`, `analysis/ANALYSIS_COMPLETE.json`, successful SSM records and `collection_verification.json` establish completion. Ten HDF5 matrices and the compressed probe cache are archived separately; their SHA-256/S3 locations are in the cache archive manifests.
- HDF5 matrices use `cells_by_genes_csr_v1` with full barcode IDs, sorted features, CSR data/indices/indptr, and original feature order. They can be reused without parsing the huge MEX files again.
