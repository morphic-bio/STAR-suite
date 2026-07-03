# MSK EmptyDrops Ambient Rescue Isolation (2026-04-02)

## Resolution / Status Update (2026-07-03)

**Verified status: RESOLVED for the stale "identified but not integrated"
claim.**

The guarded ambient policy from this investigation was integrated into the
production `libscrna` EmptyDrops path.

Landed code evidence:

- Commit `90ed49076bdf5cb5301c740f1a1770ed3e7b354e`
  (`2026-04-02`, "EmptyDrops: guarded ambient minimum, debug audit, and custom
  ambient CLI") changed `core/features/libscrna/src/OrdMagStage.cpp`.
- Both `SimpleEmptyDropsStage::runCRSimpleFilterBootstrap()` and
  `SimpleEmptyDropsStage::runCRSimpleFilter()` now compute a guarded
  `minAmbientCells` target: for `nCB >= 1000`, use
  `min(nCB/10, max(5000, nCB/50))`; for smaller inputs, use up to `100` cells.
  If the legacy rank window `[indMin, indMax)` is large enough, it is
  preserved; otherwise the fallback uses the bottom guarded pool.
- This is in the production path, not only the standalone tool:
  `core/legacy/source/SoloFeature_emptyDrops_libscrna.cpp` builds the
  `scrna_ed_config` and calls `scrna_emptydrops_run()`;
  `core/features/libscrna/src/scrna_api.cpp` calls `SimpleEmptyDropsStage` and
  builds the ambient profile from `simpleResult.ambientIndices`.
- Commit `f692f45ad41fbb035b7efc0b0795441b2bcb6d53` (`2026-03-03`) separately
  enabled bootstrap OrdMag behavior on shallow data by making the `indMin`
  skip gate apply only when bootstrap is disabled. This supports low-count /
  shallow inputs, while the ambient guard determines which ambient pool is used
  when MC rescue proceeds.

Validation / guard evidence:

- The MSK ambient-swap and guarded-policy experiments below are registered in
  `tests/ARTIFACTS.md` under `/tmp/msk_libscrna_ambientswap_*` and
  `/tmp/msk_ambient_policy_*`.
- The fresh guarded MSK full perturb rerun is tracked in
  `comparisons/msk_30polyko_full_benchmark_20260306/post_permits_20260403/`
  with run root
  `/storage/MSK-perturb-comparison/paper_bench_emptydrops_guarded_redo_20260403_214718`.
  The tracked parity report records STAR/CR cells `33,095 / 32,256`, barcode
  Jaccard `0.9742`, Gene Expression per-feature Pearson `0.994554`, and CRISPR
  set-equivalent calls `98.04%`.
- Automated tests cover related pieces (`tests/test_ordmag_bootstrap.cpp` for
  bootstrap OrdMag and `tests/emptydrops/test_simple_ed_threshold.cpp` for
  SimpleED thresholding), but there is no small unit test that directly asserts
  the guarded ambient-window selection. The strongest guard for this exact
  behavior is the tracked MSK benchmark/artifact set above.

Release-note check:

- `docs/RELEASE_NOTES_*` do not currently mention the guarded ambient fix by
  name. `v1.0.0` is the natural release-note bucket because the fix landed
  before the first production release.

## Question

Why does the legacy / CellGENI-style MSK cell list have a much higher CR9
barcode Jaccard than the current `libscrna` path, even after the OrdMag update?

The specific hypotheses tested here were:

- different recovered-cell / OrdMag behavior
- FDR vs raw p-value gating
- different ambient profile construction

UCSF was intentionally excluded. This note is MSK-only.

## Inputs

- Historical raw MEX:
  `/storage/MSK-perturb-comparison/bench_cellgeni_msk_20260401_162259/output/GeneFull/raw/`
- Instrumented STAR `soloCellFiltering` reruns:
  `/tmp/msk_ed_instrumented_20260402_215715/`
- Tool-only ambient swap experiments:
  `/tmp/msk_libscrna_ambientswap_20260402_220321/`
- CR9 reference barcodes:
  `/storage/MSK-perturb-comparison/cr_full_grna_30crispr_20260306_173247/cr9_gex_grna_30crispr/outs/per_sample_outs/cr9_gex_grna_30crispr/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz`

## Native STAR Reruns

These reruns used the same historical raw matrix and differed only in the
EmptyDrops backend.

| Run | Simple cells | Tail candidates | Ambient cells | ED rescues | Final cells |
| --- | ---: | ---: | ---: | ---: | ---: |
| legacy STAR backend | 28,527 | 4,029 | 45,000 | 3,777 | 32,304 |
| `libscrna` + legacy knee | 28,527 | 20,000 | 207,413 | 1,951 | 30,478 |

Immediate observations:

- The simple-cell set is identical: `28,527` barcodes in both runs.
- OrdMag is not the main driver on this high-cell MSK case.
- The current `libscrna` path is not testing the same surface as legacy:
  - legacy tail: `4,029` candidates
  - `libscrna` tail: `20,000` candidates
  - legacy ambient: ranks `45,000-89,999` (`45,000` cells)
  - `libscrna` ambient: `207,413` low-UMI cells from `SimpleED`

## Shared >500 UMI Tail Comparison

To avoid surface confounding, the two instrumented candidate tables were joined
by barcode on the shared `>500 UMI` tail.

- Shared `>500 UMI` tail candidates: `4,027`
- Legacy-only `>500 UMI` tail candidates: `0`
- Extra `libscrna` tail candidates with `UMI <= 500`: `15,973`

On the shared `>500 UMI` tail:

- legacy FDR pass, `libscrna` FDR pass: `1,864`
- legacy FDR pass, `libscrna` FDR fail: `1,912`
- legacy raw-p pass, `libscrna` raw-p pass: `2,030`
- legacy raw-p pass, `libscrna` raw-p fail: `1,783`

Interpretation:

- Switching `libscrna` from FDR to raw p-values would recover only `129`
  additional barcodes on the shared high-UMI tail.
- The main gap remains even after removing both the OrdMag and FDR hypotheses.

## Ambient Set Comparison

The instrumented ambient sets are completely disjoint on MSK.

- Legacy ambient set size: `45,000`
- `libscrna` ambient set size: `207,413`
- Ambient overlap: `0`

This immediately suggested that the dominant difference was ambient profile
construction rather than the Monte Carlo engine or the gate.

## Controlled Ambient Swap

To isolate that, `core/features/libscrna/tools/scrna_simpleed.cpp` was extended
with a tool-only experimental flag, `--use-legacy-rank-ambient`.

This experiment kept all of the following fixed:

- the same historical MSK raw matrix
- the same `SimpleED` threshold (`28,527` simple cells, threshold `2866`)
- the same `libscrna` Monte Carlo implementation
- the same FDR gate

It changed only the ambient pool:

- control: current `libscrna` `SimpleED` ambient set
- experiment: legacy rank window `45,000-89,999`

The tool control is very close to the current STAR `libscrna` run rather than
bit-for-bit identical:

- tool control vs STAR `libscrna` (`lib_legacyknee_auto`) Jaccard: `0.996492`
- tool-only barcodes: `20`
- STAR-only barcodes: `87`

### Results

| Run | Ambient source | ED rescues | Final cells | Jaccard vs legacy | Jaccard vs CR9 |
| --- | --- | ---: | ---: | ---: | ---: |
| `current_surface_fdr` | current `SimpleED` ambient | 1,884 | 30,411 | 0.941430 | 0.942681 |
| `legacy_rank_ambient_fdr` | legacy rank window `45k-90k` | 3,776 | 32,303 | 0.999969 | 0.997556 |
| legacy STAR backend | legacy rank window `45k-90k` | 3,777 | 32,304 | 1.000000 | 0.997525 |

This is the key result:

- swapping only the ambient pool moves `libscrna` from `30,411` to `32,303`
- that recovers `32,303 / 32,304` legacy barcodes
- it also reproduces the high CR9 Jaccard

## Residual 1-Barcode Difference

The single legacy-only barcode after the ambient swap is:

- `CTTGAGAAGTTCACTG`

In the legacy instrumented run it has:

- rank `32,554`
- `UMI=500`
- `p_value=0.000100`
- `p_adjusted=0.000109`

This barcode is consistent with a threshold-boundary effect (`500` exactly),
not a substantive backend mismatch.

## Guarded Ambient Policy Sweep

To test whether the current fallback is simply too aggressive, a tool-only
policy sweep was run on both:

- the full historical MSK raw matrix
- the MSK 100K artifact

The guarded policy keeps the legacy ambient window `[indMin, indMax)` when that
window is large enough, otherwise it falls back to the bottom `N` barcodes.

### Full MSK

Artifact root:

- `/tmp/msk_ambient_policy_20260402_221700/`

Policies tested:

- current-equivalent: `10% + abs 100`
- less aggressive fraction: `2% + abs 100`
- absolute minimum: `abs 5000`

Results:

| Policy | Ambient selection | Final cells | Jaccard vs current | Jaccard vs CR9 |
| --- | --- | ---: | ---: | ---: |
| `10% + abs 100` | bottom `207,413` cells | 30,411 | 1.000000 | 0.942681 |
| `2% + abs 100` | legacy rank window `45k-90k` | 32,303 | 0.941430 | 0.997556 |
| `abs 5000` | legacy rank window `45k-90k` | 32,303 | 0.941430 | 0.997556 |

Interpretation:

- `10% + abs 100` exactly reproduces the current behavior.
- On the full MSK dataset, both `2%` and `abs 5000` are sufficient to stop the
  fallback from overriding the `45k-90k` window.
- Once that happens, the result is identical to the legacy-window ambient
  experiment above.

### MSK 100K Artifact

Artifact roots:

- `/tmp/msk_ambient_policy_20260402_221700/`
- `10%` control recheck: `/tmp/msk_ambient_policy_recheck_20260402_221950/`

Input:

- `/storage/MSK-perturb-comparison/star_dynamic_msk30ko_larry_gex_100k_20260302_083717/Solo.out/GeneFull/raw/`

Observed nonzero barcode count after dropping zero-UMI cells:

- `29,926`

This dataset is below `indMin=45,000`, so the current STAR wrapper returns
OrdMag-only and does not run MC tail rescue.

Results:

| Policy | Ambient selection | Final cells | Jaccard vs archived 100K STAR | Jaccard vs 100K CR9 |
| --- | --- | ---: | ---: | ---: |
| current wrapper | OrdMag-only | 29,926 | 1.000000 | 0.977805 |
| `10% + abs 100` | bottom `2,992` cells | 29,926 | 1.000000 | 0.977805 |
| `2% + abs 100` | bottom `598` cells | 29,926 | 1.000000 | 0.977805 |
| `abs 5000` | bottom `5,000` cells | 29,926 | 1.000000 | 0.977805 |

Interpretation:

- On the 100K artifact, all tested guarded policies are barcode-identical to the
  archived STAR 100K result.
- This is because the dataset is shallow enough that the current wrapper stays
  in OrdMag-only mode; the ambient policy does not affect the final cell set on
  this artifact.

## Conclusion

For the MSK high-cell dataset:

- the gap is not primarily OrdMag / recovered-cells estimation
- the gap is not primarily FDR vs raw p-value gating
- the dominant driver is ambient profile construction

More specifically, the current scRNA `libscrna` path is using the `SimpleED`
ambient set, while the legacy STAR backend uses the classic rank-window ambient
pool (`45k-90k`). On MSK, changing only that ambient definition is sufficient to
recover the legacy / CellGENI-style rescue set and the corresponding CR9
barcode Jaccard.

The guarded policy sweep suggests a safer replacement for the current `10%`
fallback:

- preserve the legacy ambient window when available
- use a much smaller guard such as `2%` or an absolute minimum like `5000`
  instead of `10%`

On the full MSK dataset, both of those recover the desired ambient behavior. On
the MSK 100K artifact, neither changes the result at all.
