# UCSF Full STAR-vs-CR Comparison Audit (2026-02-20)

## Purpose
Summarize why current CR-vs-new-compatibility comparison outputs appear inconsistent.

## STAR Run Under Test
- `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_compat_on_full_forward_rerun_20260220_070445`
- Config: `--soloStrand Forward --soloCrMultimapRescue yes`

## CR Run Under Test
- `/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804`

## Comparison Modes Re-run
All runs used `scripts/run_gex_feature_parity_checks.sh`:

1. `default` (barcode normalization on, `translate_side=both`)
   - `/tmp/ucsf_full_compat_parity_20260220_073620/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
2. `no normalization`
   - `/tmp/ucsf_full_compat_parity_nonorm_20260220_075351/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
3. `translate STAR only`
   - `/tmp/ucsf_full_compat_parity_translate_star_20260220_075351/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
4. `translate CR only`
   - `/tmp/ucsf_full_compat_parity_translate_cr_20260220_075441/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`

## Key Comparison Table (CR vs New Compatibility Run)
| Mode | Raw common barcodes (GEX) | Filtered common barcodes (GEX corr block) | Filtered cell-level Pearson | Filtered gene-level Pearson (filtered genes) | Feature-call common_rows_eval | Feature-call exact match % |
|---|---:|---:|---:|---:|---:|---:|
| `default (both)` | 29292 | 527 | 0.015751 | 0.298094 | 2 | 100.0000 |
| `no normalization` | 29292 | 527 | 0.015751 | 0.298094 | 2 | 100.0000 |
| `translate STAR only` | 53083 | 1526 | 0.690069 | 0.306197 | 965 | 99.2746 |
| `translate CR only` | 53083 | 1526 | 0.690069 | 0.306197 | 965 | 99.2746 |

## Concrete Inconsistencies
1. Metrics are highly mode-sensitive:
- `filtered_vs_filtered` cell Pearson moves from `0.015751` to `0.690069` depending only on translation mode.
- Feature-call `common_rows_eval` moves from `2` to `965` with the same STAR/CR data.

2. Default mode equals no-normalization exactly:
- `default (translate both)` and `no normalization` produce identical key metrics, suggesting translation is effectively canceling out or not stabilizing domains as intended.

3. Internal report mismatch in single-side translation runs:
- In `translate STAR only` (and `translate CR only`), GEX correlation block reports `filtered_vs_filtered common_barcodes = 1526`.
- But companion filtered-overlap report still shows `common = 527`.
- This indicates mismatch between barcode normalization logic used by correlation code vs overlap code.

4. Script emits translation-domain ambiguity warning:
- `translation_warning: source/target domains overlap observed barcodes; single-side translation can shift overlap/correlation metrics.`
- Observed overlap counts are non-trivial (`translation_observed_raw_cr_in_both_domains: 33034`, `translation_observed_raw_star_in_both_domains: 37421`), confirming ambiguous mapping domain during normalization.

## Practical Conclusion
Current CR-vs-compatibility comparison outputs are not numerically stable under barcode normalization choices; therefore, current parity/correlation values should be treated as provisional.

## Recommended Fix Direction (for coding agent)
1. Unify barcode normalization into a single canonical function used by:
- MEX comparison
- Additional metrics/correlation code
- Filtered barcode overlap report
- Feature-call parity merge

2. Add explicit domain detection and enforce one canonical output domain, with:
- collision count
- number of barcodes remapped
- number dropped due ambiguity

3. Fail fast (or require explicit override) when source/target domains overlap observed barcodes enough to change overlap sets materially.

4. Add a consistency check:
- assert that `filtered_vs_filtered common_barcodes` in correlation block equals filtered-overlap common set when using same normalization options.
