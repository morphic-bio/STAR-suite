# Mixed-Chemistry Hardening (2026-03-11)

## Summary

Mixed-chemistry perturb support was already implemented and E2E-validated on the
MSK mixed-library case. The remaining work was regression hardening, not new
behavior.

This pass tightened the lightweight regression surface in two places:

1. Refreshed the namespace/filtered-barcode unit test to match the current
   strict namespace policy.
2. Added a new prepared-context regression test for the actual chemistry
   precedence logic in `PfMultiProcess.cpp`.
3. Added a small filtered-cell mixed-library smoke that exercises the full
   CR-compat merged-filtered-MEX path against real Solo filtered barcodes.

## Tests Added / Updated

### Updated

- `tests/autodetect_nxt_tru/test_namespace_unit_fixes.sh`
  - Removed the obsolete heuristic-membership fallback expectation.
  - Added assertions that missing namespace arguments hard-fail.
  - Kept deterministic overlap normalization coverage under explicit
    `sourceNamespace` / `whitelistNamespace`.

### Added

- `tests/multi_feature/test_mixed_chemistry_prepared_context.sh`
  - Builds a small harness around `buildPfMultiPreparedContext()`.
  - Covers:
    - GEX `star_chemistry` re-anchoring of `effectiveChem` and `inferredChem`
    - per-library override isolation in a multi-library config
    - distinction between empty `star_chemistry` and explicit `auto`

- `tests/multi_feature/run_mixed_chemistry_filtered_smoke.sh`
  - Builds or reuses the downsampled MSK mixed-library fixture.
  - Runs the full pf-multi / CR-compat path with filtered cell calling enabled.
  - Verifies that:
    - `outs/filtered_feature_bc_matrix/barcodes.tsv.gz` matches
      `Solo.out/GeneFull/filtered/barcodes.tsv` after suffix normalization
    - feature-library filtered barcode outputs are subsets of the merged
      filtered barcode universe after per-library chemistry normalization
    - `Log.out` records the Solo-sourced filtered-barcode path and the
      normalized GEX barcode count used for CR-compat merge filtering

## Validated Suite

The following scripts were run after a clean rebuild:

- `tests/autodetect_nxt_tru/test_autodetect.sh`
- `tests/autodetect_nxt_tru/test_star_chemistry_column.sh`
- `tests/autodetect_nxt_tru/test_two_column_wl_pairing.sh`
- `tests/autodetect_nxt_tru/test_namespace_unit_fixes.sh`
- `tests/multi_feature/test_multi_feature_config.sh`
- `tests/multi_feature/test_no_global_ref_guard.sh`
- `tests/multi_feature/test_data_driven_specs.sh`
- `tests/multi_feature/test_mixed_chemistry_prepared_context.sh`
- `tests/multi_feature/run_mixed_chemistry_filtered_smoke.sh`

## Remaining Gap

The original filtered-cell merge-path gap is now closed at smoke-test level.

What is now covered:
- strict namespace normalization behavior
- prepared-context chemistry precedence
- full mixed-library filtered-cell merge behavior on a downsampled real MSK
  fixture, including the `Solo -> filtered_feature_bc_matrix` handoff

What is still intentionally deferred:
- larger fixture or production-scale reruns beyond the smoke surface
