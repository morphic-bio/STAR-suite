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

## Remaining Gap

The original full filtered-cell merge path note is still only partially closed.

What is now covered:
- strict namespace normalization behavior
- prepared-context chemistry precedence

What is still not isolated by a small regression fixture:
- a full filtered-cell E2E mixed-library run that exercises the complete
  post-cell-calling merge path with real filtered GEX barcodes

That remaining gap is small enough to defer until there is a good reason to
add another heavier fixture-driven smoke.
