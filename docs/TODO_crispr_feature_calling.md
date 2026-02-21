# CRISPR Feature Calling TODOs

## Regression Testing Required

### `--crMinUmi` Parameter Validation

The `--crMinUmi` parameter (default: 3) needs regression testing across different assay types:

- [ ] **CRISPR Guide Capture** (general default target)
  - A375 dataset: ✓ Validated at min_umi=10 (fixture-specific override, 100% CR concordance)
  - Additional CRISPR datasets needed

- [ ] **Lineage Barcodes / Stable Features**
  - Lower threshold likely better (3 or even 2)
  - These have much less noise than CRISPR guides
  - Need to validate with real lineage barcode datasets
  - Compare calling rates at different thresholds

- [ ] **FLEX Probe Sets**
  - Current recommendation (10) may be appropriate
  - Needs validation with FLEX datasets
  - Document recommended settings

### Test Matrix

| Assay Type | Recommended `--crMinUmi` | Status |
|------------|--------------------------|--------|
| CRISPR Guide Capture (general) | 3 | Baseline default |
| A375 CR-parity fixture | 10 | ✓ Validated |
| Lineage Barcodes | 2-3 | **Needs Testing** |
| FLEX Probes | 10 (TBD) | **Needs Testing** |
| Antibody Capture | TBD | **Needs Testing** |

### Action Items

1. Collect lineage barcode test datasets
2. Run feature calling at thresholds 2, 3, 5, 10, 20
3. Compare sensitivity/specificity tradeoffs
4. Update default or provide assay-specific guidance
5. Add automated regression tests

## Notes

- The `--min-umi` default in `star_feature_call` is 3
- In CR-compat mode, `--crMinUmi` defaults to 3
- A375 parity workflows intentionally override to 10
- Users can tune lower/higher values by assay characteristics
