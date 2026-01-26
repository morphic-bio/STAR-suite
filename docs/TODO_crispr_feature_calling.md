# CRISPR Feature Calling TODOs

## Regression Testing Required

### `--crMinUmi` Parameter Validation

The new `--crMinUmi` parameter (default: 10) needs regression testing across different assay types:

- [ ] **CRISPR Guide Capture** (current default target)
  - A375 dataset: ✓ Validated at min_umi=10 (100% CR concordance)
  - Additional CRISPR datasets needed

- [ ] **Lineage Barcodes / Stable Features**
  - Lower threshold likely better (3 or even 2)
  - These have much less noise than CRISPR guides
  - Need to validate with real lineage barcode datasets
  - Compare calling rates at different thresholds

- [ ] **FLEX Probe Sets**
  - Current default (10) may be appropriate
  - Needs validation with FLEX datasets
  - Document recommended settings

### Test Matrix

| Assay Type | Recommended `--crMinUmi` | Status |
|------------|--------------------------|--------|
| CRISPR Guide Capture | 10 | ✓ Validated |
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

- The `--min-umi` default in `star_feature_call` was changed from 3 to 10 for CR compatibility
- For CR-compat mode, `--crMinUmi` defaults to 10
- Users can override with lower values for stable features like lineage barcodes
