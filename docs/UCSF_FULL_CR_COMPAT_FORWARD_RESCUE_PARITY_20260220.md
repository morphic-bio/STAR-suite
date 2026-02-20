# UCSF Full Sample: CR-Compat Forward+Rescue Parity Report (2026-02-20)

## Scope
Re-ran the UCSF full sample STAR pipeline with:
- `--soloStrand Forward`
- `--soloCrMultimapRescue yes`

Then compared STAR vs Cell Ranger using the existing parity script:
- `scripts/run_gex_feature_parity_checks.sh`

This report is a side-by-side comparison against the prior baseline STAR full run (`Unstranded`, no rescue flag).

## CR Comparison (Primary View)
Direct comparison of the new STAR run (`Forward + rescue`) against Cell Ranger.

### GEX Count-Level Parity
| Subset | CR total counts | STAR total counts | Common barcodes | Delta STAR-CR on common |
|---|---:|---:|---:|---:|
| `raw_all` | 1,957,055 | 1,116,382 | 29,292 | -705,637 |
| `raw_restricted_to_cr_filtered` | 1,350,854 | 502,437 | 6,794 | -789,953 |
| `filtered_vs_filtered` | 1,553,777 | 771,580 | 527 | 458 |

### GEX Correlations vs Cell Ranger
| Metric subset | Cell UMI Pearson | Cell UMI Spearman | Gene Pearson (all) | Gene Spearman (all) | Gene Pearson (filtered set) | Gene Spearman (filtered set) |
|---|---:|---:|---:|---:|---:|---:|
| `raw_all` | 0.057608 | 0.578707 | 0.310589 | 0.930840 | 0.280801 | 0.911062 |
| `raw_restricted_to_cr_filtered` | 0.096937 | 0.220332 | 0.298860 | 0.845007 | 0.246938 | 0.722872 |
| `filtered_vs_filtered` | 0.015751 | -0.075826 | 0.331603 | 0.874777 | 0.298094 | 0.821344 |

### Feature (CRISPR Guide Capture) Count-Level Parity
| Subset | CR total counts | STAR total counts | Common barcodes | Delta STAR-CR on common |
|---|---:|---:|---:|---:|
| `feature_raw_all` | 19,773,321 | 6,727,054 | 2,676 | -3,083,753 |
| `feature_raw_restricted_to_cr_filtered` | 16,273,843 | 3,052,450 | 823 | -2,965,581 |
| `feature_filtered_vs_filtered` | 18,865,033 | 6,304,124 | 101 | -340,754 |

### Feature Call Parity vs Cell Ranger
| Metric | Value |
|---|---:|
| CR called rows | 6119 |
| STAR called rows | 2802 |
| Common rows (all) | 344 |
| Common rows evaluated (STAR non-None) | 2 |
| Exact `num_features` match | 100.0000% |
| STAR total assigned UMIs | 6,291,228 |
| CR total assigned UMIs | 18,826,049 |
| Delta STAR-CR assigned UMIs | -12,534,821 |

## Inputs
- Cell Ranger reference run:
  - `/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804`
- Baseline STAR run:
  - `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_dynamic_off_full_20260218_203459`
- New STAR run (Forward + rescue):
  - `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_compat_on_full_forward_rerun_20260220_070445`

## Output Artifacts
- New run parity output:
  - `/tmp/ucsf_full_compat_parity_20260220_073620/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
  - `/tmp/ucsf_full_compat_parity_20260220_073620/FILTERED_BARCODE_SET_OVERLAP.txt`
- Baseline parity output:
  - `/tmp/ucsf_full_baseline_parity_20260220_073710/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
  - `/tmp/ucsf_full_baseline_parity_20260220_073710/FILTERED_BARCODE_SET_OVERLAP.txt`

## Run Confirmation
From new run logs:
- `soloStrand = Forward`
- `soloCrMultimapRescue = yes`
- Input reads: `233,828,574`

New run `Log.final.out` includes CR-compat rescue counters, confirming rescue path execution.

## STAR Mapping Summary
| Metric | Baseline (Unstranded, no rescue) | New (Forward + rescue) |
|---|---:|---:|
| Uniquely mapped reads % | 9.58% | 10.71% |
| % reads mapped to multiple loci | 2.43% | 1.30% |

## GEX Per-Barcode UMI Correlations
| Comparison subset | Baseline Pearson | New Pearson | Baseline Spearman | New Spearman |
|---|---:|---:|---:|---:|
| `raw_all` | 0.048644 | 0.057608 | 0.565255 | 0.578707 |
| `raw_restricted_to_cr_filtered` | 0.093567 | 0.096937 | 0.220249 | 0.220332 |
| `filtered_vs_filtered` | 0.007570 | 0.015751 | -0.073161 | -0.075826 |

## GEX Per-Gene Correlations (Collapsed Gene Totals)
| Comparison subset | Pearson (all genes) baseline | Pearson (all genes) new | Pearson (filtered genes) baseline | Pearson (filtered genes) new | Filtered genes baseline | Filtered genes new |
|---|---:|---:|---:|---:|---:|---:|
| `raw_all` | 0.055906 | 0.310589 | 0.109917 | 0.280801 | 206 | 219 |
| `raw_restricted_to_cr_filtered` | 0.099866 | 0.298860 | 0.201671 | 0.246938 | 303 | 337 |
| `filtered_vs_filtered` | 0.125865 | 0.331603 | 0.211446 | 0.298094 | 795 | 845 |

## Barcode Overlap (Filtered GEX Sets)
| Metric | Baseline | New |
|---|---:|---:|
| STAR filtered barcodes | 2746 | 2802 |
| CR filtered barcodes | 7325 | 7325 |
| Common barcodes | 520 | 527 |
| Jaccard | 0.054445 | 0.054896 |

## Feature (CRISPR Guide Capture) Count-Level Parity
| Comparison subset | STAR total counts baseline | STAR total counts new | Common barcodes baseline | Common barcodes new | Delta STAR-CR baseline | Delta STAR-CR new |
|---|---:|---:|---:|---:|---:|---:|
| `feature_raw_all` | 6,725,521 | 6,727,054 | 2,701 | 2,676 | -3,107,507 | -3,083,753 |
| `feature_raw_restricted_to_cr_filtered` | 3,052,452 | 3,052,450 | 829 | 823 | -2,989,451 | -2,965,581 |
| `feature_filtered_vs_filtered` | 6,298,889 | 6,304,124 | 97 | 101 | -324,663 | -340,754 |

## Feature-Call Parity Summary
| Metric | Baseline | New |
|---|---:|---:|
| `rows_cr` | 6119 | 6119 |
| `rows_star` | 2746 | 2802 |
| `rows_star_non_none` | 2022 | 2032 |
| `rows_star_none` | 724 | 770 |
| `common_rows_all` | 338 | 344 |
| `common_rows_eval` | 2 | 2 |
| `num_features_exact_match_pct` | 100.0000 | 100.0000 |
| `total_assigned_umis_star` | 6,286,261 | 6,291,228 |
| `total_assigned_umis_delta_star_minus_cr` | -12,539,788 | -12,534,821 |

## CR-Compat Rescue Counters (New Run)
From `Log.final.out`:
- Multimap reads entering rescue block: `5,685,071`
- Rescued by gene-vs-non-gene fast path: `2,341,213`
- Rescued by exonic winner (Phase 3): `298,358`
- Rescued by intronic fallback (Phase 3): `7,752`
- No rescue: multiple exonic loci: `213,322`
- No rescue: multiple intronic loci: `1,731,404`
- No rescue: all intergenic: `1,087,853`

## Bottom Line
- The new run is materially different from baseline and did execute compatibility rescue.
- On this UCSF full sample, gene-level collapsed-expression correlations improved substantially in all three barcode subsets.
- Cell-level UMI correlations changed only modestly.
- Feature parity changed slightly, with mixed direction depending on subset.
