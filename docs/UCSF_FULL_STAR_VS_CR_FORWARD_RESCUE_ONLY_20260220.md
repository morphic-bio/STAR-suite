# UCSF Full Sample: STAR (Forward+Rescue) vs Cell Ranger Only (2026-02-20)

## Inputs
- STAR run:
  - `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_compat_on_full_forward_rerun_20260220_070445`
- Cell Ranger run:
  - `/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804`
- Parity report source:
  - `/tmp/ucsf_full_compat_parity_20260220_073620/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`

## Run Configuration (STAR)
- `--soloStrand Forward`
- `--soloCrMultimapRescue yes`

## GEX Count-Level Parity vs Cell Ranger
| Subset | CR total counts | STAR total counts | Common barcodes | Delta STAR-CR on common |
|---|---:|---:|---:|---:|
| `raw_all` | 1,957,055 | 1,116,382 | 29,292 | -705,637 |
| `raw_restricted_to_cr_filtered` | 1,350,854 | 502,437 | 6,794 | -789,953 |
| `filtered_vs_filtered` | 1,553,777 | 771,580 | 527 | 458 |

## GEX Correlations vs Cell Ranger
| Metric subset | Cell UMI Pearson | Cell UMI Spearman | Gene Pearson (all) | Gene Spearman (all) | Gene Pearson (filtered set) | Gene Spearman (filtered set) |
|---|---:|---:|---:|---:|---:|---:|
| `raw_all` | 0.057608 | 0.578707 | 0.310589 | 0.930840 | 0.280801 | 0.911062 |
| `raw_restricted_to_cr_filtered` | 0.096937 | 0.220332 | 0.298860 | 0.845007 | 0.246938 | 0.722872 |
| `filtered_vs_filtered` | 0.015751 | -0.075826 | 0.331603 | 0.874777 | 0.298094 | 0.821344 |

## Feature (CRISPR Guide Capture) Count-Level Parity vs Cell Ranger
| Subset | CR total counts | STAR total counts | Common barcodes | Delta STAR-CR on common |
|---|---:|---:|---:|---:|
| `feature_raw_all` | 19,773,321 | 6,727,054 | 2,676 | -3,083,753 |
| `feature_raw_restricted_to_cr_filtered` | 16,273,843 | 3,052,450 | 823 | -2,965,581 |
| `feature_filtered_vs_filtered` | 18,865,033 | 6,304,124 | 101 | -340,754 |

## Feature-Call Parity vs Cell Ranger
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

## Notes
- This document is intentionally CR-only (no baseline STAR comparison).
- The rescue path executed in STAR (`CR-COMPAT RESCUE` counters are present in `Log.final.out`).
