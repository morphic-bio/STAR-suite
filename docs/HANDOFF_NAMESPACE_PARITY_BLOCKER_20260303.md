# Handoff: Namespace/Parity Blocker (2026-03-03)

## Scope
This handoff summarizes the work completed in this session and the current blocker for the next agent.

Primary intent was:
1. add unit/component tests before integration checks,
2. run UCSF integration parity checks (100k and 2M),
3. assess whether current namespace fixes hold end-to-end.

---

## 1) What Was Implemented

### New tests added
- `core/features/process_features/tests/test_filtered_barcode_namespace.c`
- `core/features/process_features/tests/test_pf_api_context_isolation.c`
- `tests/autodetect_nxt_tru/test_namespace_unit_fixes.sh`
- `core/features/process_features/Makefile` updated to build/run the new process_features tests.

### Unit/component test status
All new tests passed in this workspace:
- `test_filtered_barcode_namespace`: `5/5`
- `test_pf_api_context_isolation`: pass
- `test_namespace_unit_fixes.sh`: `4 passed, 0 failed`

---

## 2) Integration/Parity Runs Performed

## 2.1 100k run (current code)
STAR run:
- `/storage/ucsf-2M/star_runs/nxt_tru_namespace_fix_parity_20260302/ucsf100k_fixed_20260302_221326`

Timing:
- `real 179.34s` (from `stderr_time.log`)

Parity vs current-head 100k fixture:
- report: `/tmp/ucsf100k_fixed_vs_currenthead_20260303_041819/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
- overlap: `/tmp/ucsf100k_fixed_vs_currenthead_20260303_041819/FILTERED_BARCODE_SET_OVERLAP.txt`
- key metrics:
  - `rows_star_non_none: 163`
  - `set_equivalent_calls_pct: 77.3333`
  - `total_assigned_umis_star: 4431`
  - `total_assigned_umis_cr: 63924`
  - GEX filtered Jaccard: `0.913702`
  - Feature filtered Jaccard: `0.013000`

Interpretation: GEX barcodes are mostly aligned, but feature assignment is severely collapsed.

## 2.2 2M run (current code, TRU-configured)
STAR run:
- `/storage/ucsf-2M/star_runs/nxt_tru_namespace_fix_parity_20260302/ucsf2m_dynamic_fixed_20260302_221710`

Timing:
- `real 131.00s`

Parity vs CR baseline:
- report: `/tmp/ucsf2m_fixed_vs_cr_norm_20260303_041513/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
- overlap: `/tmp/ucsf2m_fixed_vs_cr_norm_20260303_041513/FILTERED_BARCODE_SET_OVERLAP.txt`
- key metrics:
  - `rows_star_non_none: 1818`
  - `set_equivalent_calls_pct: 99.3007`
  - `total_assigned_umis_star: 470163`
  - `total_assigned_umis_cr: 1447384`
  - filtered Jaccard: `0.048061`

Interpretation: high set-equivalence only on evaluated subset, but severe coverage collapse (too few STAR non-None rows and UMIs).

## 2.3 2M fresh rerun on current code (NXT-configured, like-for-like control chemistry)
STAR run:
- `/storage/ucsf-2M/star_runs/nxt_tru_namespace_fix_parity_20260302/ucsf2m_dynamic_fixed_nxt_20260303_042303`

Timing:
- `real 203.99s`

Parity vs CR baseline:
- report: `/tmp/ucsf2m_fixed_nxt_vs_cr_norm_20260303_042633/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
- overlap: `/tmp/ucsf2m_fixed_nxt_vs_cr_norm_20260303_042633/FILTERED_BARCODE_SET_OVERLAP.txt`
- key metrics:
  - `rows_star_non_none: 55`
  - `set_equivalent_calls_pct: 100.0000`
  - `total_assigned_umis_star: 17293`
  - `total_assigned_umis_cr: 1447384`
  - GEX filtered Jaccard: `0.616014`
  - Feature filtered Jaccard: `0.794251`

Interpretation: barcode overlap is not catastrophic, but feature calls collapse to almost all `None`.

## 2.4 Control parity proving tooling and old-good path still behave
Control (known good) STAR vs CR:
- report: `/tmp/ucsf2m_control_oldstar_vs_cr_norm_20260303_041634/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`
- overlap: `/tmp/ucsf2m_control_oldstar_vs_cr_norm_20260303_041634/FILTERED_BARCODE_SET_OVERLAP.txt`
- key metrics:
  - `rows_star_non_none: 5682`
  - `set_equivalent_calls_pct: 99.4226`
  - `total_assigned_umis_star: 1451614`
  - `total_assigned_umis_cr: 1447384`
  - filtered Jaccard: `0.631048`

Interpretation: parity script and dataset pairing are valid; regression is in current integration behavior.

---

## 3) Current Problem Statement

Current integration path shows a major collapse in feature non-`None` calls and assigned UMIs (100k and 2M), while control old-good run on the same CR baseline does not.

This is **not** a unit-test failure; it is an end-to-end integration failure.

Important caveat for reading parity reports:
- `set_equivalent_calls_pct` can appear high because evaluation excludes `STAR=None` rows (`call_parity_mode: exclude_star_none_default`).
- So this metric must be interpreted together with:
  - `rows_star_non_none`
  - `total_assigned_umis_star`
  - filtered barcode overlap.

---

## 4) Likely Failure Surface

Given observed behavior and recent namespace work, the likely surface is still around namespace/chemistry consistency across:
1. Solo filtered barcode emission namespace,
2. filtered barcode normalization before `pf_load_filtered_barcodes`,
3. assignment whitelist namespace / translation decision path,
4. final call-table join/evaluation namespace.

Key files to diff/trace first:
- `core/legacy/source/ParametersSolo.cpp`
- `core/legacy/source/PfMultiAssign.cpp`
- `core/legacy/source/PfMultiProcess.cpp`
- `core/features/process_features/src/pf_api.c`
- `core/features/process_features/src/io.c`
- `core/features/process_features/src/assignBarcodes.c`

---

## 5) Suggested Next Steps for Next Agent

1. Use the control report as hard reference (`/tmp/ucsf2m_control_oldstar_vs_cr_norm_20260303_041634/...`).
2. Diff `assignBarcodes.api_run.txt` between:
   - control-like good run,
   - failing current NXT run (`ucsf2m_dynamic_fixed_nxt_20260303_042303`).
3. Extract barcodes that are non-`None` in control but `None` in current run; trace through filtered barcode normalization and assignment lookup.
4. Verify namespace of:
   - `outs/filtered_feature_bc_matrix/barcodes.tsv.gz`,
   - `outs/crispr_analysis/protospacer_calls_per_cell.csv`,
   - any intermediate filtered barcode file fed to `pf_api`.
5. Re-run parity after fix on:
   - 100k fixture,
   - 2M fixture,
   and compare against same control metrics above.

---

## 6) Repro Command Snippets

Parity invocation used in this session (example):
```bash
bash scripts/run_gex_feature_parity_checks.sh \
  --star-run <STAR_RUN_DIR> \
  --cr-run <CR_RUN_DIR> \
  --out-dir <OUT_DIR> \
  --translation-direction left-to-right \
  --translate-side star
```

Fresh 2M NXT STAR run command used in this session:
- captured under: `/storage/ucsf-2M/star_runs/nxt_tru_namespace_fix_parity_20260302/ucsf2m_dynamic_fixed_nxt_20260303_042303/Log.out`

