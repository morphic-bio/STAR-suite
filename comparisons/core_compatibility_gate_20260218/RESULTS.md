# Core Compatibility Gate (2026-02-18)

Branch: `core-compatibility-fixes`  
HEAD at run time: `c30a106`  
Build hygiene: clean rebuild before gate (`make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR star_feature_call`)

## Gate 1: A375 GEX + Feature (GeneFull) parity

Script:
- `tests/run_a375_gex_features_cr_parity_genefull.sh`

Artifacts:
- STAR output: `/storage/A375/star_gex_features_cr_parity_genefull_20260218_013922/`
- Run log: `/storage/A375/compat_gate_a375_20260218_013922.log`
- Comparison report: `/storage/A375/star_gex_features_cr_parity_genefull_20260218_013922/comparison_report.txt`

Result summary:
- Common harmonized barcodes: `1163`
- Common genes: `38606`
- Filtered genes used for correlation: `15522`
- Pearson: `0.942599`
- Spearman: `0.951044`
- End-to-end runtime: `real 318.60s`

## Gate 2: UCSF 2M CRISPR call-only parity

Script:
- `comparisons/ucsf_ipsc2_callonly_gmm_parity_20260217/run_callonly_gmm_parity.sh`

Artifacts:
- Output base: `/storage/ucsf-2M/compat_gate_ucsf_callonly_20260218_014137/`
- Parity TSVs:
  - `/storage/ucsf-2M/compat_gate_ucsf_callonly_20260218_014137/parity_minumi10.tsv`
  - `/storage/ucsf-2M/compat_gate_ucsf_callonly_20260218_014137/parity_minumi3.tsv`
- Correlation summary:
  - `/storage/ucsf-2M/compat_gate_ucsf_callonly_20260218_014137/umi_correlation_summary.tsv`

Result summary:

`--min-umi 3`:
- Rows: CR `5767`, STAR `7346`, common `5767`
- Set-equivalent calls: `5766 / 5767` (`99.9827%`)
- CR total UMIs (common rows): `1,447,384`
- STAR total UMIs (common rows): `1,447,391` (delta `+7`, `+0.000484%`)
- Pearson: `0.9999998779`
- Spearman: `0.9999999409`
- Runtime: `real 5.74s`

`--min-umi 10`:
- Rows: CR `5767`, STAR `7346`, common `5767`
- Set-equivalent calls: `5457 / 5767` (`94.6246%`)
- CR total UMIs (common rows): `1,447,384`
- STAR total UMIs (common rows): `1,445,234` (delta `-2,150`, `-0.148544%`)
- Pearson: `0.9999627946`
- Spearman: `0.9999400299`
- Runtime: `real 11.17s`

## Compatibility status tracked in this gate

- 2-column whitelist/chemistry logic: implemented, exercised in pf-multi path.
- Chemistry precedence (`crChemistry`, `crOutputChemistry`): implemented and documented.
- Remaining explicit test gap: add dedicated automated matrix tests for 2-column parser edge cases and chemistry precedence combinations.
