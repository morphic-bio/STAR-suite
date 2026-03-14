# Docker Validation Status

Recent portability checks (Ubuntu 24.04 Docker, `master` branch):

- `docker build --no-cache --target suite-base -f docker/Dockerfile --build-arg MAKE_JOBS=8 .` passed.
- `scripts/docker/run_smokes_tier_a.sh` passed.
- `scripts/docker/run_smokes_tier_b.sh` with `STORAGE=/storage` produced:
  - pass: `run_flex_smoke`, `run_cbub_regression_test`
  - skip: `run_cr_compat_integration_smoke` (ASAN-only path), `run_cr_multi_smoke` (external script dependency), `run_unsorted_cbub_smoke_test` (missing probes fixture)
  - fail: `test_cr_compat_crispr_calling` (missing expected `outs/*` files in this containerized fixture run)
- Runtime split verified:
  - `suite-base` and `test-tier-a` have no Python runtime
  - `test-tier-b` includes Python runtime for fixture-backed tests
- `suite-base` contains: `STAR`, `star_feature_call`, `slam_requant`, `pileup_snp`, `flexfilter` (`run_flexfilter_mex` alias), and `remove_y_reads`.
