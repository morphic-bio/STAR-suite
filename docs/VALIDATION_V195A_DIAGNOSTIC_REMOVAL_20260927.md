# Corrected v1.9.5.a Validation

The optional source-derived Flex gDNA diagnostic was removed, not rederived.
See [withdrawal notice](RUNBOOK_FLEX_GDNA_DIAGNOSTIC_20260727.md).
No Cell Ranger source was consulted for this removal.

## Release Identity

The owner authorized reuse of the unfinished `v1.9.5.a` tag after canceling
publication. Its superseded tag object was
`c01b3d5c921c14cc7992b246d3e9654cc79b39bf`, pointing at
`91740391d48deb37237161f5aaef169f0a402149`.
The original `v1.9.5` tag is not changed. Historical tags/branches and previously
published development images are not erased by deleting the estimator.

Resolve the corrected release with `git rev-parse v1.9.5.a^{}` and compare
that full commit to `STAR --source-revision`. To refresh a previously fetched
copy of this specifically authorized tag:

```bash
git fetch origin +refs/tags/v1.9.5.a:refs/tags/v1.9.5.a
```

## Local Checks

Working tree: `/mnt/pikachu/STAR-suite-gdna-removal`.
Artifacts: `/tmp/star-flex-removal-20260927/`.
The clean default build used `make -C core/legacy/source clean`, then
`make -j8 core`; log: `/tmp/star-gdna-removal-build-20260927.log`.

```bash
OUT_ROOT=/tmp/star-flex-removal-20260927 \
  bash tests/run_flex_diagnostic_removal_acceptance.sh
make -C core/legacy/source -j8 test
python3 tests/test_parameters_default_generation.py
python3 tests/test_scrna_regression_report.py

TEST_WORKDIR=/tmp/star-flex-removal-20260927/jax-parity \
STAR_REF_BIN=/tmp/star-htslib-build-20260927/partial-full-v195a-clean-deps/chromap-core/source/core/legacy/source/STAR \
  bash tests/test_flex_v194_default_route.sh

BGZF_E2E_OUT_ROOT=/tmp/star-flex-removal-20260927/fused-align \
  bash tests/bgzf/test_flex_fused_align.sh

OUT_ROOT=/tmp/star-flex-removal-20260927/partial MAKE_JOBS=8 \
  bash tests/run_partial_make_regression.sh

python3 tests/run_scrna_gex_100k_regression.py \
  --outdir /tmp/star-flex-removal-20260927/pbmc100k --threads 8
```

The focused suite passed: removal/CLI guards, packed-count saturation,
cache v1/v2/v3 and khash/half-probe storage, nine HTSlib discovery tests, eight
analytical Solo count profiles, mapped-read Solo smoke, the CR-compatible
golden fixture, and the public tiny Flex end-to-end run.
Core unit tests and generated-parameter/report tests also passed.

The JAX check passed 31 assertions with zero failures: all 121
non-diagnostic outputs matched the pre-removal candidate byte-for-byte,
including raw/filtered matrices and caller outputs. New runs emitted no gDNA
reports. One spatial parameter guard was deferred because a dirty working-tree
binary has no immutable source revision; it must be checked on the committed
release build.

Legacy fused-alignment passed in all three configurations, with byte-identical
outputs across thread counts and input routes. All eight fresh partial-build
targets passed. The PBMC 100K vanilla, modern and modern-with-BAM profiles
passed their pinned golden checks, retaining 47,358/70,527 Gene/GeneFull UMIs
for vanilla and 47,203/70,270 for modern. Results are recorded in
`fused-align.log`, `partial/summary.tsv` and `pbmc100k/report.json`.
Release CI additionally builds and
validates the installed tarballs, installer and Debian packages; those gates
must pass before publication.

These are regression checks, not new paper performance benchmarks.
