# Merge Launchpad 1.9.5.b into master and 1.10 development

User authorized committing/pushing 1.9.5.b and merging the fix into master and
`dev-release-v1.10.0`. Source commit: `9a856e5`.

Master merge `4f44406` has the same Git tree as the tested 1.9.5.b source commit.
The 1.10 merge starts from `d885bce` and retains both parent histories.

## 1.10 conflict resolution

- Keep STAR version and leading Debian changelog at 1.10.0; retain the 1.9.5.b
  changelog entry below it.
- Keep the deleted native Chromap integration and built-in `morphic_multiome`
  schema deleted. The only core code change from the pre-merge 1.10 tree is the
  `--build-features` report, with `chromap_atac:false` unconditionally.
- Carry portable config/forms, installed launcher, bundled browser assets,
  input/runtime checks, managed jobs and the release test dependency.
- The unchanged pinned catalog contains the old Multiome recipe. Exercise it
  as `starsuite.official/multiome` with an explicit external 1.9.5.b runtime;
  do not present it as a standalone STAR 1.10 engine. Current Multiome binary
  ownership remains with Multiomics Suite.
- Retain both artifact inventories and clarify ownership in current docs.

## Validation before pushing

- Complete MCP/Launchpad test suite on the merged 1.10 tree: **607 passed**,
  including the installed Chromium browser run/reload test.
- Fresh isolated source build, `nice -n 10 make -j8 core-portable`: passed.
- Built `STAR --version`: `1.10.0`; `--build-features`: `chromap_atac:false`.
- Launchpad runtime check rejects that actual binary with the Multiomics Suite
  handoff message before launching a recipe.
- Official snapshot validation: **11 recipes, 10 evidence records passed**.
- JavaScript/shell syntax, whitespace and deleted-source ownership checks pass.

Ignored evidence is in `tests/launchpad_portability_output/` in the isolated
merge checkout. No biological datasets were run and no release tag is part of
this branch integration.
