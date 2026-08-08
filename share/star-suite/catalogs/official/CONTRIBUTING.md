# Contributing

Develop changes on a topic branch from `master`. Keep commits reviewable and
retain meaningful branch history; merge shared branches with `--no-ff`.

A contribution must keep recipe identity separate from site policy:

- workflow IDs use the owning catalog namespace;
- every script and schema path is relative to the catalog root;
- data locations, credentials, account names, and transfer endpoints are
  parameters, never committed defaults;
- partitioned workflows consume the public partition-manifest contract and do
  not implement or describe a private partition provider;
- execution evidence belongs in a provenance repository, not beside a recipe.

Run `tests/run_catalog_tests.sh` before proposing a merge. Recipe discovery and
help tests must not execute biological workflows.

