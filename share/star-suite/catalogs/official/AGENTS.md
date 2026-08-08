# Agent instructions

This repository owns reusable STAR Suite operational recipes and profiles, not
STAR Suite core algorithms and not project delivery state.

- Treat `catalog.yaml` as the public discovery contract.
- Keep all catalog paths relative and confined to this repository.
- Keep scripts open, reviewable, dry-run capable where practical, and free of
  credentials or private infrastructure defaults.
- Do not add raw data, large outputs, collaborator manifests, login hosts,
  transfer collection IDs, or absolute workstation paths.
- Public scatter/gather recipes may consume `biodepot.partition_manifest/v1`.
  Do not name, document, benchmark, or encode a non-public partition provider.
- Put site defaults under `profiles/`; do not bake them into portable workflows.
- A workflow change that alters output semantics belongs in STAR Suite core.
- Do not execute biological or performance benchmarks unless the user grants a
  specific execution authorization. Static validation and `--help` checks are
  safe.

