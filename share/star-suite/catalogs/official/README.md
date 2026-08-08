# STAR Suite Recipes

This public repository is the canonical home for reusable, non-core STAR Suite
recipes. It contains declarative workflows, open orchestration scripts, generic
scatter/gather contracts, and portable storage, site, and executor profiles.

STAR Suite itself owns compiled algorithms, output semantics, release tooling,
and core regression fixtures. Project-specific dataset manifests and delivery
automation belong in project catalogs such as `morphic-recipes`.

## Catalog use

The catalog follows `biodepot.recipe_catalog/v1`. A STAR Suite installation can
load a checkout directly:

```yaml
recipe_catalogs:
  - manifest: /path/to/STAR-suite-recipes/catalog.yaml
    trust: trusted
```

Use STAR Suite's catalog CLI to validate discovery without executing a recipe:

```bash
python3 -m mcp_server.recipe_cli --config /path/to/config.yaml catalog list
python3 -m mcp_server.recipe_cli --config /path/to/config.yaml recipe list
```

Released STAR Suite distributions vendor a reviewed, pinned snapshot so the
official recipes remain available offline. The Git repository is canonical for
development; the installed snapshot is canonical for a locked release run.

## Public boundary

Recipes accept data, reference, profile, and partition-manifest inputs. They do
not publish credentials, private transfer destinations, collaborator delivery
manifests, or any unpublished partition-generation implementation. Scripts in
this repository are open source under the MIT License.

See [CONTRIBUTING.md](CONTRIBUTING.md), [docs/OWNERSHIP.md](docs/OWNERSHIP.md),
and [docs/PARTITION_MANIFEST_CONTRACT.md](docs/PARTITION_MANIFEST_CONTRACT.md).

Reusable execution settings are separate from recipe identity. The profile
index includes generic Lustre, PSC Bridges-2, and Temporal/Slurm policies under
`profiles/`; project catalogs supply private accounts, roots, and delivery
settings without copying the workflows.
