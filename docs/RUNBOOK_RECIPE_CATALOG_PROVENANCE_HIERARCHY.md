# Runbook: Recipe Catalog and Provenance Hierarchies

Date: 2026-08-08

Status: catalog discovery, reconciliation, locking, and portable bundle
infrastructure implemented on `feature/recipe-catalog-stack`. Creation and
migration of the official, site, and project recipe/provenance repositories is
the next project. This work is intentionally isolated from the failed v1.7.0
release branch.

## Purpose

Make recipes a first-class STAR Suite extension surface without requiring every
recipe, dataset wrapper, site policy, or operational helper to live in the STAR
Suite source repository.

The resulting system must preserve an integrated, offline-capable STAR Suite
experience while allowing a user, laboratory, institution, or project to add
its own ordered recipe catalogs and provenance repositories. The same resolved
recipe must be usable from the command line, Workbench, and the Temporal
scheduler.

## Current Limitation

The MCP server currently has one `paths.repo_root`, one inline `workflows` list,
and one schema-loading pass. Relative schema, entry-script, and helper paths are
therefore assumed to belong to the STAR Suite checkout. An absolute path can
point at an external recipe, but its source identity and Git provenance are
then attributed through the global STAR Suite root. That is unsuitable as a
general recipe mechanism.

The first implementation tranche replaces that implicit ownership with an
explicit per-workflow origin while keeping the existing configuration fully
compatible.

## Ownership Boundaries

### STAR Suite

STAR Suite remains responsible for:

- compiled programs and installed tools that define output semantics;
- core CLI capability descriptions;
- build, release, diagnostic, and core-regression helpers;
- catalog schemas, loading, validation, locking, and discovery;
- rendering a reviewed recipe into a portable execution bundle;
- compatibility aliases for recipes moved during the transition.

### Official STAR Suite recipe catalog

A separate, public catalog will own:

- supported multi-stage operational workflows;
- generic scatter/gather orchestration;
- reusable filesystem policies such as Lustre guidance;
- public site-profile examples, including Bridges-2;
- recipe-specific orchestration helpers and smoke fixtures;
- dated validation evidence that is not core source or release evidence.

Official releases should contain a validated, pinned snapshot under
`share/star-suite/catalogs/official/`. The external repository remains
canonical, while the snapshot preserves offline integration.

### Site and project catalogs

Site catalogs provide scheduler, filesystem, container, and institutional
defaults. Project catalogs provide dataset manifests, project-specific output
selection, and handoff behavior. They extend official recipes rather than
copying them when possible.

### Temporal scheduler

Temporal executes an already resolved workflow. It does not discover catalogs,
select recipe inheritance, or interpret assay intent. Scheduler-specific
settings are supplied separately as an execution/site profile.

### Provenance repositories

Provenance repositories contain immutable records of executions. A context may
search several repositories for prior evidence but must select exactly one
writable destination for a new run.

## Catalog Stack

The target execution context is:

```yaml
schema: biodepot.execution_context/v1

suite:
  executable: /opt/star-suite/bin/STAR
  version_constraint: ">=1.7,<2"

catalogs:
  - manifest: bundle://star-suite/core/catalog.yaml
  - manifest: /opt/star-suite/share/star-suite/catalogs/official/catalog.yaml
  - manifest: /etc/star-suite/recipes/catalog.yaml
  - manifest: .star-suite/recipes/catalog.yaml

provenance:
  search:
    - id: site
      root: /shared/provenance/star-suite
    - id: project
      root: .star-suite/provenance
  write:
    id: project
    root: .star-suite/provenance

profiles:
  storage: lustre
  site: bridges2
  executor: temporal-slurm
```

Initial implementation is deliberately narrower: local catalog manifests and
local provenance roots configured through `mcp_server/config.yaml`. Bundle,
Git, and OCI sources come after local path semantics and locking are stable.

## Catalog Manifest Contract

```yaml
schema: biodepot.recipe_catalog/v1
id: lab-recipes
namespace: example.lab
version: "1"

workflows:
  - id: example.lab/bulk-rna
    logical_id: bulk-rna
    version: "2.1.0"
    title: Example bulk RNA-seq
    summary: Site specialization of the official workflow.
    kind: shell_workflow
    entry_script: scripts/run_bulk_rna.sh
    schema_file: workflows/bulk_rna.yaml
    visibility: private
    applications: [cli, workbench]
    extends: starsuite.official/bulk-rna
    replaces: [example.lab/legacy-bulk-rna]
    image: registry.example/star-suite@sha256:...
```

All relative paths are resolved against the directory containing
`catalog.yaml`, not against the STAR Suite checkout or process working
directory. External catalog paths must not escape the catalog root.

Catalog order is discovery order, not an implicit override mechanism:

- catalog IDs and namespaces must be unique;
- fully qualified workflow IDs must be unique;
- duplicate IDs are errors that name both sources;
- later inheritance uses explicit `extends` or `replaces` fields;
- a local name is never silently substituted for a built-in recipe.

## Multi-source Reconciliation

Source identity and recipe intent are separate. `id` is the unique,
source-owned workflow identity; `logical_id` groups recipes that implement the
same intent. Consumers can also filter candidates using the optional
`applications` list. An empty list means that the recipe is compatible with
all applications.

Three policies are implemented:

- `keep_separate` is the safe default. No candidate is selected when several
  sources provide the same logical recipe. Callers use a candidate workflow ID
  or the unambiguous `catalog-id::workflow-id` source reference.
- `prompt` returns `selection_required` and the complete ordered candidate
  list. Workbench uses this status to ask the user and records the chosen
  source in the lock.
- `prefer_newest` selects the highest PEP 440 recipe version. Equal versions
  are resolved by later catalog discovery order, then workflow ID. This is
  intended only for applications whose owner explicitly accepts automatic
  upgrades.

The default and application overrides are configured independently:

```yaml
recipe_resolution:
  default_policy: keep_separate
  applications:
    cli: keep_separate
    workbench: prompt
    unattended-builder: prefer_newest
```

Temporal has no catalog resolution policy. It receives a lock produced by an
upstream user, CLI, or Workbench session and executes that exact selection.

Catalog trust is assigned by the local execution context rather than by a
self-assertion in `catalog.yaml`:

```yaml
recipe_catalogs:
  - manifest: /opt/star-suite/share/star-suite/catalogs/official/catalog.yaml
    trust: trusted
  - manifest: .star-suite/recipes/catalog.yaml
    trust: untrusted
```

External catalogs default to `untrusted`; the synthesized built-in catalog is
trusted. Discovery and locks carry this state. Trust does not hide a recipe or
cause discovery-time execution; Workbench or the executor admission policy may
require confirmation for an untrusted lock.

`extends` is explicit, validated lineage in the v1 catalog contract. Every
child schema is complete and authoritative: fields are not implicitly merged
across YAML documents. This deliberately simple merge rule avoids list and
parameter ambiguity while still locking the complete parent-to-child chain.
Missing references and inheritance cycles fail catalog loading. `replaces`
records supersession but never removes a source from `keep_separate`
discovery.

## Built-in Compatibility Catalog

The existing top-level `workflows:` configuration is treated as a synthesized
catalog with:

```text
catalog id: starsuite-builtin
namespace: starsuite.core
root: paths.repo_root
```

This means existing configurations, workflow IDs, API clients, tests, and
entry-script paths continue to work without migration. New external catalogs
are additive.

After the loader is established, current workflows will be classified:

1. Core one-command capability descriptions remain built in.
2. Core smoke and parity workflows remain in a test catalog in this repository.
3. Generic multi-stage production recipes move to the official catalog.
4. Dataset- or institution-specific recipes move to their project/site catalog.
5. Existing IDs receive compatibility aliases or forwarding entries for one or
   two release cycles.

## Helper Script Policy

A helper remains in STAR Suite when it defines STAR output semantics or is
required to build, release, diagnose, or regression-test the suite. A helper
moves with a recipe when it performs orchestration, staging, output selection,
reporting, transfer, or handoff assembly.

Recipe helpers must:

- be relative to and confined within their catalog root;
- be listed in the workflow schema;
- be included by digest in a resolved execution bundle;
- never be executed merely because a catalog was discovered;
- pass trust and preflight checks before execution.

## Generic Partition Boundary

Public recipes describe a generic partition set and do not prescribe how it is
created. A scatter/gather recipe accepts a validated partition manifest whose
entries have stable ordinals and paired inputs. The contract requires complete,
non-overlapping logical input coverage, mate co-location, per-partition
completion evidence, and gather only after all required partitions succeed.

Public catalog content must not include non-public provider names, commands,
formats, transport details, or implementation-specific performance claims. A
private catalog may bind an authorized partition provider to the public
contract. Adding a provider to the official catalog requires a separate
publication and licensing review.

## Provenance Hierarchy

The configuration model separates read-only search roots from the writable
destination:

```yaml
provenance:
  search:
    - id: official-evidence
      root: /opt/star-suite/share/star-suite/evidence
    - id: site
      root: /shared/star-suite-provenance
    - id: project
      root: .star-suite/provenance
  write:
    id: project
    root: .star-suite/provenance
```

Search order may rank candidate prior runs, but records are never merged. The
chosen prior run is recorded explicitly as the provenance oracle. A resolved
run records:

- suite version, source revision, and binary digest;
- workflow ID and schema version;
- catalog ID, namespace, version, source revision, and digest;
- every inherited recipe and helper digest;
- resolved parameters and selected output profile;
- storage, site, and executor profile locks;
- input and output manifests;
- scheduler accounting, completion evidence, and validation results.

## Resolution and Execution Safety

Discovery and execution are separate trust boundaries.

1. Discovery parses and validates declarative catalog metadata only.
2. Resolution rejects duplicate IDs, path escapes, schema/config ID mismatch,
   and incompatible suite constraints.
3. Rendering produces commands and a lock; it does not execute them.
4. Workbench or CLI presents the resolved source and command for approval.
5. Preflight validates trusted roots, executability, resources, and required
   inputs.
6. Execution consumes the locked bundle. Compute workers do not re-resolve a
   mutable catalog.

Remote sources must eventually be pinned by immutable Git revision or OCI
digest. Floating branches may be used for discovery but cannot produce a
production lock without resolving to a commit.

## Workbench Integration

Workbench should display, for every recipe:

- fully qualified ID and version;
- catalog ID, namespace, and trust state;
- inheritance chain;
- compatible STAR Suite versions;
- selected provenance oracle;
- static scatter expansion and resource profile;
- helper files and digests included in the run bundle.

Workbench performs static expansion into concrete nodes before scheduler
submission. Advisory scatter candidates remain useful planning metadata, but
the resolved nodes and links are the binding execution surface.

The MCP discovery sequence is `list_recipe_candidates` -> `resolve_recipe` ->
`build_recipe_lock`. A `prompt` conflict is represented as data, not an
exception, so Workbench can render the candidates and submit the chosen
workflow/source reference. The CLI bundle command emits the same
`manifests/resolved_workflow_v1.json` path already recognized by Workbench's
Temporal/Slurm packet builder.

## Temporal Integration

Temporal receives:

- `biodepot.resolved_workflow/v1` with concrete nodes and links;
- recipe and helper locks;
- one selected site/executor profile;
- one writable provenance destination;
- retry, cancellation, cleanup, and terminal-evidence policy.

Temporal must not fetch a newer recipe during a retry. Retries use the same
locked execution bundle and partition ordinals.

The current generic bundle is a single concrete wrapper node. A future recipe
catalog may publish an application-specific static expansion into multiple
scatter/gather nodes, but that expansion must occur before scheduler
submission and must retain the recipe lock digest. Public recipes describe
only the generic partition-manifest contract; they do not name or expose the
non-public partition provider.

## Implementation Phases

### Phase 1: Local catalog stack and workflow origins

Implemented.

- Add catalog and provenance configuration models.
- Treat existing `workflows:` as the built-in compatibility catalog.
- Load ordered external `catalog.yaml` manifests.
- Track catalog root and manifest for every workflow.
- Reject duplicate catalog/workflow identities and catalog path escapes.
- Resolve scripts, schema files, rendering, and authenticated provenance from
  the owning catalog root.
- Add focused unit tests without executing biological workflows.

### Phase 2: Catalog inspection and locks

Implemented. The CLI commands are `catalog list`, `catalog describe`,
`provenance list`, and
`recipe list|conflicts|resolve|show|validate|render|lock|bundle`.

- Add `list_recipe_catalogs` and `describe_recipe_catalog` discovery surfaces.
- Surface catalog identity in workflow list/describe responses.
- Generate a deterministic recipe lock including schema and helper hashes.
- Add a CLI for `catalog list`, `recipe show`, `recipe validate`, and
  `recipe render`.

### Phase 3: Explicit inheritance and profiles

Infrastructure implemented for validated complete-schema inheritance,
supersession metadata, application compatibility, and per-application
resolution policies. Storage, site, and executor settings remain separate from
portable recipe identity and will be defined by the catalog/provenance
repositories created next.

- Add `extends` with cycle detection and deterministic merge rules.
- Add output-composition, storage, site, and executor profile schemas.
- Forbid implicit ID shadowing.
- Add compatibility aliases for migrated workflows.

### Phase 4: Resolved execution bundles

Implemented for generic wrapper recipes. A bundle contains a deterministic
recipe lock, every non-native declared schema/helper source, native dependency
digests, the portable resolved-workflow packet, the rendered command, and a
planned provenance record. Bundle creation is atomic and refuses to overwrite
an existing target.

- Copy or bind locked helpers into a self-contained run directory.
- Emit the portable resolved-workflow packet used by Workbench and Temporal.
- Emit the planned provenance record in the bundle. Writing/finalizing it is
  enabled after the provenance repositories and terminal-evidence worker are
  created in the next project.

### Phase 5: Recipe migration

Deferred by design until this catalog mechanism is merged. This is the next
project: create the repositories, author the generic recipes, validate them,
and then migrate eligible helpers.

- Create the canonical official recipe repository.
- Move generic operational workflows and their helpers.
- Retain built-in core capability and test catalogs.
- Add release-time vendoring of a pinned official catalog snapshot.
- Retire forwarding entries only after Workbench and scheduler consumers have
  migrated.

### Phase 6: Remote catalog sources

Deferred until the public repositories and licensing boundary exist. For now,
different sources are local checked-out or installed catalog roots. This keeps
discovery offline and makes Git/OCI acquisition policy independent of recipe
resolution.

- Add pinned Git and OCI catalog sources.
- Add signature/trust policy and an offline cache.
- Never fetch or execute remote content during scheduler retries.

## Phase 1 Acceptance Tests

- Existing configuration loads with no catalog settings and returns the same
  workflows and rendered paths.
- Two external catalogs load in declared order.
- External schema, entry, and helper paths resolve from the owning catalog.
- Public discovery exposes safe catalog identity but not host paths or Git
  metadata.
- Authenticated discovery reports the owning catalog root and Git provenance.
- Duplicate catalog IDs, namespaces, and workflow IDs fail closed.
- Missing manifests and schema files fail with source-specific errors.
- External absolute paths and `..` escapes are rejected.
- A failed reload restores the previous schemas, workflow configurations, and
  origins atomically.
- Provenance configuration rejects duplicate search IDs and an ambiguous write
  destination.

## Implemented Reconciliation and Bundle Tests

- all three policies return deterministic, structured results;
- application compatibility filters source variants before reconciliation;
- an explicit user source overrides a prompt without changing global policy;
- `prefer_newest` uses PEP 440 ordering and deterministic equal-version ties;
- missing lineage references, invalid versions, and cycles fail closed;
- recipe locks are stable and change when a declared helper changes;
- public catalog discovery redacts host roots and Git metadata;
- public provenance discovery exposes ordered IDs and the write selection but
  redacts host roots;
- bundles are created atomically, never overwrite a target, and contain a
  valid `biodepot.resolved_workflow/v1` packet at the Workbench handoff path;
- the CLI and MCP tools use the same loader and resolution functions.

## Definition of Done

This effort is complete when a fresh STAR Suite installation can discover its
built-in and pinned official recipes offline; a user can add site and project
catalogs without editing STAR Suite; Workbench resolves the same catalog stack;
Temporal executes only a locked resolved bundle; and every completed run points
unambiguously to its suite, recipe, catalog, helper, profile, and provenance
identities.
