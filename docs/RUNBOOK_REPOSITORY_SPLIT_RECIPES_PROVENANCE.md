# Repository Split Runbook: Core, Recipes, Provenance

Date: 2026-05-22
Status: Phase 1 cutover implemented locally and pushed to the split repos.

## Goal

Keep STAR-suite focused on core processing while moving production recipes and
run provenance into repositories with lifecycles that match operational use.
During the transition, duplication is allowed, but every duplicated script or
document must have an explicit canonical owner and a drift policy.

The target end state:

1. `STAR-suite`: core STAR/Flex/Solo processing, compiled tools, core tests,
   release packaging, and generic workflow machinery.
2. `morphic-recipes`: dataset and project recipes, launch
   wrappers, downstream h5ad/QC/CellBender/celltyping/Globus helpers, and
   workflow schemas for production runs.
3. `morphic-provenance`: immutable run records, manifests,
   checksums, rendered commands, environment pins, and handoff provenance.

The boundaries below are the source of truth.

## Initial Local Repos

Initial local split targets were scaffolded on 2026-05-22:

```text
/mnt/pikachu/morphic-recipes
/mnt/pikachu/morphic-provenance
```

GitHub remotes:

```text
git@github.com:morphic-bio/morphic-recipes.git
git@github.com:morphic-bio/morphic-provenance.git
```

Initial commits:

```text
morphic-recipes:     682f554 Initial Morphic recipe mirror
morphic-provenance:  19f0fe9 Initial Morphic provenance scaffold
```

No STAR-suite files have been removed. The first recipes repo commit is a mirror
of selected production scripts, runbooks, manifests, and workflow YAMLs copied
from STAR-suite commit `43a5853af0c627925f827ab576814b770d1874c1`.

Phase 1 cutover is now active:

- `morphic-recipes` is the canonical home for new production recipe work.
- Recipe launchers that need core binaries use an external STAR-suite checkout
  through `STAR_SUITE_ROOT` instead of assuming `core/` exists in the recipes
  repo.
- STAR-suite keeps compatibility launchers for active MSK/JAX/Multiome
  production entrypoints and remote downstream helpers. Those launchers delegate
  to `MORPHIC_RECIPES_ROOT` and should not receive new recipe logic.

## Definitions

- Core processing: code that changes STAR-suite outputs for a general class of
  inputs, such as barcode parsing, feature assignment, STARsolo counting,
  Velocyto MEX writing, OCM materialization, Chromap integration, or shared
  EmptyDrops behavior.
- Recipe: an operational route for running core processing on real datasets,
  including dataset manifests, wrapper defaults, remote staging, downstream
  conversion, CellBender, celltyping, packet assembly, and transfers.
- Provenance record: the immutable record of one executed run: exact code
  revisions, rendered commands, parameters, environment, inputs, outputs,
  checksums, logs, and handoff notes.
- Handoff packet: the deliverable directory sent to collaborators. It may copy
  or summarize provenance, but it is not the canonical run registry.

## Repository Ownership Policy

| Artifact | Canonical home | Notes |
| --- | --- | --- |
| STAR core source, feature overlays, Flex, SLAM, compiled tools | `STAR-suite` | Includes `core/`, `flex/`, `slam/`, release build logic, and toolchain tests. |
| Core algorithm docs and regression runbooks | `STAR-suite` | Keep docs that explain behavior or validated core contracts. |
| Generic public demos and core smoke fixtures | `STAR-suite` | Small public fixtures and CI-relevant scripts stay with the code they test. |
| Dataset-specific production wrappers | recipes repo | Examples: MSK/JAX/UCSF production launchers, sample-specific manifests, remote worker orchestration. |
| Downstream h5ad/QC/CellBender/celltyping scripts | recipes repo | Move scripts whose main purpose is post-MEX operational analysis or handoff generation. |
| Packet builders and Globus upload helpers | recipes repo | The recipe repo owns deliverable layout policy. |
| Workflow schemas for production datasets | recipes repo | During transition they may remain duplicated under `mcp_server/workflows/`; final source should be recipes. |
| MCP server implementation | `STAR-suite` unless split later | The server engine can remain core tooling. It should eventually load external recipe workflow directories. |
| Executed run manifests, checksums, command JSON, logs, environment pins | provenance repo | Do not rewrite after handoff except to add an explicit correction note. |
| Dockerfiles for STAR-suite release images | `STAR-suite` | These define the core distribution. |
| Dockerfiles or lockfiles for downstream recipes | recipes repo | The provenance repo records the resolved image digest and exact package locks used for each run. |
| Large data outputs | outside git | Reference by URI/path/checksum from provenance. Do not commit large h5ad/BAM/FASTQ outputs. |

## Lookup Policy For Users And Agents

When trying to reproduce a STAR-suite processing pipeline, look in this order:

1. Provenance repo:
   `runs/<project>/<run_id>/README.md` and `run.json`.
   This is the authority for what actually ran.
2. Recipes repo:
   workflow schema, dataset manifest, launch script, downstream recipe, and
   packet-builder used by the provenance record.
3. STAR-suite:
   core commit, core docs, compiled binary provenance, tests, and release notes
   referenced by the run record.
4. Handoff packet:
   collaborator-facing README and checksums. Use this as a convenience view,
   not as the complete execution record.

If the provenance record and a recipe default disagree, trust the provenance
record for that run. If a recipe and STAR-suite docs disagree about core
behavior, verify the STAR-suite commit recorded in provenance.

## Transition Phases

### Phase 0: Inventory And Mark Canonical Owners

Create an inventory of `scripts/`, `docs/RUNBOOK_*`, `mcp_server/workflows/`,
and packet-generation helpers. Classify each item as:

- `keep-core`
- `move-recipe`
- `move-provenance`
- `duplicate-transition`
- `retire`

For duplicated files, add a short header or adjacent note with:

```text
Canonical: <repo:path>
Mirror: <repo:path>
Last synced from: <commit>
Drift policy: change canonical first, then sync mirror in the same task
```

### Phase 1: Duplicate Recipes Without Breaking Existing Agents

Copy production recipes into the recipes repo while leaving compatibility copies
or thin wrappers in STAR-suite. STAR-suite wrappers should either call the
recipe repo path or print the canonical location clearly before running.

Allowed duplication window: one or two release cycles, or until all active
agents and runbooks point to the recipes repo.

During this phase:

- New production launches should be authored in recipes first.
- STAR-suite copies should receive only compatibility fixes.
- Any change to a duplicated recipe must update the canonical copy and record
  the source commit in the mirror note.

### Phase 2: Make Recipes Canonical

Move dataset-specific operational docs and workflow schemas to the recipes repo.
STAR-suite should keep only:

- short migration notes,
- links to recipe docs,
- core behavior documentation,
- regression docs needed to validate STAR-suite itself.

The MCP server should support an external workflow directory, for example:

```text
STAR_SUITE_WORKFLOW_DIRS=/path/to/star-recipes/workflows:/path/to/local/workflows
```

Until that exists, keep workflow schema mirrors in `mcp_server/workflows/` with
canonical-owner notes.

### Phase 3: Make Provenance The Run Registry

Every production or handoff run gets an immutable provenance folder before
handoff. The folder is created at launch time and finalized after checksums and
transfer status are known.

Required layout:

```text
runs/<project>/<run_id>/
  README.md
  run.json
  inputs/
    manifest.tsv
    checksums.tsv
  commands/
    00_preflight.argv.json
    01_star.argv.json
    02_downstream.argv.json
    03_transfer.argv.json
  environment/
    image-digest.txt
    Dockerfile
    conda-lock.yml
    pip-freeze.txt
    nvidia-smi.txt
    host.txt
  outputs/
    manifest.tsv
    checksums.tsv
  logs/
    README.md
  handoff/
    README.md
    transfer-status.json
```

Large logs may be stored outside git if necessary, but the provenance record
must contain their path, size, checksum, and retention policy.

### Phase 4: Retire Old Locations

After active recipes have moved, remove or archive STAR-suite scripts that are
not needed for core tests. Leave a short tombstone document only when old
runbooks or collaborators are likely to search for the previous path.

## Provenance Schema Minimum

`run.json` must be machine-readable and include at least:

```json
{
  "schema_version": "0.1",
  "run_id": "20260522T000000Z_msk40ko_prod",
  "project": "MSK-40KO",
  "created_utc": "2026-05-22T00:00:00Z",
  "operator": "human-or-agent-name",
  "status": "complete",
  "star_suite": {
    "repo": "STAR-suite",
    "commit": "<sha>",
    "dirty": false,
    "binary": "<path-or-image>",
    "binary_sha256": "<sha256>"
  },
  "recipes": {
    "repo": "star-recipes",
    "commit": "<sha>",
    "workflow_id": "<workflow-id>",
    "workflow_schema_version": "<version>"
  },
  "environment": {
    "container_image": "<registry/name@sha256:digest>",
    "cuda_required": true,
    "cellbender_cuda": true
  },
  "inputs": [
    {
      "name": "<sample-or-file>",
      "path": "<path-or-uri>",
      "sha256": "<sha256>",
      "bytes": 0
    }
  ],
  "outputs": [
    {
      "name": "<output>",
      "path": "<path-or-uri>",
      "sha256": "<sha256>",
      "bytes": 0
    }
  ]
}
```

The schema should be append-only for completed runs. Corrections should add a
new `corrections` entry with author, timestamp, reason, and changed fields.

## Environment Policy

Dockerfiles are not sufficient provenance by themselves. Every run must record
the resolved container image digest and package versions actually used.

For GPU downstream work:

- Record `nvidia-smi` output before or during the run.
- Record CUDA, driver, PyTorch, CellBender, Scanpy, AnnData, scikit-learn, and
  scvi-tools versions when present.
- CellBender production runs must use CUDA. The rendered command must contain
  the recipe-level GPU flag and the CellBender-level `--cuda` flag, and Docker
  launches must include GPU access.

## Initial Migration Candidates

Likely `move-recipe` or `duplicate-transition`:

- `scripts/run_msk_40ko_pipeline_from_manifest.py`
- `scripts/run_jax_scrnaseq02_ocm_production_batch.sh`
- `scripts/run_remote_scrna_downstream_rsync.sh`
- `scripts/run_remote_cellbender_*.sh`
- `scripts/run_scrna_downstream_gene_full_velocyto.sh`
- `scripts/build_gene_full_velocyto_h5ad.py`
- `scripts/postprocess_downstream_filters.py`
- `scripts/integrate_feature_library.py`
- packet builders and Globus handoff helpers under `plans/artifacts/` once
  stabilized into real scripts.

Likely `keep-core`:

- `core/`, `flex/`, `slam/`
- `scripts/release/`
- `scripts/docker/` for STAR-suite release images
- `tests/`
- core parity and regression tools used by STAR-suite CI or paper validation
- `mcp_server/` engine code

Likely `move-provenance`:

- final handoff README/provenance files,
- run manifests for MSK/JAX/UCSF production runs,
- rendered command records,
- checksums and transfer status files,
- environment lockfiles for completed runs.

## Authoring Rules After The Split

Use this decision tree for new work:

1. Does it change STAR/Flex/Solo behavior or a compiled tool?
   Put it in STAR-suite.
2. Does it launch a real dataset, choose operational defaults, run downstream
   analysis, or assemble a handoff packet?
   Put it in recipes.
3. Does it describe one completed run or make that run reproducible?
   Put it in provenance.
4. Does it validate a core invariant for all future code?
   Put it in STAR-suite tests or docs.
5. Does it validate one dataset recipe?
   Put it in recipes, and record executed results in provenance.

Do not add new dataset-specific production scripts to STAR-suite unless they
are temporary mirrors with a canonical recipes path.

## Handoff Packet Policy

Handoff packets should include collaborator-facing summaries and checksums, but
they should point back to the provenance run id. Recommended packet files:

```text
README.md
PROVENANCE.md
MANIFEST.tsv
CHECKSUMS.sha256
```

`PROVENANCE.md` should include:

- provenance repo path and run id,
- STAR-suite commit and binary checksum,
- recipes commit and workflow id,
- container image digest,
- input manifest checksum,
- output manifest checksum,
- transfer task id and completion status when applicable.

## Acceptance Criteria For The Split

The migration is complete when:

- a new agent can start from a handoff packet, find the provenance run id, and
  identify the exact STAR-suite and recipes commits without asking a human;
- production recipe changes no longer land primarily in STAR-suite;
- STAR-suite can run its core build and smoke tests without private production
  datasets;
- recipes can render a full command without reading STAR-suite-local runbooks;
- provenance records contain enough information to reconstruct the run
  environment and input/output file inventories without relying on memory or
  chat history.
