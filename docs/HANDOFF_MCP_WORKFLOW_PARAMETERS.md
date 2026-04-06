# Handoff: MCP Workflow Parameters Implementation

> **Date**: 2026-04-06
> **Status**: Ready to implement
> **Runbook**: `plans/mcp_workflow_parameters_runbook.md`
> **Repo**: `/mnt/pikachu/STAR-suite`
> **Branch**: `master` (create feature branch `feature-mcp-workflow-params`)

## What This Is

Extend the existing STAR-suite MCP server so it can serve structured workflow
parameter metadata and validate/render workflow requests. The first supported
workflow is the UCSF STAR-suite production wrapper family.

Today the MCP server is a discovery / preflight / execution server with opaque
`args: list[str]`. This work adds a parallel structured workflow contract so
agents can discover what parameters a workflow accepts, validate a request, and
render it into a command — without parsing Bash help text.

## What To Read First

1. **The runbook** — `plans/mcp_workflow_parameters_runbook.md` — complete
   specification with acceptance criteria, implementation order, test plan,
   and non-goals. This is the source of truth. This handoff document provides
   context and architecture details to accelerate implementation.

2. **The existing MCP server** — `mcp_server/` — understand the patterns
   before adding to them.

## Existing MCP Server Architecture

### Framework and Patterns

- **Framework**: FastMCP (Starlette + Uvicorn)
- **Config**: YAML (`mcp_server/config.yaml`) validated by Pydantic models
  (`mcp_server/schemas/config.py`), loaded via singleton in `mcp_server/config.py`
- **Tool registration**: `@mcp.tool()` decorator in `mcp_server/app.py`
- **Response models**: Pydantic `BaseModel` subclasses in
  `mcp_server/schemas/responses.py`, returned as `.model_dump()` dicts
- **Security**: All paths validated against `trusted_roots` via
  `mcp_server/tools/utils.py:is_path_allowed()`
- **Auth**: Optional token in config; checked via `check_auth()` helper
- **Tests**: pytest with fixtures in `mcp_server/tests/conftest.py`

### Key Files

```
mcp_server/
├── app.py                 # FastMCP server, tool registration (600 lines)
├── config.py              # Config loading singleton (129 lines)
├── config.yaml            # Production config
├── schemas/
│   ├── config.py          # MCPConfig, ScriptConfig, etc. Pydantic models
│   ├── run_config.py      # RunConfig model (script, args, module, ...)
│   └── responses.py       # All response models (PreflightResponse, etc.)
├── tools/
│   ├── preflight.py       # 9 check functions + run_preflight() orchestrator
│   ├── executor.py        # ThreadPoolExecutor job queue + subprocess runner
│   ├── discovery.py       # list_datasets, list_test_suites, find_docs, find_tests
│   ├── build.py           # Source hash change detection, make build
│   ├── reload.py          # Hot config reload
│   └── utils.py           # is_path_allowed, validate_path, find_binary, etc.
└── tests/
    ├── conftest.py        # temp_dir, sample_config_dict, loaded_config fixtures
    ├── test_config.py
    ├── test_preflight.py
    ├── test_executor.py
    ├── test_discovery.py
    ├── test_utils.py
    ├── test_auth.py
    └── test_integration.py
```

### How To Add a New Tool (pattern to follow)

1. **Create response models** in `mcp_server/schemas/responses.py`:
   ```python
   class WorkflowInfo(BaseModel):
       id: str
       title: str
       summary: str
       ...
   ```

2. **Create handler module** in `mcp_server/tools/workflows.py`:
   ```python
   from ..config import get_config
   from ..schemas.responses import WorkflowInfo
   
   def list_workflows() -> ListWorkflowsResponse:
       config = get_config()
       ...
   ```

3. **Register in `app.py`**:
   ```python
   @mcp.tool()
   def list_workflows(auth_token: Optional[str] = None) -> dict:
       """Return supported workflow templates."""
       auth_error = check_auth(auth_token, is_discovery=True)
       if auth_error:
           return auth_error.model_dump()
       try:
           result = _list_workflows()
           return result.model_dump()
       except Exception as e:
           return ErrorResponse(code="...", message=str(e)).model_dump()
   ```

4. **Write tests** following the existing pattern with `loaded_config` fixture.

### Existing Config Structure (relevant sections)

```yaml
# config.yaml (existing)
server:
  host: "0.0.0.0"
  port: 8765
  auth_token: "${MCP_AUTH_TOKEN}"
  transport: "streamable-http"
  public_discovery: true

paths:
  repo_root: "/mnt/pikachu/STAR-suite"
  artifact_log_root: "/mnt/pikachu/STAR-suite/mcp_server/artifacts"
  temp_root: "/tmp/star_suite_mcp"

trusted_roots:
  - "/mnt/pikachu/STAR-suite"
  - "/mnt/pikachu"
  - "/storage"
  - "/home/lhhung"
  - "/tmp/star_suite_mcp"

scripts:
  - name: "ucsf_corrected_production"
    path: "scripts/paper/run_ucsf_corrected_production_workflow.sh"
    ...
```

The new `workflows:` section goes alongside `scripts:` (not replacing it).

## UCSF Workflow Parameter Schema (Complete Reference)

The following is the complete parameter inventory extracted from the actual
wrapper scripts. Use this to build the YAML schema file.

### Top-Level: `scripts/paper/run_ucsf_corrected_production_workflow.sh`

This is the entry script for workflow id `ucsf_star_suite_production`.

#### Selection Group

| Flag | Type | Required | Default | Notes |
|------|------|----------|---------|-------|
| `--samples` | string (CSV) | * | `""` | Comma-separated sample IDs |
| `--all-samples` | bool | * | `0` | Run every corrected sample |

**Constraint**: Mutually exclusive. If neither provided, defaults to
`all_samples=true` (lines 102–104 of script).

#### Inputs Group

| Flag | Type | Required | Default | Notes |
|------|------|----------|---------|-------|
| `--dataset-root` | directory | No | `/mnt/pikachu/ucsf-perturb-seq-corrected` | path_must_exist |
| `--feature-ref` | file | No | `/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv` | path_must_exist |
| `--genome-dir` | directory | No | `/storage/autoindex_110_44/bulk_index` | path_must_exist |
| `--solo-cb-whitelist` | file | No | `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt` | path_must_exist |
| `--cr-whitelist` | file | No | `/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt` | path_must_exist |
| `--star-bin` | file | No | `${REPO_ROOT}/core/legacy/source/STAR.release` | must be executable |

#### Performance Group

| Flag | Type | Required | Default | Notes |
|------|------|----------|---------|-------|
| `--threads` | int | No | `24` | positive |
| `--cellbender-cpu-cores` | int | No | `${threads}` | positive |
| `--downsample-reads` | int | No | `0` | non-negative; 0 = no downsampling |
| `--downsample-seed` | int | No | `1` | exported as env var |

#### Execution Mode Group

| Flag | Type | Required | Default | Notes |
|------|------|----------|---------|-------|
| `--cellbender-gpu` | bool | No | `0` | only effective when not `star_only` |
| `--trim-qc` | bool | No | `0` | emit read-level trim-QC |
| `--star-only` | bool | No | `0` | skip downstream + CellBender |
| `--dry-run` | bool | No | `0` | prepare manifests only |

**Dependencies**:
- If `star_only=true`, downstream + CellBender stages are skipped entirely
- `cellbender_gpu` only matters when `star_only=false`
- When not `star_only`, the script always adds `--run-downstream`,
  `--adaptive-filter`, `--run-cellbender` to the batch runner invocation

#### Transfer Group

| Flag | Type | Required | Default | Notes |
|------|------|----------|---------|-------|
| `--globus-src-endpoint` | string (UUID) | No | `""` | |
| `--globus-dst-endpoint` | string (UUID) | No | `""` | |
| `--globus-dst-root` | string | No | `""` | |
| `--globus-poll-seconds` | int | No | `30` | positive |

**Constraint**: If any one Globus param is supplied, all three endpoint/root
params should be present (warn/error if incomplete set).

#### Output Group

| Flag | Type | Required | Default | Notes |
|------|------|----------|---------|-------|
| `--out-root` | directory | No | `/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_$(date)` | auto-timestamped |

### Sub-Scripts (for reference, not first-pass schema)

The batch runner (`run_ucsf_perturb_yremove_batch.sh`) adds:
- `--sample-map`, `--cr-config`, `--cr-sample-id`
- `--run-downstream`, `--adaptive-filter`, `--run-cellbender`
- `--downstream-output-name`, `--trim-qc-max-reads`

The downstream wrapper (`run_scrna_downstream_gene_full_velocyto.sh`) adds:
- `--run-dir` (required), `--output-dir`
- `--min-genes`, `--max-genes`, `--mt-pct-cutoff`, `--n-mad`
- `--docker-image`, `--feature-gather-image`, `--cellbender-image`
- `--reuse-cellbender`, `--cellbender-layer`, `--cellbender-flags`

The remote CellBender wrapper (`run_remote_cellbender_rsync.sh`) adds:
- `--downstream-dir` (required), `--remote-host` (required),
  `--remote-root` (required)
- `--cellbender-gpu-device`, `--local-log`, `--no-sync-image`, `--keep-remote`

For the first implementation, only the top-level production workflow parameters
need to be in the schema. The sub-script parameters are documented here for
completeness and can be added in a later pass.

### Cross-Script Parameter Propagation

```
run_ucsf_corrected_production_workflow.sh
  ├─ --samples / --all-samples        → batch runner
  ├─ --downsample-reads               → batch runner
  ├─ --downsample-seed                → exported env var DOWNSAMPLE_SEED
  ├─ --threads                        → batch runner
  ├─ --star-bin                       → batch runner
  ├─ --cellbender-cpu-cores           → batch runner (gated by !star_only)
  ├─ --cellbender-gpu                 → batch runner (gated by !star_only)
  ├─ --out-root                       → batch runner
  ├─ --dataset-root                   → batch runner
  ├─ --feature-ref                    → batch runner
  ├─ --genome-dir                     → batch runner
  ├─ --solo-cb-whitelist              → batch runner
  ├─ --cr-whitelist                   → batch runner
  ├─ --globus-*                       → batch runner
  ├─ --trim-qc                        → batch runner
  ├─ --star-only                      → gates whether downstream flags are added
  └─ --dry-run                        → batch runner
```

## Rendering Rules

The renderer must produce an argv array for the entry script. Key rules:

1. **Booleans**: Emit flag only when value is `true` (e.g., `--star-only`
   appears or is omitted, never `--star-only 0`)
2. **Strings/ints/files/dirs**: Emit `--flag value`
3. **CSV list** (`--samples`): Emit `--samples "EBs1_1,EBs1_2,..."`
4. **Omit absent optionals**: Don't emit flags for unset optional params
5. **Deterministic order**: Emit flags in a fixed order (by group)
6. **`out_root` auto-generation**: If not provided, note that the script
   generates a timestamped default — the renderer can either leave it out
   (let the script default) or generate one

## Implementation Order

From the runbook §16:

1. Add workflow config and schema Pydantic models
   → `mcp_server/schemas/workflow.py`
2. Add workflow schema YAML loader
   → extend `mcp_server/config.py`
3. Add `list_workflows` and `describe_workflow`
   → `mcp_server/tools/workflows.py` + register in `app.py`
4. Add `get_workflow_parameter_schema`
   → same module
5. Add `validate_workflow_parameters`
   → same module
6. Add `render_workflow_command`
   → same module
7. Add unit tests
   → `mcp_server/tests/test_workflow_*.py`
8. Add the UCSF E2E dry-run contract test
   → `mcp_server/tests/test_ucsf_workflow_e2e.py`
9. Update `mcp_server/README.md`

## New Files To Create

```
mcp_server/
├── schemas/
│   └── workflow.py                              # NEW: Pydantic models for workflow metadata
├── tools/
│   └── workflows.py                             # NEW: workflow tool handlers
├── workflows/
│   └── ucsf_star_suite_production.yaml          # NEW: UCSF workflow schema
├── config.yaml                                  # MODIFY: add workflows: section
├── config.py                                    # MODIFY: load workflow schemas
├── schemas/config.py                            # MODIFY: add WorkflowConfig model
├── schemas/responses.py                         # MODIFY: add workflow response models
├── app.py                                       # MODIFY: register 5 new tools
├── README.md                                    # MODIFY: document new tools
└── tests/
    ├── test_workflow_config.py                   # NEW
    ├── test_workflow_discovery.py                # NEW
    ├── test_workflow_validation.py               # NEW
    ├── test_workflow_render.py                   # NEW
    ├── test_workflow_integration.py              # NEW
    └── test_ucsf_workflow_e2e.py                 # NEW
```

## Config YAML Addition

Add to `mcp_server/config.yaml`:

```yaml
workflows:
  - id: "ucsf_star_suite_production"
    title: "UCSF STAR-suite Production Workflow"
    summary: "Corrected UCSF GEX+guide workflow with optional downstream and remote CellBender stages."
    entry_script: "scripts/paper/run_ucsf_corrected_production_workflow.sh"
    kind: "shell_workflow"
    schema_file: "mcp_server/workflows/ucsf_star_suite_production.yaml"
```

## Acceptance Criteria

From runbook §15:

1. MCP server can list supported workflows separately from test scripts
2. MCP server can return machine-readable parameter schema for the UCSF workflow
3. MCP server can validate structured params for that workflow
4. MCP server can render those params into the exact UCSF wrapper command
5. Unit tests cover config loading, discovery, validation, and rendering
6. E2E-style UCSF dry-run contract test passes
7. Existing discovery / preflight / executor tests still pass

## Key Constraints

- **Do not replace** the existing `scripts:` config or `run_script` flow
- **Do not parse Bash** — hand-author the schema from the reference above
- **Do not block on STAR core param export** from `Parameters.cpp`
- **Argv arrays are canonical** — not shell string assembly
- **Security**: All path params must go through `is_path_allowed()` during validation
- **Existing tests must keep passing** — no breaking changes

## Testing Notes

- Use `loaded_config` fixture pattern from `conftest.py`
- Create temp dirs/files for path params in validation tests
- Patch `trusted_roots` to include temp fixture root
- Use `dry_run=true` in E2E tests — never execute the real workflow
- Compare rendered output against argv arrays, not shell strings
- The E2E test flow is: schema load → validate → render → preflight (reuse
  existing `run_preflight` with the rendered script + args)
