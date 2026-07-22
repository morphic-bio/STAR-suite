# Review: MCP Workflow Parameters Implementation

> **Date**: 2026-04-06
> **Runbook**: `plans/mcp_workflow_parameters_runbook.md`
> **Handoff**: `docs/HANDOFF_MCP_WORKFLOW_PARAMETERS.md`
> **Branch**: `master` (not yet committed)

## What Was Done

Extended the existing STAR-suite MCP server with a structured workflow parameter
layer. The server can now serve machine-readable workflow metadata, validate
structured parameter requests, and render them into shell commands -- all without
parsing Bash scripts or executing anything.

The first supported workflow is the UCSF STAR-suite production wrapper
(`scripts/paper/run_ucsf_corrected_production_workflow.sh`), with 21 parameters
across 6 groups, 5 constraint rules, and deterministic argv rendering.

## New MCP Tools

| Tool | Auth | Purpose |
|------|------|---------|
| `list_workflows` | discovery | Return all configured workflow templates |
| `describe_workflow` | discovery | Full metadata: stages, parameter groups, caveats |
| `get_workflow_parameter_schema` | discovery | Machine-readable param definitions, types, defaults, constraints |
| `validate_workflow_parameters` | required | Validate structured params (types, paths, mutual exclusion, etc.) |
| `render_workflow_command` | required | Convert validated params into argv array + shell preview |

The workflow tools live alongside (not replacing) the existing `scripts:`
allowlist and `run_script` / `preflight` tools. An agent can compose them:
`list_workflows` -> `get_workflow_parameter_schema` -> `validate_workflow_parameters`
-> `render_workflow_command` -> optionally `preflight` / `run_script`.

## Extensibility: Adding New Workflows

Adding a new workflow requires **zero code changes**. The process is:

1. Create a YAML schema file in `mcp_server/workflows/<workflow_id>.yaml`
2. Add a config entry under `workflows:` in `mcp_server/config.yaml`
3. Call `reload_config` (or restart the server)

The YAML schema is the complete contract. All rendering behavior is driven by
schema-level fields, not hardcoded logic in the renderer:

| Schema field | Purpose | Example |
|-------------|---------|---------|
| `env_var` | Export param as named env var in render response | `env_var: "DOWNSAMPLE_SEED"` |
| `skip_when_default` | Omit CLI flag when value equals schema default | `skip_when_default: true` on `downsample_reads` (default 0 = disabled) |
| `is_output_root` | Report this param's value as the workflow output root | `is_output_root: true` on `out_root` |
| `path_must_exist` | Validate file/dir existence at validation time | `path_must_exist: true` on input paths |

Full authoring spec: `mcp_server/workflows/AUTHORING.md`

### Minimal new workflow (no code changes needed)

```yaml
# mcp_server/workflows/my_pipeline.yaml
id: "my_pipeline"
title: "My Pipeline"
summary: "Runs the analysis pipeline."
entry_script: "scripts/run_my_pipeline.sh"

parameters:
  - name: "input_dir"
    cli_flag: "--input-dir"
    type: "directory"
    required: true
    path_must_exist: true
  - name: "threads"
    cli_flag: "--threads"
    type: "int"
    default: 8
  - name: "output_dir"
    cli_flag: "--output-dir"
    type: "directory"
    is_output_root: true

constraints:
  - kind: "positive"
    params: ["threads"]

rendering:
  flag_order: ["input_dir", "threads", "output_dir"]
```

```yaml
# In config.yaml, add:
workflows:
  - id: "my_pipeline"
    title: "My Pipeline"
    summary: "Runs the analysis pipeline."
    entry_script: "scripts/run_my_pipeline.sh"
    schema_file: "mcp_server/workflows/my_pipeline.yaml"
```

## Files Created

| File | Purpose |
|------|---------|
| `mcp_server/schemas/workflow.py` | Pydantic models: `WorkflowSchema`, `WorkflowParameterDef`, `WorkflowConstraint`, `WorkflowStage`, `WorkflowParameterGroup`, `WorkflowRenderingRule` |
| `mcp_server/tools/workflows.py` | Tool handlers for all 5 workflow tools |
| `mcp_server/workflows/ucsf_star_suite_production.yaml` | UCSF workflow parameter schema (hand-authored from the actual shell script) |
| `mcp_server/workflows/AUTHORING.md` | Authoring spec: how to add a new workflow, full field reference, checklist |
| `mcp_server/tests/test_workflow_config.py` | Config loading, schema resolution, error cases (6 tests) |
| `mcp_server/tests/test_workflow_discovery.py` | list/describe/schema tools, field coverage (10 tests) |
| `mcp_server/tests/test_workflow_validation.py` | Type, path, enum, constraint, normalization checks (23 tests) |
| `mcp_server/tests/test_workflow_render.py` | Bool flags, CSV, int, declarative rendering hints, deterministic order (24 tests) |
| `mcp_server/tests/test_workflow_integration.py` | Validate-then-render pipeline composition (4 tests) |
| `mcp_server/tests/test_ucsf_workflow_e2e.py` | E2E contract test with real UCSF schema + temp fixtures (15 tests) |

## Files Modified

| File | Change |
|------|--------|
| `mcp_server/schemas/config.py` | Added `WorkflowConfig` model, `get_workflow()` helper on `MCPConfig`, `workflows` field |
| `mcp_server/config.py` | Added `_load_workflow_schemas()`, `get_workflow_schema()`, `get_workflow_schemas()` accessors; workflow schemas loaded alongside config |
| `mcp_server/schemas/responses.py` | Added workflow response models: `WorkflowInfo`, `ListWorkflowsResponse`, `WorkflowStageInfo`, `DescribeWorkflowResponse`, `ParameterInfo`, `ConstraintInfo`, `ParameterGroupInfo`, `WorkflowParameterSchemaResponse`, `ValidateWorkflowResponse`, `RenderWorkflowResponse` |
| `mcp_server/app.py` | Registered 5 new `@mcp.tool()` functions with auth, error handling |
| `mcp_server/config.yaml` | Added `workflows:` section with UCSF entry |
| `mcp_server/tests/conftest.py` | Reset `_workflow_schemas` in test teardown |
| `mcp_server/README.md` | Documented workflow tools, updated project structure, added Phase 6 |

## UCSF Workflow Schema Coverage

The schema covers all top-level parameters from
`scripts/paper/run_ucsf_corrected_production_workflow.sh`:

| Group | Parameters |
|-------|-----------|
| Selection | `samples` (string_list), `all_samples` (bool) |
| Inputs | `dataset_root`, `feature_ref`, `genome_dir`, `solo_cb_whitelist`, `cr_whitelist`, `star_bin` |
| Performance | `threads`, `cellbender_cpu_cores`, `downsample_reads` (skip_when_default), `downsample_seed` (env_var) |
| Execution Mode | `cellbender_gpu`, `trim_qc`, `star_only`, `dry_run` |
| Transfer | `globus_src_endpoint`, `globus_dst_endpoint`, `globus_dst_root`, `globus_poll_seconds` |
| Output | `out_root` (is_output_root) |

### Constraints encoded

- **Mutual exclusion**: `samples` vs `all_samples`
- **Group required**: all three Globus params must be supplied together (warns if partial)
- **Dependency**: `cellbender_gpu` has no effect when `star_only` is set
- **Positive**: `threads`, `cellbender_cpu_cores`, `globus_poll_seconds`
- **Non-negative**: `downsample_reads`

### Rendering rules (all schema-driven, no hardcoded param names)

- Booleans emit flag-only when true, omitted when false
- `skip_when_default: true` omits flag when value equals schema default (e.g. `downsample_reads=0`)
- `env_var: "DOWNSAMPLE_SEED"` exports value as env override (matches script's `export`)
- `is_output_root: true` reports the param value as `output_root` in the render response
- Flags emitted in deterministic order defined in `rendering.flag_order`
- `omit_absent_optionals` and `csv_style` configurable per workflow schema

## Validation Coverage

The validator checks:

- Missing required parameters
- Unknown parameters
- Type correctness (int, float, bool, enum, string, string_list, file, directory)
- Enum membership
- Path existence and trusted-root enforcement for `file`/`directory` params with `path_must_exist`
- Mutual exclusion constraints
- Group-required constraints (warn on partial set)
- Dependency constraints (warn when gated param is irrelevant)
- Positive integer constraints
- Non-negative integer constraints

All checks produce structured errors/warnings in the response. `check_paths`
can be set to `false` to skip filesystem checks during dry planning.

## Test Results

```
189 passed in 11.61s
```

- 82 new workflow-specific tests
- 107 existing tests all still pass (auth, config, discovery, executor, preflight, utils)
- 0 failures, 0 warnings

### Test categories

| Category | Tests | What's covered |
|----------|-------|----------------|
| Config loading | 6 | YAML load, schema resolution, missing file error, invalid config rejection |
| Discovery tools | 10 | list/describe/schema responses, parameter metadata fields, group typing, constraint details |
| Validation | 23 | Required, unknown, type, enum, path existence, trusted roots, mutual exclusion, group required, positive/non-negative, normalization |
| Rendering | 24 | Bool flags, CSV join, int rendering, declarative `env_var`/`skip_when_default`/`is_output_root`, empty-string omission, shell preview, deterministic order |
| Integration | 4 | Validate-then-render pipeline, defaults flow-through, preflight compatibility |
| UCSF E2E | 15 | Real schema, temp fixtures, declarative field assertions, full contract flow (schema -> validate -> render -> assert argv) |

## Design Decisions

1. **Curated schema, not Bash parsing.** The UCSF schema was hand-authored from
   the actual script source. This is deliberate -- a brittle partial parser would
   be worse than a correct curated schema.

2. **Workflows parallel to scripts.** The `workflows:` config section sits
   alongside `scripts:`, not replacing it. Existing `run_script` / `preflight`
   paths are unchanged.

3. **Argv arrays are canonical.** The renderer produces `list[str]`, not shell
   strings. The `shell_preview` field is derived from argv via `shlex.quote`.

4. **Validation is separate from rendering.** `render_workflow_command` does not
   call validation internally. The intended flow is validate -> render, giving
   agents explicit control over the pipeline.

5. **Top-level params only.** The schema covers only the production wrapper's
   flags. Sub-script parameters (batch runner, downstream, remote CellBender)
   are documented in the handoff for a later pass.

6. **No hardcoded param names in the renderer.** All workflow-specific rendering
   behavior (env var export, skip-when-default, output root designation) is
   driven by declarative fields on the parameter schema, not by `if p_name ==`
   checks. This makes the system genuinely extensible: new workflows with
   different env vars, different skip rules, and different output parameters
   work without touching the renderer code.

## Acceptance Criteria (from runbook section 15)

| # | Criterion | Status |
|---|-----------|--------|
| 1 | MCP server can list supported workflows separately from test scripts | Done |
| 2 | MCP server can return machine-readable parameter schema for the UCSF workflow | Done |
| 3 | MCP server can validate structured params for that workflow | Done |
| 4 | MCP server can render those params into the exact UCSF wrapper command | Done |
| 5 | Unit tests cover config loading, discovery, validation, and rendering | Done (63 tests) |
| 6 | E2E-style UCSF dry-run contract test passes | Done (15 tests) |
| 7 | Existing discovery / preflight / executor tests still pass | Done (107 tests unchanged) |

## What Is Not Included

Per runbook section 14 (non-goals):

- No generic shell script parser
- No STAR core parameter registry export from `Parameters.cpp`
- No replacement of existing `scripts:` config
- No direct `bwb-nextflow-utils` integration in this repo
- No real UCSF production execution in tests
- No `run_workflow` tool (can be added later; current flow composes existing tools)
