# Workflow Schema Authoring Guide

How to add a new workflow to the STAR-suite MCP server.

## Overview

A workflow is a structured parameter contract for a shell script. Adding one
requires two things: a YAML schema file and a config entry. No code changes are
needed.

## Step 1: Create the schema YAML

Create a file in `mcp_server/workflows/<workflow_id>.yaml`.

### Minimal example

```yaml
id: "my_workflow"
title: "My Workflow"
summary: "One-line description of what this workflow does."
entry_script: "scripts/my_workflow.sh"

parameters:
  - name: "input_dir"
    cli_flag: "--input-dir"
    type: "directory"
    required: true
    description: "Input data directory."
    path_must_exist: true
    category: "inputs"

  - name: "threads"
    cli_flag: "--threads"
    type: "int"
    default: 8
    description: "Thread count."
    category: "performance"

  - name: "dry_run"
    cli_flag: "--dry-run"
    type: "bool"
    default: false
    description: "Dry run mode."
    category: "mode"

constraints:
  - kind: "positive"
    params: ["threads"]
    message: "Must be a positive integer."

rendering:
  flag_order: ["input_dir", "threads", "dry_run"]
```

### Full field reference

#### Top-level fields

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `id` | yes | | Unique workflow identifier (snake_case) |
| `title` | yes | | Human-readable display name |
| `summary` | yes | | One-line description |
| `kind` | no | `"shell_workflow"` | Workflow kind |
| `entry_script` | yes | | Path to entry script (relative to repo root) |
| `supported_modes` | no | `["local", "dry-run-capable"]` | Execution modes |
| `caveats` | no | `[]` | Known caveats (shown to agents) |
| `default_output_layout` | no | `""` | Description of output directory structure |
| `stages` | no | `[]` | Semantic stages within the workflow |
| `parameter_groups` | no | `[]` | Ordered groups for display |
| `parameters` | no | `[]` | Parameter definitions |
| `constraints` | no | `[]` | Cross-parameter constraints |
| `rendering` | no | defaults | Rendering configuration |

#### Parameter fields

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `name` | yes | | Identifier (snake_case, must be unique within workflow) |
| `cli_flag` | yes | | CLI flag string (e.g. `--my-flag`) |
| `type` | yes | | One of: `string`, `int`, `float`, `bool`, `enum`, `file`, `directory`, `string_list` |
| `required` | no | `false` | Whether the parameter must be supplied |
| `default` | no | `null` | Default value; `null` means no default |
| `description` | no | `""` | Human-readable description |
| `choices` | no | `null` | Allowed values (only for `enum` type) |
| `repeatable` | no | `false` | Whether the flag can appear multiple times |
| `path_must_exist` | no | `false` | For `file`/`directory` types: validate existence at validation time |
| `category` | no | `"general"` | Logical group (matches `parameter_groups.name`) |
| `stage` | no | `"top_level"` | Which workflow stage owns this parameter |
| `source` | no | `"workflow_wrapper"` | Which script layer defines this parameter |
| `env_var` | no | `null` | If set, export the value as this env var in the rendered command |
| `skip_when_default` | no | `false` | If true, omit the CLI flag when the value equals the schema default |
| `is_output_root` | no | `false` | If true, report this parameter's value as the workflow `output_root` |

#### Parameter types

| Type | CLI rendering | Validation |
|------|--------------|------------|
| `string` | `--flag value` | Non-empty string |
| `int` | `--flag 42` | Must parse as integer |
| `float` | `--flag 0.5` | Must parse as float |
| `bool` | `--flag` (when true, omitted when false) | Boolean-like value |
| `enum` | `--flag value` | Must be in `choices` list |
| `file` | `--flag /path/to/file` | Optionally checked for existence |
| `directory` | `--flag /path/to/dir` | Optionally checked for existence |
| `string_list` | `--flag "a,b,c"` | CSV or list |

#### Rendering hints

These three fields on parameters control how the renderer handles special cases
without hardcoding workflow-specific logic:

**`env_var`**: When set, the parameter's value is added to `env_overrides` in the
render response under the given name. Use this for parameters that the entry
script exports as environment variables rather than passing as CLI flags.

```yaml
  - name: "my_seed"
    cli_flag: "--seed"
    type: "int"
    default: 1
    env_var: "MY_SEED"  # Will appear in env_overrides as MY_SEED=<value>
```

**`skip_when_default`**: When true, the CLI flag is omitted from the rendered
argv when the value equals the schema default. Use this for parameters where
the default value means "disabled" and the script already handles the default
internally.

```yaml
  - name: "downsample_reads"
    cli_flag: "--downsample-reads"
    type: "int"
    default: 0              # 0 means no downsampling
    skip_when_default: true  # Don't emit --downsample-reads 0
```

**`is_output_root`**: Marks which parameter holds the workflow's output
directory. Its value is returned as `output_root` in the render response so
agents can find where results will land.

```yaml
  - name: "out_dir"
    cli_flag: "--out-dir"
    type: "directory"
    is_output_root: true
```

#### Constraint kinds

| Kind | Behavior |
|------|----------|
| `mutual_exclusion` | Error if more than one of `params` is set |
| `group_required` | Warn if some but not all of `params` are set |
| `dependency` | Warn if `params[0]` is set but `params[1]` makes it irrelevant |
| `positive` | Error if any of `params` is <= 0 |
| `non_negative` | Error if any of `params` is < 0 |

#### Rendering configuration

| Field | Default | Description |
|-------|---------|-------------|
| `bool_style` | `"flag_only"` | How booleans render (currently only `flag_only` is supported) |
| `flag_order` | `[]` | Deterministic emission order (by param name); if empty, uses declaration order |
| `omit_absent_optionals` | `true` | Skip flags for optional params with no value |
| `csv_style` | `"quoted"` | How `string_list` values are joined |

#### Stages

Stages are informational metadata describing the workflow's semantic phases.
They don't affect validation or rendering but help agents understand the
workflow structure.

```yaml
stages:
  - name: "alignment"
    title: "STAR Alignment"
    script: "scripts/align.sh"
    description: "Run STAR alignment."
  - name: "qc"
    title: "Quality Control"
    script: "scripts/qc.sh"
    gated_by: "skip_qc"  # This stage is skipped when skip_qc=true
```

#### Parameter groups

Groups control display ordering. Each group has a name, title, and a list of
parameter names. Parameters can only belong to one group.

```yaml
parameter_groups:
  - name: "inputs"
    title: "Input Paths"
    parameters: ["input_dir", "reference"]
  - name: "performance"
    title: "Performance"
    parameters: ["threads", "memory"]
```

## Step 2: Register in config.yaml

Add an entry under the `workflows:` section:

```yaml
workflows:
  - id: "my_workflow"
    title: "My Workflow"
    summary: "One-line description."
    entry_script: "scripts/my_workflow.sh"
    kind: "shell_workflow"
    schema_file: "mcp_server/workflows/my_workflow.yaml"
```

That's it. The MCP server will automatically pick up the new workflow on the
next config load (or call `reload_config`).

## Step 3: Verify

```python
# List workflows -- your new one should appear
client.call_tool("list_workflows", {})

# Get the parameter schema
client.call_tool("get_workflow_parameter_schema", {
    "workflow_id": "my_workflow"
})

# Validate params
client.call_tool("validate_workflow_parameters", {
    "workflow_id": "my_workflow",
    "params": {"input_dir": "/path/to/data", "threads": 4},
    "check_paths": True,
})

# Render command
client.call_tool("render_workflow_command", {
    "workflow_id": "my_workflow",
    "params": {"input_dir": "/path/to/data", "threads": 4, "dry_run": True},
})
```

## Step 4: Add tests (recommended)

Create `mcp_server/tests/test_<workflow_id>_e2e.py` following the pattern in
`test_ucsf_workflow_e2e.py`:

1. Use the real schema YAML with temporary directories and stub executables
2. Test the full contract: schema load -> validate -> render -> assert argv
3. Use `dry_run=true` -- never execute the real workflow in tests

## Checklist

- [ ] Schema YAML created in `mcp_server/workflows/`
- [ ] Config entry added to `config.yaml` under `workflows:`
- [ ] `list_workflows` returns the new workflow
- [ ] `get_workflow_parameter_schema` returns correct params, groups, constraints
- [ ] `validate_workflow_parameters` catches missing required params, bad types, constraint violations
- [ ] `render_workflow_command` produces correct argv with deterministic flag order
- [ ] E2E test passes with temp fixtures and dry-run mode
