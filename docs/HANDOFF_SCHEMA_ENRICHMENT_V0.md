# Handoff: STAR Schema Enrichment v0

**Branch:** `star-schema-enrichment-v0`
**Date:** 2026-04-12
**Status:** Implementation complete, not yet committed

## Goal

Enrich STAR MCP parameter schemas so planners and browser forms can render
better structured UX and make safer parameter decisions. Changes are additive
and non-breaking. No product-routing, browser session, or cross-step recipe
logic was added to STAR-suite.

## What changed

### 1. New optional fields on WorkflowParameterDef / ParameterInfo

All fields default to `None` or `False`. Existing consumers that only read
`name`, `type`, `required`, `default`, `description` are unaffected.

| Field | Type | Purpose |
|-------|------|---------|
| `label` | `str?` | Stable human-readable label for form headers (e.g. "Reference FASTA") |
| `help` | `str?` | Extended guidance beyond the short `description` |
| `example` | `str?` | Example value for placeholders and documentation |
| `widget_hint` | `str?` | UI rendering hint: `text`, `number`, `path`, `readonly`, `checkbox`, `select`, `textarea`, `hidden` |
| `aliases` | `list[str]?` | Canonical alternate names (e.g. camelCase STAR flags: `genomeDir` for `genome_dir`) |
| `advanced` | `bool` | UIs may collapse into an advanced section |
| `display_order` | `int?` | Explicit sort key within group |
| `min_value` | `float?` | Minimum numeric value hint |
| `max_value` | `float?` | Maximum numeric value hint |

All are advisory. Clients may ignore them.

### 2. Group description

`WorkflowParameterGroup` / `ParameterGroupInfo` gained a `description` field
(default `""`). Used for section-level help text in forms.

### 3. Structured field_errors on validation

`ValidateWorkflowResponse` gained `field_errors: list[FieldValidationError]`
alongside the existing flat `errors` and `warnings` lists.

Each `FieldValidationError` has:
- `field` -- parameter name the error applies to
- `message` -- same string as the flat error/warning
- `kind` -- `"error"` or `"warning"`

This lets form UIs attach validation messages directly to fields instead of
regex-matching against flat strings. The flat lists remain unchanged for
backward compatibility.

### 4. Alias resolution in validation

New internal function `_resolve_aliases()` normalizes aliased parameter names
before validation. If a consumer passes `{"genomeDir": "/tmp/idx"}`, it is
resolved to `{"genome_dir": "/tmp/idx"}` using the `aliases` metadata from
the schema. `normalized_params` always uses canonical names.

**Collision detection**: if both an alias and its canonical name appear in the
same parameter bag (e.g. `genomeDir=/tmp/a` + `genome_dir=/tmp/b`), validation
returns an error instead of silently picking one.

**Multi-operand safety**: aliases are only assigned to parameters that own a
CLI flag exclusively. Parameters sharing a flag via `operand_group` (e.g.
`read_files_mate1`/`read_files_mate2` both use `--readFilesIn`, and
`out_sam_kind`/`out_sam_sort` both use `--outSAMtype`) do **not** carry
aliases, because a single alias key cannot resolve to two canonical parameters.

### 5. Enriched YAML schemas

Two workflows were enriched:

**star_genome_generate** (4 parameters):
- All params have `label`, `help`, `example`, `widget_hint`, `aliases`, `display_order`
- `run_thread_n` has `min_value: 1`, `max_value: 128`, `advanced: true`
- `run_mode` has `widget_hint: readonly`
- Group `core` has a `description`
- Aliases: `runMode`, `genomeDir`, `genomeFastaFiles`, `runThreadN`

**star_bulk_pe_batch** (10 parameters):
- All params enriched with the same field set
- Groups reorganized from one flat `core` group into three: `inputs`, `output`, `performance` -- each with a `description`
- Aliases: `runMode`, `genomeDir`, `outFileNamePrefix`, `outFileNamePrefixAuto`, `batchMode`, `runThreadN`
- No aliases on multi-operand params: `read_files_mate1`/`read_files_mate2` (`--readFilesIn`), `out_sam_kind`/`out_sam_sort` (`--outSAMtype`)
- `display_order` is sequential 0-9 across all params

## Files modified

| File | Change |
|------|--------|
| `mcp_server/schemas/workflow.py` | 9 new optional fields on `WorkflowParameterDef`; `description` on `WorkflowParameterGroup` |
| `mcp_server/schemas/responses.py` | Matching fields on `ParameterInfo` and `ParameterGroupInfo`; new `FieldValidationError` model; `field_errors` on `ValidateWorkflowResponse` |
| `mcp_server/tools/workflows.py` | Surface new fields in `get_workflow_parameter_schema` and `describe_workflow`; `_resolve_aliases()` + structured `field_errors` in `validate_workflow_parameters` |
| `mcp_server/workflows/star_genome_generate.yaml` | Enriched all params and group |
| `mcp_server/workflows/star_bulk_pe_batch.yaml` | Enriched all params; reorganized groups |
| `mcp_server/tests/test_schema_enrichment.py` | New -- 30 tests |

## Files not modified

- `mcp_server/app.py` -- no changes needed; it already passes through all fields
- Other workflow YAMLs -- only `star_genome_generate` and `star_bulk_pe_batch` were in scope
- No changes to the Launchpad, browser, or any product-routing code

## Test results

- **31 new tests** in `test_schema_enrichment.py` -- all pass
- **504 existing tests** (including the 31 new) -- all pass
- **2 pre-existing failures** (`test_validate_render_preflight_composition` in two files) -- confirmed to fail on `master` due to low disk space on the test host; not related to this branch

## Canonical vs advisory fields

**Canonical** (breaking to remove, depended on by existing consumers):
- `name`, `cli_flag`, `type`, `required`, `default`, `description`, `choices`, `category`, `stage`, `source`
- `valid`, `normalized_params`, `warnings`, `errors` on validation response

**Advisory** (additive, clients may ignore):
- `label`, `help`, `example`, `widget_hint`, `aliases`, `advanced`, `display_order`, `min_value`, `max_value`
- `description` on parameter groups
- `field_errors` on validation response

## What the thin client can now do

1. Use `label` instead of raw `name` in form headers
2. Use `widget_hint` to select input widgets (`readonly` for `runMode`, `path` for directories, `number` for integers)
3. Use `help` for tooltip/popover text
4. Use `example` as `placeholder` on text inputs
5. Use `aliases` for parameter name resolution from camelCase STAR CLI args
6. Use `advanced` to collapse tuning params by default
7. Use `display_order` for deterministic field ordering
8. Use `group.description` as section-level help text
9. Use `field_errors` to render inline validation messages per field
10. Use `min_value`/`max_value` for HTML `min`/`max` attributes

## What was explicitly not done

- No browser/product/coordinator logic added to STAR-suite
- No recipe construction semantics
- No cross-workflow composition
- No approval/review semantics
- No changes to other workflow YAMLs beyond the two in scope
- No speculative aliases (only real STAR camelCase flag names)

## Suggested next steps

1. Commit and merge this branch
2. Enrich remaining workflow YAMLs (`star_scrna_solo_droplet`, `star_perturb_cr_compat`, `star_flex_fixed_rna`, `ucsf_star_suite_production`) with the same metadata pattern
3. Consume the enriched schema in the thin client form rendering
4. Consume `field_errors` in the thin client validation display
5. Consider adding `aliases` resolution to the thin client's parameter submission path
