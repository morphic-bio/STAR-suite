# STAR Launchpad — Composition tab (Phase C) — shareable summary

**Date:** 2026-04-20  
**Repo:** `STAR-suite` (coordination context: BioDepot stack; canonical composition logic remains in `bwb-nextflow-utils`)

## One-line outcome

Launchpad now has a third tab, **Composition**, that lists constructor profiles, shows backend-computed recipe inputs, builds execution-recipe **drafts**, and generates **review artifacts** (validate / normalize / preview / script / manifest) without duplicating profile or constructor semantics in `STAR-suite`.

## What was built

### HTTP API (under `/launchpad/api/composition/`)

| Method & path | Role |
|---------------|------|
| `GET .../profiles` | Lists profiles (`list_constructor_profiles`). |
| `GET .../profiles/{profile_id}` | Describes one profile (`describe_constructor_profile`). |
| `POST .../recipe-inputs` | Recipe input schema for a profile (`describe_constructor_profile_recipe_inputs`). Optional `workflow_ir_by_id` / `component_binding`. Bulk STAR defaults: server merges **STAR-suite** IR from `mcp_server/launchpad/composition_ir/*.json` + standard `mkindex` / `align` binding (see `composition_ir/README.md`). |
| `POST .../draft` | `build_headless_recipe_draft` + optional `recipe_input_values` (literals wired structurally, not profile-hard-coded). |
| `POST .../artifacts` | Same draft path, then validate → normalize → preview → script → manifest (side-effect free; no execution). |

**Capabilities:** `GET /launchpad/api/capabilities` includes `composition_utils_ready` (same readiness as Script Lane: `bwb-nextflow-utils` import via `BWB_NEXTFLOW_UTILS_ROOT` or sibling checkout).

### UI

- New **Composition** tab (Workflows and Script Lane unchanged).
- Profile selection, archetype / `helper_ids_by_role` / notes, required vs optional inputs from backend `kind` (with optional/advanced section and JSON overrides for power users).
- Optional/defaulted inputs stay out of the draft until the user explicitly sets them; they are not auto-injected just because the backend advertises a default.
- Invalid override JSON is rejected consistently on reload / draft / artifacts (no silent fallback).
- **Build draft** and **Generate artifacts**; results shown in-page (draft JSON, normalized recipe, preview, generated script, manifest).

### Bridge / ownership

- **`mcp_server/launchpad/composition_bridge.py`**: locates `bwb-nextflow-utils`, applies bulk STAR IR defaults from **`mcp_server/launchpad/composition_ir/`** (STAR-owned, aligned with `mcp_server/workflows/star_*.yaml`), delegates constructor calls to `mcp_workflow_services`.
- **`mcp_server/launchpad/script_lane_bridge.py`**: adds `get_mcp_workflow_services()` shared by Script Lane and Composition.

## First-wave profiles (verified in tests)

At minimum, end-to-end support aligns with:

- Bulk: `star_bulk_de_deseq2_default_v0`
- scRNA: `ucsf_scrna_default_v0`, `ucsf_scrna_mex_default_v0`, `ucsf_scrna_full_default_v0`, `ucsf_scrna_mex_full_default_v0`

The profile list comes from the backend; additional profiles appear as `bwb-nextflow-utils` adds them.

## Intentionally not in this slice

- Executing composed recipes from Launchpad.
- Redesign of the Workflows tab or Script Lane.
- Planner/chat integration, generic arbitrary composition authoring, new package/helper contracts in `bwb-nextflow-utils` (unless a blocking bug).

## Files touched

- `mcp_server/launchpad/composition_bridge.py` (new)
- `mcp_server/launchpad/composition_ir/*.json`, `README.md` (STAR bulk composition IR)
- `mcp_server/launchpad/script_lane_bridge.py`
- `mcp_server/launchpad/api.py`
- `mcp_server/launchpad/static/app.js`
- `mcp_server/launchpad/static/index.html`
- `mcp_server/tests/test_launchpad.py`
- `docs/HANDOFF_LAUNCHPAD_SPA.md` (Composition section)
- This file: `docs/launchpad_composition_phase_c_summary.md`

## How to verify

```bash
cd /path/to/STAR-suite
python -m pytest mcp_server/tests/test_launchpad.py -q --tb=short
```

Composition tests need a usable `bwb-nextflow-utils` checkout (`BWB_NEXTFLOW_UTILS_ROOT` or `../bwb-nextflow-utils` next to `STAR-suite`).

## Remaining gap (product)

To make composition **runnable** from Launchpad: connect a trusted execution path (e.g. Temporal / serial recipe runner / MCP run tools) with explicit policy for `recipe_input_values` and host paths—out of scope for this Phase C UI/API slice.
