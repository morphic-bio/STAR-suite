# Handoff: STAR-suite Launchpad SPA

**Branch:** `launchpad`
**Date:** 2026-04-08
**Status:** Branch created, no implementation yet

## Goal

Build a form-based single-page application ("Launchpad") served from the same
HTTP server that hosts the MCP protocol endpoints. The SPA lets users select a
STAR-suite workflow, fill in parameters via a guided form, and receive a
validated, copy-pasteable shell command they can run on their own data.

This is **not** a remote execution UI. It generates commands; users run them.

## Why

- Workflow schemas already exist in `mcp_server/workflows/*.yaml` with typed
  parameters, defaults, constraints, and rendering rules.
- The MCP server already validates and renders commands via
  `validate_workflow_parameters` and `render_workflow_command`.
- A form-based frontend over this data eliminates the need for users to read
  script source to figure out flags and defaults.

## Architecture

### Server side

The MCP server (`mcp_server/app.py`) builds a Starlette app in
`build_http_app()` (line 88). It currently mounts MCP streamable-HTTP at `/`
and SSE at `/sse` + `/messages`. The Launchpad SPA should be mounted alongside
these on the same Starlette app.

**Mounting strategy:**

```python
# In build_http_app(), add Launchpad routes:
from starlette.routing import Route, Mount
from starlette.staticfiles import StaticFiles

# API endpoints for the SPA (JSON, no auth required for public workflows)
Route("/launchpad/api/workflows", endpoint=launchpad_list_workflows),
Route("/launchpad/api/workflows/{workflow_id}/schema", endpoint=launchpad_get_schema),
Route("/launchpad/api/workflows/{workflow_id}/validate", endpoint=launchpad_validate, methods=["POST"]),
Route("/launchpad/api/workflows/{workflow_id}/render", endpoint=launchpad_render, methods=["POST"]),

# Static assets + SPA fallback
Mount("/launchpad", app=StaticFiles(directory="mcp_server/launchpad/static", html=True)),
```

**Key constraint:** The MCP streamable-HTTP transport occupies `POST /`. The
Launchpad must not conflict. All Launchpad routes live under `/launchpad/`.

### API layer

Create `mcp_server/launchpad/api.py` with thin Starlette endpoint functions.
These call the existing workflow tool functions directly (the same ones MCP
tools call):

| Endpoint | Method | Calls | Auth |
|----------|--------|-------|------|
| `/launchpad/api/workflows` | GET | `list_workflows(authenticated=False)` | None |
| `/launchpad/api/workflows/{id}/schema` | GET | `get_workflow_parameter_schema(id)` | None |
| `/launchpad/api/workflows/{id}/validate` | POST | `validate_workflow_parameters(id, params, check_paths=False)` | None |
| `/launchpad/api/workflows/{id}/render` | POST | `render_workflow_command(id, params)` | None |

**Important:** `check_paths=False` for validate, since users are filling in
paths for their own machines — the server cannot verify them. Render does not
execute anything.

**Security:** These endpoints expose only public workflows (visibility=public).
They do not execute scripts, touch the filesystem, or return host-specific
paths. The rendered command uses relative paths from the workflow schema; the
user supplies their own absolute paths. No auth token is needed.

### Frontend

A self-contained SPA in `mcp_server/launchpad/static/`. No build step, no
npm, no bundler. Use vanilla HTML + CSS + JS (or a single-file Vue/Alpine.js
if you prefer, loaded from CDN).

**File structure:**

```
mcp_server/launchpad/
    __init__.py
    api.py              # Starlette endpoints
    static/
        index.html      # SPA shell
        app.js          # Form logic
        style.css       # Styling
```

**Page flow:**

1. **Workflow picker** — GET `/launchpad/api/workflows`, show cards/list of
   public workflows with title + summary.
2. **Parameter form** — GET `/launchpad/api/workflows/{id}/schema`, render a
   form with:
   - Grouped sections (from `parameter_groups`)
   - Typed inputs: text for string, number for int/float, checkbox for bool,
     select for enum, file/directory text inputs with appropriate labels
   - Required markers, default values pre-filled, descriptions as help text
   - `skip_when_default` params shown but visually de-emphasized
3. **Validate** — POST params to `/launchpad/api/workflows/{id}/validate`.
   Show errors inline, warnings as notices. On success, enable "Generate
   Command".
4. **Render** — POST params to `/launchpad/api/workflows/{id}/render`. Display:
   - The `shell_preview` in a copyable code block
   - `env_overrides` as `export VAR=val` lines above the command
   - `output_root` if available
   - A "Copy to clipboard" button
5. **Constraint display** — Show `constraints` from the schema as
   inline warnings (e.g., mutual exclusion between `--samples` and
   `--all-samples`).

### Form generation from schema

The parameter schema has everything needed to generate form fields:

```javascript
// Pseudo-code for field rendering
function renderField(param) {
    switch (param.type) {
        case "bool":    return checkbox(param.name, param.default);
        case "int":     return numberInput(param.name, param.default, { step: 1 });
        case "float":   return numberInput(param.name, param.default, { step: 0.01 });
        case "enum":    return select(param.name, param.choices, param.default);
        case "string":  return textInput(param.name, param.default);
        case "file":    return textInput(param.name, param.default, { placeholder: "/path/to/file" });
        case "directory": return textInput(param.name, param.default, { placeholder: "/path/to/dir" });
        case "string_list": return textInput(param.name, param.default, { placeholder: "val1,val2,..." });
    }
}
```

Group fields by `parameter_groups[].parameters` ordering. Parameters not in
any group go in an "Advanced" section.

## Files to modify

| File | Change |
|------|--------|
| `mcp_server/app.py` | Import and mount Launchpad routes in `build_http_app()` |
| `mcp_server/launchpad/__init__.py` | New (empty) |
| `mcp_server/launchpad/api.py` | New — 4 Starlette endpoint functions |
| `mcp_server/launchpad/static/index.html` | New — SPA HTML shell |
| `mcp_server/launchpad/static/app.js` | New — form generation + API calls |
| `mcp_server/launchpad/static/style.css` | New — styling |
| `mcp_server/tests/test_launchpad.py` | New — API endpoint tests |
| `mcp_server/README.md` | Add Launchpad section |

## Existing code to reuse (do not duplicate)

All workflow logic already exists. The Launchpad API layer should be a thin
HTTP wrapper around these functions from `mcp_server/tools/workflows.py`:

- `list_workflows(authenticated=False)` → returns `ListWorkflowsResponse`
- `get_workflow_parameter_schema(workflow_id, authenticated=False)` → returns `WorkflowParameterSchemaResponse`
- `validate_workflow_parameters(workflow_id, params, check_paths=False)` → returns `ValidateWorkflowResponse`
- `render_workflow_command(workflow_id, params)` → returns `RenderWorkflowResponse`

Response models are in `mcp_server/schemas/responses.py`. All have
`.model_dump()` for JSON serialization.

## Reference: existing workflow schemas

These are the workflows the form will expose (public visibility only):

- `ucsf_star_suite_production` — 23 parameters, 5 stages, full production
  workflow with STAR + downstream + CellBender

Private workflows (slam, adapter_clip, etc.) are hidden from unauthenticated
callers by `_is_visible()` and should not appear in Launchpad.

## Schema structure reference

See `mcp_server/schemas/workflow.py` for models. Key fields the form needs:

**WorkflowParameterDef:** `name`, `cli_flag`, `type`, `required`, `default`,
`description`, `choices`, `category`, `skip_when_default`, `env_var`,
`is_output_root`

**WorkflowParameterGroup:** `name`, `title`, `parameters` (ordered list of
param names)

**WorkflowConstraint:** `kind` (mutual_exclusion, dependency, group_required,
positive, non_negative), `params`, `message`

**WorkflowRenderingRule:** `bool_style`, `flag_order`, `omit_absent_optionals`

## Constraints

- No npm / node / build step. Static files only.
- No new Python dependencies beyond what is already installed (starlette,
  jinja2, pydantic, uvicorn are all available).
- The SPA must not conflict with MCP endpoints (`POST /`, `GET /sse`,
  `POST /messages`).
- Do not add authentication to Launchpad endpoints — they expose only public
  workflow metadata and command rendering (no execution, no host paths).
- Do not duplicate workflow logic — call existing functions.
- Keep the frontend simple and functional. No framework overhead.

## Testing

- API tests: `mcp_server/tests/test_launchpad.py` using Starlette's
  `TestClient`. Test all 4 endpoints, including validation errors and unknown
  workflow handling.
- Manual: start server, visit `http://localhost:8765/launchpad/`, select a
  workflow, fill the form, verify the generated command.

## Out of scope

- Remote execution (Launchpad generates commands, not runs them)
- Authentication for Launchpad endpoints
- Private workflow access from Launchpad
- File browser / path autocomplete
- Saving/loading parameter sets (future)
