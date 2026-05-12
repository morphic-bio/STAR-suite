# STAR-suite MCP Server

A Model Context Protocol (MCP) server for agents working with the STAR-suite repository.
Provides discovery, preflight validation, and script execution capabilities.

**Terminology**

| Term | Meaning |
|------|---------|
| **STAR-suite** | This repository and its tooling. |
| **STAR Server** | The Python process (`python -m mcp_server.app`): serves MCP over HTTP, SSE, and the browser UI on one port. |
| **STAR MCP** | The agent-facing MCP tool surface (`POST /` streamable-HTTP, `GET /sse` + `POST /messages` SSE). |
| **Shared core** | `mcp_server/tools/workflows.py` and workflow YAML under `mcp_server/workflows/` — validation, rendering, schemas (used by both STAR MCP and STAR Launchpad). |
| **STAR Launchpad** | Static SPA at `/launchpad/` with JSON API under `/launchpad/api/`. **STAR workflows** tab: validates and renders recipe commands; **Load/Save parameters** uses a JSON file in the browser (client-side). **Run in shell** (loopback) starts the rendered argv on the **server host**. **Script Lane** tab: annotated Bash to Workflow IR, simple-local viability, and local execute via `bwb-nextflow-utils` (`POST /launchpad/api/script-lane/*`); separate from **Run in shell**. Discovery: set `BWB_NEXTFLOW_UTILS_ROOT` or place `bwb-nextflow-utils` beside this repo. The **Include test & other recipes** checkbox (off by default) limits the recipe list to **`star_*`** workflows; turn it on for UCSF production, SLAM PE smoke/production, E2E tests, and private entries on localhost. See `plans/star_launchpad_v1_runbook.md`. |

## Quick Start

### Installation

```bash
cd mcp_server
pip install -r requirements.txt
```

### Running the Server

**HTTP/SSE mode** (default, for networked access):

```bash
# Set auth token (optional but recommended)
export MCP_AUTH_TOKEN="your-secret-token"

# Run server
python -m mcp_server.app
```

Open `http://<host>:<port>/launchpad/` in a browser for **STAR Launchpad** (recipe builder). Remote browsers typically see **public** workflows only; on **loopback**, authenticated discovery can list **private** workflows too. The UI defaults to **`star_*`** recipes only; enable **Include test & other recipes** to show the full list, including the private SLAM PE recipes. MCP clients continue to use `POST /` (streamable-HTTP) or `GET /sse` + `POST /messages` (SSE).

#### Launchpad quick start / stop (recommended)

Start in the background (writes a pidfile + log under `plans/artifacts/`):

```bash
bash scripts/launchpad_server.sh up
```

Stop it:

```bash
bash scripts/launchpad_server.sh down
```

Show status / tail logs:

```bash
bash scripts/launchpad_server.sh status
bash scripts/launchpad_server.sh logs
```

**stdio mode** (for local development):

```bash
python -m mcp_server.app --transport stdio
```

### Configuration

Edit `config.yaml` to customize:

- **Server settings**: host, port, auth token
- **Trusted paths**: allowed filesystem roots
- **Datasets**: data directories to expose
- **Scripts**: allowlisted test/build scripts
- **Test suites**: grouped tests by module

Environment variables in config use `${VAR_NAME}` syntax.

## Available Tools

### Discovery (Phase 1)

| Tool | Description |
|------|-------------|
| `list_datasets(auth_token?)` | List configured datasets with metadata |
| `list_test_suites(module?, auth_token?)` | List test suites with fixture availability |
| `find_docs(topic, auth_token?)` | Search documentation by keyword |
| `find_tests(tag, auth_token?)` | Search test scripts by keyword |
| `reload_config(auth_token?)` | Hot-reload config.yaml without restart |

### Validation (Phase 2)

| Tool | Description |
|------|-------------|
| `preflight(script, args?, module?, dataset_id?, out_dir?, env_overrides?, auth_token?)` | Validate run configuration before execution |

**Preflight checks**:
- Script is in allowlist and marked runnable
- Script file exists and is within trusted roots
- Script working_dir is within trusted roots (if specified)
- Dataset exists and is within trusted roots (if specified)
- Output directory is writable (warns if deep chain needs creation)
- All paths are within trusted roots
- Required binaries (STAR, samtools) present
- Sufficient disk space (per-suite thresholds)
- Required fixtures present and within trusted roots

### Build Tools

| Tool | Description |
|------|-------------|
| `build_star(target?, clean?, force?, auth_token?)` | Build STAR-suite binaries |
| `check_build_status(target?, auth_token?)` | Check if rebuild is needed |
| `ensure_fresh_build(target?, auth_token?)` | Clean build (recommended before tests) |

**Targets**: `core` (STAR binary), `flex`, `slam`, `feature-tools`

**IMPORTANT**: Always use `ensure_fresh_build()` or `build_star(clean=True)` before running test suites to prevent stale binary issues that can cause segfaults.

### Workflow Tools

| Tool | Description |
|------|-------------|
| `list_workflows(auth_token?)` | List supported workflow templates |
| `describe_workflow(workflow_id, auth_token?)` | Full workflow metadata with stages and parameter groups |
| `get_workflow_scripts(workflow_id, auth_token?)` | Scripts composing a workflow with provenance (for script-backed encoders) |
| `get_workflow_parameter_schema(workflow_id, auth_token?)` | Machine-readable parameter schema with types, defaults, and constraints |
| `validate_workflow_parameters(workflow_id, params, check_paths?, auth_token?)` | Validate structured params without executing |
| `render_workflow_command(workflow_id, params, auth_token?)` | Render params into argv array and shell preview |

**Workflow vs Script tools**: Workflows provide structured parameter contracts for agent consumption. Scripts provide raw allowlisted execution. They are complementary:

Current local/private SLAM workflows:

| Workflow ID | Purpose | Default render |
|-------------|---------|----------------|
| `slam_pe_100k_smoke` | Reproduce the 100K R1-only vs R1/R2 smoke with fixed noSU-derived trims. | Runs the two-sample ARID1A smoke unless `dry_run` is set. |
| `slam_pe_production` | Full PE panel runner with TranscriptVB tximport counts, GrandSLAM, cB, Y/noY outputs, and Globus cleanup. | `--pilot --dry-run` for safety. |
| `slam_deseq2_container` | Build and verify the pinned R/Bioconductor/DESeq2/tximport container. | Docker build with pinned versions via environment overrides. |

Do not launch the SLAM production or DESeq2 container workflows while another
STAR-SLAM production/benchmark run is active on the same host.

```python
# Agent workflow: discover -> inspect scripts -> validate -> render -> (optionally execute)
wf = client.call_tool("list_workflows", {})
schema = client.call_tool("get_workflow_parameter_schema", {
    "workflow_id": "ucsf_star_suite_production"
})

# Get scripts composing the workflow (entry + helpers + provenance)
scripts = client.call_tool("get_workflow_scripts", {
    "workflow_id": "ucsf_star_suite_production"
})
# scripts["scripts"] lists each script with role, path, description, language, exists
# Public callers: absolute_path is null; provenance has only workflow_schema
# Authenticated callers: absolute_path populated; provenance adds repo_root,
#   git_commit, git_remote

# Validate with path checks (default) or skip path checks for dry planning
val = client.call_tool("validate_workflow_parameters", {
    "workflow_id": "ucsf_star_suite_production",
    "params": {
        "all_samples": True,
        "threads": 16,
        "dry_run": True,
        "dataset_root": "/mnt/pikachu/ucsf-perturb-seq-corrected",
        "genome_dir": "/storage/autoindex_110_44/bulk_index",
    },
    "check_paths": True,  # set False to skip file/dir existence checks
})

# Render normalized params into argv + shell preview
cmd = client.call_tool("render_workflow_command", {
    "workflow_id": "ucsf_star_suite_production",
    "params": val["normalized_params"]
})
# cmd["argv"] is the canonical argv array
# cmd["shell_preview"] is a shell-safe joined string
# cmd["env_overrides"] contains e.g. {"DOWNSAMPLE_SEED": "1"}
```

### Execution (Phase 3)

| Tool | Description |
|------|-------------|
| `run_script(script, args?, module?, dataset_id?, out_dir?, env_overrides?, auth_token?)` | Execute allowlisted script (runs preflight first) |
| `collect_outputs(run_id, auth_token?)` | Retrieve outputs and logs for a run |
| `get_run_status(run_id, auth_token?)` | Get run status with log tail |

**Execution features**:
- Preflight validation before execution
- Job queue (1 concurrent, up to 10 queued)
- Configurable timeout (per-script or default 2 hours)
- Process group cleanup on timeout
- Log files: `stdout.log`, `stderr.log`, `combined.log`
- Run metadata: `run.json`, `summary.json`

## Authentication

Authentication is enforced per-tool via the `auth_token` parameter.

**Discovery tools are public by default** (`public_discovery: true` in config):
- `list_datasets`, `list_test_suites`, `find_docs`, `find_tests` work without auth
- Set `public_discovery: false` to require auth for all tools

**Sensitive tools always require auth** (when `auth_token` is configured):
- `reload_config`, `preflight`, `run_script`, `collect_outputs`

**To enable auth**:
1. Set `MCP_AUTH_TOKEN` environment variable
2. For sensitive tools, pass the token via the `auth_token` parameter

**Example**:
```python
# Discovery tools work without auth (public_discovery: true)
response = client.call_tool("list_datasets", {})

# Sensitive tools require auth
response = client.call_tool("preflight", {
    "script": "cbub_regression",
    "auth_token": "your-secret-token"
})
```

To disable all auth, leave `auth_token` empty in `config.yaml`.

## Example Usage

```python
# Using an MCP client
response = client.call_tool("list_test_suites", {
    "module": "core",
    "auth_token": "your-token"
})
print(response["suites"])

response = client.call_tool("find_docs", {
    "topic": "cbub",
    "auth_token": "your-token"
})
for doc in response["matches"]:
    print(f"{doc['path']}: {doc['title']}")
```

## Build Workflow (Important!)

**ALWAYS ensure a clean build before running tests** to prevent stale binary issues:

```python
# Recommended workflow for test execution
# Step 1: Ensure fresh build
response = client.call_tool("ensure_fresh_build", {
    "target": "core",
    "auth_token": "your-token"
})
if not response.get("success"):
    print(f"Build failed: {response.get('error')}")
    sys.exit(1)

# Step 2: Run preflight
response = client.call_tool("preflight", {
    "script": "cbub_regression",
    "auth_token": "your-token"
})

# Step 3: Run test
response = client.call_tool("run_script", {
    "script": "cbub_regression",
    "auth_token": "your-token"
})
```

**Why this matters**: Stale object files from incremental builds can cause hard-to-debug segfaults that only appear with `-O3` optimization. The `ensure_fresh_build` tool tracks source file changes and always runs `make clean` before building.

**For incremental development** (not tests):
```python
# Quick incremental build (only rebuilds if sources changed)
response = client.call_tool("build_star", {
    "target": "core",
    "auth_token": "your-token"
})
```

## Endpoints and Remote Access

The server now exposes **both** MCP HTTP transports on the same port:

- **Streamable-HTTP** (Codex CLI / agent runtime): `POST /`
- **SSE** (Cursor MCP client): `GET /sse` + `POST /messages`

This is required because Codex and Cursor use different MCP transports and
endpoints.

### Remote Servers (SSH)

If the MCP server runs on a remote host, Codex agents cannot see your local SSH
tunnel (e.g., `localhost:8765`). Use a **public tunnel** so agents can reach the
server directly.

Recommended (no firewall changes):

```bash
cloudflared tunnel --url http://localhost:8765
```

This prints a `https://*.trycloudflare.com` URL that can be used by agents.

### Codex MCP config (streamable-http)

```toml
[mcp_servers.star-suite]
type = "streamable-http"
url = "https://your-public-url"
```

### Cursor MCP config (SSE)

```toml
[mcp_servers.star-suite]
type = "sse"
url = "https://your-public-url/sse"
```

### Sanity check (streamable-http)

```bash
curl -i -X POST https://your-public-url/ \
  -H 'content-type: application/json' \
  -d '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"clientInfo":{"name":"probe","version":"0.0"},"protocolVersion":"2024-11-05","capabilities":{}}}'
```

## Testing

```bash
pytest mcp_server/tests/ -v
```

## Project Structure

```
mcp_server/
├── app.py              # FastMCP server entry point
├── config.py           # Config loading & management
├── config.yaml         # Server configuration
├── tools/
│   ├── discovery.py    # list_datasets, list_test_suites, find_docs, find_tests
│   ├── preflight.py    # preflight validation
│   ├── executor.py     # run_script, collect_outputs
│   ├── workflows.py    # workflow discovery, validation, rendering
│   ├── reload.py       # reload_config
│   └── utils.py        # Path validation, run ID generation
├── schemas/
│   ├── config.py       # Config models (incl. WorkflowConfig)
│   ├── workflow.py     # Workflow schema models (WorkflowSchema, params, constraints)
│   ├── run_config.py   # Run request models
│   └── responses.py    # Response models
├── workflows/
│   └── ucsf_star_suite_production.yaml  # UCSF workflow parameter + script schema
└── tests/
    ├── conftest.py
    ├── test_config.py
    ├── test_discovery.py
    ├── test_preflight.py
    ├── test_executor.py
    ├── test_auth.py
    ├── test_utils.py
    ├── test_workflow_config.py
    ├── test_workflow_discovery.py
    ├── test_workflow_validation.py
    ├── test_workflow_render.py
    ├── test_workflow_integration.py
    └── test_ucsf_workflow_e2e.py
```

## Implementation Status

- [x] **Phase 1**: Discovery tools + Auth + Config reload
- [x] **Phase 2**: Preflight validation
- [x] **Phase 3**: Script execution (queue, timeout, logs, true concurrency)
- [x] **Phase 4**: Stabilization + Docs (107 tests, AGENTS.md updated)
- [x] **Phase 5**: Containerization (Dockerfile, docker-compose.yml)
- [x] **Phase 6**: Workflow parameter service (structured schemas, validation, rendering)

## Docker Deployment

### Quick Start

```bash
# Build and run with docker-compose
cd mcp_server
export MCP_AUTH_TOKEN="your-secret-token"
docker-compose up -d

# Or build manually
docker build -t star-mcp-server .
docker run -p 8765:8765 \
  -e MCP_AUTH_TOKEN="your-token" \
  -v /mnt/pikachu/STAR-suite:/repo:ro \
  -v /storage:/storage:ro \
  star-mcp-server
```

### Volume Mounts

| Container Path | Purpose | Mode |
|----------------|---------|------|
| `/repo` | STAR-suite repository | read-only |
| `/storage` | Data storage | read-only |
| `/tmp/mcp_runs` | Temporary outputs | read-write |
| `/app/artifacts` | Run logs | read-write |

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `MCP_AUTH_TOKEN` | (empty) | Authentication token |
| `MCP_HOST` | `0.0.0.0` | Server bind address |
| `MCP_PORT` | `8765` | Server port |
| `MCP_TRANSPORT` | `http` | Transport mode |
| `STAR_SUITE_PATH` | `..` | Host path to STAR-suite |
| `STORAGE_PATH` | `/storage` | Host path to data storage |

### Custom Configuration

Mount your own config file:

```bash
docker run -p 8765:8765 \
  -v ./my-config.yaml:/app/config.yaml:ro \
  star-mcp-server
```

## Security

All paths are validated against trusted roots configured in `config.yaml`:
- Script paths, working directories, output directories, datasets, and fixtures
- Binaries found via PATH are also validated against trusted roots
- Deep directory chains (>3 levels) generate warnings in preflight

## Concurrency

The executor uses a `ThreadPoolExecutor` for true parallel job execution:
- `max_concurrent_jobs` controls the thread pool size (default: 1)
- Jobs are dispatched from queue to pool as slots become available
- Settings can be updated via `reload_config()` (new pool created for changed worker count)
