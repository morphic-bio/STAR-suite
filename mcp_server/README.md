# STAR-suite MCP Server

A Model Context Protocol (MCP) server for agents working with the STAR-suite repository.
Provides discovery, preflight validation, and script execution capabilities.

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
│   ├── reload.py       # reload_config
│   └── utils.py        # Path validation, run ID generation
├── schemas/
│   ├── config.py       # Config models
│   ├── run_config.py   # Run request models
│   └── responses.py    # Response models
└── tests/
    ├── conftest.py
    ├── test_config.py
    ├── test_discovery.py
    ├── test_preflight.py
    ├── test_executor.py
    ├── test_auth.py
    └── test_utils.py
```

## Implementation Status

- [x] **Phase 1**: Discovery tools + Auth + Config reload
- [x] **Phase 2**: Preflight validation
- [x] **Phase 3**: Script execution (queue, timeout, logs, true concurrency)
- [x] **Phase 4**: Stabilization + Docs (107 tests, AGENTS.md updated)
- [x] **Phase 5**: Containerization (Dockerfile, docker-compose.yml)

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
