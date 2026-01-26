# MCP Server Plan (Internal Agent Development)

Goal: provide a lightweight, stable MCP service for agents to discover data,
run preflight checks, and execute standard test/run scripts. HPC/Temporal
executor integration is deferred to a later stage.

---

## Design Decisions (Finalized)

### Transport & Clients
- **Primary**: HTTP/SSE for networked access (Cursor, Codex, other MCP clients)
- **Secondary**: stdio for local dev/testing only
- **Auth**: Bearer token (configurable in `config.yaml`)

### Security Model
- **Trusted roots** (all paths validated against these):
  - Repo: `/mnt/pikachu/STAR-suite`
  - Data: `/storage/*` (explicit allowlist per dataset)
  - Temp: `/tmp/*`
  - Tools: `/usr/bin` (read-only, for binary discovery)
- **Script execution**: Hard allowlist only—no arbitrary commands
- **Path validation**: Reject all paths outside trusted roots

### Execution Constraints
- **Timeout**: 2 hours default, per-script override in config
- **Concurrency**: 1 job at a time (queue subsequent requests)
- **Output**: Stream logs during execution + final artifact summary

### State & Logging
- **Run IDs**: `YYYYMMDD_HHMMSS_<rand4>` (e.g., `20260126_143052_a7f3`)
- **Logs**: Filesystem only (`plans/artifacts/mcp_runs_YYYYMMDD/`)
- **SQLite**: Deferred (add if search/history needed)
- **Config reload**: Hot reload supported (manual tool + optional file watcher)

### Dependencies
- **Python**: 3.10+
- **Packages**: `fastmcp`, `pydantic`, `pyyaml`, `aiofiles`

---

## File Structure

```
mcp_server/
├── __init__.py
├── app.py                 # FastMCP server entry point
├── config.py              # Config loading & validation
├── config.yaml            # Server configuration
├── tools/
│   ├── __init__.py
│   ├── discovery.py       # list_datasets, list_test_suites, find_docs, find_tests
│   ├── preflight.py       # preflight validation
│   ├── executor.py        # run_script, collect_outputs
│   ├── reload.py          # reload_config
│   └── utils.py           # path validation, run ID generation
├── schemas/
│   ├── __init__.py
│   ├── config.py          # Pydantic models for config
│   ├── run_config.py      # Pydantic models for run requests
│   └── responses.py       # Pydantic models for tool responses
├── tests/
│   ├── __init__.py
│   ├── test_preflight.py
│   ├── test_discovery.py
│   ├── test_executor.py
│   └── conftest.py
├── README.md
└── requirements.txt
```

---

## Configuration Schema (`config.yaml`)

```yaml
server:
  host: "0.0.0.0"
  port: 8765
  auth_token: "${MCP_AUTH_TOKEN}"  # env var expansion
  transport: "http"  # or "stdio"

paths:
  repo_root: "/mnt/pikachu/STAR-suite"
  artifact_log_root: "/mnt/pikachu/STAR-suite/plans/artifacts"
  temp_root: "/tmp"

trusted_roots:
  - "/mnt/pikachu/STAR-suite"
  - "/storage"
  - "/tmp"
  - "/usr/bin"  # read-only for tool discovery

datasets:
  - id: "a375"
    path: "/storage/A375"
    description: "A375 cell line data (GEX + features)"
    metadata_file: "dataset.yaml"  # optional, relative to path
  # Add more as discovered

execution:
  default_timeout_seconds: 7200  # 2 hours
  max_concurrent_jobs: 1
  queue_max_size: 10

# Disk space thresholds (per-suite configurable)
disk_space:
  default_min_gb: 20           # smoke tests
  full_depth_min_gb: 100       # full A375/GeneFull runs
  # Per-suite overrides below in test_suites

# Allowlisted scripts with metadata
scripts:
  # Core / Solo
  - name: "cbub_regression"
    path: "tests/run_cbub_regression_test.sh"
    module: "core"
    timeout_seconds: 3600
    description: "CB/UB tag regression test"
    
  - name: "sorted_bam_cbub_nonflex"
    path: "tests/run_sorted_bam_cbub_nonflex_test.sh"
    module: "core"
    description: "Sorted BAM CB/UB injection (non-Flex)"
    
  - name: "unsorted_cbub_smoke"
    path: "tests/run_unsorted_cbub_smoke_test.sh"
    module: "core"
    description: "Unsorted BAM CB/UB smoke test"

  # Flex
  - name: "flex_smoke"
    path: "tests/flex_smoke/run_flex_smoke.sh"
    module: "flex"
    description: "Flex basic smoke test"
    fixtures:
      - "/storage/flex_fixtures"  # example
    
  - name: "flex_cbub_validation"
    path: "tests/run_flex_cbub_validation_test.sh"
    module: "flex"
    description: "Flex CB/UB validation"
    
  - name: "flex_multisample"
    path: "tests/run_flex_multisample_test.sh"
    module: "flex"
    description: "Flex multi-sample test"

  # CR-compat / CRISPR
  - name: "crispr_calling"
    path: "tests/test_cr_compat_crispr_calling.sh"
    module: "cr-compat"
    description: "CR-compat CRISPR feature calling"
    
  - name: "a375_gex_features_cr_parity"
    path: "tests/run_a375_gex_features_cr_parity_genefull.sh"
    module: "cr-compat"
    description: "A375 GEX+features CR parity test (GeneFull)"
    fixtures:
      - "/storage/A375"

  # Feature tools
  - name: "assignbarcodes_regression"
    path: "core/features/process_features/tests/test_assignbarcodes_regression.sh"
    module: "feature-tools"
    description: "assignBarcodes regression test"

  # SLAM (optional - fixtures not yet confirmed; paths are placeholders)
  - name: "slam_smoke"
    path: "slam/tests/run_slam_smoke.sh"
    module: "slam"
    description: "SLAM-seq basic smoke test"
    runnable: false  # fixtures not confirmed
    fixtures:
      - "/storage/slam_fixtures"  # placeholder
    
  - name: "slam_qc"
    path: "slam/tests/run_slam_qc_test.sh"
    module: "slam"
    description: "SLAM-seq QC validation"
    runnable: false  # fixtures not confirmed
    fixtures:
      - "/storage/slam_fixtures"  # placeholder

  # Build commands (special category)
  - name: "build_core"
    path: "make"
    args: ["core"]
    module: "build"
    working_dir: "."
    description: "Build STAR core binary"
    
  - name: "build_flex"
    path: "make"
    args: ["flex"]
    module: "build"
    description: "Build Flex tools"
    
  - name: "build_slam"
    path: "make"
    args: ["slam"]
    module: "build"
    description: "Build SLAM tools"

# Test suites (grouped by module, with disk space requirements)
test_suites:
  - module: "core"
    description: "Core STAR / Solo tests (non-Flex)"
    scripts: ["cbub_regression", "sorted_bam_cbub_nonflex", "unsorted_cbub_smoke"]
    disk_space_min_gb: 20  # smoke tests
    
  - module: "flex"
    description: "STAR-Flex tests"
    scripts: ["flex_smoke", "flex_cbub_validation", "flex_multisample"]
    disk_space_min_gb: 20  # smoke tests
    
  - module: "cr-compat"
    description: "CellRanger compatibility tests"
    scripts: ["crispr_calling", "a375_gex_features_cr_parity"]
    disk_space_min_gb: 100  # full A375/GeneFull runs
    
  - module: "feature-tools"
    description: "Feature barcode tools tests"
    scripts: ["assignbarcodes_regression"]
    disk_space_min_gb: 20  # smoke tests
    
  - module: "slam"
    description: "SLAM-seq tests (fixtures not yet confirmed)"
    scripts: ["slam_smoke", "slam_qc"]
    disk_space_min_gb: 20  # TBD once fixtures confirmed
    runnable: false  # promote to true once fixtures available

# Required binaries for preflight checks
required_binaries:
  - name: "STAR"
    paths:
      - "core/legacy/source/STAR"
      - "/usr/local/bin/STAR"
  - name: "samtools"
    paths:
      - "/usr/bin/samtools"
      - "/usr/local/bin/samtools"
```

---

## MCP Tool Specifications

### 1. `list_datasets`

**Purpose**: Return known dataset roots with metadata.

**Parameters**: None (or optional `module` filter)

**Response**:
```json
{
  "datasets": [
    {
      "id": "a375",
      "path": "/storage/A375",
      "description": "A375 cell line data (GEX + features)",
      "size_bytes": 15234567890,
      "mtime": "2026-01-15T10:30:00Z",
      "metadata": { ... }  // from dataset.yaml if present
    }
  ]
}
```

---

### 2. `list_test_suites`

**Purpose**: Return module-aligned test suites with fixture availability.

**Parameters**:
- `module` (optional): Filter by module name

**Response**:
```json
{
  "suites": [
    {
      "module": "core",
      "description": "Core STAR / Solo tests (non-Flex)",
      "scripts": [
        {
          "name": "cbub_regression",
          "path": "tests/run_cbub_regression_test.sh",
          "description": "CB/UB tag regression test",
          "runnable": true,
          "missing_fixtures": []
        },
        {
          "name": "flex_smoke",
          "path": "tests/run_flex_smoke.sh",
          "description": "Flex basic smoke test",
          "runnable": false,
          "missing_fixtures": ["/storage/flex_fixtures"]
        }
      ]
    }
  ]
}
```

---

### 3. `find_docs`

**Purpose**: Search documentation by topic/keyword.

**Parameters**:
- `topic` (required): Search term

**Response**:
```json
{
  "matches": [
    {
      "path": "docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md",
      "title": "CRISPR Feature Calling Implementation",
      "snippet": "CR-compat mode runs GMM-based CRISPR calling..."
    }
  ]
}
```

---

### 4. `find_tests`

**Purpose**: Search test scripts by tag/keyword.

**Parameters**:
- `tag` (required): Search term (e.g., "cbub", "flex", "crispr")

**Response**:
```json
{
  "matches": [
    {
      "name": "cbub_regression",
      "path": "tests/run_cbub_regression_test.sh",
      "module": "core",
      "description": "CB/UB tag regression test"
    }
  ]
}
```

---

### 5. `preflight`

**Purpose**: Validate a run configuration before execution.

**Parameters** (`run_config`):
```json
{
  "module": "core",
  "script": "cbub_regression",
  "args": ["--verbose"],
  "dataset_id": "a375",           // optional
  "out_dir": "/tmp/test_output",  // optional, defaults to temp
  "env_overrides": {              // optional
    "STAR_PATH": "/custom/path/STAR"
  }
}
```

**Checks performed**:
1. Script is in allowlist
2. All paths exist and are readable/writable as needed
3. Dataset exists (if specified)
4. Output directory is writable (or can be created)
5. Required binaries are present (STAR, samtools)
6. Sufficient disk space (per-suite: 20GB smoke, 100GB full-depth)
7. Required fixtures exist for the script

**Response**:
```json
{
  "valid": true,
  "checks": [
    {"name": "script_allowed", "passed": true},
    {"name": "paths_valid", "passed": true},
    {"name": "binaries_present", "passed": true, "details": {"STAR": "/mnt/pikachu/STAR-suite/core/legacy/source/STAR"}},
    {"name": "disk_space", "passed": true, "details": {"available_gb": 150}},
    {"name": "fixtures_present", "passed": true}
  ],
  "warnings": [],
  "errors": []
}
```

---

### 6. `run_script`

**Purpose**: Execute an allowlisted script.

**Parameters**:
```json
{
  "script": "cbub_regression",
  "args": ["--verbose"],
  "dataset_id": "a375",
  "out_dir": "/tmp/test_output",
  "env_overrides": {},
  "stream_logs": true
}
```

**Behavior**:
1. Generate run ID: `20260126_143052_a7f3`
2. Create log directory: `plans/artifacts/mcp_runs_20260126/20260126_143052_a7f3/`
3. Run preflight (fail fast if invalid)
4. Execute script with timeout
5. Stream stdout/stderr if `stream_logs: true`
6. On completion, write summary JSON

**Response** (immediate):
```json
{
  "run_id": "20260126_143052_a7f3",
  "status": "running",
  "log_path": "plans/artifacts/mcp_runs_20260126/20260126_143052_a7f3/",
  "stream_url": "/runs/20260126_143052_a7f3/stream"  // if HTTP
}
```

**Response** (on completion via SSE or poll):
```json
{
  "run_id": "20260126_143052_a7f3",
  "status": "completed",  // or "failed", "timeout"
  "exit_code": 0,
  "duration_seconds": 127,
  "log_path": "plans/artifacts/mcp_runs_20260126/20260126_143052_a7f3/",
  "log_tail": "... last 50 lines ...",
  "outputs": [
    "/tmp/test_output/results.tsv",
    "/tmp/test_output/summary.json"
  ]
}
```

---

### 7. `collect_outputs`

**Purpose**: Retrieve outputs and logs for a completed run.

**Parameters**:
- `run_id` (required): Run identifier

**Response**:
```json
{
  "run_id": "20260126_143052_a7f3",
  "status": "completed",
  "exit_code": 0,
  "duration_seconds": 127,
  "started_at": "2026-01-26T14:30:52Z",
  "completed_at": "2026-01-26T14:32:59Z",
  "script": "cbub_regression",
  "args": ["--verbose"],
  "log_files": {
    "stdout": "plans/artifacts/mcp_runs_20260126/20260126_143052_a7f3/stdout.log",
    "stderr": "plans/artifacts/mcp_runs_20260126/20260126_143052_a7f3/stderr.log",
    "combined": "plans/artifacts/mcp_runs_20260126/20260126_143052_a7f3/combined.log"
  },
  "outputs": [
    {
      "path": "/tmp/test_output/results.tsv",
      "size_bytes": 12345,
      "mtime": "2026-01-26T14:32:58Z"
    }
  ]
}
```

---

### 8. `reload_config`

**Purpose**: Reload `config.yaml` without restarting the server.

**Parameters**: None

**Response**:
```json
{
  "reloaded": true,
  "config_path": "mcp_server/config.yaml",
  "loaded_at": "2026-01-26T15:02:11Z"
}
```

---

## Implementation Phases

### Phase 1: Foundation (Transport + Auth + Discovery)

**Deliverables**:
- `app.py` with FastMCP HTTP/SSE server
- `config.py` with Pydantic config loading
- `config.yaml` with initial datasets/scripts
- Token-based auth middleware
- `list_datasets` tool
- `list_test_suites` tool (with fixture checks)
- `find_docs` tool
- `find_tests` tool
- `reload_config` tool (manual hot reload)

**Tests**:
- Config validation
- Auth rejection for bad tokens
- Discovery returns expected structure

**Exit criteria**: Agent can connect, authenticate, and discover available tests/datasets.

---

### Phase 2: Preflight Validation

**Deliverables**:
- `preflight.py` with all checks
- Path validation utilities
- Binary discovery logic
- Disk space checker
- Fixture presence checker

**Tests**:
- Preflight passes for valid configs
- Preflight fails with clear errors for:
  - Unknown script
  - Missing paths
  - Missing binaries
  - Insufficient disk space
  - Missing fixtures

**Exit criteria**: Agent can validate any run config before execution.

---

### Phase 3: Script Execution

**Deliverables**:
- `executor.py` with subprocess management
- Run ID generation
- Log directory creation
- Stdout/stderr capture
- Timeout handling
- Job queue (1 concurrent, queue up to 10)
- `run_script` tool
- `collect_outputs` tool
- SSE streaming for logs

**Tests**:
- Successful script execution
- Timeout handling
- Queue behavior
- Log file creation
- Output collection

**Exit criteria**: Agent can run allowlisted scripts and retrieve results.

---

### Phase 4: Stabilization + Docs

**Deliverables**:
- `mcp_server/README.md` (setup, auth, examples)
- Update `AGENTS.md` with MCP server reference
- Integration test suite
- Error message improvements
- Optional file watcher to auto-reload config on change (behind flag)

**Tests**:
- End-to-end: discover → preflight → run → collect

**Exit criteria**: Another developer/agent can use the server with docs alone.

---

### Phase 5: Containerization (Optional)

**Deliverables**:
- `Dockerfile`
- `docker-compose.yml`
- Bind mount configuration
- Environment-based config

**Tests**:
- Container builds
- Container runs with mounted volumes
- Auth works in container

---

### Phase 6 (Deferred): HPC/Temporal Integration

**Requirements** (for future):
- `submit_hpc_job(run_config)` → Temporal workflow
- `job_status(job_id)`
- `collect_outputs(job_id)` (from HPC output location)
- Audit trail for all submissions

---

## Run Log Structure

```
plans/artifacts/
└── mcp_runs_20260126/
    └── 20260126_143052_a7f3/
        ├── run.json          # Full run config + metadata
        ├── stdout.log
        ├── stderr.log
        ├── combined.log
        └── summary.json      # Exit code, duration, outputs
```

**`run.json`**:
```json
{
  "run_id": "20260126_143052_a7f3",
  "script": "cbub_regression",
  "args": ["--verbose"],
  "dataset_id": "a375",
  "out_dir": "/tmp/test_output",
  "env_overrides": {},
  "started_at": "2026-01-26T14:30:52Z",
  "preflight_result": { ... }
}
```

**`summary.json`**:
```json
{
  "run_id": "20260126_143052_a7f3",
  "status": "completed",
  "exit_code": 0,
  "started_at": "2026-01-26T14:30:52Z",
  "completed_at": "2026-01-26T14:32:59Z",
  "duration_seconds": 127,
  "outputs": ["/tmp/test_output/results.tsv"]
}
```

---

## Error Handling

All tools return structured errors:

```json
{
  "error": true,
  "code": "SCRIPT_NOT_ALLOWED",
  "message": "Script 'unknown_script' is not in the allowlist",
  "details": {
    "requested": "unknown_script",
    "allowed": ["cbub_regression", "flex_smoke", ...]
  }
}
```

Error codes:
- `AUTH_FAILED` - Invalid or missing token
- `SCRIPT_NOT_ALLOWED` - Script not in allowlist
- `PATH_OUTSIDE_ROOTS` - Path validation failed
- `PREFLIGHT_FAILED` - One or more preflight checks failed
- `BINARY_NOT_FOUND` - Required binary missing
- `FIXTURE_MISSING` - Required fixture not present
- `DISK_SPACE_LOW` - Insufficient disk space
- `TIMEOUT` - Script exceeded timeout
- `QUEUE_FULL` - Job queue at capacity
- `RUN_NOT_FOUND` - Unknown run ID

---

## Testing Strategy

### Unit Tests

- Config loading (valid/invalid YAML, missing fields, env expansion)
- Path validation (allowed vs disallowed paths)
- Dataset discovery (metadata file parsing, size/mtime collection)
- Test suite discovery (fixture detection, runnable flag)
- Preflight checks (missing binaries, missing fixtures, disk space)
- Executor timeout and queue handling (mocked subprocess)
- Hot reload (invalid config rejected, old config retained)

### Integration Tests

- Auth required for HTTP endpoints
- End-to-end: discover → preflight → run → collect
- SSE log streaming behavior
- Artifact log writing and output discovery

### Test Fixtures

- Use temp directories under `/tmp` for fake outputs/logs
- Mock dataset roots for discovery and preflight
- Keep tests hermetic (no `/storage` dependency)

### CI / Local

- `pytest -q` for unit tests
- Integration tests optional (flagged) if no datasets present

---

## Open Items (Track as Discovered)

- [ ] **Pending**: Dataset inventory from `/storage/` scan (user to provide)
- [ ] Confirm fixture paths for Flex tests
- [x] ~~Decide if SLAM tests should be included~~ → Added as optional/not-runnable
- [x] ~~Determine disk space threshold~~ → 20GB smoke, 100GB full-depth, per-suite
- [ ] Decide on log retention policy (auto-prune after N days?)
- [ ] Confirm SLAM fixture paths once available
