"""FastMCP server entry point for STAR-suite MCP server."""

import argparse
import sys
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Optional

from fastmcp import FastMCP
from starlette.applications import Starlette
from starlette.middleware import Middleware
import uvicorn

from .config import get_config, load_config
from .schemas.responses import (
    CollectOutputsResponse,
    ErrorResponse,
    FindDocsResponse,
    FindTestsResponse,
    ListDatasetsResponse,
    ListTestSuitesResponse,
    PreflightResponse,
    ReloadConfigResponse,
    RunScriptResponse,
)
from .schemas.run_config import RunConfig
from .tools.discovery import (
    find_docs as _find_docs,
    find_tests as _find_tests,
    list_datasets as _list_datasets,
    list_test_suites as _list_test_suites,
)
from .tools.executor import (
    collect_outputs as _collect_outputs,
    get_log_tail as _get_log_tail,
    run_script as _run_script,
)
from .tools.preflight import run_preflight as _run_preflight
from .tools.reload import reload_config as _reload_config
from .tools.build import (
    build_target as _build_target,
    needs_rebuild as _needs_rebuild,
    load_build_state as _load_build_state,
)
from .tools.workflows import (
    list_workflows as _list_workflows,
    describe_workflow as _describe_workflow,
    get_workflow_scripts as _get_workflow_scripts,
    get_workflow_parameter_schema as _get_workflow_parameter_schema,
    validate_workflow_parameters as _validate_workflow_parameters,
    render_workflow_command as _render_workflow_command,
)
from .tools.scaffold import (
    scaffold_workflow_schema as _scaffold_workflow_schema,
    validate_draft_workflow_schema as _validate_draft_workflow_schema,
)
from .launchpad.api import get_launchpad_routes


# Create the MCP server
mcp = FastMCP("star-suite")

# Set to True when running in stdio transport (local to the host).
# Stdio callers are implicitly trusted — they are pikachu-local processes.
_trust_local: bool = False


class AcceptHeaderMiddleware:
    """Ensure clients accept both JSON and event-stream responses."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope.get("type") == "http":
            headers = dict(scope.get("headers") or [])
            accept = headers.get(b"accept", b"")
            if b"application/json" not in accept or b"text/event-stream" not in accept:
                parts = [p.strip() for p in accept.split(b",") if p.strip()]
                if b"application/json" not in accept:
                    parts.append(b"application/json")
                if b"text/event-stream" not in accept:
                    parts.append(b"text/event-stream")
                headers[b"accept"] = b", ".join(parts)
                scope = {**scope, "headers": list(headers.items())}
        await self.app(scope, receive, send)


def build_http_app() -> Starlette:
    """Build an HTTP app that supports both SSE and streamable-http clients."""
    sse_app = mcp.http_app(transport="sse")
    stream_app = mcp.http_app(transport="streamable-http", path="/")

    @asynccontextmanager
    async def lifespan(app: Starlette):
        async with sse_app.lifespan(sse_app), stream_app.lifespan(stream_app):
            yield

    middleware = [Middleware(AcceptHeaderMiddleware)]
    middleware += stream_app.user_middleware + sse_app.user_middleware

    # Launchpad API + static SPA under /launchpad/ (must precede generic /launchpad static catch-all)
    launchpad = get_launchpad_routes()

    return Starlette(
        routes=launchpad + stream_app.routes + sse_app.routes,
        middleware=middleware,
        lifespan=lifespan,
    )


def check_auth(token: Optional[str], is_discovery: bool = False) -> Optional[ErrorResponse]:
    """Check if the provided token is valid.

    Args:
        token: Bearer token from request.
        is_discovery: If True, this is a discovery tool that may be public.

    Returns:
        ErrorResponse if auth fails, None if auth passes.
    """
    # Stdio transport is local to the host — implicitly trusted.
    if _trust_local:
        return None

    config = get_config()
    expected_token = config.server.auth_token

    # If no token configured, skip auth
    if not expected_token:
        return None

    # Discovery tools can be public (no auth required)
    if is_discovery and config.server.public_discovery:
        return None

    if not token:
        return ErrorResponse(
            code="AUTH_FAILED",
            message="Authentication required. Provide auth_token parameter.",
        )

    if token != expected_token:
        return ErrorResponse(
            code="AUTH_FAILED",
            message="Invalid authentication token.",
        )

    return None


def is_authenticated(token: Optional[str]) -> bool:
    """Return True if *token* is valid (or if auth is not configured)."""
    if _trust_local:
        return True
    config = get_config()
    expected_token = config.server.auth_token
    if not expected_token:
        return True
    return bool(token and token == expected_token)


def get_workflow_visibility(workflow_id: str) -> str:
    """Return the visibility setting for a workflow ('public' or 'private')."""
    config = get_config()
    wf_cfg = config.get_workflow(workflow_id)
    if wf_cfg is None:
        return "public"
    return wf_cfg.visibility


# --- MCP Tools ---


@mcp.tool()
def list_datasets(auth_token: Optional[str] = None) -> dict:
    """List all configured datasets with metadata.

    Args:
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns information about available datasets including path,
    description, size, and modification time.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _list_datasets()
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="INTERNAL_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def list_test_suites(
    module: Optional[str] = None,
    auth_token: Optional[str] = None,
) -> dict:
    """List test suites with fixture availability.

    Args:
        module: Optional filter by module name (e.g., "core", "flex").
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns test suites grouped by module, with information about
    which scripts are runnable and any missing fixtures.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _list_test_suites(module=module)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="INTERNAL_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def find_docs(topic: str, auth_token: Optional[str] = None) -> dict:
    """Search documentation by topic or keyword.

    Args:
        topic: Search term to look for in documentation.
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns matching documents with title and snippet.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _find_docs(topic=topic)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="INTERNAL_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def find_tests(tag: str, auth_token: Optional[str] = None) -> dict:
    """Search test scripts by tag or keyword.

    Args:
        tag: Search term (e.g., "cbub", "flex", "crispr").
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns matching test scripts with module and description.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _find_tests(tag=tag)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="INTERNAL_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def reload_config(auth_token: Optional[str] = None) -> dict:
    """Reload configuration without restarting the server.

    Args:
        auth_token: Authentication token (required if server has auth configured).

    Use this after modifying config.yaml to apply changes.
    If the new config is invalid, the old config is retained.
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _reload_config()
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="INTERNAL_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def preflight(
    script: str,
    args: Optional[list[str]] = None,
    module: Optional[str] = None,
    dataset_id: Optional[str] = None,
    out_dir: Optional[str] = None,
    env_overrides: Optional[dict[str, str]] = None,
    auth_token: Optional[str] = None,
) -> dict:
    """Validate a run configuration before execution.

    Performs preflight checks to ensure a script can be run successfully:
    - Script is in the allowlist and marked as runnable
    - Script file exists
    - Dataset exists (if specified)
    - Output directory is writable
    - Required binaries (STAR, samtools) are present
    - Sufficient disk space available
    - Required fixtures are present

    Args:
        script: Name of the script to run (must be in allowlist).
        args: Optional additional arguments to pass to the script.
        module: Optional module context for validation.
        dataset_id: Optional dataset to use.
        out_dir: Optional output directory (defaults to temp).
        env_overrides: Optional environment variable overrides.
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        PreflightResponse with validation results and any errors.
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        run_config = RunConfig(
            script=script,
            args=args,
            module=module,
            dataset_id=dataset_id,
            out_dir=out_dir,
            env_overrides=env_overrides,
        )
        result = _run_preflight(run_config)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="PREFLIGHT_FAILED",
            message=str(e),
        ).model_dump()


@mcp.tool()
def run_script(
    script: str,
    args: Optional[list[str]] = None,
    module: Optional[str] = None,
    dataset_id: Optional[str] = None,
    out_dir: Optional[str] = None,
    env_overrides: Optional[dict[str, str]] = None,
    auth_token: Optional[str] = None,
) -> dict:
    """Execute an allowlisted script.

    Runs preflight validation first, then queues the script for execution.
    Returns immediately with a run_id that can be used to check status.

    Args:
        script: Name of the script to run (must be in allowlist).
        args: Optional additional arguments to pass to the script.
        module: Optional module context for validation.
        dataset_id: Optional dataset to use.
        out_dir: Optional output directory (defaults to temp).
        env_overrides: Optional environment variable overrides.
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        RunScriptResponse with run_id and status, or ErrorResponse on failure.
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        run_config = RunConfig(
            script=script,
            args=args,
            module=module,
            dataset_id=dataset_id,
            out_dir=out_dir,
            env_overrides=env_overrides,
        )
        result = _run_script(run_config)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="EXECUTION_FAILED",
            message=str(e),
        ).model_dump()


@mcp.tool()
def collect_outputs(run_id: str, auth_token: Optional[str] = None) -> dict:
    """Retrieve outputs and logs for a completed run.

    Args:
        run_id: Run identifier returned by run_script.
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        CollectOutputsResponse with run details, logs, and outputs.
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _collect_outputs(run_id)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="COLLECTION_FAILED",
            message=str(e),
        ).model_dump()


@mcp.tool()
def get_run_status(run_id: str, auth_token: Optional[str] = None) -> dict:
    """Get the current status of a run.

    Args:
        run_id: Run identifier returned by run_script.
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        Status information including run state and log tail.
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _collect_outputs(run_id)
        response = result.model_dump()

        # Add log tail for convenience
        log_tail = _get_log_tail(run_id, lines=50)
        if log_tail:
            response["log_tail"] = log_tail

        return response
    except Exception as e:
        return ErrorResponse(
            code="STATUS_FAILED",
            message=str(e),
        ).model_dump()


# --- Build Tools ---


@mcp.tool()
def build_star(
    target: str = "core",
    clean: bool = False,
    force: bool = False,
    auth_token: Optional[str] = None,
) -> dict:
    """Build STAR-suite binaries with proper state tracking.

    IMPORTANT: Always use clean=True before running test suites to prevent
    stale binary issues that can cause segfaults.

    Args:
        target: Build target - 'core' (STAR binary), 'flex', 'slam', or 'feature-tools'.
        clean: If True, run 'make clean' before building (recommended for tests).
        force: If True, rebuild even if sources haven't changed.
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        Build result with success status, duration, and any error messages.

    Example:
        # Clean build before running tests (recommended)
        build_star(target="core", clean=True)

        # Quick incremental build during development
        build_star(target="core")
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        config = get_config()
        repo_root = Path(config.paths.repo_root)
        result = _build_target(repo_root, target, clean=clean, force=force)
        return result
    except Exception as e:
        return ErrorResponse(
            code="BUILD_FAILED",
            message=str(e),
        ).model_dump()


@mcp.tool()
def check_build_status(
    target: str = "core",
    auth_token: Optional[str] = None,
) -> dict:
    """Check if a build target needs rebuilding.

    Uses source file hashing to detect changes since last build.
    This helps avoid stale binary issues.

    Args:
        target: Build target to check ('core', 'flex', 'slam', 'feature-tools').
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        Status including:
        - needs_rebuild: Whether sources have changed
        - reason: Explanation of why rebuild is/isn't needed
        - last_build: Timestamp of last successful build
        - was_clean_build: Whether last build was a clean build
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        config = get_config()
        repo_root = Path(config.paths.repo_root)

        needs_build, reason = _needs_rebuild(repo_root, target)
        state = _load_build_state(repo_root)
        target_state = state.get(target, {})

        return {
            "target": target,
            "needs_rebuild": needs_build,
            "reason": reason,
            "last_build": target_state.get("build_time"),
            "was_clean_build": target_state.get("clean_build"),
        }
    except Exception as e:
        return ErrorResponse(
            code="STATUS_CHECK_FAILED",
            message=str(e),
        ).model_dump()


@mcp.tool()
def ensure_fresh_build(
    target: str = "core",
    auth_token: Optional[str] = None,
) -> dict:
    """Ensure a fresh, clean build exists (recommended before test suites).

    This ALWAYS runs 'make clean' first to prevent stale binary issues
    that can cause hard-to-debug segfaults.

    Use this tool before running any test suite to ensure consistent results.

    Args:
        target: Build target - 'core' (STAR), 'flex', 'slam', or 'feature-tools'.
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        Build result with success status and timing.
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        config = get_config()
        repo_root = Path(config.paths.repo_root)
        result = _build_target(repo_root, target, clean=True, force=True)
        return result
    except Exception as e:
        return ErrorResponse(
            code="BUILD_FAILED",
            message=str(e),
        ).model_dump()


# --- Workflow Tools ---


@mcp.tool()
def list_workflows(auth_token: Optional[str] = None) -> dict:
    """List supported workflow templates.

    Returns structured workflows with metadata, separate from the scripts allowlist.

    Args:
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns:
        ListWorkflowsResponse with workflow summaries.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _list_workflows(authenticated=is_authenticated(auth_token))
        payload = result.model_dump()
        payload["agent_protocol"] = get_config().agent_protocol
        return payload
    except Exception as e:
        return ErrorResponse(
            code="WORKFLOW_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def describe_workflow(
    workflow_id: str,
    auth_token: Optional[str] = None,
) -> dict:
    """Return full workflow metadata including stages and parameter groups.

    Args:
        workflow_id: Workflow identifier (e.g. "ucsf_star_suite_production").
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns:
        DescribeWorkflowResponse with stages, parameter groups, and caveats.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _describe_workflow(workflow_id, authenticated=is_authenticated(auth_token))
        payload = result.model_dump()
        payload["agent_protocol"] = get_config().agent_protocol
        return payload
    except Exception as e:
        return ErrorResponse(
            code="WORKFLOW_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def get_workflow_scripts(
    workflow_id: str,
    auth_token: Optional[str] = None,
) -> dict:
    """Return the scripts composing a workflow with provenance metadata.

    Provides entry script path, helper script paths, existence checks,
    and repo provenance (git commit, remote). Intended for downstream
    script-backed encoders that need to know which scripts a workflow uses.

    Args:
        workflow_id: Workflow identifier (e.g. "ucsf_star_suite_production").
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns:
        GetWorkflowScriptsResponse with script details and provenance.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _get_workflow_scripts(
            workflow_id, authenticated=is_authenticated(auth_token)
        )
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="WORKFLOW_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def get_workflow_parameter_schema(
    workflow_id: str,
    auth_token: Optional[str] = None,
) -> dict:
    """Return the machine-readable parameter schema for a workflow.

    Includes parameter definitions, ordered groups, constraints,
    and mutual exclusion / dependency rules.

    Args:
        workflow_id: Workflow identifier.
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns:
        WorkflowParameterSchemaResponse with parameters, groups, and constraints.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _get_workflow_parameter_schema(
            workflow_id, authenticated=is_authenticated(auth_token)
        )
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="WORKFLOW_ERROR",
            message=str(e),
        ).model_dump()


@mcp.tool()
def validate_workflow_parameters(
    workflow_id: str,
    params: dict,
    check_paths: bool = True,
    auth_token: Optional[str] = None,
) -> dict:
    """Validate structured params for a workflow without executing anything.

    Checks required params, types, enum values, path existence,
    mutual exclusion, and dependency constraints.

    Args:
        workflow_id: Workflow identifier.
        params: Parameter dict (e.g. {"samples": "EBs1_1,EBs1_2", "threads": 8}).
        check_paths: Whether to verify file/directory existence (default true).
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        ValidateWorkflowResponse with valid flag, normalized params, warnings, errors.
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _validate_workflow_parameters(workflow_id, params, check_paths=check_paths)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="WORKFLOW_VALIDATION_FAILED",
            message=str(e),
        ).model_dump()


@mcp.tool()
def render_workflow_command(
    workflow_id: str,
    params: dict,
    auth_token: Optional[str] = None,
) -> dict:
    """Render validated params into the actual command invocation.

    Does not execute anything. Returns the argv array, a shell preview,
    environment overrides, and the output root if derivable.

    Args:
        workflow_id: Workflow identifier.
        params: Parameter dict (should be validated first).
        auth_token: Authentication token (required if server has auth configured).

    Returns:
        RenderWorkflowResponse with argv, shell_preview, env_overrides.
    """
    auth_error = check_auth(auth_token)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _render_workflow_command(workflow_id, params)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="WORKFLOW_RENDER_FAILED",
            message=str(e),
        ).model_dump()


@mcp.tool()
def scaffold_workflow_schema(
    script_path: str,
    auth_token: Optional[str] = None,
) -> dict:
    """Parse a shell script and propose a draft workflow schema YAML.

    Best-effort extraction from while/case flag-parsing patterns.
    Returns a draft YAML string that should be reviewed, refined,
    and committed to the repo before it becomes active.

    Does NOT modify the running server configuration.

    Args:
        script_path: Path to the shell script (absolute or relative to repo root).
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns:
        ScaffoldWorkflowResponse with draft YAML, parameter count, and review notes.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _scaffold_workflow_schema(script_path)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="SCAFFOLD_FAILED",
            message=str(e),
        ).model_dump()


@mcp.tool()
def validate_draft_workflow_schema(
    draft_yaml: str,
    auth_token: Optional[str] = None,
) -> dict:
    """Validate a draft workflow schema YAML against the WorkflowSchema model.

    Checks structural validity, semantic consistency (flag_order references,
    constraint params, duplicate names), and flags TODO placeholders.

    Does NOT load the schema into the running server configuration.

    Args:
        draft_yaml: YAML string to validate.
        auth_token: Authentication token (optional if public_discovery is enabled).

    Returns:
        ValidateDraftWorkflowResponse with valid flag, errors, warnings,
        and parameter count.
    """
    auth_error = check_auth(auth_token, is_discovery=True)
    if auth_error:
        return auth_error.model_dump()

    try:
        result = _validate_draft_workflow_schema(draft_yaml)
        return result.model_dump()
    except Exception as e:
        return ErrorResponse(
            code="VALIDATION_FAILED",
            message=str(e),
        ).model_dump()


def main():
    """Main entry point for the MCP server."""
    parser = argparse.ArgumentParser(description="STAR-suite MCP Server")
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Path to config.yaml (default: mcp_server/config.yaml)",
    )
    parser.add_argument(
        "--transport",
        choices=["stdio", "http"],
        default=None,
        help="Transport mode (overrides config)",
    )
    parser.add_argument(
        "--host",
        default=None,
        help="Host to bind to (overrides config, HTTP only)",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=None,
        help="Port to bind to (overrides config, HTTP only)",
    )

    args = parser.parse_args()

    # Load configuration
    try:
        config = load_config(args.config)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading config: {e}", file=sys.stderr)
        sys.exit(1)

    # Apply CLI overrides
    transport = args.transport or config.server.transport
    host = args.host or config.server.host
    port = args.port or config.server.port

    # Run the server
    global _trust_local
    if transport == "stdio":
        _trust_local = True
        mcp.run(transport="stdio")
    else:
        # HTTP mode (serve both streamable-http at "/" and SSE at "/sse"/"/messages")
        app = build_http_app()
        uvicorn.run(app, host=host, port=port, log_level="info")


if __name__ == "__main__":
    main()
