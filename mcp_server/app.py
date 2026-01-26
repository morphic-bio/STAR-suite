"""FastMCP server entry point for STAR-suite MCP server."""

import argparse
import sys
from pathlib import Path
from typing import Optional

from fastmcp import FastMCP

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


# Create the MCP server
mcp = FastMCP("star-suite")


def check_auth(token: Optional[str], is_discovery: bool = False) -> Optional[ErrorResponse]:
    """Check if the provided token is valid.

    Args:
        token: Bearer token from request.
        is_discovery: If True, this is a discovery tool that may be public.

    Returns:
        ErrorResponse if auth fails, None if auth passes.
    """
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
    if transport == "stdio":
        mcp.run(transport="stdio")
    else:
        # HTTP/SSE mode
        mcp.run(transport="sse", host=host, port=port)


if __name__ == "__main__":
    main()
