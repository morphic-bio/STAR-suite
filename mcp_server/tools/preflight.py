"""Preflight validation for MCP server."""

from pathlib import Path
from typing import Any, Optional

from ..config import get_config
from ..schemas.config import ScriptConfig
from ..schemas.responses import PreflightCheck, PreflightResponse
from ..schemas.run_config import RunConfig
from .utils import (
    find_binary,
    get_disk_space_gb,
    is_path_allowed,
    validate_path,
)


def check_script_allowed(run_config: RunConfig) -> PreflightCheck:
    """Check if the requested script is in the allowlist.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()
    script = config.get_script(run_config.script)

    if script is None:
        allowed_scripts = [s.name for s in config.scripts]
        return PreflightCheck(
            name="script_allowed",
            passed=False,
            message=f"Script '{run_config.script}' is not in the allowlist",
            details={
                "requested": run_config.script,
                "allowed": allowed_scripts,
            },
        )

    # Check if script is marked as runnable
    if not script.runnable:
        return PreflightCheck(
            name="script_allowed",
            passed=False,
            message=f"Script '{run_config.script}' is not currently runnable",
            details={"reason": "Script marked as not runnable in config"},
        )

    return PreflightCheck(
        name="script_allowed",
        passed=True,
        details={"script": script.name, "module": script.module},
    )


def check_script_path(run_config: RunConfig) -> PreflightCheck:
    """Check if the script file exists and is within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()
    script = config.get_script(run_config.script)

    if script is None:
        return PreflightCheck(
            name="script_path",
            passed=False,
            message="Script not found in config",
        )

    # Resolve script path
    script_path = Path(script.path)
    if not script_path.is_absolute():
        script_path = config.paths.repo_root / script_path

    # For build commands (like 'make'), check if the command exists
    if script.path in ("make", "cmake", "ninja"):
        # find_binary enforces trusted roots by default
        binary = find_binary(script.path, enforce_trusted_roots=True)
        if binary:
            return PreflightCheck(
                name="script_path",
                passed=True,
                details={"path": str(binary), "type": "build_command"},
            )
        return PreflightCheck(
            name="script_path",
            passed=False,
            message=f"Build command '{script.path}' not found in trusted roots",
        )

    # Check script path is within trusted roots
    if not is_path_allowed(script_path):
        return PreflightCheck(
            name="script_path",
            passed=False,
            message=f"Script path outside trusted roots: {script_path}",
            details={"path": str(script_path)},
        )

    # Check script file exists
    if not script_path.exists():
        return PreflightCheck(
            name="script_path",
            passed=False,
            message=f"Script file not found: {script_path}",
            details={"path": str(script_path)},
        )

    # Check script is a file, not a directory
    if not script_path.is_file():
        return PreflightCheck(
            name="script_path",
            passed=False,
            message=f"Script path is not a file: {script_path}",
            details={"path": str(script_path)},
        )

    return PreflightCheck(
        name="script_path",
        passed=True,
        details={"path": str(script_path)},
    )


def check_working_dir(run_config: RunConfig) -> PreflightCheck:
    """Check if the script's working_dir is within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()
    script = config.get_script(run_config.script)

    if script is None:
        return PreflightCheck(
            name="working_dir",
            passed=False,
            message="Script not found in config",
        )

    # If no working_dir specified, defaults to repo_root which is trusted
    if not script.working_dir:
        return PreflightCheck(
            name="working_dir",
            passed=True,
            message="No working_dir specified (defaults to repo_root)",
            details={"path": str(config.paths.repo_root)},
        )

    # Resolve working_dir path
    working_dir = Path(script.working_dir)
    if not working_dir.is_absolute():
        working_dir = config.paths.repo_root / working_dir

    # Check working_dir is within trusted roots
    if not is_path_allowed(working_dir):
        return PreflightCheck(
            name="working_dir",
            passed=False,
            message=f"Script working_dir outside trusted roots: {working_dir}",
            details={"path": str(working_dir), "script": script.name},
        )

    # Check working_dir exists
    if not working_dir.exists():
        return PreflightCheck(
            name="working_dir",
            passed=False,
            message=f"Script working_dir does not exist: {working_dir}",
            details={"path": str(working_dir), "script": script.name},
        )

    # Check working_dir is a directory
    if not working_dir.is_dir():
        return PreflightCheck(
            name="working_dir",
            passed=False,
            message=f"Script working_dir is not a directory: {working_dir}",
            details={"path": str(working_dir), "script": script.name},
        )

    return PreflightCheck(
        name="working_dir",
        passed=True,
        details={"path": str(working_dir)},
    )


def check_dataset(run_config: RunConfig) -> PreflightCheck:
    """Check if the requested dataset exists and is within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    if not run_config.dataset_id:
        return PreflightCheck(
            name="dataset",
            passed=True,
            message="No dataset specified",
        )

    config = get_config()
    dataset = config.get_dataset(run_config.dataset_id)

    if dataset is None:
        available = [d.id for d in config.datasets]
        return PreflightCheck(
            name="dataset",
            passed=False,
            message=f"Dataset '{run_config.dataset_id}' not found",
            details={"requested": run_config.dataset_id, "available": available},
        )

    dataset_path = Path(dataset.path)

    # Check dataset path is within trusted roots
    if not is_path_allowed(dataset_path):
        return PreflightCheck(
            name="dataset",
            passed=False,
            message=f"Dataset path outside trusted roots: {dataset.path}",
            details={"dataset_id": dataset.id, "path": str(dataset.path)},
        )

    # Check dataset path exists
    if not dataset_path.exists():
        return PreflightCheck(
            name="dataset",
            passed=False,
            message=f"Dataset path does not exist: {dataset.path}",
            details={"dataset_id": dataset.id, "path": str(dataset.path)},
        )

    return PreflightCheck(
        name="dataset",
        passed=True,
        details={"dataset_id": dataset.id, "path": str(dataset.path)},
    )


def check_output_directory(run_config: RunConfig) -> PreflightCheck:
    """Check if the output directory is writable.

    This check validates that:
    1. The path is within trusted roots
    2. If the directory exists, it's writable
    3. If it doesn't exist, a writable parent exists in the chain

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()

    # Determine output directory
    if run_config.out_dir:
        out_dir = Path(run_config.out_dir)
    else:
        out_dir = config.paths.temp_root / "mcp_runs"

    # Check path is within trusted roots
    if not is_path_allowed(out_dir):
        return PreflightCheck(
            name="output_directory",
            passed=False,
            message=f"Output directory outside trusted roots: {out_dir}",
            details={"path": str(out_dir)},
        )

    # Check if directory exists or can be created
    valid, error = validate_path(out_dir, must_be_writable=True)
    if not valid:
        return PreflightCheck(
            name="output_directory",
            passed=False,
            message=error or f"Output directory not writable: {out_dir}",
            details={"path": str(out_dir)},
        )

    # Additional check: verify we can actually create missing parent chain
    # validate_path finds an existing ancestor, but we should report how deep that is
    if not out_dir.exists():
        ancestor = out_dir.parent
        depth = 1
        while not ancestor.exists() and ancestor.parent != ancestor:
            ancestor = ancestor.parent
            depth += 1

        if depth > 3:
            # Warn if the directory chain is deep - may indicate a typo
            return PreflightCheck(
                name="output_directory",
                passed=True,
                message=f"Output directory will be created ({depth} levels deep)",
                details={
                    "path": str(out_dir),
                    "existing_ancestor": str(ancestor),
                    "depth": depth,
                    "warning": "Deep directory chain - verify path is correct",
                },
            )

    return PreflightCheck(
        name="output_directory",
        passed=True,
        details={"path": str(out_dir), "exists": out_dir.exists()},
    )


def check_binaries(run_config: RunConfig) -> PreflightCheck:
    """Check if required binaries are present.

    Skips this check for scripts in the "build" module since those scripts
    are meant to build the required binaries.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()

    # Get the script to check its module
    script = config.get_script(run_config.script)

    # Skip binary checks for build scripts - they create the binaries
    if script and script.module == "build":
        return PreflightCheck(
            name="binaries_present",
            passed=True,
            message="Skipped for build scripts",
            details={"reason": "Build scripts create binaries, not require them"},
        )

    missing: list[str] = []
    found: dict[str, str] = {}

    for binary_config in config.required_binaries:
        binary_path = find_binary(binary_config.name)
        if binary_path:
            found[binary_config.name] = str(binary_path)
        else:
            missing.append(binary_config.name)

    if missing:
        return PreflightCheck(
            name="binaries_present",
            passed=False,
            message=f"Required binaries not found: {', '.join(missing)}",
            details={"missing": missing, "found": found},
        )

    return PreflightCheck(
        name="binaries_present",
        passed=True,
        details=found,
    )


def check_disk_space(run_config: RunConfig) -> PreflightCheck:
    """Check if there is sufficient disk space.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()

    # Determine output directory for space check
    if run_config.out_dir:
        check_path = Path(run_config.out_dir)
    else:
        check_path = config.paths.temp_root

    # Get available space
    available_gb = get_disk_space_gb(check_path)

    # Determine required space based on script's test suite
    script = config.get_script(run_config.script)
    required_gb = config.disk_space.default_min_gb

    if script:
        # Find the test suite this script belongs to
        for suite in config.test_suites:
            if script.name in suite.scripts:
                if suite.disk_space_min_gb is not None:
                    required_gb = suite.disk_space_min_gb
                break

    if available_gb < required_gb:
        return PreflightCheck(
            name="disk_space",
            passed=False,
            message=f"Insufficient disk space: {available_gb:.1f}GB available, {required_gb}GB required",
            details={
                "available_gb": round(available_gb, 1),
                "required_gb": required_gb,
                "check_path": str(check_path),
            },
        )

    return PreflightCheck(
        name="disk_space",
        passed=True,
        details={
            "available_gb": round(available_gb, 1),
            "required_gb": required_gb,
        },
    )


def check_fixtures(run_config: RunConfig) -> PreflightCheck:
    """Check if required fixtures are present and within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    config = get_config()
    script = config.get_script(run_config.script)

    if script is None or not script.fixtures:
        return PreflightCheck(
            name="fixtures_present",
            passed=True,
            message="No fixtures required",
        )

    missing: list[str] = []
    found: list[str] = []
    outside_roots: list[str] = []

    for fixture_path in script.fixtures:
        path = Path(fixture_path)

        # Check fixture is within trusted roots first
        if not is_path_allowed(path):
            outside_roots.append(fixture_path)
        elif path.exists():
            found.append(fixture_path)
        else:
            missing.append(fixture_path)

    # Report paths outside trusted roots as a security error
    if outside_roots:
        return PreflightCheck(
            name="fixtures_present",
            passed=False,
            message=f"Fixture paths outside trusted roots: {len(outside_roots)}",
            details={"outside_roots": outside_roots, "found": found, "missing": missing},
        )

    if missing:
        return PreflightCheck(
            name="fixtures_present",
            passed=False,
            message=f"Required fixtures not found: {len(missing)} missing",
            details={"missing": missing, "found": found},
        )

    return PreflightCheck(
        name="fixtures_present",
        passed=True,
        details={"fixtures": found},
    )


def check_paths_valid(run_config: RunConfig) -> PreflightCheck:
    """Check that all specified paths are within trusted roots.

    Args:
        run_config: Run configuration.

    Returns:
        PreflightCheck result.
    """
    invalid_paths: list[str] = []

    # Check out_dir if specified
    if run_config.out_dir:
        if not is_path_allowed(run_config.out_dir):
            invalid_paths.append(run_config.out_dir)

    if invalid_paths:
        return PreflightCheck(
            name="paths_valid",
            passed=False,
            message="Paths outside trusted roots",
            details={"invalid_paths": invalid_paths},
        )

    return PreflightCheck(
        name="paths_valid",
        passed=True,
    )


def run_preflight(run_config: RunConfig) -> PreflightResponse:
    """Run all preflight checks.

    Args:
        run_config: Run configuration to validate.

    Returns:
        PreflightResponse with all check results.
    """
    checks: list[PreflightCheck] = []
    warnings: list[str] = []
    errors: list[str] = []

    # Run all checks
    check_functions = [
        check_script_allowed,
        check_script_path,
        check_working_dir,  # Validate script's working_dir is within trusted roots
        check_dataset,
        check_output_directory,
        check_binaries,
        check_disk_space,
        check_fixtures,
        check_paths_valid,
    ]

    for check_fn in check_functions:
        try:
            result = check_fn(run_config)
            checks.append(result)
            if not result.passed:
                errors.append(result.message or f"{result.name} check failed")
        except Exception as e:
            checks.append(
                PreflightCheck(
                    name=check_fn.__name__.replace("check_", ""),
                    passed=False,
                    message=f"Check error: {str(e)}",
                )
            )
            errors.append(f"{check_fn.__name__}: {str(e)}")

    # Determine overall validity
    valid = all(check.passed for check in checks)

    return PreflightResponse(
        valid=valid,
        checks=checks,
        warnings=warnings,
        errors=errors,
    )
