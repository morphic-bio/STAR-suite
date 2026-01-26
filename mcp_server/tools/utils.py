"""Utility functions for MCP tools."""

import os
import random
import shutil
import string
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from ..config import get_config


def generate_run_id() -> str:
    """Generate a unique run ID.

    Format: YYYYMMDD_HHMMSS_<rand4>
    Example: 20260126_143052_a7f3
    """
    now = datetime.now(timezone.utc)
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    rand_suffix = "".join(random.choices(string.ascii_lowercase + string.digits, k=4))
    return f"{timestamp}_{rand_suffix}"


def is_path_allowed(path: Path | str, trusted_roots: Optional[list[str]] = None) -> bool:
    """Check if a path is within trusted roots.

    Args:
        path: Path to check.
        trusted_roots: List of trusted root paths. If None, uses config.

    Returns:
        True if path is within a trusted root, False otherwise.
    """
    if trusted_roots is None:
        config = get_config()
        trusted_roots = config.trusted_roots

    # Resolve to absolute path
    path = Path(path).resolve()

    for root in trusted_roots:
        root_path = Path(root).resolve()
        try:
            path.relative_to(root_path)
            return True
        except ValueError:
            continue

    return False


def validate_path(
    path: Path | str,
    must_exist: bool = False,
    must_be_file: bool = False,
    must_be_dir: bool = False,
    must_be_readable: bool = False,
    must_be_writable: bool = False,
) -> tuple[bool, Optional[str]]:
    """Validate a path against various criteria.

    Args:
        path: Path to validate.
        must_exist: Path must exist.
        must_be_file: Path must be a file.
        must_be_dir: Path must be a directory.
        must_be_readable: Path must be readable.
        must_be_writable: Path must be writable.

    Returns:
        Tuple of (is_valid, error_message).
    """
    path = Path(path)

    # Check trusted roots first
    if not is_path_allowed(path):
        return False, f"Path outside trusted roots: {path}"

    if must_exist and not path.exists():
        return False, f"Path does not exist: {path}"

    if must_be_file and path.exists() and not path.is_file():
        return False, f"Path is not a file: {path}"

    if must_be_dir and path.exists() and not path.is_dir():
        return False, f"Path is not a directory: {path}"

    if must_be_readable and path.exists() and not os.access(path, os.R_OK):
        return False, f"Path is not readable: {path}"

    if must_be_writable:
        if path.exists():
            if not os.access(path, os.W_OK):
                return False, f"Path is not writable: {path}"
        else:
            # Find the first existing ancestor and check if it's writable
            ancestor = path.parent
            while not ancestor.exists() and ancestor.parent != ancestor:
                ancestor = ancestor.parent

            if not ancestor.exists():
                return False, f"No existing ancestor directory found for: {path}"

            if not os.access(ancestor, os.W_OK):
                return False, f"Ancestor directory is not writable: {ancestor}"

    return True, None


def get_disk_space_gb(path: Path | str) -> float:
    """Get available disk space in GB for a path.

    Args:
        path: Path to check (uses the mount point of this path).

    Returns:
        Available space in GB.
    """
    path = Path(path)

    # Find an existing parent to check
    check_path = path
    while not check_path.exists() and check_path.parent != check_path:
        check_path = check_path.parent

    if not check_path.exists():
        check_path = Path("/")

    usage = shutil.disk_usage(check_path)
    return usage.free / (1024**3)


def get_directory_size(path: Path | str) -> int:
    """Get total size of a directory in bytes.

    Args:
        path: Directory path.

    Returns:
        Total size in bytes.
    """
    path = Path(path)
    if not path.is_dir():
        return 0

    total = 0
    for entry in path.rglob("*"):
        if entry.is_file():
            try:
                total += entry.stat().st_size
            except (OSError, PermissionError):
                pass
    return total


def find_binary(
    name: str,
    search_paths: Optional[list[str]] = None,
    enforce_trusted_roots: bool = True,
) -> Optional[Path]:
    """Find a binary in search paths or PATH.

    Args:
        name: Binary name.
        search_paths: Additional paths to search. If None, uses config.
        enforce_trusted_roots: If True, only return binaries within trusted roots.

    Returns:
        Path to binary if found (and within trusted roots if enforced), None otherwise.
    """
    config = get_config()

    # Build search list
    paths_to_check: list[Path] = []

    # Check configured paths for this binary first
    for binary_config in config.required_binaries:
        if binary_config.name == name:
            for p in binary_config.paths:
                # Handle relative paths (relative to repo root)
                if not Path(p).is_absolute():
                    p = str(config.paths.repo_root / p)
                paths_to_check.append(Path(p))

    # Add custom search paths
    if search_paths:
        paths_to_check.extend(Path(p) for p in search_paths)

    # Check each configured path (these are explicitly trusted via config)
    for path in paths_to_check:
        if path.exists() and os.access(path, os.X_OK):
            # Configured paths must still be within trusted roots
            if enforce_trusted_roots and not is_path_allowed(path):
                continue
            return path

    # Fall back to which/PATH - but enforce trusted roots
    result = shutil.which(name)
    if result:
        result_path = Path(result)
        if enforce_trusted_roots:
            # PATH-resolved binaries must be within trusted roots
            if is_path_allowed(result_path):
                return result_path
            # Binary found but outside trusted roots - don't return it
            return None
        return result_path

    return None


def get_run_log_dir(run_id: str) -> Path:
    """Get the log directory for a run.

    Args:
        run_id: Run identifier.

    Returns:
        Path to log directory.
    """
    config = get_config()
    date_str = run_id.split("_")[0]  # Extract YYYYMMDD
    return config.paths.artifact_log_root / f"mcp_runs_{date_str}" / run_id


def ensure_run_log_dir(run_id: str) -> Path:
    """Ensure the log directory for a run exists.

    Args:
        run_id: Run identifier.

    Returns:
        Path to log directory.
    """
    log_dir = get_run_log_dir(run_id)
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir
