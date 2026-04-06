"""Configuration loading and management for MCP server."""

import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import yaml

from .schemas.config import MCPConfig
from .schemas.workflow import WorkflowSchema


# Global config instance
_config: Optional[MCPConfig] = None
_config_path: Optional[Path] = None
_config_loaded_at: Optional[datetime] = None
_workflow_schemas: dict[str, WorkflowSchema] = {}


def expand_env_vars(value: Any) -> Any:
    """Recursively expand environment variables in config values.

    Supports ${VAR_NAME} syntax. Missing vars expand to empty string.
    """
    if isinstance(value, str):
        # Match ${VAR_NAME} pattern
        pattern = re.compile(r"\$\{([^}]+)\}")

        def replacer(match: re.Match) -> str:
            var_name = match.group(1)
            return os.environ.get(var_name, "")

        return pattern.sub(replacer, value)
    elif isinstance(value, dict):
        return {k: expand_env_vars(v) for k, v in value.items()}
    elif isinstance(value, list):
        return [expand_env_vars(item) for item in value]
    else:
        return value


def _load_workflow_schemas(config: MCPConfig, repo_root: Path) -> dict[str, WorkflowSchema]:
    """Load all workflow schema YAML files referenced in config.

    Args:
        config: The loaded MCPConfig.
        repo_root: Repository root for resolving relative paths.

    Returns:
        Dict mapping workflow id to WorkflowSchema.
    """
    schemas: dict[str, WorkflowSchema] = {}
    for wf_cfg in config.workflows:
        schema_path = Path(wf_cfg.schema_file)
        if not schema_path.is_absolute():
            schema_path = repo_root / schema_path

        if not schema_path.exists():
            raise FileNotFoundError(
                f"Workflow schema file not found for '{wf_cfg.id}': {schema_path}"
            )

        with open(schema_path) as f:
            raw_schema = yaml.safe_load(f)

        schemas[wf_cfg.id] = WorkflowSchema(**raw_schema)

    return schemas


def load_config(config_path: Optional[Path] = None) -> MCPConfig:
    """Load configuration from YAML file.

    Args:
        config_path: Path to config file. If None, uses default location.

    Returns:
        MCPConfig instance.

    Raises:
        FileNotFoundError: If config file doesn't exist.
        ValueError: If config is invalid.
    """
    global _config, _config_path, _config_loaded_at, _workflow_schemas

    if config_path is None:
        # Default to config.yaml in the same directory as this file
        config_path = Path(__file__).parent / "config.yaml"

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        raw_config = yaml.safe_load(f)

    # Expand environment variables
    expanded_config = expand_env_vars(raw_config)

    # Parse into Pydantic model
    config = MCPConfig(**expanded_config)

    # Load workflow schemas
    repo_root = config.paths.repo_root
    _workflow_schemas = _load_workflow_schemas(config, repo_root)

    # Store globally
    _config = config
    _config_path = config_path
    _config_loaded_at = datetime.now(timezone.utc)

    return config


def get_config() -> MCPConfig:
    """Get the current configuration.

    Returns:
        MCPConfig instance.

    Raises:
        RuntimeError: If config hasn't been loaded yet.
    """
    if _config is None:
        raise RuntimeError("Configuration not loaded. Call load_config() first.")
    return _config


def get_workflow_schemas() -> dict[str, WorkflowSchema]:
    """Get all loaded workflow schemas."""
    return _workflow_schemas


def get_workflow_schema(workflow_id: str) -> Optional[WorkflowSchema]:
    """Get a workflow schema by ID."""
    return _workflow_schemas.get(workflow_id)


def get_config_path() -> Optional[Path]:
    """Get the path to the loaded config file."""
    return _config_path


def get_config_loaded_at() -> Optional[datetime]:
    """Get when the config was last loaded."""
    return _config_loaded_at


def reload_config() -> tuple[MCPConfig, datetime]:
    """Reload configuration from the same file.

    Returns:
        Tuple of (new config, load time).

    Raises:
        RuntimeError: If no config was previously loaded.
        ValueError: If new config is invalid (old config is retained).
    """
    global _config, _config_loaded_at, _workflow_schemas

    if _config_path is None:
        raise RuntimeError("No config previously loaded. Call load_config() first.")

    # Load into a new variable first to validate
    old_config = _config
    old_schemas = _workflow_schemas
    try:
        new_config = load_config(_config_path)
        return new_config, _config_loaded_at  # type: ignore
    except Exception:
        # Restore old config on failure
        _config = old_config
        _workflow_schemas = old_schemas
        raise
